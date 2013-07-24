#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include <assert.h>
#include "kvec.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };

typedef kvec_t(int) ivec_t;

typedef struct {
	ivec_t list;
} sites_t;

int binary_search_abs(int A[], int key, int imin, int imax)
{
	while (imin < imax) {
		int imid = (imin + imax) >> 1;
		if (abs(A[imid]) < key) imin = imid + 1;
		else imax = imid;
	}
	return imax == imin && abs(A[imin]) == key? imin : -1;
}

int main_sites(int argc, char *argv[])
{
	int c;
	kvec_t(sites_t) sites = {0, 0, 0};

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: alktk sites <in.alk> <in.srt.list>\n\n");
		fprintf(stderr, "Note: Each line of <in.srt.list> consists of 1-based chr index, position and\n");
		fprintf(stderr, "      a zero or one indicating whether the reference allele is ancestral.\n\n");
		return 1;
	}
	{ // read the sites
		kstring_t str = {0, 0, 0};
		int dret;
		gzFile fp;
		kstream_t *ks;

		fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p, *q, *r;
			int n = 0, chr = -1, pos = -1, anc = -1, het = -1, flip = -1, is_digit = -1, tmp;
			p = str.s;
			for (q = str.s + 1;; ++q) {
				if (isspace(*q) || *q == 0) {
					if (n == 0) chr = strtol(p, &r, 10);
					else if (n == 1) pos = strtol(p, &r, 10);
					else if (n == 2) {
						anc = seq_nt16_table[(uint8_t)*p];
						is_digit = isdigit(*p);
						flip = is_digit? (*p != '0') : (bitcnt_table[anc] == 2);
					} else if (n == 3) het = seq_nt16_table[(uint8_t)*p];
					++n; p = q + 1;
					if (*q == 0) break;
				}
			}
			if (n < 3 || chr < 0 || pos < 0 || is_digit < 0) continue;
			tmp = bitcnt_table[anc];
			--chr;
			if (!is_digit && bitcnt_table[anc] > 2) continue;
			if (n >= 4 && !is_digit && bitcnt_table[het|anc] != 2) continue;
			if (chr >= sites.m) {
				int i, oldm = sites.m;
				kv_resize(sites_t, sites, chr + 1);
				for (i = oldm; i < sites.m; ++i)
					memset(&sites.a[i].list, 0, sizeof(kvec_t(int)));
			}
			if (chr >= sites.n) sites.n = chr + 1;
			tmp = flip? -pos : pos;
			kv_push(int, sites.a[chr].list, tmp);
			//fprintf(stderr, "%s | %d,%d,%d,%d | %d\n", str.s, chr, pos, anc, het, flip);
		}
		ks_destroy(ks);
		gzclose(fp);
		free(str.s);
	}
	{ // read alk
		gzFile fp, fpout;
		int32_t M, size, ret;
		uint8_t *rec;

		fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
		fpout = gzdopen(fileno(stdout), "w1");
		gzread(fp, &M, 4);
		gzwrite(fpout, &M, 4);
		size = 4 + 4 + sizeof(double) + (M + 1) * sizeof(float);
		rec = malloc(size);
		while ((ret = gzread(fp, rec, size)) == size) {
			int ret, chr, pos;
			ivec_t *p;

			chr = *(int32_t*)rec;
			pos = *(int32_t*)(rec+4) + 1;
			assert(chr < sites.n);
			p = &sites.a[chr].list;
			ret = binary_search_abs(p->a, pos, 0, (int)p->n - 1);
			if (ret < 0) continue;
			if (p->a[ret] < 0) {
				double *pd = (double*)(rec+8);
				float tmp, *pf = (float*)(rec + 8 + sizeof(double));
				int i, end = (M + 1) >> 1;
				*pd = 1. - *pd;
				for (i = 0; i < end; ++i)
					tmp = pf[i], pf[i] = pf[M-i], pf[M-i] = tmp;
			}
			gzwrite(fpout, rec, size);
		}
		free(rec);
		gzclose(fpout);
		gzclose(fp);
	}
	
	return 0;
}

int main_afs(int argc, char *argv[])
{
	int32_t c, M, i, size, ret;
	kvec_t(double) phi = {0, 0, 0};
	gzFile fp;
	kstream_t *ks;
	uint8_t *rec;
	long double *afs, *tmp;
	int L;

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: alktk afs <in.alk> [init_afs.txt]\n");
		return 1;
	}
	if (optind + 1 < argc) { // read the afs
		kstring_t str = {0, 0, 0};
		int dret;
		long double sum;
		fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p, *q;
			double x;
			for (p = str.s; *p && !isspace(*p); ++p); // find the first space/tab
			if (*p == 0 || *(p+1) == 0) continue;
			x = strtod(p+1, &q);
			kv_push(double, phi, x);
		}
		ks_destroy(ks);
		gzclose(fp);
		for (i = 0, sum = 0.; i <= phi.n; ++i) sum += phi.a[i];
		for (i = 0; i <= phi.n; ++i) phi.a[i] /= sum;
	}
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	gzread(fp, &M, 4);
	if (phi.n) { // then check consistency
		if (phi.n != M + 1) {
			fprintf(stderr, "[E::%s] inconsistent number of chromosomes: %d != %ld\n", __func__, M+1, phi.n);
			free(phi.a);
			gzclose(fp);
			return 1;
		}
	} else { // then initialize a flat prior
		kv_resize(double, phi, M+1);
		phi.n = M + 1;
		for (i = 0; i <= M; ++i)
			phi.a[i] = 1. / (M + 1);
	}
	size = 4 + 4 + sizeof(double) + (M + 1) * sizeof(float);
	rec = malloc(size);
	afs = calloc(M + 1, sizeof(long double));
	tmp = calloc(M + 1, sizeof(long double));
	L = 0;
	while ((ret = gzread(fp, rec, size)) == size) {
		long double sum = 0;
		float *pf = (float*)(rec + 8 + sizeof(double));
		for (i = 0; i <= M; ++i)
			sum += (tmp[i] = phi.a[i] * pf[i]);
		if (sum == 0.0) {
			fprintf(stderr, "%d\t%d\t%f\n", *(int32_t*)rec, *(int32_t*)(rec+4), *(double*)(rec+8));
			continue;
		}
		sum = 1. / sum;
		for (i = 0; i <= M; ++i)
			afs[i] += tmp[i] * sum;
		++L;
	}
	for (i = 0, tmp[0] = 0.; i <= M; ++i) {
		tmp[0] += afs[i] / L;
		printf("%d\t%.9Lf\n", i, afs[i] / L);
	}
	fprintf(stderr, "[M::%s] sum of AFS: %Lf\n", __func__, tmp[0]);
	free(tmp); free(afs); free(rec);
	gzclose(fp);
	return 0;
}

/* a record is packed as follows:
 *
 * +----+----+--------+----+----+-   -+----+
 * |CHR |POS |REF-AF  |L(0)|L(1)| ... |L(M)|
 * +----+----+--------+----+----+-   -+----+
 *
 * where CHR is the index of chromsome in the 1000g reference file and
 * L(k)=Pr{d|AC=k} is the reference allele likelihood. Each "-" represents a
 * byte. The length of a record is "4+4+sizeof(double)+sizeof(float)*(M+1)".
 * On almost all machines, it equals "16+4*(M+1)".
 */
int main_cat(int argc, char *argv[])
{
	gzFile fp, fpout = 0;
	int32_t M = 0, c, k, is_text = 0, is_pos_only = 0;

	while ((c = getopt(argc, argv, "tp")) >= 0)
		if (c == 't') is_text = 1;
		else if (c == 'p') is_pos_only = 1;
	if (argc == optind) {
		fprintf(stderr, "Usage: alktk cat [-tp] <in1.alk> <in2.alk> [...]\n");
		return 1;
	}
	if (!is_text) fpout = gzdopen(fileno(stdout), "w");

	for (k = optind; k < argc; ++k) {
		uint8_t *rec;
		int32_t ret, size, i;
		fp = gzopen(argv[k], "r");
		if (gzread(fp, &size, 4) != 4) continue; // read the number of chromosomes
		if (k > optind) {
			if (M != size) {
				fprintf(stderr, "[E::%s] `%s' has different number of chromosomes\n", __func__, argv[k]);
				return 1;
			}
		} else {
			M = size;
			if (!is_text) gzwrite(fpout, &M, 4);
		}
		size = 4 + 4 + sizeof(double) + (M + 1) * sizeof(float); // the size of a record
		rec = malloc(size);
		while ((ret = gzread(fp, rec, size)) == size) { // read a record
			if (!fpout) {
				// unpack the record with explicit pointer conversions
				printf("%d\t%d\t%f", *(int32_t*)rec, *(int32_t*)(rec+4), *(double*)(rec+8));
				if (!is_pos_only) {
					float *pf;
					pf = (float*)(rec + 8 + sizeof(double)); // alias to the L(k) array
					for (i = 0; i <= M; ++i) printf("\t%f", pf[i]);
				}
				putchar('\n');
			} else gzwrite(fpout, rec, size);
		}
		if (ret && ret != size)
			fprintf(stderr, "[W::%s] `%s' is truncated\n", __func__, argv[k]);
		free(rec);
		gzclose(fp);
	}
	gzclose(fpout);
	return 0;
}

int main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: alktk (Toolkit for processing .alk files)\n");
		fprintf(stderr, "Version: 0.0-alpha\n\n");
		fprintf(stderr, "Usage:   alktk <command> [arguments]\n\n");
		fprintf(stderr, "Command: cat       display/concatenate .alk file(s)\n");
		fprintf(stderr, "         sites     extract given a list of sites and their ancestral alleles\n");
		fprintf(stderr, "         afs       one AFS iteration\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "cat") == 0) main_cat(argc-1, argv+1);
	else if (strcmp(argv[1], "sites") == 0) main_sites(argc-1, argv+1);
	else if (strcmp(argv[1], "afs") == 0) main_afs(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized commad '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
