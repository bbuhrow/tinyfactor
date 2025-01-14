/*
Copyright (c) 2024, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "gmp.h"
#include "ytools.h"
#include "tinyfactor.h"
#include "arith.h"
#include "cmdOptions.h"

#define MAX_FACTORS 128

enum test_algorithm {
	LEHMAN,
	SIQS,
	UECM,
	UECM_PAR,
	TECM,
	TECM_PAR,
	RHO,
	SQUFOF,
	SQUFOF_PAR,
	SQUFOF_RDS,
	SQUFOF_ALPERN,
	FERMAT,
	UNKNOWN_ALG
};

enum test_status {
	success,
	fail,
	warnings
};

typedef struct
{
	// input file and some stats on the inputs
	char* testfile;
	uint32_t min_bits;
	uint32_t max_bits;
	double avg_bits;

	// test to run
	enum test_algorithm alg;
	uint32_t num_to_test;
	uint64_t lcg_state;

	// results
	uint32_t num_ff;
	uint32_t num_pf;
	double parse_time;
	double elasped_time;
} test_t;

enum factor_type {
	FTYPE_PRIME,
	FTYPE_LOGPRIME,
	FTYPE_COMPOSITE,
	FTYPE_RANDOM
};

typedef struct
{
	enum factor_type type;
	int minbits;
	int maxbits;
} gen_t;

int run_test(test_t* test);
int factor_tiny(mpz_t in, mpz_t* out,
	uint64_t* primes, uint64_t nump, uint64_t* prng);

int main(int argc, char** argv)
{
	FILE* in;
	uint64_t* comp, f64;
	uint32_t* f1;
	uint32_t* f2, bits, totBits, minBits, maxBits;
	uint64_t* fout64;
	double avgBits;
	double t_time;
	double f_time;
	int i, j, k, numtest, num = 100000, correct;
	mpz_t gmptmp;
	struct timeval gstart;
	struct timeval gstop;
	int nf;
	int num_files;
	int num_big_files;
	char filenames[30][80];
	char bigfilenames[30][80];
	uint64_t lcg_state;
	options_t *options;
	gmp_randstate_t gmp_randstate;
	uint64_t* primes;
	uint64_t nump;
	uint32_t seed1, seed2;

	get_random_seeds(&seed1, &seed2);
	mpz_init(gmptmp);
	gmp_randinit_default(gmp_randstate);
	mpz_set_ui(gmptmp, (uint64_t)seed2 << 32 + (uint64_t)seed1);
	gmp_randseed(gmp_randstate, gmptmp);
	mpz_urandomb(gmptmp, gmp_randstate, 64);
	lcg_state = mpz_get_ui(gmptmp);

	options = initOpt();
	processOpts(argc, argv, options);

	primes = tiny_soe(10000, &nump);

	if (options->generate)
	{
		FILE* out;
		char tokens[8][32];
		gen_t generators[8];
		int i;
		mpz_t n, tmp1, tmp2, tmp3;
		uint32_t num = options->num;
		
		if (strlen(options->outFile) == 0)
		{
			printf("no output file name specified, defaulting to test_data.dat");
			strcpy(options->outFile, "test_data.dat");
		}

		if (options->num == 0)
		{
			printf("number of tests to generate == 0.  Use -n to set number to generate\n");
			exit(0);
		}

		if (strlen(options->genCmd) == 0)
		{
			printf("Test generation creates lists of input integers to factor by "
				"specifying their consituient factors.\n");
			printf("specify a factor using a type and size, e.g., p50\n");
			printf("factor types: p (uniform random prime), q (log random prime 1/log(N)),\n"
				"c (composite), r (random: prime or composite)\n");
			printf("factor sizes: in bits with optional range e.g., 50 or 50-70\n");
			printf("The number of primes, their types and sizes comprising each test case\n"
				"can be specified in a comma-delimited list.\n");
			printf("e.g., ""p30,p30"" specifies a 60-bit semiprime test case\n");
			printf("and ""p26,p50-70"" specifies a semiprime comprised of one 26-bit prime\n");
			printf("and one 50 to 70-bit prime (uniform random size in range 50 to 70 bits).\n");
			printf("To use 1/log(N) as the distribution instead of uniform, use type 'q'\n");
			printf("Generation will abort if any generated number exceeds 128 bits.\n");
			printf("Complete the test specification by including a number of cases to generate and output file.\n");
			printf("Examples:\n");
			printf("\ttfactor -g ""p30,p30"" -n 10000 -o hard_semiprime_60b.dat\n");
			printf("\ttfactor -g ""c80-100"" -n 10000 -o composite_80b_100b.dat\n");
			printf("\ttfactor -g ""q27-42,q27-42,q27-42"" -n 10000 -o tlp_81b_126b_logweighted.dat\n");

			exit(0);
		}

		mpz_init(n);
		mpz_init(tmp1);
		mpz_init(tmp2);
		mpz_init(tmp3);

		// parse the string.  First break into tokens delimited by ','.
		char *pch = strtok(options->genCmd, ",");
		int numtok = 0;
		while (pch != NULL)
		{
			if (numtok >= 8)
			{
				printf("too many specifiers (8 max)\n");
				exit(0);
			}
			if (strlen(pch) > 32)
			{
				printf("specifier too long (32 character max)\n");
				exit(0);
			}
			strncpy(tokens[numtok++], pch, 32);
			pch = strtok(NULL, ",");
		}

		// create a generator for each token
		for (i = 0; i < numtok; i++)
		{
			switch (tokens[i][0])
			{
			case 'p':
				generators[i].type = FTYPE_PRIME;
				break;
			case 'q':
				generators[i].type = FTYPE_LOGPRIME;
				break;
			case 'c':
				generators[i].type = FTYPE_COMPOSITE;
				break;
			case 'r':
				generators[i].type = FTYPE_RANDOM;
				break;
			default:
				printf("unrecognized factor type '%c', valid types are p (prime), "
					"q (logprime), c (composite), r (random)\n", tokens[i][0]);
				exit(0);
				break;
			}

			int numscanned = sscanf(&tokens[i][1], "%d-%d", 
				&generators[i].minbits, &generators[i].maxbits);
			if (numscanned == 1)
			{
				generators[i].maxbits = generators[i].minbits;
			}
			if (numscanned == 0)
			{
				printf("formatting issue with token %s: expected an integer or a range e.g., 50-60\n", tokens[i]);
				exit(0);
			}

			printf("from token '%s': generated test specification type:%d minbits:%d maxbits:%d\n",
				tokens[i], generators[i].type, generators[i].minbits, generators[i].maxbits);
		}

		// generate the test file
		out = fopen(options->outFile, "w");
		if (out == NULL)
		{
			printf("couldn't open %s for writing\n", options->outFile);
			exit(0);
		}

		char buf[1024];
		char buf2[1024];

		printf("generating %u test cases\n", options->num);
		gettimeofday(&gstart, NULL);

		totBits = 0;
		minBits = 999;
		maxBits = 0;
		for (i = 0; i < num; i++)
		{
			int g;
			int nf = 0;
			mpz_set_ui(n, 1);

			strcpy(buf, "");
			for (g = 0; g < numtok; g++)
			{
				int range = generators[g].maxbits - generators[g].minbits;
				if (range < 0)
				{
					printf("maxbits should be greater than minbits\n");
					exit(0);
				}

				if (generators[g].type == FTYPE_LOGPRIME)
				{

					nf++;
				}
				else
				{
					// all other size distributions are uniform
					int sz = generators[g].minbits;

					if (range > 0)
					{
						sz += (int)floor((double)range * lcg_rand_d(&lcg_state) + 0.5);
					}

					int done = 0;
					do
					{
						mpz_urandomb(tmp2, gmp_randstate, sz);
						mpz_setbit(tmp2, sz - 1);

						if (generators[g].type == FTYPE_PRIME)
						{
							mpz_nextprime(tmp1, tmp2);
							gmp_sprintf(buf2, "%s,%Zd", buf,tmp1);
							strncpy(buf, buf2, 1024);
							nf++;
							done = 1;
						}
						else 
						{
							if (mpz_probab_prime_p(tmp2, 1))
							{
								if (generators[g].type == FTYPE_COMPOSITE)
								{
									done = 0;
								}
								else if (generators[g].type == FTYPE_RANDOM)
								{
									mpz_set(tmp1, tmp2);
									gmp_sprintf(buf2, "%s,%Zd", buf, tmp1);
									strncpy(buf, buf2, 1024);
									nf++;
									done = 1;
								}
							}
							else
							{
								// need to factor this composite to determine factors to list.
								mpz_t factors[MAX_FACTORS];
								mpz_set(tmp1, tmp2);

								int numf = factor_tiny(tmp2, factors, primes, nump, &lcg_state);
								if (numf == 0)
								{
									gmp_printf("error processing composite number %Zd\n", tmp2);
									exit(0);
								}
								
								for (j = 0; j < numf; j++)
								{
									gmp_sprintf(buf2, "%s,%Zd", buf, factors[j]);
									strncpy(buf, buf2, 1024);
								}
								nf += numf;
								done = 1;
							}
						}
					} while (!done);
				}

				
				mpz_mul(n, n, tmp1);
			}

			j = mpz_sizeinbase(n, 2);
			totBits += j;
			if ((uint32_t)j > maxBits)
				maxBits = j;
			if ((uint32_t)j < minBits)
				minBits = j;

			gmp_fprintf(out, "%Zd,%u%s\n", n, nf, buf);
		}
		fclose(out);

		gettimeofday(&gstop, NULL);
		t_time = ytools_difftime(&gstart, &gstop);

		printf("generated %d test cases in %1.2lf seconds\n", num, t_time);
		printf("test case statistics:\n");
		printf("\tmin bits: %u\n", minBits);
		printf("\tmax bits: %u\n", maxBits);
		printf("\tavg bits: %1.2lf\n", (double)totBits / num);
		mpz_clear(tmp1);
		mpz_clear(tmp2);
		mpz_clear(tmp3);
		mpz_clear(n);
	}
	else if (options->test)
	{
		if (strlen(options->testPlan) == 0)
		{
			printf("specify a test plan to run with option -i or --infile\n");
			printf("documentation about how to construct a test plan and \n");
			printf("example plans should be contained in the repository at\n");
			printf("https://github.com/bbuhrow/tinyfactor\n");
			exit(0);
		}

		in = fopen(options->testPlan, "r");
		if (in == NULL)
		{
			printf("could not open %s to read\n", options->testPlan);
		}

		while (~feof(in))
		{
			char buf[1024];
			char testFile[1024];
			char* ptr;
			char* algname[80];
			test_t test;

			test.alg = UNKNOWN_ALG;

			ptr = fgets(buf, 1024, in);
			if (ptr == NULL)
				break;

			if (strlen(buf) < 3)
				continue;

			if (ptr[0] == '#')
				continue;

			ptr = strtok(buf, ",");
			if (strcmp(ptr, "lehman") == 0) test.alg = LEHMAN;
			if (strcmp(ptr, "squfof") == 0) test.alg = SQUFOF;
			if (strcmp(ptr, "squfof-par") == 0) test.alg = SQUFOF_PAR;
			if (strcmp(ptr, "rho") == 0) test.alg = RHO;
			if (strcmp(ptr, "fermat") == 0) test.alg = FERMAT;
			if (strcmp(ptr, "uecm") == 0) test.alg = UECM;
			if (strcmp(ptr, "uecm-par") == 0) test.alg = UECM_PAR;
			if (strcmp(ptr, "tecm") == 0) test.alg = TECM;
			if (strcmp(ptr, "tecm-par") == 0) test.alg = TECM_PAR;
			if (strcmp(ptr, "siqs") == 0) test.alg = SIQS;
			strncpy(algname, ptr, 80);

			if (test.alg == UNKNOWN_ALG)
			{
				printf("skipping unknown algorithm %s specified in %s\n", ptr, options->testPlan);
				continue;
			}

			ptr = strtok(NULL, ",");
			strncpy(testFile, ptr, 1024);

			ptr = strtok(NULL, ",");
			test.num_to_test = atoi(ptr);

			test.testfile = testFile;
			test.lcg_state = lcg_state;

			//printf("===============================================================\n");
			//printf("running %u test cases from file %s using algorithm %d\n",
			//	test.num_to_test, test.testfile, test.alg);
			//printf("===============================================================\n");
			
			run_test(&test);
			lcg_state = test.lcg_state;

			printf("alg %s: (%d, %d, %1.1f) bits; "
				"%d ff, %d pf, %2.2f sec (%1.2f us per testcase)\n",
				algname, test.min_bits, test.max_bits, test.avg_bits, test.num_ff, test.num_pf, 
				test.elasped_time, 1000000 * test.elasped_time / (double)test.num_to_test);
		}

		fclose(in);

		mpz_clear(gmptmp);
		free(primes);
		return 0;

#if 0
		int do_squfof_p;
		int do_tinyecm;

		if (do_squfof_p)
		{
			printf("=================== Squfof Method (Parallel) ==================\n");
			numtest = num;
			for (nf = 0; nf < num_files; nf++)
			{
				in = fopen(filenames[nf], "r");

				gettimeofday(&gstart, NULL);
				i = 0;
				totBits = 0;
				minBits = 999;
				maxBits = 0;
				while (!feof(in))
				{
					fscanf(in, "%" PRIu64 ",%u,%u", comp + i, f1 + i, f2 + i);
					mpz_set_ui(gmptmp, comp[i]);
					j = mpz_sizeinbase(gmptmp, 2);
					if (j < 32) continue;
					totBits += j;
					if ((uint32_t)j > maxBits)
						maxBits = j;
					if ((uint32_t)j < minBits && j != 0)
						minBits = j;
					i++;
				}
				fclose(in);
				gettimeofday(&gstop, NULL);

				f_time = ytools_difftime(&gstart, &gstop);
				avgBits = (double)totBits / (double)i;

				if (avgBits > 62.5)
					continue;

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;

				par_shanks_loop(comp, fout64, numtest);

				for (i = 0; i < numtest; i++)
				{
					if (((uint32_t)fout64[i] == f1[i]) || ((uint32_t)fout64[i] == f2[i]))
						correct++;
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("squfof-parallel setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

		if (do_tinyecm)
		{
			mpz_t fact1, fact2, gmp_comp;

			mpz_init(fact1);
			mpz_init(fact2);
			mpz_init(gmp_comp);

			printf("====================== tinyecm Method =========================\n");
			numtest = 10000;
			for (nf = 0; nf < num_big_files; nf++)
			{
				uint64_t known1, known2;
				char buf[1024];
				in = fopen(bigfilenames[nf], "r");

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;
				i = 0;
				totBits = 0;
				minBits = 999;
				maxBits = 0;
				for (i = 0; i < numtest; i++)
				{
					fgets(buf, 768, in);

#ifdef _MSC_VER
					gmp_sscanf(buf, "%Zd, %llu, %llu",
						gmp_comp, &known1, &known2);
#else
					gmp_sscanf(buf, "%Zd,%lu,%lu", gmp_comp, &known1, &known2);
#endif

					j = mpz_sizeinbase(gmp_comp, 2);
					if (j < 32) continue;
					totBits += j;
					if ((uint32_t)j > maxBits)
						maxBits = j;
					if ((uint32_t)j < minBits && j != 0)
						minBits = j;

					getfactor_tecm(gmp_comp, fact1, 0, &lcg_state);
					f64 = mpz_get_64(fact1);

					if ((f64 == known1) || (f64 == known2))
					{
						correct++;
					}
				}

				fclose(in);
				avgBits = (double)totBits / (double)numtest;

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("tinyecm setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);
			}

			mpz_clear(gmp_comp);
			mpz_clear(fact1);
			mpz_clear(fact2);
		}
#endif
	}

	mpz_clear(gmptmp);
	free(primes);
	return 0;
}

int run_test(test_t* test)
{
	int status = success;
	uint32_t minBits;
	uint32_t maxBits;
	uint32_t totBits;
	double avgBits;
	int i;
	int j;
	char buf[1024];
	char* fptr;
	FILE* in;
	mpz_t *composites, **factors;
	uint64_t* composites64, **factors64;
	uint64_t* fout64;
	uint32_t *num_factors;
	struct timeval gstart;
	struct timeval gstop;
	uint32_t num_read;
	mpz_t f1, f2;

	in = fopen(test->testfile, "r");
	if (in == NULL)
	{
		printf("failed to open input file %s\n", test->testfile);
		exit(0);
	}

	mpz_init(f1);
	mpz_init(f2);

	switch (test->alg)
	{
		case SIQS:
			// special case inputs/outputs:
			// takes an mpz_t as an input and returns 2 factors as mpz_t's
		{
			
			gettimeofday(&gstart, NULL);

			composites = (mpz_t*)xmalloc(test->num_to_test * sizeof(mpz_t));
			factors = (mpz_t**)xmalloc(test->num_to_test * sizeof(mpz_t *));
			num_factors = (uint32_t *)xmalloc(test->num_to_test * sizeof(uint32_t));

			totBits = 0;
			minBits = 999;
			maxBits = 0;
			num_read = 0;
			for (i = 0; i < test->num_to_test && !feof(in); i++)
			{
				fptr = fgets(buf, 1024, in);

				if (fptr == NULL)
					break;

				mpz_init(composites[i]);

				fptr = strtok(buf, ",");
				gmp_sscanf(fptr, "%Zd", composites[i]);
				fptr = strtok(NULL, ",");
				gmp_sscanf(fptr, "%u", &num_factors[i]);

				factors[i] = (mpz_t*)xmalloc(num_factors[i] * sizeof(mpz_t));
				
				for (j = 0; j < num_factors[i]; j++)
				{
					fptr = strtok(NULL, ",");
					mpz_init(factors[i][j]);
					gmp_sscanf(fptr, "%Zd", factors[i][j]);
				}

				j = mpz_sizeinbase(composites[i], 2);

				totBits += j;
				if ((uint32_t)j > maxBits)
					maxBits = j;
				if ((uint32_t)j < minBits && j != 0)
					minBits = j;

				num_read++;
			}

			fclose(in);
			gettimeofday(&gstop, NULL);
			test->parse_time = ytools_difftime(&gstart, &gstop);
			test->avg_bits = (double)totBits / (double)i;
			test->min_bits = minBits;
			test->max_bits = maxBits;

			if (num_read != test->num_to_test)
			{
				printf("warning: did not read %d requested test cases in file %s\n"
					"proceeding with %u test cases\n", 
					test->num_to_test, test->testfile, num_read);
				status = warnings;
			}

			// now run the test
			gettimeofday(&gstart, NULL);

			test->num_ff = 0;
			test->num_pf = 0;

			tiny_qs_params* params;
			params = init_tinyqs();

			for (i = 0; i < num_read; i++)
			{
				int k;

				k = tinyqs(params, composites[i], f1, f2);

				if ((mpz_cmp(f1, factors[i][0]) == 0) || (mpz_cmp(f1, factors[i][1]) == 0))
				{
					test->num_ff++;
				}
			}

			gettimeofday(&gstop, NULL);
			test->elasped_time = ytools_difftime(&gstart, &gstop);
			
			for (i = 0; i < num_read; i++)
			{
				for (j = 0; j < num_factors[i]; j++)
				{
					mpz_clear(factors[i][j]);
				}
				mpz_clear(composites[i]);
				free(factors[i]);
			}

			free(composites);
			free(factors);
			free(num_factors);
		}
		
		mpz_init(f1);
		mpz_init(f2);
		return status;
		break;

		case TECM_PAR:
			// special case inputs/outputs:
			// takes a list of 104-bit inputs formatted as 52 bits in each
			// of 2 contiguous uint64_t's.  Output factors are the same.
		{
			fout64 = (uint64_t*)xmalloc(test->num_to_test * 2 * sizeof(uint64_t));
			composites64 = (uint64_t*)xmalloc(test->num_to_test * 2 * sizeof(uint64_t));
			factors64 = (uint64_t**)xmalloc(test->num_to_test * sizeof(uint64_t*));
			num_factors = (uint32_t*)xmalloc(test->num_to_test * sizeof(uint32_t));

			totBits = 0;
			minBits = 999;
			maxBits = 0;
			num_read = 0;
			for (i = 0; i < test->num_to_test && !feof(in); i++)
			{
				fptr = fgets(buf, 1024, in);

				if (fptr == NULL)
					break;

				fptr = strtok(buf, ",");
				gmp_sscanf(fptr, "%Zd", f1);
				fptr = strtok(NULL, ",");
				gmp_sscanf(fptr, "%u", &num_factors[i]);

				factors64[i] = (uint64_t*)xmalloc(num_factors[i] * sizeof(uint64_t));

				j = mpz_sizeinbase(f1, 2);
				totBits += j;
				if ((uint32_t)j > maxBits)
					maxBits = j;
				if ((uint32_t)j < minBits && j != 0)
					minBits = j;

				composites64[2 * i + 0] = mpz_get_ui(f1) & 0x000fffffffffffffull;
				mpz_tdiv_q_2exp(f1, f1, 52);
				composites64[2 * i + 1] = mpz_get_ui(f1) & 0x000fffffffffffffull;

				for (j = 0; j < num_factors[i]; j++)
				{
					fptr = strtok(NULL, ",");
					
					gmp_sscanf(fptr, "%Zd", f1);
					factors64[i][j] = mpz_get_ui(f1);
				}

				num_read++;
			}

			fclose(in);
			gettimeofday(&gstop, NULL);
			test->parse_time = ytools_difftime(&gstart, &gstop);
			test->avg_bits = (double)totBits / (double)i;
			test->min_bits = minBits;
			test->max_bits = maxBits;

			if (num_read != test->num_to_test)
			{
				printf("warning: did not read %d requested test cases in file %s\n"
					"proceeding with %u test cases\n",
					test->num_to_test, test->testfile, num_read);
				status = warnings;
			}

			gettimeofday(&gstart, NULL);

			test->num_ff = 0;
			test->num_pf = 0;
			
			getfactor_tecm_x8_list(composites64, fout64, maxBits / 2, num_read, &test->lcg_state);

			for (i = 0; i < num_read; i++)
			{
				uint64_t f64 = (fout64[2 * i + 1] << 52) | fout64[2 * i + 0];

				if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
				{
					test->num_ff++;
				}
			}

			gettimeofday(&gstop, NULL);
			test->elasped_time = ytools_difftime(&gstart, &gstop);

			free(num_factors);
			free(fout64);
			free(composites64);
			
			for (i = 0; i < num_read; i++)
			{
				free(factors64[i]);
			}
			free(factors64);
		}

		mpz_init(f1);
		mpz_init(f2);
		return status;
		break;

		case TECM:
			// special case inputs/outputs:
			// takes an mpz_t as an input and returns 2 factors as mpz_t's
		{

			gettimeofday(&gstart, NULL);

			composites = (mpz_t*)xmalloc(test->num_to_test * sizeof(mpz_t));
			factors = (mpz_t**)xmalloc(test->num_to_test * sizeof(mpz_t*));
			num_factors = (uint32_t*)xmalloc(test->num_to_test * sizeof(uint32_t));

			totBits = 0;
			minBits = 999;
			maxBits = 0;
			num_read = 0;
			for (i = 0; i < test->num_to_test && !feof(in); i++)
			{
				fptr = fgets(buf, 1024, in);

				if (fptr == NULL)
					break;

				mpz_init(composites[i]);

				fptr = strtok(buf, ",");
				gmp_sscanf(fptr, "%Zd", composites[i]);
				fptr = strtok(NULL, ",");
				gmp_sscanf(fptr, "%u", &num_factors[i]);

				factors[i] = (mpz_t*)xmalloc(num_factors[i] * sizeof(mpz_t));

				for (j = 0; j < num_factors[i]; j++)
				{
					fptr = strtok(NULL, ",");
					mpz_init(factors[i][j]);
					gmp_sscanf(fptr, "%Zd", factors[i][j]);
				}

				j = mpz_sizeinbase(composites[i], 2);

				totBits += j;
				if ((uint32_t)j > maxBits)
					maxBits = j;
				if ((uint32_t)j < minBits && j != 0)
					minBits = j;

				num_read++;
			}

			fclose(in);
			gettimeofday(&gstop, NULL);
			test->parse_time = ytools_difftime(&gstart, &gstop);
			test->avg_bits = (double)totBits / (double)i;
			test->min_bits = minBits;
			test->max_bits = maxBits;

			if (num_read != test->num_to_test)
			{
				printf("warning: did not read %d requested test cases in file %s\n"
					"proceeding with %u test cases\n",
					test->num_to_test, test->testfile, num_read);
				status = warnings;
			}

			// now run the test
			gettimeofday(&gstart, NULL);

			test->num_ff = 0;
			test->num_pf = 0;

			for (i = 0; i < num_read; i++)
			{
				getfactor_tecm(composites[i], f1, 0, &test->lcg_state);

				if ((mpz_cmp(f1, factors[i][0]) == 0) || (mpz_cmp(f1, factors[i][1]) == 0))
				{
					test->num_ff++;
				}
			}

			gettimeofday(&gstop, NULL);
			test->elasped_time = ytools_difftime(&gstart, &gstop);

			for (i = 0; i < num_read; i++)
			{
				for (j = 0; j < num_factors[i]; j++)
				{
					mpz_clear(factors[i][j]);
				}
				mpz_clear(composites[i]);
				free(factors[i]);
			}

			free(composites);
			free(factors);
			free(num_factors);
		}

		mpz_init(f1);
		mpz_init(f2);
		return status;
		break;

		case FERMAT:
		case RHO:
		case LEHMAN:
		case UECM:
		case UECM_PAR:
		case SQUFOF:
		case SQUFOF_PAR:
		case SQUFOF_RDS:
		case SQUFOF_ALPERN:

			composites64 = (uint64_t*)xmalloc(test->num_to_test * sizeof(uint64_t));
			factors64 = (uint64_t**)xmalloc(test->num_to_test * sizeof(uint64_t*));
			num_factors = (uint32_t*)xmalloc(test->num_to_test * sizeof(uint32_t));
			fout64 = (uint64_t*)xmalloc(test->num_to_test * sizeof(uint64_t));

			totBits = 0;
			minBits = 999;
			maxBits = 0;
			num_read = 0;
			for (i = 0; i < test->num_to_test && !feof(in); i++)
			{
				fptr = fgets(buf, 1024, in);

				if (fptr == NULL)
					break;

				fptr = strtok(buf, ",");
				gmp_sscanf(fptr, "%Zd", f1);
				fptr = strtok(NULL, ",");
				gmp_sscanf(fptr, "%u", &num_factors[i]);

				factors64[i] = (uint64_t*)xmalloc(num_factors[i] * sizeof(uint64_t));

				j = mpz_sizeinbase(f1, 2);
				totBits += j;
				if ((uint32_t)j > maxBits)
					maxBits = j;
				if ((uint32_t)j < minBits && j != 0)
					minBits = j;

				if ((test->alg == LEHMAN) && (j > 60))
				{
					printf("invalid test case for algorithm LEHMAN: inputs must be <= 60 bits\n");
					exit(0);
				}
				if ((test->alg == SQUFOF) && (j > 62))
				{
					printf("invalid test case for algorithm SQUFOF: inputs must be <= 62 bits\n");
					exit(0);
				}
				if ((test->alg == SQUFOF_PAR) && (j > 62))
				{
					printf("invalid test case for algorithm SQUFOF_PAR: inputs must be <= 62 bits\n");
					exit(0);
				}
				if ((test->alg == SQUFOF_RDS) && (j > 62))
				{
					printf("invalid test case for algorithm SQUFOF_RDS: inputs must be <= 62 bits\n");
					exit(0);
				}
				if ((test->alg == SQUFOF_ALPERN) && (j > 62))
				{
					printf("invalid test case for algorithm SQUFOF_ALPERN: inputs must be <= 62 bits\n");
					exit(0);
				}
				if ((test->alg == UECM_PAR) && (j > 52))
				{
					printf("invalid test case for algorithm UECM_PAR: inputs must be <= 52 bits\n");
					exit(0);
				}

				composites64[i] = mpz_get_ui(f1);

				for (j = 0; j < num_factors[i]; j++)
				{
					fptr = strtok(NULL, ",");

					gmp_sscanf(fptr, "%Zd", f1);
					factors64[i][j] = mpz_get_ui(f1);
				}

				num_read++;
			}

			fclose(in);
			gettimeofday(&gstop, NULL);
			test->parse_time = ytools_difftime(&gstart, &gstop);
			test->avg_bits = (double)totBits / (double)i;
			test->min_bits = minBits;
			test->max_bits = maxBits;

			if (num_read != test->num_to_test)
			{
				printf("warning: did not read %d requested test cases in file %s\n"
					"proceeding with %u test cases\n",
					test->num_to_test, test->testfile, num_read);
				status = warnings;
			}

			break;

		default:

			break;

	}

	uint64_t f64;
	uint32_t fmt_iterations = 1000000;
	
	gettimeofday(&gstart, NULL);

	test->num_ff = 0;
	test->num_pf = 0;

	switch (test->alg)
	{

	case FERMAT:
		for (i = 0; i < num_read; i++)
		{
			f64 = spfermat(fmt_iterations, 1, composites64[i]);
			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
				continue;
			}
			f64 = spfermat(fmt_iterations, 3, composites64[i]);
			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
				continue;
			}
			f64 = spfermat(fmt_iterations, 5, composites64[i]);
			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
				continue;
			}
			f64 = spfermat(fmt_iterations, 7, composites64[i]);
			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
				continue;
			}
			f64 = spfermat(fmt_iterations, 11, composites64[i]);
			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
				continue;
			}
			f64 = spfermat(fmt_iterations, 15, composites64[i]);
			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
				continue;
			}
		}
		break;
	case RHO:
		for (i = 0; i < num_read; i++)
		{
			f64 = spbrent64(composites64[i], 8192);

			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	case LEHMAN:
		for (i = 0; i < num_read; i++)
		{
			f64 = LehmanFactor(composites64[i], 3.5, 0, 0.1);

			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	case UECM:
		for (i = 0; i < num_read; i++)
		{
			f64 = getfactor_uecm(composites64[i], 0, &test->lcg_state);

			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	case UECM_PAR:
		getfactor_uecm_x8_list(composites64, fout64, num_read, &test->lcg_state);

		for (i = 0; i < num_read; i++)
		{
			if ((fout64[i] == factors64[i][0]) || (fout64[i] == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	case SQUFOF:
		for (i = 0; i < num_read; i++)
		{
			mpz_set_64(f1, composites64[i]);
			f64 = sp_shanks_loop(f1);

			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	case SQUFOF_PAR:
		par_shanks_loop(composites64, fout64, num_read);

		for (i = 0; i < num_read; i++)
		{
			if ((fout64[i] == factors64[i][0]) || (fout64[i] == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	case SQUFOF_RDS:
		for (i = 0; i < num_read; i++)
		{
			int fac1, fac2, numfound;
			numfound = squfof_rds(composites64[i], &fac1, &fac2);

			if ((fac1 == factors64[i][0]) || (fac1 == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	case SQUFOF_ALPERN:
		for (i = 0; i < num_read; i++)
		{
			f64 = alpern_SQUFOF(composites64[i]);

			if ((f64 == factors64[i][0]) || (f64 == factors64[i][1]))
			{
				test->num_ff++;
			}
		}
		break;
	}

	gettimeofday(&gstop, NULL);
	test->elasped_time = ytools_difftime(&gstart, &gstop);

	mpz_init(f1);
	mpz_init(f2);
	free(num_factors);
	free(composites64);
	free(fout64);
	for (i = 0; i < num_read; i++)
	{
		free(factors64[i]);
	}
	free(factors64);


	return status;
}




int factor_tiny(mpz_t in, mpz_t* out,
	uint64_t* primes, uint64_t nump, uint64_t* prng)
{
	// factor input 'in', which is assumed to be <= 128 bits in size.
	// avoid the overhead associated with the main factor() routine.
	// utilize an input list of primes for trial division.
	// also accept a 64-bit PRNG seed for LCG-RNG
	mpz_t gmpf;
	mpz_init(gmpf);

	// first a bit of trial division.
	int k = 0;
	int numout = 0;
	while ((mpz_cmp_ui(in, 1) > 0) && (primes[k] < 10000) && (k < nump))
	{
		uint64_t q = primes[k];
		uint64_t r = mpz_tdiv_ui(in, q);

		if (r != 0)
		{
			k++;
		}
		else
		{
			mpz_tdiv_q_ui(in, in, q);
			mpz_init(out[numout]);
			mpz_set_64(out[numout], q);
			numout++;
		}
	}

	// survived TD, proceed to ECM
	// this is the lasieve5 cofactorization strategy, modified
	// to handle an arbitrary number of small factors.
	while (mpz_cmp_ui(in, 1) > 0)
	{
		if (mpz_probab_prime_p(in, 1) > 0)
		{
			// prime residue
			mpz_init(out[numout]);
			mpz_set(out[numout], in);
			numout++;
			break;
		}

		if (mpz_sizeinbase(in, 2) <= 64) {
			uint64_t n64 = mpz_get_ui(in);
			uint64_t f = getfactor_uecm(n64, 1, &prng);
			if (f > 1)
			{
				if (prp_uecm(f) == 0)
				{
					// found a composite factor.  try P-1 and rho on the factor.
					uint64_t f1 = getfactor_upm1(f, 33);
					if (f1 > 1) {
						if (prp_uecm(f1) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f1);
							numout++;
							mpz_tdiv_q_ui(in, in, f1);
							if (prp_uecm(f / f1) == 1)
							{
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f / f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f / f1);
							}
							continue;
						}
					}
					f1 = getfactor_upm1(f, 100);
					if (f1 > 1) {
						if (prp_uecm(f1) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f1);
							numout++;
							mpz_tdiv_q_ui(in, in, f1);
							if (prp_uecm(f / f1) == 1)
							{
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f / f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f / f1);
							}
							continue;
						}
					}
					f1 = getfactor_upm1(f, 333);
					if (f1 > 1) {
						if (prp_uecm(f1) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f1);
							numout++;
							mpz_tdiv_q_ui(in, in, f1);
							if (prp_uecm(f / f1) == 1)
							{
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f / f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f / f1);
							}
							continue;
						}
					}
					int imax = 64;
					int found = 0;
					for (; imax < 8192; imax *= 2)
					{
						f1 = spbrent64(f, imax);
						if (f1 > 1) {
							if (prp_uecm(f1) == 1)
							{
								found = 1;
								mpz_init(out[numout]);
								mpz_set_ui(out[numout], f1);
								numout++;
								mpz_tdiv_q_ui(in, in, f1);
								if (prp_uecm(f / f1) == 1)
								{
									mpz_init(out[numout]);
									mpz_set_ui(out[numout], f / f1);
									numout++;
									mpz_tdiv_q_ui(in, in, f / f1);
								}
								break;
							}
						}
					}
					if (!found)
					{
						printf("failed to split composite factor %lu of input %lu\n", f, n64);
						break;
					}
				}
				else
				{
					mpz_init(out[numout]);
					mpz_set_ui(out[numout], f);
					numout++;
					mpz_tdiv_q_ui(in, in, f);
				}
			}
			else
			{
				// uecm failed. try PM1, rho, then MPQS.
				f = getfactor_upm1(n64, 33);
				if (f > 1) {
					if (prp_uecm(f) == 1)
					{
						mpz_init(out[numout]);
						mpz_set_ui(out[numout], f);
						numout++;
						mpz_tdiv_q_ui(in, in, f);
						continue;
					}
				}
				f = getfactor_upm1(n64, 100);
				if (f > 1) {
					if (prp_uecm(f) == 1)
					{
						mpz_init(out[numout]);
						mpz_set_ui(out[numout], f);
						numout++;
						mpz_tdiv_q_ui(in, in, f);
						continue;
					}
				}
				f = getfactor_upm1(n64, 333);
				if (f > 1) {
					if (prp_uecm(f) == 1)
					{
						mpz_init(out[numout]);
						mpz_set_ui(out[numout], f);
						numout++;
						mpz_tdiv_q_ui(in, in, f);
						continue;
					}
				}
				int imax = 64;
				for (; imax < 8192; imax *= 2)
				{
					f = spbrent64(n64, imax);
					if (f > 1) {
						if (prp_uecm(f) == 1)
						{
							mpz_init(out[numout]);
							mpz_set_ui(out[numout], f);
							numout++;
							mpz_tdiv_q_ui(in, in, f);
							break;
						}
					}
				}
				printf("failed to find factor of %lu\n", n64);
				break;
			}
		}
		else
		{
#if 0
			if (getfactor_tecm(in, gmpf,
				mpz_sizeinbase(in, 2) / 3 - 2, &prng) > 0)
			{
				if (mpz_sizeinbase(gmpf, 2) <= max_primebits[s])
				{
					mpz_tdiv_q(fac[1], large_factors[s], fac[0]);

					// if the remaining residue is obviously too big, we're done.
					if (mpz_sizeinbase(fac[1], 2) > ((max_primebits[s] * 2)))
					{
						nf = 0;
						goto done;
					}

					// check if the residue is prime.  could again use
					// a cheaper method.
					if (mpz_probab_prime_p(fac[1], 1) > 0)
					{
						if (mpz_sizeinbase(fac[1], 2) <= max_primebits[s])
						{
							// we just completed a DLP factorization involving
							// 2 primes whos product was > 64 bits.
							nf = 2;
							goto done;
						}
						nf = 0;
						goto done;
					}

					// ok, so we have extracted one suitable factor, and the 
					// cofactor is not prime and a suitable size.  Do more work to 
					// split the cofactor.
					// todo: target this better based on expected factor size.
					uint64_t q64;
					uint64_t f64;
					if (mpz_sizeinbase(fac[1], 2) <= 64)
					{
						q64 = mpz_get_ui(fac[1]);
						f64 = getfactor_uecm(q64, 0, &pran);
						mpz_set_ui(fac[2], f64);
					}
					else
					{
						// we have a composite residue > 64 bits.  
						// use ecm first with high effort.
						getfactor_tecm(fac[1], fac[2], 32, &pran);
					}
					f64 = mpz_get_ui(fac[2]);

					if (f64 > 1)
					{
						mpz_tdiv_q_ui(fac[1], fac[1], f64);
						nf = 3;

						if (mpz_sizeinbase(fac[1], 2) > max_primebits[s]) {
							nf = 0;
						}
						if (mpz_sizeinbase(fac[2], 2) > max_primebits[s]) {
							nf = 0;
						}
						if (mpz_probab_prime_p(fac[0], 1) == 0)
						{
							nf = 0;
						}
						if (mpz_probab_prime_p(fac[1], 1) == 0)
						{
							nf = 0;
						}
						if (mpz_probab_prime_p(fac[2], 1) == 0)
						{
							nf = 0;
						}
					}
					else
					{
						// uecm/tecm failed, which does sometimes happen
						nf = mpqs_factor(fac[1], max_primebits[s], &fac);
						if (nf == 2)
						{
							// fac is now set to mpqs's statically allocated
							// set of mpz_t's.  copy in the one we found by ecm.
							nf = 3;
							mpz_set(fac[2], uecm_factors[0]);
						}
						else
						{
							nf = 0;
						}
					}
				}
				else
				{
					// check if the factor is prime.  could again use
					// a cheaper method.
					if (mpz_probab_prime_p(fac[0], 1) > 0)
					{
						// if the factor is obviously too big, give up.  This isn't a
						// failure since we haven't expended much effort yet.
						nf = 0;
					}
					else
					{
						// tecm found a composite first factor.
						// if it is obviously too big, we're done.
						if (mpz_sizeinbase(fac[0], 2) > ((max_primebits[s] * 2)))
						{
							nf = 0;
							goto done;
						}

						// isolate the 2nd smaller factor, and check its size.
						mpz_tdiv_q(fac[1], large_factors[s], fac[0]);

						if (mpz_sizeinbase(fac[1], 2) > (max_primebits[s]))
						{
							nf = 0;
							goto done;
						}

						// todo: target this better based on expected factor size.
						uint64_t q64;
						uint64_t f64;
						if (mpz_sizeinbase(fac[0], 2) <= 64)
						{
							q64 = mpz_get_ui(fac[0]);
							f64 = getfactor_uecm(q64, 0, &pran);
							mpz_set_ui(fac[2], f64);
						}
						else
						{
							// we have a composite residue > 64 bits.  
							// use ecm with high effort first
							getfactor_tecm(fac[0], fac[2], 32, &pran);
						}
						f64 = mpz_get_ui(fac[2]);

						if (f64 > 1)
						{
							mpz_tdiv_q_ui(fac[0], fac[0], f64);
							nf = 3;

							if (mpz_sizeinbase(fac[0], 2) > max_primebits[s]) {
								nf = 0;
							}
							if (mpz_sizeinbase(fac[2], 2) > max_primebits[s]) {
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[0], 1) == 0)
							{
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[1], 1) == 0)
							{
								nf = 0;
							}
							if (mpz_probab_prime_p(fac[2], 1) == 0)
							{
								nf = 0;
							}

						}
						else
						{
							// uecm/tecm failed, which does sometimes happen
							nf = mpqs_factor(fac[0], max_primebits[s], &fac);
							if (nf == 2)
							{
								// fac is now set to mpqs's statically allocated
								// set of mpz_t's.  copy in the one we found by ecm.
								nf = 3;
								mpz_set(fac[2], uecm_factors[1]);
							}
							else
							{
								nf = 0;
							}
						}
					}
				}
			}
			else
			{
				// if ecm can't find a factor, give up.  
				// unless this is a DLP with lpbr/a > 32... i.e., if the
				// large factor size is greater than 64 bits but less than
				// lpbr/a * 2.  In that case run mpqs... or tecm with
				// greater effort.


#if 0
				if (mpz_sizeinbase(large_factors[s1], 2) <= (max_primebits[s1] * 2))
				{
					if (getfactor_tecm(large_factors[s1], factor1, 33, &pran) > 0)
					{
						if (mpz_sizeinbase(factor1, 2) <= max_primebits[s1])
						{
							mpz_tdiv_q(factor2, large_factors[s1], factor1);

							// check if the residue is prime.  could again use
							// a cheaper method.
							if (mpz_probab_prime_p(factor2, 1) > 0)
							{
								if (mpz_sizeinbase(factor2, 2) <= max_primebits[s1])
								{
									// we just completed a DLP factorization involving
									// 2 primes whos product was > 64 bits.
									mpz_set(large_primes[s1][0], factor1);
									mpz_set(large_primes[s1][1], factor2);
									nlp[s1] = 2;
								}
								else
									break;
							}
							else
								break;
						}
						else
							break;
					}
					else
						break;
				}
				else
					break;
#else

				if (mpz_sizeinbase(large_factors[s], 2) <= (max_primebits[s] * 2))
				{
					nf = mpqs_factor(large_factors[s], max_primebits[s], &fac);
				}
				else
				{
#if 0
					// try for a lucky p-1 hit on the 3LP before we go?
					// testing on an input with LPB=33 and 3LP enabled
					// saw that p-1 finds lots of factors but the residues
					// are all (99.9%) large primes.  I.e., exactly the
					// kind of inputs we want to not waste time on.
					if (getfactor_tpm1(large_factors[s], fac[0], 333))
					{
						mpz_tdiv_q(fac[1], large_factors[s], fac[0]);
						if (mpz_sizeinbase(fac[1], 2) <= max_primebits[s])
						{
							gmp_printf("P-1 Success! %Zd = %Zd * %Zd\n",
								large_factors[s], fac[0], fac[1]);
						}
						else if (mpz_probab_prime_p(fac[1], 1) == 0)
						{
							gmp_printf("Residue %Zd with %d bits is composite\n",
								fac[1], mpz_sizeinbase(fac[1], 2));
							gmp_printf("3LP = ");

							mpz_set(fac[2], fac[0]);
							nf = 1 + mpqs_factor(fac[2], max_primebits[s], &fac);

							for (i = 0; i < nf; i++)
								gmp_printf("%Zd ", fac[i]);
							printf("\n");
						}
					}
#else
					nf = 0;
#endif
				}
#endif
			}
#endif
		}
	}


	mpz_clear(gmpf);
	return numout;
}









