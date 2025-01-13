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
	lehman,
	siqs,
	microecm,
	microecm_par,
	tinyecm,
	tinyecm_par,
	rho,
	squfof,
	squfof_par,
	squfof_rds,
	squfof_alpern,
	fermat
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

	// results
	uint32_t num_ff;
	uint32_t num_pf;
	double parse_time;
	double elasped_time;
} test_t;

int run_test(test_t* test);

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
	uint64_t lcg_state = 0xdeadbeef0badcafe;
	options_t *options;
	gmp_randstate_t gmp_randstate;

	gmp_randinit_default(gmp_randstate);

	options = initOpt();
	processOpts(argc, argv, options);

	if (options->generate)
	{
		FILE* out;
		int i;
		mpz_t tmp1, tmp2, tmp3;
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
			printf("Test generation creates lists of input integers to factor.\n");
			printf("The number of primes and their sizes comprising each test case\n"
				"can be specified as a comma-delimited list of bit-sizes or bit ranges,\n");
			printf("e.g., ""30,30"" specifies a 60-bit semiprime test case\n");
			printf("and ""26,50-70"" specifies a semiprime comprised of one 26-bit prime\n");
			printf("and one 50 to 70-bit prime (uniform random size in range 50 to 70 bits).\n");
			printf("The sum of the sizes should not exceed 128.\n");
			
			exit(0);
		}

		mpz_init(tmp1);
		mpz_init(tmp2);
		mpz_init(tmp3);

		out = fopen(options->outFile, "w");
		if (out == NULL)
		{
			printf("couldn't open %s for writing\n", options->outFile);
			exit(0);
		}

		bits = 64;

		for (i = 0; i < num; i++)
		{
			mpz_urandomb(tmp3, gmp_randstate, bits / 2);
			mpz_setbit(tmp3, bits / 2 - 1);
			mpz_nextprime(tmp1, tmp3);
			if (bits & 1)
			{
				mpz_urandomb(tmp3, gmp_randstate, bits / 2 + 1);
				mpz_setbit(tmp3, bits / 2);
			}
			else
			{
				mpz_urandomb(tmp3, gmp_randstate, bits / 2);
				mpz_setbit(tmp3, bits / 2 - 1);
			}
			mpz_nextprime(tmp2, tmp3);
			mpz_mul(tmp3, tmp2, tmp1);
			gmp_fprintf(out, "%Zd,%Zd,%Zd\n",
				tmp3, tmp1, tmp2);
		}
		fclose(out);

		printf("generated %d semiprimes in file semiprimes.dat\n", num);
		mpz_clear(tmp1);
		mpz_clear(tmp2);
		mpz_clear(tmp3);
	}
	else if (options->test)
	{
		if (strlen(options->testPlan) == 0)
		{
			printf("specify a test plan to run with option -i or --infile\n");
			printf("documentation about how to construct a test plan and \n");
			printf("example plans should be contained in the repository at\n");
			printf("https://github.com/bbuhrow/tinyfactor\n");
			//exit(0);
		}

		int do_lehman = 0;
		int do_siqs = 0;
		int do_microecm = 0;
		int do_microecm_p = 0;
		int do_tinyecm = 0;
		int do_tinyecm_p = 1;
		int do_rho = 0;
		int do_squfof_s = 0;
		int do_squfof_p = 0;
		int do_squfof_rds = 0;
		int do_squfof_alpern = 0;
		int do_fermat = 0;
		int do_rds_qs = 0;

		printf("============================================================\n");
		printf("64-bit tests												\n");
		printf("============================================================\n\n");

		mpz_init(gmptmp);
		comp = (uint64_t*)malloc(2000000 * sizeof(uint64_t));
		f1 = (uint32_t*)malloc(2000000 * sizeof(uint32_t));
		f2 = (uint32_t*)malloc(2000000 * sizeof(uint32_t));
		fout64 = (uint64_t*)malloc(2000000 * sizeof(uint64_t));

		i = 0;
		//strcpy(filenames[i++], "semiprimes_32bit.dat");
		strcpy(filenames[i++], "semiprimes_34bit.dat");
		strcpy(filenames[i++], "semiprimes_36bit.dat");
		strcpy(filenames[i++], "semiprimes_38bit.dat");
		strcpy(filenames[i++], "semiprimes_40bit.dat");
		strcpy(filenames[i++], "semiprimes_42bit.dat");
		strcpy(filenames[i++], "semiprimes_44bit.dat");
		strcpy(filenames[i++], "semiprimes_46bit.dat");
		strcpy(filenames[i++], "semiprimes_48bit.dat");
		strcpy(filenames[i++], "semiprimes_50bit.dat");
		strcpy(filenames[i++], "semiprimes_52bit.dat");
		strcpy(filenames[i++], "semiprimes_54bit.dat");
		strcpy(filenames[i++], "semiprimes_56bit.dat");
		strcpy(filenames[i++], "semiprimes_58bit.dat");
		strcpy(filenames[i++], "semiprimes_60bit.dat");
		strcpy(filenames[i++], "semiprimes_62bit.dat");
		strcpy(filenames[i++], "semiprimes_64bit.dat");
		num_files = i;

		i = 0;
		strcpy(bigfilenames[i++], "semiprimes_50bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_52bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_54bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_56bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_58bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_60bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_62bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_64bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_70bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_80bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_90bit.dat");
		strcpy(bigfilenames[i++], "semiprimes_100bit.dat");
		num_big_files = i;


		if (do_lehman)
		{
			printf("=================== Lehman Method ============================\n");

			numtest = num;

			for (nf = 0; nf < num_files; nf++)
			{
				in = fopen(filenames[nf], "r");

				i = 0;
				totBits = 0;
				minBits = 999;
				maxBits = 0;
				gettimeofday(&gstart, NULL);

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

				if (avgBits > 60.5)
					continue;

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;
				int j;

				for (i = 0; i < numtest; i++)
				{
					f64 = LehmanFactor(comp[i], 3.5, 0, 0.1);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
					}
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("Lehman setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

		if (do_rho)
		{
			printf("====================== Rho Method ============================\n");
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

				gettimeofday(&gstart, NULL);

				correct = 0;
				for (i = 0; i < numtest; i++)
				{
					f64 = spbrent64(comp[i], 8192);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
					}
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("rho setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

		if (do_squfof_s)
		{
			printf("=================== Squfof Method (Serial)=====================\n");
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
				for (i = 0; i < numtest; i++)
				{
					mpz_set_64(gmptmp, comp[i]);
					f64 = sp_shanks_loop(gmptmp);

					if (((uint32_t)f64 == f1[i]) || ((uint32_t)f64 == f2[i]))
						correct++;
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("squfof-serial setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

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

		if (do_squfof_rds)
		{
			printf("=================== Squfof Method (RDS) ========================\n");
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
				for (i = 0; i < numtest; i++)
				{
					int fac1, fac2;

					f64 = squfof_rds(comp[i], &fac1, &fac2);

					if (((uint32_t)fac1 == f1[i]) || ((uint32_t)fac1 == f2[i]))
						correct++;
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("squfof-rds setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

		if (do_squfof_alpern)
		{
			printf("=================== Squfof Method (Alpern) =====================\n");
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

				if (avgBits > 58.5)
					continue;

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;
				for (i = 0; i < numtest; i++)
				{
					f64 = alpern_SQUFOF(comp[i]);

					if (((uint32_t)f64 == f1[i]) || ((uint32_t)f64 == f2[i]))
						correct++;
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("squfof-alpern setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

		if (do_siqs)
		{
			tiny_qs_params* params;
			mpz_t fact1, fact2, gmp_comp;

			mpz_init(fact1);
			mpz_init(fact2);
			mpz_init(gmp_comp);
			params = init_tinyqs();

			printf("====================== SIQS Method ============================\n");
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

					k = tinyqs(params, gmp_comp, fact1, fact2);
					f64 = mpz_get_64(fact1);

					//printf("known factors: %lu, %lu\n", known1, known2);
					if ((f64 == known1) || (f64 == known2))
					{
						correct++;
					}
				}

				fclose(in);
				avgBits = (double)totBits / (double)numtest;

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("tinyqs setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);
			}

			params = free_tinyqs(params);
			mpz_clear(gmp_comp);
			mpz_clear(fact1);
			mpz_clear(fact2);
		}

		if (do_fermat)
		{
			printf("=================== Fermat Method ============================\n");
			numtest = num;
			for (nf = 0; nf < num_files; nf++)
			{
				uint32_t iterations = 1000000;

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

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;
				for (i = 0; i < numtest; i++)
				{
					f64 = spfermat(iterations, 1, comp[i]);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
						continue;
					}
					f64 = spfermat(iterations, 3, comp[i]);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
						continue;
					}
					f64 = spfermat(iterations, 5, comp[i]);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
						continue;
					}
					f64 = spfermat(iterations, 7, comp[i]);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
						continue;
					}
					f64 = spfermat(iterations, 11, comp[i]);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
						continue;
					}
					f64 = spfermat(iterations, 15, comp[i]);
					if ((f64 == f1[i]) || (f64 == f2[i]))
					{
						correct++;
						continue;
					}
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("fermat setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

		if (do_microecm)
		{
			printf("====================== ECM Method ============================\n");
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

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;
				for (i = 0; i < numtest; i++)
				{
					f64 = getfactor_uecm(comp[i], 0, &lcg_state);

					if (((uint32_t)f64 == f1[i]) || ((uint32_t)f64 == f2[i]))
						correct++;
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("uecm setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

		if (do_microecm_p)
		{
			printf("==================== ECM-Parallel Method =====================\n");
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

				if ((maxBits > 52) || (avgBits > 52.0))
					continue;

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;

				// for (i = 0; i < numtest; i++)
				{
					getfactor_uecm_x8_list(comp, fout64, numtest, &lcg_state);

					for (i = 0; i < numtest; i++)
					{
						if (((uint32_t)fout64[i] == f1[i]) || ((uint32_t)fout64[i] == f2[i]))
							correct++;
					}
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("uecm-parallel setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}

#ifdef INCLUDE_RDS_QS
		if (do_rds_qs)
		{
			printf("====================== RDS QS Method ==========================\n");
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

				gettimeofday(&gstart, NULL);

				correct = 0;
				k = 0;
				for (i = 0; i < numtest; i++)
				{
					unsigned int fac1, fac2;

					f64 = small_qs64((uint64)comp[i], &fac1, &fac2);

					if (((uint32_t)f64 == f1[i]) || ((uint32_t)f64 == f2[i]))
						correct++;
				}

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("RDS qs setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);

				if (t_time > 10.0)
					numtest /= 2;
			}
		}
#endif

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

		if (do_tinyecm_p)
		{
			mpz_t fact1, fact2, gmp_comp;

			mpz_init(fact1);
			mpz_init(fact2);
			mpz_init(gmp_comp);

			printf("=================== tinyecm-parallel Method ===================\n");
			numtest = 10000;
			for (nf = 0; nf < num_big_files; nf++)
			{
				uint64_t known1[20000], known2[20000];
				char buf[1024];
				in = fopen(bigfilenames[nf], "r");

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
					gmp_sscanf(buf, "%Zd,%lu,%lu", gmp_comp, &known1[i], &known2[i]);
#endif

					j = mpz_sizeinbase(gmp_comp, 2);
					if (j < 32) continue;
					totBits += j;
					if ((uint32_t)j > maxBits)
						maxBits = j;
					if ((uint32_t)j < minBits && j != 0)
						minBits = j;

					//gmp_printf("n[%d] = %Zx", i, gmp_comp);

					comp[2 * i + 0] = mpz_get_ui(gmp_comp) & 0x000fffffffffffffull;
					mpz_tdiv_q_2exp(gmp_comp, gmp_comp, 52);
					comp[2 * i + 1] = mpz_get_ui(gmp_comp) & 0x000fffffffffffffull;

					//printf(" = %013lx%013lx\n", comp[2 * i + 1], comp[2 * i + 0]);

				}

				fclose(in);
				avgBits = (double)totBits / (double)numtest;

				gettimeofday(&gstart, NULL);

				{
					getfactor_tecm_x8_list(comp, fout64, maxBits / 2, numtest, &lcg_state);

					for (i = 0; i < numtest; i++)
					{
						f64 = (fout64[2 * i + 1] << 52) + fout64[2 * i + 0];
						if ((f64 == known1[i]) || (f64 == known2[i]))
						{
							correct++;
						}
					}
				}

				

				gettimeofday(&gstop, NULL);
				t_time = ytools_difftime(&gstart, &gstop);

				printf("tinyecm-parallel setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
					avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);
			}

			mpz_clear(gmp_comp);
			mpz_clear(fact1);
			mpz_clear(fact2);
		}

		free(f1);
		free(f2);
		free(fout64);
		free(comp);
		mpz_clear(gmptmp);

	}

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
	uint32_t *num_factors;
	struct timeval gstart;
	struct timeval gstop;
	uint32_t num_read;

	if (in = fopen(test->testfile, "r") == NULL)
	{
		printf("failed to open input file %s\n", test->testfile);
		exit(0);
	}

	switch (test->alg)
	{
		case siqs:
		case tinyecm_par:
		case microecm_par:

		{
			// these are > 64-bit algorithms.  Some take
			// mpz_t's as input so we just use gmp to read
			// and process the file.
			gettimeofday(&gstart, NULL);

			composites = (mpz_t*)xmalloc(test->num_to_test * sizeof(mpz_t));
			factors = (mpz_t**)xmalloc(test->num_to_test * sizeof(mpz_t *));
			for (i = 0; i < test->num_to_test; i++)
			{
				factors[i] = (mpz_t*)xmalloc(MAX_FACTORS * sizeof(mpz_t));
				for (j = 0; j < MAX_FACTORS; j++)
				{
					mpz_init(factors[i][j]);
				}
			}

			totBits = 0;
			minBits = 999;
			maxBits = 0;
			num_read = 0;
			for (i = 0; i < test->num_to_test && !feof(in); i++)
			{
				fptr = gmp_fscanf(in, "%Zd,%u,", composites[i], &num_factors[i]);

				if (fptr == NULL)
					break;
				
				for (j = 0; j < num_factors[i] - 1; j++)
				{
					gmp_fscanf(in, "%Zd,", factors[i][j]);
				}
				gmp_fscanf(in, "%Zd", factors[i][j]);

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
				printf("warning: did not read %d requested test cases in file %s\n",
					"proceeding with %u test cases\n", 
					test->num_to_test, test->testfile, num_read);
				status = warnings;
			}

			// now run the test
			gettimeofday(&gstart, NULL);


			



			
			for (i = 0; i < test->num_to_test; i++)
			{
				for (j = 0; j < MAX_FACTORS; j++)
				{
					mpz_clear(factors[i][j]);
				}
				free(factors[i]);
			}

			free(composites);
			free(factors);
		}
		
		break;

		case fermat:
		case rho:
		case lehman:
		case microecm:
		case squfof:
		case squfof_rds:
		case squfof_alpern:

			// take uint64_t's as input
			totBits = 0;
			minBits = 999;
			maxBits = 0;
			gettimeofday(&gstart, NULL);

			while (!feof(in))
			{
				fscanf(in, "%" PRIu64 ",%u,", composites64 + i, num_factors[i]);
				for (j = 0; j < num_factors - 1; j++)
				{
					gmp_fscanf(in, "%Zd,", factors[i][j]);
				}
				gmp_fscanf(in, "%Zd", factors[i][j]);

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

			if (avgBits > 60.5)
				continue;


			break;

		default:

	}


	


	return status;
}













