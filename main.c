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

	int do_lehman = 0;
	int do_siqs = 0;
	int do_microecm = 0;
	int do_tinyecm = 1;
	int do_rho = 0;
	int do_squfof_s = 0;
	int do_squfof_p = 0;
	int do_squfof_rds = 0;
	int do_squfof_alpern = 0;
	int do_fermat = 0;

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

	free(f1);
	free(f2);
	free(fout64);
	free(comp);
	mpz_clear(gmptmp);
	return;
}
