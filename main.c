/*
Copyright (c) 2014, Ben Buhrow
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
	double t_time;
	int i, j, k, numtest, num = 100000, correct;
	mpz_t gmptmp;
	struct timeval gstart;
	struct timeval gstop;
	int nf;
	int num_files;
	char filenames[30][80];
	uint64_t lcg_state = 0xdeadbeef0badcafe;

	int do_lehman = 0;
	int do_siqs = 0;
	int do_microecm = 1;
	int do_tinyecm = 1;
	int do_rho = 1;
	int do_squfof = 1;
	int do_fermat = 1;

	printf("============================================================\n");
	printf("64-bit tests												\n");
	printf("============================================================\n\n");

	mpz_init(gmptmp);
	comp = (uint64_t*)malloc(2000000 * sizeof(uint64_t));
	f1 = (uint32_t*)malloc(2000000 * sizeof(uint32_t));
	f2 = (uint32_t*)malloc(2000000 * sizeof(uint32_t));

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
			t_time = ytools_difftime(&gstart, &gstop);

			printf("data read in %2.4f sec\n", t_time);
			printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
			printf("minimum bits of input numbers = %d\n", minBits);
			printf("maximum bits of input numbers = %d\n", maxBits);

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

			printf("Lehman got %d of %d correct in %2.2f sec\n", correct, numtest, t_time);
			printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)numtest);
			printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)numtest);

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
			t_time = ytools_difftime(&gstart, &gstop);

			printf("data read in %2.4f sec\n", t_time);
			printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
			printf("minimum bits of input numbers = %d\n", minBits);
			printf("maximum bits of input numbers = %d\n", maxBits);

			gettimeofday(&gstart, NULL);

			correct = 0;
			for (i = 0; i < numtest; i++)
			{
				f64 = spbrent64(comp[i], 8192);
				if ((f64 == f1[i]) || (f64 == f2[i]))
				{
					correct++;
					break;
				}
			}

			gettimeofday(&gstop, NULL);
			t_time = ytools_difftime(&gstart, &gstop);

			printf("rho got %d of %d correct in %2.2f sec \n", correct, numtest, t_time);
			printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)numtest);
			printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)numtest);

			if (t_time > 10.0)
				numtest /= 2;
		}
	}

	if (do_squfof)
	{
		printf("=================== Squfof Method ============================\n");
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
			t_time = ytools_difftime(&gstart, &gstop);

			if (((double)totBits / (double)i) > 62.5)
				continue;

			printf("data read in %2.4f sec\n", t_time);
			printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
			printf("minimum bits of input numbers = %d\n", minBits);
			printf("maximum bits of input numbers = %d\n", maxBits);

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

			printf("squfof got %d of %d correct in %2.2f sec\n", correct, numtest, t_time);
			printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)numtest);
			printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)numtest);

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

		for (nf = 0; nf < num_files; nf++)
		{
			uint64_t known1, known2;
			char buf[1024];
			in = fopen(filenames[nf], "r");

			gettimeofday(&gstart, NULL);

			correct = 0;
			k = 0;
			for (i = 0; i < num; i++)
			{
				fgets(buf, 768, in);

#ifdef _MSC_VER
				gmp_sscanf(buf, "%Zd, %llu, %llu",
					gmp_comp, &known1, &known2);
#else
				gmp_sscanf(buf, "%Zd,%u,%u", gmp_comp, &known1, &known2);
#endif

				k = tinyqs(params, gmp_comp, fact1, fact2);

				if (k > 0)
				{
					correct++;
				}
			}

			fclose(in);

			gettimeofday(&gstop, NULL);
			t_time = ytools_difftime(&gstart, &gstop);

			printf("tinyqs got %d of %d correct in %2.2f sec\n", correct, num, t_time);
			printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)num);
			printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);
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
			t_time = ytools_difftime(&gstart, &gstop);

			printf("data read in %2.4f sec\n", t_time);
			printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
			printf("minimum bits of input numbers = %d\n", minBits);
			printf("maximum bits of input numbers = %d\n", maxBits);

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

			printf("Fermat got %d of %d correct in %2.2f sec\n", correct, numtest, t_time);
			printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)numtest);
			printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)numtest);

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
			t_time = ytools_difftime(&gstart, &gstop);

			printf("data read in %2.4f sec\n", t_time);
			printf("average bits of input numbers = %.2f\n", (double)totBits / (double)i);
			printf("minimum bits of input numbers = %d\n", minBits);
			printf("maximum bits of input numbers = %d\n", maxBits);

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

			printf("uecm got %d of %d correct in %2.2f sec\n", correct, numtest, t_time);
			printf("percent correct = %.2f\n", 100.0 * (double)correct / (double)numtest);
			printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)numtest);

			if (t_time > 10.0)
				numtest /= 2;
		}
	}

	free(f1);
	free(f2);
	free(comp);
	mpz_clear(gmptmp);
	return;
}
