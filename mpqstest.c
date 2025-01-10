/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>	//for gettimeofday using gcc
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include "gmp.h"
#include "ggnfs_mpqs/siever-config.h"
#include "if.h"
#include "gmp-aux.h"
#include <string.h>
#include <stdint.h>

int iter=0;
u64_t stat_asm_eval=0, stat_asm_td=0;
u64_t stat_td_cand=0,stat_td_surv=0;
u64_t stat_ff=0,stat_pf=0,stat_comb=0;
u32_t stat_asm_div=0, stat_final_mulmod=0;
u32_t stat_counter0=0, stat_counter1=0, stat_retry=0;
u32_t stat_size[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double _difftime(struct timeval* start, struct timeval* end)
{
	double secs;
	double usecs;

	if (start->tv_sec == end->tv_sec) {
		secs = 0;
		usecs = end->tv_usec - start->tv_usec;
	}
	else {
		usecs = 1000000 - start->tv_usec;
		secs = end->tv_sec - (start->tv_sec + 1);
		usecs += end->tv_usec;
		if (usecs >= 1000000) {
			usecs -= 1000000;
			secs += 1;
		}
	}

	return secs + usecs / 1000000.;
}

int main(int argc, char **argv)
{
  mpz_t N, *f;
  int n, i;

  setbuf(stdout,NULL);
  initzeit(55); zeita(0);
  mpz_ull_init();
//  init_montgomery_multiplication();
  mpz_init(N);

{
	  FILE* in;
	  uint64_t* comp;
	  uint32_t bits, totBits, minBits, maxBits;
	  double avgBits;
	  double t_time;
	  double f_time;
	  int i, j, k, numtest, num = 100000, correct;
	  struct timeval gstart;
	  struct timeval gstop;
	  int nf;
  int num_big_files;
  char bigfilenames[30][80];

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
  num_big_files = i;

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

		  gmp_sscanf(buf, "%Zd,%lu,%lu", N, &known1, &known2);

		  j = mpz_sizeinbase(N, 2);
		  if (j < 32) continue;
		  totBits += j;
		  if ((uint32_t)j > maxBits)
			  maxBits = j;
		  if ((uint32_t)j < minBits && j != 0)
			  minBits = j;

		  n = mpqs_factor(N, j, &f);
		  if (n >= 1) correct++;
	  }

	  fclose(in);
	  avgBits = (double)totBits / (double)numtest;

	  gettimeofday(&gstop, NULL);
	  t_time = _difftime(&gstart, &gstop);

	  printf("tinyqs setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, avg)\n",
		  avgBits, maxBits, correct, numtest, t_time, 1000000 * t_time / (double)numtest);
  }

}

  zeitb(0);
  if (stat_retry) printf("Warning: %u retries\n",stat_retry);
  printf("Stat: sieves %Lu, td %Lu->%Lu->%Lu\n",
	 stat_asm_eval,stat_td_cand,stat_asm_td,stat_td_surv);
  printf("  ff %llu, pf %llu->%llu\n",stat_ff,stat_pf,stat_comb);
  for (i=0; i<14; i++) if (stat_size[i])
    printf("size %d: %u  ",i,stat_size[i]); printf("\n");

  printf("timing:  total: "); printzeit(0); printf("\n");
  printf("  init       : "); printzeit(9); printf("\n");
  printf("  init fb    : "); printzeit(1); printf("\n");
  printf("  nextpol    : "); printzeit(2); printf("\n");
  printf("  sieve      : "); printzeit(3); printf("\n");
  printf("  eval       : "); printzeit(4); printf("\n");
  printf("  decompose  : "); printzeit(5); printf("\n");
  printf("  final      :"); printzeit(10); printf("\n");
  printf("    matrix   : "); printzeit(8); printf("\n");
printzeit(20); printzeit(21); printzeit(22); printzeit(23); printzeit(24); printzeit(25); printzeit(26); printf("\n");
}

