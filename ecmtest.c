/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "ggnfs_mpqs/siever-config.h"
#include "if.h"
#include <string.h>
#include <stdint.h>

#include "ggnfs_mpqs/montgomery_mul.h"

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

int main(int argc, char *argv[])
{
  mpz_t N;
  u32_t n, i, necms=0, iter=0, B1, B2;
  FILE *fi;
  mpz_t *f;
  char *input_line=NULL;
  size_t input_line_alloc=0;
  u32_t imax;
  int nc;
  
  setbuf(stdout,NULL);
  if (argc<4) complain("Usage: ecmtest B1 filename ncurves [ B2 ]\n");
  B1=(u32_t)atoi(argv[1]);
  initzeit(13); zeita(0);
  mpz_init(N); n=0;

  imax=strtoul(argv[3],NULL,10);
  if (argc>4) B2=strtoul(argv[4],NULL,10); else B2=80*B1;
  
  
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
  //for (nf = 0; nf < num_big_files; nf++)
  {
	  uint64_t known1, known2;
	  char buf[1024];
	  in = fopen(argv[2], "r");

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

		  nc=ecm_factor(N,B1,B2,&f,imax);
		  if (nc > 0) { necms+=(u32_t)nc; correct++; /* gmp_printf("%Zd  ",*f); */ } else necms+=imax;
	  }

	  fclose(in);
	  avgBits = (double)totBits / (double)numtest;

	  gettimeofday(&gstop, NULL);
	  t_time = _difftime(&gstart, &gstop);

	  printf("ggnfs-ecm setup: %1.1f avg bits, %d max; results: %d of %d %2.2f sec (%1.2f us, %1.2f curves avg)\n",
		  avgBits, maxBits, correct, numtest, t_time, 
		  1000000 * t_time / (double)numtest, (double)necms / (double)numtest);
  }

}

  
  exit(0);
}

