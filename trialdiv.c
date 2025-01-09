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



#include "ytools.h"
#include <math.h>

#define setbit(a,b) (((a)[(b) >> 3]) |= (nmasks[(b) & 7])) 
#define getbit(a,b) (((a)[(b) >> 3]) & (nmasks[(b) & 7])) 


int sptestsqr(uint64_t n)
{
	uint64_t t;
	t = n & 31;
	if (t == 0 || t == 1 || t == 4 ||
		t == 9 || t == 16 || t == 17 || t == 25)
	{
		t = (uint64_t)sqrt((int64_t)n);
		if (n == t * t)
			return 1;
	}
	return 0;
}

static const uint32_t smM = 2 * 2 * 2 * 2 * 3 * 3; //72
static const uint32_t smM1 = 7 * 11; //77
static const uint32_t smM2 = 5 * 13; //65

static uint8_t* smsqr, * smsqr1, * smsqr2, * smmod, * smmod1, * smmod2;
static uint16_t* smskip;
static int sm_fermat_initialized = 0;

uint64_t spfermat(uint64_t limit, uint32_t mult, uint64_t n)
{
    // Fermat's factorization method with a sieve-based improvement
    // provided by 'neonsignal'
    uint64_t a, b2, tmp, multN, a2;
    int i;
    uint64_t count;
    uint64_t i64;
    uint32_t m, mmn, s, d;
    uint32_t iM = 0, iM1 = 0, iM2 = 0;
    uint8_t *thismod, * thismod1, * thismod2;
    int sz = (smM / 8 + 1);
    int sz1 = (smM1 / 8 + 1);
    int sz2 = (smM2 / 8 + 1);
    uint8_t masks[8] = { 0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f };
    uint8_t nmasks[8];

    if (sptestsqr(n))
    {
        return sqrt(n);
    }

    a = 0;
    b2 = 0;
    tmp = 0;
    multN = 0;
    a2 = 0;

    // apply the user supplied multiplier
    multN = n * mult;

    // compute ceil(sqrt(multN))
    a = ceil(sqrt((int64_t)multN));

    // form b^2
    b2 = a * a;
    b2 = b2 - multN;

    for (i = 0; i < 8; i++)
        nmasks[i] = ~masks[i];

    // test successive 'a' values using a sieve-based approach.
    // the idea is that not all 'a' values allow a^2 or b^2 to be square.  
    // we pre-compute allowable 'a' values modulo various smooth numbers and 
    // build tables to allow us to quickly iterate over 'a' values that are 
    // more likely to produce squares.
    // init sieve structures
    if (!sm_fermat_initialized)
    {
        smsqr = (uint8_t*)calloc((smM / 8 + 1), sizeof(uint8_t));
        smsqr1 = (uint8_t*)calloc((smM1 / 8 + 1), sizeof(uint8_t));
        smsqr2 = (uint8_t*)calloc((smM2 / 8 + 1), sizeof(uint8_t));
        smmod = (uint8_t*)calloc((smM / 8 + 1) * smM * smM, sizeof(uint8_t));
        smmod1 = (uint8_t*)calloc((smM1 / 8 + 1) * smM1 * smM1, sizeof(uint8_t));
        smmod2 = (uint8_t*)calloc((smM2 / 8 + 1) * smM2 * smM2, sizeof(uint8_t));
        smskip = (uint16_t*)malloc(smM * sizeof(uint16_t));

        

        // marks locations where squares can occur mod M, M1, M2
        for (i64 = 0; i64 < smM; ++i64)
            setbit(smsqr, (i64 * i64) % smM);

        for (i64 = 0; i64 < smM1; ++i64)
            setbit(smsqr1, (i64 * i64) % smM1);

        for (i64 = 0; i64 < smM2; ++i64)
            setbit(smsqr2, (i64 * i64) % smM2);

        // for the modular sequence of b*b = a*a - n values 
        // (where b2_2 = b2_1 * 2a + 1), mark locations where
        // b^2 can be a square
        int j, k;
        for (j = 0; j < smM; j++)
        {
            for (k = 0; k < smM; k++)
            {
                m = j;
                mmn = k;
                for (i = 0; i < smM; ++i)
                {
                    if (getbit(smsqr, mmn)) setbit(smmod + j * smM * sz + k * sz, i);
                    mmn = (mmn + m + m + 1) % smM;
                    m = (m + 1) % smM;
                }
            }
        }

        // for the modular sequence of b*b = a*a - n values 
        // (where b2_2 = b2_1 * 2a + 1), mark locations where the
        // modular sequence can be a square mod M1.  These will
        // generally differ from the sequence mod M.
        for (j = 0; j < smM1; j++)
        {
            for (k = 0; k < smM1; k++)
            {
                m = j;
                mmn = k;
                for (i = 0; i < smM1; ++i)
                {
                    if (getbit(smsqr1, mmn)) setbit(smmod1 + j * smM1 * sz1 + k * sz1, i);
                    mmn = (mmn + m + m + 1) % smM1;
                    m = (m + 1) % smM1;
                }
            }
        }

        // for the modular sequence of b*b = a*a - n values 
        // (where b2_2 = b2_1 * 2a + 1), mark locations where the
        // modular sequence can be a square mod M2.  These will
        // generally differ from the sequence mod M or M1.
        for (j = 0; j < smM2; j++)
        {
            for (k = 0; k < smM2; k++)
            {
                m = j;
                mmn = k;
                for (i = 0; i < smM2; ++i)
                {
                    if (getbit(smsqr2, mmn)) setbit(smmod2 + j * smM2 * sz2 + k * sz2, i);
                    mmn = (mmn + m + m + 1) % smM2;
                    m = (m + 1) % smM2;
                }
            }
        }

        sm_fermat_initialized = 1;
    }

    // test it.  This will be good enough if |u*p-v*q| < 2 * N^(1/4), where
    // mult = u*v
    count = 0;
    if (sptestsqr(b2))
    {
        return sqrt(b2);
    }

    // select locations where b^2 can be a square mod M from precomputed list
    m = a % smM;
    mmn = b2 % smM;
    thismod = smmod + m * smM * sz + mmn * sz;

    // we only consider locations where the modular sequence mod M can
    // be square, so compute the distance to the next square location
    // at each possible value of i mod M.
    s = 0;
    d = 0;
    for (i = 0; !getbit(thismod, i); ++i)
        ++s;
    for (i = smM; i > 0;)
    {
        --i;
        ++s;
        smskip[i] = s;
        if (s > d) d = s;
        if (getbit(thismod, i)) s = 0;
    }

    // select locations where b^2 can be a square mod M1 from precomputed list
    m = a % smM1;
    mmn = b2 % smM1;
    thismod1 = smmod1 + m * smM1 * sz1 + mmn * sz1;

    // select locations where b^2 can be a square mod M2 from precomputed list
    m = a % smM2;
    mmn = b2 % smM2;
    thismod2 = smmod2 + m * smM2 * sz2 + mmn * sz2;

    // loop, checking for perfect squares
    //mpz_mul_2exp(a2, a, 1);
    a2 = a << 1;
    count = 0;

    do
    {
        d = 0;
        i64 = 0;
        do
        {
            // skip to the next possible square residue of b*b mod M
            s = smskip[iM];

            // remember how far we skipped
            d += s;

            // update the other residue indices
            if ((iM1 += s) >= smM1) iM1 -= smM1;
            if ((iM2 += s) >= smM2) iM2 -= smM2;
            if ((iM += s) >= smM) iM -= smM;

            // some multpliers can lead to infinite loops.  bail out if so.
            if (++i64 > smM) return 0;

            // continue if either of the other residues indicates non-square.
        } while (!getbit(thismod1, iM1) || !getbit(thismod2, iM2));

        // form b^2 by incrementing by many factors of 2*a+1
        tmp = a2 + d;
        tmp = tmp * d;
        b2 = b2 + tmp;

        // accumulate so that we can reset d 
        // (and thus keep it single precision)
        a2 = a2 + d * 2;

        count += d;
        if (count > limit)
            break;

    } while (!sptestsqr(b2));

    // 'count' is how far we had to scan 'a' to find a square b
    a = a + count;

    if ((b2 > 0) && sptestsqr(b2))
    {
        tmp = (uint64_t)sqrt((int64_t)b2);
        tmp = a + tmp;
        tmp = spGCD(n, tmp);
        return tmp;
    }

    return 0;
}




