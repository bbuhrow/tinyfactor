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



#include "arith.h"
#include "ytools.h"
#include "monty.h"

uint64_t spbrent(uint64_t N, uint64_t c, int imax)
{


    /*
    run pollard's rho algorithm on n with Brent's modification,
    returning the first factor found in f, or else 0 for failure.
    use f(x) = x^2 + c
    see, for example, bressoud's book.
    */
    uint64_t x, y, q, g, ys, t1, f = 0, nhat;
    uint32_t i = 0, k, r, m;
    int it;
    
    // start out checking gcd fairly often
    r = 1;

    // under 48 bits, don't defer gcd quite as long
    i = _trail_zcnt64(N);
    if (i > 20)
        m = 32;
    else if (i > 16)
        m = 160;
    else if (i > 3)
        m = 256;
    else
        m = 384;

    it = 0;
    q = 1;
    g = 1;

    x = (((N + 2) & 4) << 1) + N; // here x*a==1 mod 2**4
    x *= 2 - N * x;               // here x*a==1 mod 2**8
    x *= 2 - N * x;               // here x*a==1 mod 2**16
    x *= 2 - N * x;               // here x*a==1 mod 2**32         
    x *= 2 - N * x;               // here x*a==1 mod 2**64
    nhat = (uint64_t)0 - x;

    // Montgomery representation of c
    c = u64div(c, N);
    y = c;

    do
    {
        x = y;
        for (i = 0; i <= r; i++)
        {
            y = mulredc63(y, y + c, N, nhat);
        }

        k = 0;
        do
        {
            ys = y;
            for (i = 1; i <= MIN(m, r - k); i++)
            {
                y = mulredc63(y, y + c, N, nhat);
                t1 = x > y ? y - x + N : y - x;
                q = mulredc63(q, t1, N, nhat);
            }

            g = bingcd64(N, q);
            k += m;
            it++;

            if (it>imax)
            {
                f = 0;
                goto done;
            }

        } while ((k<r) && (g == 1));
        r *= 2;
    } while (g == 1);

    if (g == N)
    {
        //back track
        do
        {
            ys = mulredc63(ys, ys + c, N, nhat);
            t1 = x > ys ? ys - x + N : ys - x;
            g = bingcd64(N, t1);
        } while (g == 1);

        if (g == N)
        {
            f = 0;
        }
        else
        {
            f = g;
        }
    }
    else
    {
        f = g;
    }

done:

    return f;
}

uint64_t spbrent64(uint64_t N, int imax)
{
    /*
	run pollard's rho algorithm on n with Brent's modification,
	returning the first factor found in f, or else 0 for failure.
	use f(x) = x^2 + c
	see, for example, bressoud's book.
	*/
	uint64_t x, y, q, g, ys, t1, f = 0, nhat;
	uint64_t c = 1;
	uint32_t i = 0, k, r, m;
	int it;

	// start out checking gcd fairly often
	r = 1;

	// under 48 bits, don't defer gcd quite as long
	i = _lead_zcnt64(N);
	if (i > 20)
		m = 32;
	else if (i > 16)
		m = 160;
	else if (i > 3)
		m = 256;
	else
		m = 384;

	it = 0;
	q = 1;
	g = 1;

	x = (((N + 2) & 4) << 1) + N; // here x*a==1 mod 2**4
	x *= 2 - N * x;               // here x*a==1 mod 2**8
	x *= 2 - N * x;               // here x*a==1 mod 2**16
	x *= 2 - N * x;               // here x*a==1 mod 2**32         
	x *= 2 - N * x;               // here x*a==1 mod 2**64
	nhat = (uint64_t)0 - x;

	// Montgomery representation of c
	c = u64div(c, N);
	y = c;

	do
	{
		x = y;
		for (i = 0; i <= r; i++)
		{
			y = mulredc(y, y + c, N, nhat);
		}

		k = 0;
		do
		{
			ys = y;
			for (i = 1; i <= MIN(m, r - k); i++)
			{
				y = mulredc(y, y + c, N, nhat);
				t1 = x > y ? y - x + N : y - x;
				q = mulredc(q, t1, N, nhat);
			}

			g = bingcd64(N, q);
			k += m;
			it++;

			if (it > imax)
			{
				f = 0;
				goto done;
			}

		} while ((k < r) && (g == 1));
		r *= 2;
	} while (g == 1);

	if (g == N)
	{
		//back track
		do
		{
			ys = mulredc(ys, ys + c, N, nhat);
			t1 = x > ys ? ys - x + N : ys - x;
			g = bingcd64(N, t1);
		} while (g == 1);

		if (g == N)
		{
			f = 0;
		}
		else
		{
			f = g;
		}
	}
	else
	{
		f = g;
	}

done:

	return f;
}
