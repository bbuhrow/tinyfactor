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

// Code originally written by Bob Silverman
// Modified for standalone compilation by Ben Buhrow
// Available: https://www.mersenneforum.org/node/9657?p=282510#post282510


#include <stdlib.h>
#include <stdint.h>
#include <math.h>

int qqueue[100];
int qpoint;

void enqu(int q, int* iter)
{
    qqueue[qpoint] = q;
    if (++qpoint > 100) *iter = -1;
}


int squfof_rds(int64_t n, int* fac1, int* fac2)
{   /* start of squfof: factor n as fac1*fac2  faster in FP?????*/
    int64_t temp, temp1;
    register int iq, ll, l2, p, pnext, q, qlast, r, s, t, i;
    int jter, iter;

    qlast = 1;
    s = (int)sqrt(n);

    p = s;
    temp1 = s * s;
    temp = n - temp1;                 /* temp = n - floor(sqrt(n))^2   */
    if (temp == 0)
    {                                   /* Here n is a square            */
        *fac1 = s;
        *fac2 = s;
        return(1);
    }

    q = (int)temp;              /* q = excess of n over next smaller square */
    ll = 1 + 2 * (int)sqrt((double)(p + p));
    l2 = ll / 2;
    qpoint = 0;

    /*   In the loop below, we need to check if q is a square right before   */
    /*  the end of the loop.  Is there a faster way? The current way is      */
    /*   EXPENSIVE! (many branches and double prec sqrt)                     */

    for (jter = 0; jter < 800000; jter++)      /* I see no way to speed this   */
    {                                     /*  main loop                   */
        iq = (s + p) / q;
        pnext = iq * q - p;
        if (q <= ll)
        {
            if ((q & 1) == 0) enqu(q / 2, &jter);
            else if (q <= l2) enqu(q, &jter);
            if (jter < 0)
            {
                return 0;
            }
        }
        t = qlast + iq * (p - pnext);
        qlast = q;
        q = t;
        p = pnext;                          /* check for square; even iter   */
        if (jter & 1) continue;             /* jter is odd:omit square test  */
        r = (int)sqrt((double)q);                 /* r = floor(sqrt(q))      */
        if (q != r * r) continue;
        if (qpoint == 0) goto gotit;
        for (i = 0; i < qpoint - 1; i += 2)      /* treat queue as list for simplicity*/
        {
            if (r == qqueue[i]) goto contin;
            if (r == qqueue[i + 1]) goto contin;
        }
        if (r == qqueue[qpoint - 1]) continue;
        goto gotit;
    contin:;
    }   /* end of main loop */

gotit:;
    qlast = r;
    p = p + r * ((s - p) / r);
    temp = (int64_t)p * (int64_t)p;
    temp = n - temp;
    temp1 = temp / qlast;
    q = (int)temp1;					/* q = (n - p*p)/qlast (div is exact)*/
    for (iter = 0; iter < 40000; iter++)
    {                              /* begin second main loop            */
        iq = (s + p) / q;                /* unroll it, of course              */
        pnext = iq * q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq * (p - pnext);
        qlast = q;
        q = t;
        p = pnext;
        iq = (s + p) / q;
        pnext = iq * q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq * (p - pnext);
        qlast = q;
        q = t;
        p = pnext;
        iq = (s + p) / q;
        pnext = iq * q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq * (p - pnext);
        qlast = q;
        q = t;
        p = pnext;
        iq = (s + p) / q;
        pnext = iq * q - p;
        if (p == pnext) goto gotfac;
        t = qlast + iq * (p - pnext);
        qlast = q;
        q = t;
        p = pnext;
    }


    return(0);                               /* this shouldn't happen      */

gotfac:; if ((q & 1) == 0) q /= 2;      /* q was factor or 2*factor   */
    *fac1 = q;
    temp = n / q;
    *fac2 = (int)temp;
    return(1);
}