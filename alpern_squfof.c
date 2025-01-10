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

// Code originally written by Dario Alpern
// Available: https://www.mersenneforum.org/node/9657#post9657


#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/* Implementation of algorithm explained in Gower and Wagstaff paper */
int alpern_SQUFOF(int64_t N)
{
    int queue[1000];
    double dsqrt;
    int B, Q, Q1, P, P1, L, S;
    int i, r, s, t, q, u;
    int queueHead, queueTail, queueIndex;
    /* Step 1: Initialize */
    if ((N & 3) == 1)
    {
        N <<= 1;
    }
    dsqrt = sqrt(N);
    S = (int)dsqrt;
    if ((long)(S + 1) * (long)(S + 1) <= N)
    {
        S++;
    }
    if ((long)S * (long)S > N)
    {
        S--;
    }
    if ((long)S * (long)S == N)
    {
        return S;
    }
    Q1 = 1;
    P = S;
    Q = (int)N - P * P;
    L = (int)(2 * sqrt(2 * dsqrt));
    B = L << 1;
    queueHead = 0;
    queueTail = 0;
    /* Step 2: Cycle forward to find a proper square form */
    for (i = 0; i <= B; i++)
    {
        q = (S + P) / Q;
        P1 = q * Q - P;
        if (Q <= L)
        {
            if ((Q & 1) == 0)
            {
                queue[queueHead++] = Q >> 1;
                queue[queueHead++] = P % (Q >> 1);
                if (queueHead == 100)
                {
                    queueHead = 0;
                }
            }
            else if (Q + Q <= L)
            {
                queue[queueHead++] = Q;
                queue[queueHead++] = P % Q;
                if (queueHead == 100)
                {
                    queueHead = 0;
                }
            }
        }
        t = Q1 + q * (P - P1);
        Q1 = Q;
        Q = t;
        P = P1;
        if ((i & 1) == 0 && ((Q & 7) < 2 || (Q & 7) == 4))
        {
            r = (int)sqrt(Q);
            if (r * r == Q)
            {
                queueIndex = queueTail;
                for (;;)
                {
                    if (queueIndex == queueHead)
                    {
                        /* Step 3: Compute inverse square root of the square form */
                        Q1 = r;
                        u = (S - P) % r;
                        u += (u >> 31) & r;
                        P = S - u;
                        Q = (int)((N - (long)P * (long)P) / Q1);
                        /* Step 4: Cycle in the reverse direction to find a factor of N */
                        for (;;)
                        {
                            q = (S + P) / Q;
                            P1 = q * Q - P;
                            if (P == P1)
                            {
                                break;
                            }
                            t = Q1 + q * (P - P1);
                            Q1 = Q;
                            Q = t;
                            P = P1;
                        }
                        /* Step 5: Get the factor of N */
                        if ((Q & 1) == 0)
                        {
                            return Q >> 1;
                        }
                        return Q;
                    }
                    s = queue[queueIndex++];
                    t = queue[queueIndex++];
                    if (queueIndex == 100)
                    {
                        queueIndex = 0;
                    }
                    if ((P - t) % s == 0)
                    {
                        break;
                    }
                }
                if (r > 1)
                {
                    queueTail = queueIndex;
                }
                if (r == 1)
                {
                    queueIndex = queueTail;
                    for (;;)
                    {
                        if (queueIndex == queueHead)
                        {
                            break;
                        }
                        if (queue[queueIndex] == 1)
                        {
                            return 0;
                        }
                        queueIndex += 2;
                        if (queueIndex == 100)
                        {
                            queueIndex = 0;
                        }
                    }
                }
            }
        }
    }
    return 0;
}