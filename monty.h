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



#ifndef MONTY_H
#define MONTY_H

#include "arith.h"
#include "common.h"

/********************* arbitrary-precision Montgomery arith **********************/
typedef struct
{
    mpz_t rhat;
    mpz_t nhat;
    mpz_t n;
    mpz_t r;
    mpz_t tmp;

	mpz_t x;
	mpz_t y;
	mpz_t c; 
	mpz_t q; 
	mpz_t g; 
	mpz_t ys; 
	mpz_t t1;

} monty_t;

monty_t * monty_alloc();
void monty_init(mpz_t n, monty_t *mdata);
void to_monty(monty_t *mdata, mpz_t x);
void monty_add(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_mul(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_sub(monty_t *mdata, mpz_t x, mpz_t y, mpz_t res);
void monty_free(monty_t *mdata);
void monty_redc(monty_t *mdata, mpz_t x);

/********************* 128-bit Montgomery arith **********************/
typedef struct
{
	uint64_t r[2];
	uint64_t n[2];
	uint64_t np[2];
	uint64_t nhat[2];
	uint64_t rhat[2];
	uint64_t rmask[2];
	uint64_t one[2];
	uint64_t mtmp1[2];
	uint64_t mtmp2[2];
	uint64_t mtmp3[2];
	uint64_t mtmp4[2];
	uint64_t rho;
} monty128_t;


void to_monty128(monty128_t *mdata, uint64_t * x);
void monty128_init(monty128_t * mdata, uint64_t * n);
void mulmod128(uint64_t * u, uint64_t * v, uint64_t * w, monty128_t *mdata);
void sqrmod128(uint64_t * u, uint64_t * w, monty128_t *mdata);
void addmod128(uint64_t * u, uint64_t * v, uint64_t * w, uint64_t * n);
void submod128(uint64_t * u, uint64_t * v, uint64_t * w, uint64_t * n);


/********************* 64-bit Montgomery arith **********************/

#if (defined(GCC_ASM64X) || defined(__MINGW64__))


__inline uint64_t _umul128(uint64_t x, uint64_t y, uint64_t* hi);

#if defined(USE_AVX512F) || defined(USE_BMI2)
__inline uint64_t mulx64(uint64_t x, uint64_t y, uint64_t* hi) {
    __asm__(
        "mulx %3, %0, %1	\n\t"
        : "=&d"(x), "=&a"(y)
        : "0"(x), "1"(y)
    );

    *hi = y;
    return x;
}
#endif
__inline uint64_t mul64(uint64_t x, uint64_t y, uint64_t* hi) {
    __asm__(
        "mulq %3	\n\t"
        : "=&a"(x), "=&d"(y)
        : "0"(x), "1"(y)
        : "cc"
    );

    *hi = y;
    return x;
}



__inline uint64_t submod(uint64_t a, uint64_t b, uint64_t n)
{
    __asm__(
        "xorq %%r8, %%r8 \n\t"
        "subq %1, %0 \n\t"
        "cmovc %2, %%r8 \n\t"
        "addq %%r8, %0 \n\t"
        : "+r"(a)
        : "r"(b), "r"(n)
        : "r8", "cc");

    return a;
}

__inline uint64_t addmod(uint64_t x, uint64_t y, uint64_t n)
{
    uint64_t t = x - n;
    x += y;
    __asm__("add %2, %1\n\t"
        "cmovc %1, %0\n\t"
        :"+r" (x), "+&r" (t)
        : "r" (y)
        : "cc"
        );
    return x;
}

__inline uint64_t u64div(uint64_t c, uint64_t n)
{
#if 1
    __asm__("divq %4"
        : "=a"(c), "=d"(n)
        : "1"(c), "0"(0ULL), "r"(n));
#else
// this should work if the above won't compile (e.g. on clang)
    uint64_t tmp = 0;
    __asm__("divq %4"
        : "=a"(tmp), "=d"(n)
        : "1"(c), "0"(tmp), "r"(n));
#endif
    return n;
}


#if defined(USE_AVX512F) || defined(USE_BMI2)

__inline uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulx %2, %%r10, %%r11	\n\t"
            "movq %%r10, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "xorq %%r12, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=&d"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulx %2, %%r10, %%r11	\n\t"
            "movq %3, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "mulq %%r10 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=d"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulx %1, %%r10, %%r11	\n\t"
            "movq %%r10, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "xorq %%r12, %%r12 \n\t"
            "mulq %2 \n\t"
            "mulq %3 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %3, %%r12 \n\t"
            "subq %3, %%rdx \n\t"
            "cmovc %%r12, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=&d"(x)
            : "0"(x), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulx %1, %%r10, %%r11	\n\t"
            "movq %2, %%rax		\n\t"
            "xorq %%r8, %%r8 \n\t"
            "mulq %%r10 \n\t"
            "mulq %3 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "subq %3, %%rdx \n\t"
            "cmovc %3, %%r8 \n\t"
            "addq %%r8, %%rdx \n\t"
            : "=d"(x)
            : "0"(x), "r"(nhat), "r"(n)
            : "rax", "r8", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulx %2, %%r10, %%r11	\n\t"
        "movq %3, %%rax		\n\t"
        "xorq %%r8, %%r8 \n\t"
        "mulq %%r10 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%r8 \n\t"
        "addq %%r8, %%rdx \n\t"
        : "=d"(x)
        : "0"(x), "r"(y), "r"(nhat), "r"(n)
        : "rax", "r8", "r10", "r11", "cc");

    return x;
}

__inline uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulx %1, %%r10, %%r11	\n\t"
        "movq %2, %%rax		\n\t"
        "xorq %%r8, %%r8 \n\t"
        "mulq %%r10 \n\t"
        "mulq %3 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "subq %3, %%rdx \n\t"
        "cmovc %3, %%r8 \n\t"
        "addq %%r8, %%rdx \n\t"
        : "=d"(x)
        : "0"(x), "r"(nhat), "r"(n)
        : "rax", "r8", "r10", "r11", "cc");

    return x;
}
#else


__inline uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "movq $0, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "xorq %%rax, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "movq $0, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulq %2	\n\t"
        "movq %%rax, %%r10		\n\t"
        "movq %%rdx, %%r11		\n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(x)
        : "0"(x), "r"(y), "r"(nhat), "r"(n)
        : "rdx", "r10", "r11", "cc");

    return x;
}

__inline uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "movq $0, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "xorq %%rax, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(x), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "movq $0, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(x), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "cc");

    }
    return x;
}

__inline uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
{
    __asm__(
        "mulq %2	\n\t"
        "movq %%rax, %%r10		\n\t"
        "movq %%rdx, %%r11		\n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(x)
        : "0"(x), "r"(x), "r"(nhat), "r"(n)
        : "rdx", "r10", "r11", "cc");

    return x;
}

#endif


#elif _MSC_VER

// TODO: need something portable to replace 64-bit assembler versions
// of modular multiplication.  This is getting closer, but for now these things
// are spread out locally where they are needed instead of being gathered here,
// apologies for the uglyness.

#include <immintrin.h>
#include <intrin.h>

__inline uint64_t u64div(uint64_t c, uint64_t n)
{
    uint64_t r;
    //mpz_t a;
    //mpz_init(a);
    //mpz_set_ui(a, c);
    //mpz_mul_2exp(a, a, 64);
    //r = mpz_tdiv_ui(a, n);
    //mpz_clear(a);
    // first available in Visual Studio 2019
    _udiv128(c, 0, n, &r);

    return r;
}

__inline uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

__inline uint64_t mulredc63(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    return x;
}

__inline uint64_t sqrredc(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}

__inline uint64_t sqrredc63(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    return x;
}

__inline uint64_t submod(uint64_t a, uint64_t b, uint64_t n)
{
    uint64_t r0;
    if (_subborrow_u64(0, a, b, &r0))
        r0 += n;
    return r0;
}

__inline uint32_t submod32(uint32_t a, uint32_t b, uint32_t n)
{
    uint32_t r0;
    if (_subborrow_u32(0, a, b, &r0))
        r0 += n;
    return r0;
}

__inline uint64_t addmod(uint64_t x, uint64_t y, uint64_t n)
{
#if 0
    uint64_t r;
    uint64_t tmp = x - n;
    uint8_t c = _addcarry_u64(0, tmp, y, &r);
    return (c) ? r : x + y;
#else
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x>=n-y) ? x-(n-y) : x+y;
#endif
}

__inline uint32_t addmod32(uint32_t x, uint32_t y, uint32_t n)
{
    // FYI: The clause above often compiles with a branch in MSVC.
    // The statement below often compiles without a branch (uses cmov) in MSVC.
    return (x >= n - y) ? x - (n - y) : x + y;
}



// good to 60 bit inputs
__inline uint64_t sqrredc60(uint64_t x, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    uint8_t c;
    tl = _umul128(x, x, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// good to 60 bit inputs
__inline uint64_t mulredc60(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    uint8_t c;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    c = _addcarry_u64(0, al, tl, &al);
    _addcarry_u64(c, th, ah, &x);
    return x;
}


// this works if inputs are 62 bits or less
#define addmod60(x, y, n) ((x) + (y))

#endif

#endif
