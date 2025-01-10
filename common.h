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


#ifndef _COMMON_H_
#define _COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#define _(x) #x
#define STRING(x) _(x)

	/* this byzantine complexity sets up the correct assembly
   language syntax based on the compiler, OS and word size 
   
   Where an inline assembler segment is provided in both 
   GCC and MSC format (i.e. alternative sections), the 
   Intel compiler is configured using guards with an A
   suffix to prefer the native version (GCC on Linux/Unix, 
   MSC on Windows). Where an inline assembler segment is 
   only provided in GCC or MSC format but not both (i.e. 
   exclusive sections) the guards have an X suffix.

   The Intel compiler on Windows appears to have some
   bugs in its processing of GCC inline assembler code.
   These are escaped with the _ICL_WIN_ define
   */

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)

	#define ASM_G __asm__
	#define ASM_M __asm

	/* for inline assembler on Unix/Linux */
	#if defined(__unix__)
		#if defined(__x86_64__)
			#define GCC_ASM64A
			#define GCC_ASM64X
			#define MSC_ASM64X
		#elif defined(__i386__)
			#define GCC_ASM32A
			#define GCC_ASM32X
			#define MSC_ASM32X
		#endif
	#endif

	/* for inline assembler on Windows */
	#if defined(_WIN32)
		#define _ICL_WIN_
		#if defined(_M_X64)
			#define MSC_ASM64A
			#define MSC_ASM64X
			#define GCC_ASM64X
		#elif defined(_M_IX86)
			#define MSC_ASM32A
			#define MSC_ASM32X
			#define GCC_ASM32X
		#endif
	#endif

#elif defined(__GNUC__)

	#define ASM_G __asm__

	#if defined(__x86_64__) 
		#define GCC_ASM64A
		#define GCC_ASM64X
	#elif defined(__i386__)
		#define GCC_ASM32A
		#define GCC_ASM32X
	#endif

#elif defined(_MSC_VER)

	#define ASM_M __asm

	#if defined(_M_IX86) && !defined(_WIN64)
		#define MSC_ASM32A
		#define MSC_ASM32X
	#elif defined(_WIN64)	
		

	#endif
#endif



/* loop alignment directives need to know whether
   we're using MSVC */

#ifndef _MSC_VER
	#define ALIGN_LOOP   ".p2align 4,,7 \n\t" 
#else
	#define ALIGN_LOOP /* nothing */
#endif


#ifdef __cplusplus
}
#endif

#endif /* _COMMON_H_ */

