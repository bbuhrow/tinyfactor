
CC = gcc
CFLAGS = -fno-common -g -m64 -std=gnu99 
WARN_FLAGS = -Wall -Wconversion
OPT_FLAGS = -O2 -DNDEBUG 
OBJ_EXT = .o


# ===================== path options =============================

# standard search directories for headers/libraries within yafu.
# These should normally not be modified.
INC = -I. 
LIBS = -L. 

# we require additional search directories for gmp and gmp-ecm, 
# for libraries and headers.  Change
# these if your installation locations differ.
INC += -I../gmp_install/gmp_6.2.0/include
LIBS += -L../gmp_install/gmp_6.2.0/lib



# ===================== compiler options =========================

ifeq ($(COMPILER),icc)
	CC = icc
	INC += -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4
	CFLAGS += -qopt-report=5
endif



# ===================== architecture features =========================
# several functions in yafu can take advantage of advanced processor
# features (instruction sets).  Specify on the command line, e.g., 
# USE_AVX2=1

ifeq ($(USE_SSE41),1)
	CFLAGS += -DUSE_SSE41 -msse4.1
endif

ifeq ($(USE_AVX2),1)
	USE_SSE41=1
	CFLAGS += -DUSE_AVX2 -DUSE_SSE41 

ifeq ($(COMPILER),icc)
	CFLAGS += -march=core-avx2  
else
	CFLAGS += -mavx2 
endif

endif

ifeq ($(ICELAKE),1)
	CFLAGS += -DUSE_BMI2 -DUSE_AVX2 -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -DIFMA -march=icelake-client
	SKYLAKEX = 1
else

ifeq ($(SKYLAKEX),1)
	CFLAGS += -DUSE_BMI2 -DUSE_AVX2 -DUSE_AVX512F -DUSE_AVX512BW -DSKYLAKEX -march=skylake-avx512 
endif
	
endif

ifeq ($(USE_BMI2),1)
# -mbmi enables _blsr_u64 and -mbmi2 enables _pdep_u64 when using gcc
  CFLAGS += -mbmi2 -mbmi -DUSE_BMI2
endif


# make sure we get the correct libgmp linked by using an absolute path
LIBS += -lgmp

ifeq ($(SKYLAKEX),1)
    # define KNL now for skylakex, after handling an actual command line KNL
    KNL=1
endif

# attempt to get static builds to work... unsuccessful so far
ifeq ($(STATIC),1)
# https://software.intel.com/en-us/articles/error-ld-cannot-find-lm
	CFLAGS += -static-intel -static
	LIBS += -L/usr/lib/x86_64-redhat-linux6E/lib64/ -lpthread -lm
else
	LIBS += -lpthread -lm
endif

ifeq ($(MINGW),1)
# not needed with mingw
#	-ldl
else
	LIBS += -ldl
endif

ifeq ($(COMPILER),icc)
	LIBS +=  -lsvml
endif

CFLAGS += $(OPT_FLAGS) $(WARN_FLAGS) $(INC)


#---------------------------YAFU file lists -------------------------
SRCS = \
	arith.c \
	cofactorize_siqs.c \
	main.c \
	microecm.c \
	micropm1.c \
	rho.c \
	squfof.c \
	tinyecm.c \
	trialdiv.c \
	util.c \
	monty.c
	
HDR = \
	microecm.h \
	cofactorize.h \
	ytools.h \
	tinyfactor.h \
	monty.h \
	arith.h \
	common.h
	
	
OBJS = $(SRCS:.c=.o)


#---------------------------Make Targets -------------------------


all: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o tfactor $(LIBS)
	
clean:
	rm -f $(OBJS) 
#---------------------------Build Rules -------------------------

%.o: %.c $(COMMON_HDR)
	$(CC) $(CFLAGS) -c -o $@ $<