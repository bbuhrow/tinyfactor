
CC = gcc
CFLAGS = -fno-common -g -m64 -std=gnu99 -fPIE
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
INC += -I../gmp-install/6.2.1/include
LIBS += -L../gmp-install/6.2.1/lib



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
	monty.c \
	rds_squfof.c \
	alpern_squfof.c \
	cmdOptions.c
	
HDR = \
	microecm.h \
	cofactorize.h \
	ytools.h \
	tinyfactor.h \
	monty.h \
	arith.h \
	common.h \
	cmdOptions.h
	
	
OBJS = $(SRCS:.c=.o)


#---------------------------Make Targets -------------------------

gmp-aux.o: gmp-aux.c if.h ggnfs_mpqs/siever-config.h
	$(CC) $(CFLAGS)  -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<
	
mpz-ull.o: mpz-ull.c if.h ggnfs_mpqs/siever-config.h
	$(CC) $(CFLAGS)  -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<
	
libgmp-aux.a: gmp-aux.o mpz-ull.o
	$(AR) rcs $@ $^
	
mpqs.o: mpqs.c ggnfs_mpqs/mpqs-config.h ggnfs_mpqs/siever-config.h
#	gcc -c -O2 -o $@ $<
	$(CC) $(CFLAGS) -c -o $@ $<

mpqs3.o: mpqs3.c ggnfs_mpqs/mpqs-config.h ggnfs_mpqs/siever-config.h
#	gcc -c -O2 -o $@ $<
	$(CC) $(CFLAGS) -c -o $@ $<
	
mpqsz.o: mpqs.c ggnfs_mpqs/mpqs-config.h
	$(CC) $(CFLAGS)  -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<
#	gcc -g -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<
	
mpqstest.o: mpqstest.c if.h ggnfs_mpqs/siever-config.h
	$(CC) $(CFLAGS)  -DMPQS_STAT -DMPQS_ZEIT -c -o $@ $<

mpqstest: mpqstest.o mpqsz.o if.o mpz-ull.o ggnfs_mpqs/liblasieve.a
	$(CC) $(CFLAGS) -L../gmp-install/6.2.1/lib -o $@ $^ -lgmp -lm

mpqs3z.o: mpqs3.c ggnfs_mpqs/mpqs-config.h
	$(CC) $(CFLAGS)  -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<
#	gcc -g -DMPQS3_STAT -DMPQS3_ZEIT -c -o $@ $<

mpqs3test: mpqs3test.o mpqs3z.o if.o libgmp-aux.a ggnfs_mpqs/liblasieve.a
	$(CC) $(CFLAGS) -L../gmp-install/6.2.1/lib -o $@ $^ -lgmp -lm

	
all: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o tfactor $(LIBS)
	
clean:
	rm -f $(OBJS) 
#---------------------------Build Rules -------------------------

%.o: %.c $(COMMON_HDR)
	$(CC) $(CFLAGS) -c -o $@ $<
	
