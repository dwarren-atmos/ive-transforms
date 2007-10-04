LD = /bin/ld
F77 = f77
CC = cc
CPP = /usr/lib/cpp

machine := ${shell /bin/uname}
ifndef machine
machine := ${shell /usr/bin/arch}
endif
ifndef machine
machine := ${shell /usr/bin/uname}
endif
ifndef machine
machine := ${shell /usr/ucb/arch}
endif

ifeq (${machine},sun4)
   SUN = -pic
   LSUN = -assert pure-text
   LIBS = 
   OPT = -fast
   AUTOMATIC = 
else 
ifeq (${machine},Linux)
   F77 = ifort
   LD = ifort
   DEC =   
   LDEC = -shared
   LIBS = -L/usr/local/lib -lnetcdf
   OPT = -O
   INCLUDES = -I/usr/local/include
   AUTOMATIC = 
else
ifeq (${machine},alpha)
   CC = c89
   DEC =  -common -DSIXTY_FOUR_bptr 
   LDEC = -shared 
   OPT = -O2
   AUTOMATIC = -automatic
   INCLUDES = -I /usr/local/include
endif
endif
endif
#library for intermediate level cdf-calls
LIBCDFIO = -L./cdfio -lcdfio
#include path for netcdf.inc
INCCDF = '/usr/include'
#
#
RM = -/bin/rm -f
#
DEBUG = 
CFLAGS = ${OPT} ${SUN} ${DEC} -D_NO_PROTO \
	-I/usr/local/include/udposix -I/usr/include/mit 
FFLAGS = ${OPT} ${SUN} ${DEC} ${DEBUG} -u

OBJS_EM    = trafos.o calc_field.o datelib.o drawing.o derivat.o filter.o \
	     heading.o means.o privat.o thermodyn.o utilities.o compare.o

.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${F77} -c ${FFLAGS} $*.f -o $*.o
.F.o:  
	${F77} ${INCLUDES} -c ${FFLAGS} $*.F

all: cdfio TRANS_EM 


TRANS_EM: ${OBJS_EM}
	${LD} ${LDEC} ${LSUN} ${OBJS_EM} -o trans_em.o \
	${LIBCDFIO}


cdfio: cdfio/libcdfio.o
	${MAKE} -C ./cdfio F77=${F77} FFLAGS="${FFLAGS}" 
clean:
	${RM} *.o cdfio/libcdfio.o cdfio/libcdfio.a

${OBJS_EM}    : attributes.icl constants.icl rotpol.icl cdfconsts.icl text.icl




