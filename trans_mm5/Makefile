LD = ld
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
version := ${shell /usr/bin/uname -r}
ifeq ($(patsubst 5.%,5,${version}),5)
   SUN = -KPIC
   LSUN = -G -z text
   LIBS =
   OPT = -fast
   AUTOMATIC =
   INCLUDES =  -I/opt/local/include
else
   SUN = -pic
   LSUN = -assert pure-text
   LIBS =
   OPT = -fast
   AUTOMATIC =
   INCLUDES =  -I/usr/local/include/udposix
endif
else
ifeq (${machine},mips)
   DEC =   -G 0
   LDEC = -r -G 0
   OPT = -O
   AUTOMATIC = -automatic
   INCLUDES = -I/usr/local/include/udposix -I/usr/include/mit
else
ifeq (${machine},alpha)
   CC = c89
   DEC =  -common -DSIXTY_FOUR_bptr
   LDEC = -shared
   OPT = -O2
   AUTOMATIC = -automatic
   INCLUDES = -I /usr/local/include
else
ifeq (${machine}, IRIX)
   AUTOMATIC =
   LDEC = -shared
   INCLUDES = -I/usr/local/include
else
ifeq (${machine}, Linux)
   F77 = ifort
   SUN = -fpic
   LIBS =
   OPT = -O3
   LDEC = -shared
   INCLUDES = -I/usr/local/include
endif
endif
endif
endif
endif
#

#
#
RM = -/bin/rm -f
#
DEBUG = 
CFLAGS = ${OPT} ${SUN} ${DEC} -D_NO_PROTO ${INCLUDES}
FFLAGS = ${OPT} ${SUN} ${DEC} ${DEBUG} -u

.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${F77} -c ${FFLAGS} $*.f -o $*.o
.F.o:  
	${F77} -c  ${FFLAGS} $*.F -o $*.o


all: trans_mm5.o

trans_mm5.o: default.sigma.o default.ght.o default.icl calc_field.o
	${LD} ${LDEC} ${LSUN} default.sigma.o calc_field.o \
	 -o trans_mm5_sigma.o \
	${LIBS}
	${LD} ${LDEC} ${LSUN} default.ght.o calc_field.o \
	 -o trans_mm5_ght.o \
	${LIBS}

clean:
	${RM} default*.o calc_field.o