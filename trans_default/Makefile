LD = ld
FC = f77
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
   SUN = -K PIC
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
   INCLUDES = -I/usr/local/include
endif
else 
ifeq (${machine},Linux)
   F77 = ifort
   LD = ifort
   CC = icc
   DEC =  
   LDEC = -shared
   OPT = -O
   AUTOMATIC = -automatic
   INCLUDES = -I/usr/local/include
else
ifeq (${machine},alpha)
   CC = c89
   FC = f90
   DEC =  -common -DSIXTY_FOUR_bptr 
   LDEC = -shared 
   OPT = -O2
   AUTOMATIC = -automatic
   INCLUDES = -I /usr/local/include
endif
endif
endif
#
#
RM = /bin/rm -f

#
DEBUG = 
CFLAGS = ${OPT} ${SUN} ${DEC} -D_NO_PROTO ${INCLUDES}
FFLAGS = ${OPT} ${SUN} ${DEC} ${DEBUG} -u

.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${FC} -c ${FFLAGS} $*.f -o $*.o
.F.o:
	${FC} -c ${FFLAGS} $*.F -o $*.o

all: trans_default.o

trans_default.o: default.o default.icl
	${LD} ${LDEC} ${LSUN} default.o -o trans_default.o \
	${LIBS}

clean:
	${RM} *.o *.mod
