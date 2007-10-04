# $Id: Makefile,v 1.4 1996/10/25 19:21:38 harry Exp $
#
# $Log: Makefile,v $
# Revision 1.4  1996/10/25 19:21:38  harry
# Add support for SGI.
#

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
   INCLUDES =  -I/usr/local/include/udposix 
endif
else 
ifeq (${machine},Linux)
   CC = icc
   F77 = ifort
   LD = ifort
   DEC =  
   LDEC = -shared
   OPT = -O
   AUTOMATIC =
   INCLUDES = -I/usr/local/include/udposix -I/usr/include/mit 
else
ifeq (${machine},alpha)
   CC = cc -migrate
   DEC =  -common -DSIXTY_FOUR_bptr 
   LDEC = -shared 
   OPT = -O
   AUTOMATIC = -automatic
   INCLUDES = -I/usr/local/include -I/usr/include
else
ifeq (${machine}, IRIX)	
   LDEC = -shared
   INCLUDES = -I/usr/local/include
endif
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

OBJS_NPS   = nps.o

.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${F77} -c ${FFLAGS} $*.f -o $*.o
.F.f:  
	${CPP} -P $*.F > $*.f
	${F77} -c  ${FFLAGS} $*.f -o $*.o
	${RM} $*.f
	${RCSRM} $*.F


all: TRANS_NPS

TRANS_NPS: ${OBJS_NPS}
	${LD} ${LDEC} ${LSUN} ${OBJS_NPS} -o trans_nps.o \
	${LIBS}
	/bin/rm ${OBJS_NPS}

clean:
	${RM} *.o *.mod 
 