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
RM = /bin/rm -f

#
DEBUG = 
CFLAGS = ${OPT} ${SUN} ${DEC} -D_NO_PROTO ${INCLUDES}
FFLAGS = ${OPT} ${SUN} ${DEC} ${DEBUG} -u

OBJS_UARS  = uars.o 

.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${F77} -c ${FFLAGS} $*.f -o $*.o
.F.o:  
	${F77} -c  ${FFLAGS} $*.F -o $*.o


all: TRANS_UARS

TRANS_UARS: ${OBJS_UARS}
	${LD} ${LDEC} ${LSUN} ${OBJS_UARS} -o trans_uars.o \
	${LIBS}
	/bin/rm ${OBJS_UARS}

clean:
	${RM} ${OBJS_UARS}

${OBJS_UARS}  : default.icl 

