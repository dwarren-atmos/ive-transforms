# $Id: Makefile,v 1.4 1997/01/07 00:40:52 cepi Exp $
#
# $Log: Makefile,v $
# Revision 1.4  1997/01/07 00:40:52  cepi
# Added new routines to object list.
#
# Revision 1.3  1996/10/25 19:53:07  harry
# Add support for SGI IRIX.
#

LD = ld
F77 = f77
CC = cc
CPP = /usr/lib/cpp

machine := ${shell /usr/bin/arch}
ifndef machine
machine := ${shell /usr/bin/uname}
endif
ifndef machine
machine := ${shell /bin/uname}
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
   F77=ifort
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
RM = -/bin/rm -f
RCSRM = -/usr/local/bin/rcsclean
#
DEBUG = 
CFLAGS = ${OPT} ${SUN} ${DEC} -D_NO_PROTO ${INCLUDES}
FFLAGS = ${OPT} ${SUN} ${DEC} ${DEBUG} -u

OBJS_TER  = default.o calc_field.o deriv.o diff.o buoy_freq.o \
            def_cmp.o mix_tpts.o mix_wpts.o rho.o press_drag.o \
            per_vol.o pert.o vol_int.o get_args.o run_int.o \
            area_int_z.o area_int_x.o bar.o units_smpfy.o grid.o \
            area_int_y.o list.o

.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${F77} -c ${FFLAGS} $*.f -o $*.o
.F.o:  
	${F77} -c ${FFLAGS} $*.F -o $*.o


all: TRANS_TER

TRANS_TER: ${OBJS_TER}
	${LD} ${LDEC} ${LSUN} ${OBJS_TER} -o trans_ter.o \
	${LIBS}
	/bin/rm ${OBJS_TER}

clean:
	${RM} ${OBJS_TER}

${OBJS_TER}  : default.icl 

