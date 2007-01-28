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
F90 = f90
CC = cc
CPP = /usr/lib/cpp

machine := ${shell /usr/bin/arch}
ifeq (${machine},alpha)
   CC = c89
   DEC =  -common -DSIXTY_FOUR_bptr 
   LDEC = -shared 
   LIBS = -L/usr/local/matlab4.2/extern/lib/alpha -lmat -lfftpack
   OPT = -O2
   AUTOMATIC = -automatic
   INCLUDES = -I /usr/local/include
else
   LD=icc
   F90=ifort
   FC=ifort
   CC = icc
   DEC = 
   LDEC = -shared
   LIBS = -lfftpack 
   OPT = -O2
   AUTOMATIC = -automatic
   INCLUDES = -I/usr/include -I/usr/local/include
endif
#
RM = -/bin/rm -f
RCSRM = -/usr/local/bin/rcsclean
#
DEBUG = 
CFLAGS = ${OPT} ${SUN} ${DEC} -D_NO_PROTO ${INCLUDES}
FFLAGS = ${OPT} ${SUN} ${DEC} ${AUTOMATIC} ${DEBUG} ${INCLUDES} -u
#
OBJS_TER =	\
	handle_err.o	\
	get_args.o		\
 	z2zint.o	\
	s2pint.o		\
	units_smpfy.o	\
    ij2ll.o		\
	ll2ij.o		\
	hm2uv.o		\
	grid.o		 	\
	grdij.o		\
	map_factors.o	\
	default.o

#		 calc_field_coamps.o   \
	     buoy_freq.o		\
	     deriv.o		\
	     deriv_nsew.o		\
	     geopot.o		\
	     geowind.o		\
	     ive_funcs.o        \
	     list.o             \
	     mkfield.o		\
	     metrics.o		\
	     misc.o			\
	     move_parcel.o	\
	     omega.o		\
	     uvrot.o		\
	     uvrot2.o		\
	     vert_grid.o		\
	     x_grid.o		\
	     y_grid.o		\
	     massflx_lat.o 	\
	     massflx_max.o 	\
	     massflx_top.o	\
	     rho.o			\
	     interp.o		\
	     interp_field.o	\
	     stagering_calc.o	\
	     test.o			\
	     uvg2uv.o		\
	     windspeed.o	\
	     reldiff.o	\
	     time_diff.o	\
	     variance.o	\
# 
.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${F90} -c ${FFLAGS} $*.f -o $*.o
.F.f:  
	${F90} -c  ${FFLAGS} $*.F -o $*.o

all: TRANS_TER

TRANS_TER: ${OBJS_TER}
	${LD} ${LDEC} ${LSUN} ${OBJS_TER} -o trans_coamps_prl.o ${LIBS}

clean:
	${RM} ${OBJS_TER} trans_coamps_prl.o

${OBJS_TER}  : default.icl 

