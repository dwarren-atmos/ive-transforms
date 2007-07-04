# $Id: Makefile,v 1.4 1997/01/07 00:40:52 cepi Exp $
#
# $Log: Makefile,v $
# Revision 1.4  1997/01/07 00:40:52  cepi
# Added new routines to object list.
#
# Revision 1.3  1996/10/25 19:53:07  harry
# Add support for SGI IRIX.
#

MODLIB	= ./

machine := ${shell /usr/bin/arch}
ifeq (${machine},alpha)
   F90 = f90
   CC = cc
   CPP = /usr/lib/cpp
   LD = ld
   CC = c89
   CPP = /usr/lib/cpp
   DEC =  -common -DSIXTY_FOUR_bptr 
   LDEC = -shared 
   LIBS = -L/usr/local/matlab4.2/extern/lib/alpha -lmat -lfftpack
   OPT = -O2
   AUTOMATIC = -automatic
   INCLUDES = -I /usr/local/include
else
   LD=icc
   F90=ifort
   CPP = /usr/bin/cpp
   FC=ifort
   CC = icc
   CPP = ${F90} -E
   DEC = 
   LDEC = -shared
   LIBS = -lfftpack -lnetcdf
   OPT = -O2
   AUTOMATIC = -automatic
   INCLUDES = -I/usr/include -I/usr/local/include
endif

RM = -/bin/rm -f
RCSRM = -/usr/local/bin/rcsclean

DEBUG = 
CFLAGS = ${OPT} ${DEC} -D_NO_PROTO ${INCLUDES}
FFLAGS = ${OPT} ${DEC} ${AUTOMATIC} ${DEBUG} ${INCLUDES} -u
#
#
MODS	=	default_module.o

OBJS_F = 						\
		 default.o				\
 		 calc_field_coamps.o	\
	     buoy_freq.o			\
	     deriv.o				\
	     get_args.o				\
	     geopot.o				\
	     geowind.o				\
	     ive_funcs.o			\
	     list.o     			\
	     mkfield.o				\
	     metrics.o				\
	     destagger.o			\
	     move_parcel.o			\
	     omega.o				\
	     units_smpfy.o			\
	     uvrot2.o				\
	     vert_grid.o			\
	     xy_grid.o				\
	     massflx_lat.o			\
	     massflx_max.o			\
	     massflx_top.o			\
	     rho.o					\
	     interp_pres.o			\
	     interp_height.o		\
	     vert_interp_field.o	\
	     horiz_interp_field.o	\
	     test.o					\
	     windspeed.o			\
	     reldiff.o				\
	     time_diff.o			\
	     time_mean.o			\
	     the2temp.o				\
	     variance.o				\
		 max_min.o 				\
		 sound_out.o 			\
		 point_val.o			\
		 kgain.o				\
		 ens_cov.o

OBJS_f = grdij.o grid.o hm2uv.o ij2ll.o ll2ij.o map_factors.o s2pint.o uvg2uv.o z2zint.o

all: TRANS_TER

MODS: ${MODS}
OBJS_F: ${OBJS_F}
OBJS_f: ${OBJS_f}

TRANS_TER: ${MODS} ${OBJS_f} ${OBJS_F}
	${LD} ${LDEC} ${OBJS_f} ${OBJS_F} ${MODS} -o trans_coamps.o ${LIBS}
	cp trans_coamps.o ../

clean:
	${RM} ${OBJS_F} ${OBJS_f} ${MODS} *.mod trans_coamps.o

${OBJS_F}  : pointer.icl ${MODS}
${MODS}  : pointer.icl

.f.o:
	$(F90) $(FFLAGS) -module $(MODLIB) -c $<
.F.o:
	$(F90) -module $(MODLIB) $(FFLAGS) -c $<
.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
