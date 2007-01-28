# $Id: Makefile,v 1.4 1997/01/07 00:40:52 cepi Exp $
#
# $Log: Makefile,v $
# Revision 1.4  1997/01/07 00:40:52  cepi
# Added new routines to object list.
#
# Revision 1.3  1996/10/25 19:53:07  harry
# Add support for SGI IRIX.
#


machine := ${shell /usr/bin/arch}
ifeq (${machine},alpha)
   LD = ld
   F90 = f90
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
   FC=ifort
   CC = icc
   CPP = /usr/lib/cpp
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
OBJS_TER = default.o metrics.o calc_field.o deriv.o rho.o run_int.o \
           mix.o press_drag.o press_drag_vp.o area.o area_int_z.o \
           area_int_z_vp.o area_int_x.o area_int_x_vp.o area_int_y.o \
           area_int_y_vp.o per_vol.o pert.o vol_int.o vol_int_vp.o \
           get_args.o units_smpfy.o grid.o list.o buoy_freq.o \
           ref.o abs.o diff_area_int.o diff_int_vp.o mkfield.o \
           z_ave.o xy_ave.o aves.o misc.o div.o momflx.o \
           vort.o ufcv.o z_coord.o tangent.o thres_vol.o spectra.o \
           area_int_xz0_y.o area_int_yz0_x.o vol_int_xyz0.o area_int.o \
           vort_spec.o ive_funcs.o bernoulli.o util.o massflx.o sound_out.o
# 
#		   cdf_out.o 
.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
.f.o:
	${F90} -c ${FFLAGS} $*.f -o $*.o
.F.f:  
#	${CPP} -P $*.F > $*.f
	${F90} -c  ${FFLAGS} $*.F -o $*.o

all: TRANS_TER

TRANS_TER: ${OBJS_TER}
	${LD} ${LDEC} ${LSUN} ${OBJS_TER} -o trans_ter.o ${LIBS}

clean:
	${RM} ${OBJS_TER} trans_ter.o

${OBJS_TER}  : default.icl 

