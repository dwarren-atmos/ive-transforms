# $Id: Makefile,v 1.4 1997/01/07 00:40:52 cepi Exp $
#
# $Log: Makefile,v $
# Revision 1.4  1997/01/07 00:40:52  cepi
# Added new routines to object list.
#
# Revision 1.3  1996/10/25 19:53:07  harry
# Add support for SGI IRIX.
#

MODLIB	= ./mods

machine := ${shell /usr/bin/arch}
MACH := ${shell uname -m}

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
   LD=ifort
   F90=ifort
   CPP = /usr/bin/cpp
   FC=ifort
   CC = icc
   CPP = ${F90} -E
ifeq (${MACH}, x86_64)
   DEC = -fPIC -DSIXTY_FOUR_bptr
   LDEC = -shared 
else
   DEC =
   LDEC = -shared
endif
   LIBS = -L/usr/lib -L/usr/local/lib -lnetcdf -lnetcdff #-lfftpack
   OPT = -O2
   AUTOMATIC = -automatic -recursive
   INCLUDES = -I/usr/include -I/usr/local/include
endif

RM = -/bin/rm -f
RCSRM = -/usr/local/bin/rcsclean

DEBUG = 
CFLAGS = ${OPT} ${DEC} -D_NO_PROTO ${INCLUDES}
FFLAGS = ${OPT} ${DEC} ${AUTOMATIC} ${DEBUG} ${INCLUDES} -module $(MODLIB) -u
#
#
#			module_numerics.o
MODS	=	default_module.o   	\
            module_mapping.o   	\
			misc_module.o

OBJS_F =  default.o rho.o calc_field.o regrid.o

#calc_field.o run_sub.o \
#	   metrics.o omega.o \
#	  destagger.o units_smpfy.o
#	  buoy_freq.o  ive_funcs.o  move_parcel.o\
#	  mkfield.o \
#	  vert_grid.o xy_grid.o rho.o interp_pres.o interp_height.o \
#	  vert_interp_field.o horiz_interp_field.o the2temp.o \
#	  max_min.o point_val.o vort.o filt9.o

OBJS_f = grdij.o grid.o hm2uv.o ij2ll.o ll2ij.o map_factors.o uvg2uv.o z2zint.o s2pint.o \
	deriv.o list.o massflux.o get_args.o

all: TRANS_WRF

OBJS_F: ${OBJS_F} ${MODS}
OBJS_f: ${OBJS_f}

TRANS_WRF: ${MODS} ${OBJS_f} ${OBJS_F}
	${LD} ${LDEC} ${OBJS_f} ${OBJS_F} ${MODS} -o trans_wrf.o ${LIBS}
#	cp trans_wrf.o ../

clean:
	${RM} ${OBJS_F} ${OBJS_f} ${MODS} ${MODLIB}/*.mod trans_wrf.o

${OBJS_F}  : ${MODS}

.f.o:
	$(F90) $(FFLAGS) -c $<
.F.o:
	$(F90) $(FFLAGS) -c $<
.c.o:
	${CC}  -c ${CFLAGS} $*.c -o $*.o
