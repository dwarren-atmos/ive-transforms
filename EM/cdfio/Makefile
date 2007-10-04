# simple makefile
#
CC = cc
F77 = f77
INCLUDEDIR = /usr/local/include
LIBRARY = libcdfio.a
FFLAGS = 
.F.o:
	${F77} -c -I${INCLUDEDIR} ${FFLAGS} $*.F

${LIBRARY}: libcdfio.o
	ar rv libcdfio.a libcdfio.o
	ranlib libcdfio.a

