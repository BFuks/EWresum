### Compilers r###
CC=g++## for C++
FF=gfortran## for Fortran

###  Directories ###
SRCDIR=${PWD}/src
INCDIR=${PWD}/src
LIBDIR=${PWD}/lib
BLDDIR=${PWD}/bld

### Flags ###
CFLAG=-O3 -Wall --std=c++14 -I${INCDIR}
FFLAG=-O3 -I${INCDIR}

### Paths ###
VPATH=${SRCDIR}

### Source files ###
CFILES=$(wildcard ${SRCDIR}/*.cpp)
FFILES=$(wildcard ${SRCDIR}/*.f)

### Object files ###
COBJS=$(subst .cpp,.o,$(subst ${SRCDIR},${BLDDIR},${CFILES}))
FOBJS=$(subst .f,.o,$(subst ${SRCDIR},${BLDDIR},${FFILES}))

### Libraries ###
LIB=${LIBDIR}/libresum.a
GSLIB= -lgsl -lgslcblas
STLIB= -L/usr/local/gfortran/lib -lm -lgfortran 
LHAPDF= -L/usr/local/lib -lLHAPDF

### Commands ###
all: RUN
lib: ${LIB}

RUN: main.cpp ${LIB}
	${CC} ${CFLAG} -o $@ main.cpp ${LIB} ${GSLIB} ${STLIB} ${LHAPDF}

${LIB}:	${COBJS} ${FOBJS}
	ar -ruc $@ ${BLDDIR}/*.o

${BLDDIR}%.o: ${SRCDIR}%.f
	cd ${BLDDIR}; ${FF} -c ${FFLAG} $?;

${BLDDIR}%.o: ${SRCDIR}%.cpp
	cd ${BLDDIR}; ${CC} -c ${CFLAG} $<;
	
### Cleaning ###
clean:
	rm -f RUN ${LIBDIR}/*.a ${BLDDIR}/*.o *.log output/*;
