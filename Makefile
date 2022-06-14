### Compilers r###
CC=g++## for C++
FF=gfortran## for Fortran

###  Directories ###
SRCDIR=${PWD}/src
INCDIR=${PWD}/src
INCHEP=-I/ada1/lpthe/fuks/hubble/benj_local/include
LIBDIR=${PWD}/lib
BLDDIR=${PWD}/bld

### Flags ###
CFLAG=-O3 -Wall --std=c++14 -I${INCDIR} ${INCHEP}
FFLAG=-O3 -I${INCDIR} -I${INCLT}

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
STLIB= -lm -lgfortran 
HEPTLS= -L/ada1/lpthe/fuks/hubble/benj_local/lib -L/ada1/lpthe/fuks/hubble/benj_local/lib64 -lLHAPDF -looptools

### Commands ###
all: RUN
lib: ${LIB}

RUN: main.cpp ${LIB}
	${CC} ${CFLAG} -o $@ main.cpp ${LIB} ${GSLIB} ${HEPTLS} ${STLIB} 

${LIB}:	${COBJS} ${FOBJS}
	ar -ruc $@ ${BLDDIR}/*.o

${BLDDIR}%.o: ${SRCDIR}%.f
	cd ${BLDDIR}; ${FF} -c ${FFLAG} $?;

${BLDDIR}%.o: ${SRCDIR}%.cpp
	cd ${BLDDIR}; ${CC} -c ${CFLAG} $<;
	
### Cleaning ###
clean:
	rm -f RUN ${LIBDIR}/*.a ${BLDDIR}/*.o *.log output/*;

