CC=gcc
CXX=g++
CFLAGS=-fopenmp -g -O2
CXXFLAGS=-fopenmp -g -O2

all:
	make mxv
	make vpv
	make cgsolve

mxv:
	${CC} ${CFLAGS} -o mxv mxv_openmp.c

vpv:
	${CC} ${CFLAGS} -o vpv vpv_openmp.c

cgsolve: ./cg-solve/generate_matrix.hpp
	${CXX} ${CXXFLAGS} -o cgsolve ./cg-solve/cgsolve_omp.cpp

.PHONY: clean

clean:
	rm mxv
	rm vpv
	rm cgsolve
