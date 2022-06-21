CC=gcc
CFLAGS=-fopenmp -g -O2

all:
	make mxv
	make vpv

mxv:
	${CC} ${CFLAGS} -o mxv mxv_openmp.c

vpv:
	${CC} ${CFLAGS} -o vpv vpv_openmp.c

.PHONY: clean

clean:
	rm mxv
	rm vpv
