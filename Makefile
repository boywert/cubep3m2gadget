all: convert read
convert: cubep3m2gadget.f90
	mpif90 -fpp cubep3m2gadget.f90 -o convert
read: read_catalogues.c
	mpicc read_catalogues.c -o read