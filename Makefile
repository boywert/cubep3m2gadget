all: convert read
convert: cubep3m2gadget.f90 intbitmask.c
	mpicc -c intbitmask.c
	mpif90  intbitmask.o cubep3m2gadget.f90 -o convert
read: read_catalogues.c
	mpicc read_catalogues.c -o read