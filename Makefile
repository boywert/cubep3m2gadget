all: convert read
convert: cubep3m2gadget.f90 intbitmask.c
	mpicc -c intbitmask.c
	mpif90 -c -free -fpp cubep3m2gadget.f90
	mpif90 -fpp intbitmask.o cubep3m2gadget.o -o convert
read: read_catalogues.c
	mpicc read_catalogues.c -o read