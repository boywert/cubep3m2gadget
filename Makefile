all: cubep3m2gadget read
cubep3m2gadget: cubep3m2gadget.f90 
	mpif90 -fpp cubep3m2gadget.f90 -o cubep3m2gadget
read: read_catalogues.c
	mpicc read_catalogues.c -o read