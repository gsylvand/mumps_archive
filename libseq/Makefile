# This file is part of MUMPS VERSION 4.3.2
# This Version was built on Wed Nov 12 16:57:09 2003
# COPYRIGHT (C) 1996-2003 P. R. Amestoy, I. S. Duff, J. Koster, J.-Y. L Excellent
#
#
all: libmpiseq.a

include ../Makefile.inc

libmpiseq.a: mpi.o mpic.o elapse.o
	$(AR) libmpiseq.a mpi.o mpic.o elapse.o
.f.o:
	$(FC) $(OPTF) -c $*.f
.c.o:
	$(CC) $(OPTC) -I. -c $*.c

clean:
	$(RM) *.o *.a
