# This file is part of MUMPS VERSION 4.5.4
# This Version was built on Wed Sep 28 15:22:28 2005
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
