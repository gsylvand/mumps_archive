# This file is part of MUMPS VERSION 4.4.2
# This Version was built on Mon May 23 09:41:55 2005
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
