# This file is part of MUMPS VERSION 4.7.2
# This Version was built on Mon Apr 16 13:49:29 2007
#
#
.PHONY: libmpiseq

all: libmpiseq

include ../Makefile.inc

libmpiseq: libmpiseq$(PLAT).a

libmpiseq$(PLAT).a: mpi.o mpic.o elapse.o
	$(AR) libmpiseq$(PLAT).a mpi.o mpic.o elapse.o
.f.o:
	$(FC) $(OPTF) -c $*.f
.c.o:
	$(CC) $(OPTC) -I. -c $*.c

clean:
	$(RM) *.o *.a
