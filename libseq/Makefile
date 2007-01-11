# This file is part of MUMPS VERSION 4.6.4
# This Version was built on Thu Jan 11 13:32:35 2007
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
