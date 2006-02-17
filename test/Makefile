# This file is part of MUMPS VERSION 4.6.1
# This Version was built on Fri Feb 17 14:27:51 2006
#
#


default:double

all: simple double cmplx cmplx16

simple : ssimpletest
double : dsimpletest c_example
cmplx  : csimpletest
cmplx16: zsimpletest

include ../Makefile.inc

.SUFFIXES: .F .o
.F.o:
	$(FC) $(OPTF) $(INC) -I. -I../include -c $*.F

.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) -I. -I../include -c $*.c

clean:
	$(RM) *.o ssimpletest dsimpletest csimpletest zsimpletest c_example

dsimpletest: ../lib/libdmumps.a dsimpletest.o
	$(FL) -o dsimpletest $(OPTL) dsimpletest.o ../lib/libdmumps.a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)
c_example: ../lib/libdmumps.a c_example.o
	$(FL) -o c_example $(OPTL) c_example.o ../lib/libdmumps.a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)

ssimpletest: ../lib/libsmumps.a ssimpletest.o
	$(FL) -o ssimpletest $(OPTL) ssimpletest.o ../lib/libsmumps.a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)

zsimpletest: ../lib/libzmumps.a zsimpletest.o
	$(FL) -o zsimpletest $(OPTL) zsimpletest.o ../lib/libzmumps.a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)

csimpletest: ../lib/libcmumps.a csimpletest.o
	$(FL) -o csimpletest $(OPTL) csimpletest.o ../lib/libcmumps.a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)


../lib/libsmumps.a:
	@ echo 'Error: you should build the library ../lib/libsmumps.a first'
	exit -1
../lib/libdmumps.a:
	@ echo 'Error: you should build the library ../lib/libdmumps.a first'
	exit -1
../lib/libcmumps.a:
	@ echo 'Error: you should build the library ../lib/libcmumps.a first'
	exit -1
../lib/libzmumps.a:
	@ echo 'Error: you should build the library ../lib/libzmumps.a first'
	exit -1


