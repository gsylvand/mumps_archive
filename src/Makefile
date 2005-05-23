# This file is part of MUMPS VERSION 4.4.2
# This Version was built on Mon May 23 09:41:55 2005
#
#
.PHONY: double simple cmplx cmplx16

double:../lib/libdmumps.a
simple:../lib/libsmumps.a
cmplx16:../lib/libzmumps.a
cmplx:../lib/libcmumps.a

include ../Makefile.inc

.SUFFIXES: .F .o
.F.o:
	$(FC) $(OPTF) $(INC) $(ORDERINGSF) -I. -I../include -c $*.F
.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) $(CDEFS) $(IORDERINGS) $(ORDERINGSC) -I../include -c $*.c
clean:
	$(RM) *.o
################################################################################
#                            Double precision
################################################################################
../lib/libdmumps.a:dmumps_part1.o dmumps_part2.o dmumps_part3.o dmumps_part4.o dmumps_part5.o dmumps_part6.o dmumps_part7.o dmumps_part8.o dmumps_part9.o dmumps_buffer.o dmumps_bloc2.o dmumps_cv_load.o dmumps_cv.o dmumps_distribution.o dmumps_c.o dmumps_orderings.o
	$(AR) $@ $?
	$(RANLIB) $@

dmumps_cv_load.o  : dmumps_buffer.o
dmumps_buffer.o: dmumps_bloc2.o
dmumps_part1.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o
dmumps_part2.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o
dmumps_part3.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o
dmumps_part4.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o
dmumps_part5.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_cv.o
dmumps_part6.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o
dmumps_part7.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o
dmumps_part8.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o
dmumps_part9.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o


################################################################################
#                            Simple precision
################################################################################
../lib/libsmumps.a:smumps_part1.o smumps_part2.o smumps_part3.o smumps_part4.o smumps_part5.o smumps_part6.o smumps_part7.o smumps_part8.o smumps_part9.o smumps_buffer.o smumps_bloc2.o smumps_cv_load.o smumps_cv.o smumps_distribution.o smumps_c.o smumps_orderings.o
	$(AR) $@ $?
	$(RANLIB) $@

smumps_cv_load.o  : smumps_buffer.o
smumps_buffer.o: smumps_bloc2.o
smumps_part1.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o
smumps_part2.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o
smumps_part3.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o
smumps_part4.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o
smumps_part5.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_cv.o
smumps_part6.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o
smumps_part7.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o
smumps_part8.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o
smumps_part9.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o


################################################################################
#                            Complex precision
################################################################################
../lib/libcmumps.a:cmumps_part1.o cmumps_part2.o cmumps_part3.o cmumps_part4.o cmumps_part5.o cmumps_part6.o cmumps_part7.o cmumps_part8.o cmumps_part9.o cmumps_buffer.o cmumps_bloc2.o cmumps_cv_load.o cmumps_cv.o cmumps_distribution.o cmumps_c.o cmumps_orderings.o
	$(AR) $@ $?
	$(RANLIB) $@

cmumps_cv_load.o  : cmumps_buffer.o
cmumps_buffer.o: cmumps_bloc2.o
cmumps_part1.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o
cmumps_part2.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o
cmumps_part3.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o
cmumps_part4.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o
cmumps_part5.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_cv.o
cmumps_part6.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o
cmumps_part7.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o
cmumps_part8.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o
cmumps_part9.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o


################################################################################
#                            Double complex precision
################################################################################
../lib/libzmumps.a:zmumps_part1.o zmumps_part2.o zmumps_part3.o zmumps_part4.o zmumps_part5.o zmumps_part6.o zmumps_part7.o zmumps_part8.o zmumps_part9.o zmumps_buffer.o zmumps_bloc2.o zmumps_cv_load.o zmumps_cv.o zmumps_distribution.o zmumps_c.o zmumps_orderings.o
	$(AR) $@ $?
	$(RANLIB) $@

zmumps_cv_load.o  : zmumps_buffer.o
zmumps_buffer.o: zmumps_bloc2.o
zmumps_part1.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o
zmumps_part2.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o
zmumps_part3.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o
zmumps_part4.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o
zmumps_part5.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_cv.o
zmumps_part6.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o
zmumps_part7.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o
zmumps_part8.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o
zmumps_part9.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o


