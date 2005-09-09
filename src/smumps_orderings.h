/*

   THIS FILE IS PART OF MUMPS VERSION 4.5.2
   This Version was built on Fri Sep  9 16:03:33 2005

  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999). 
  It has also been supported by CERFACS, ENSEEIHT-IRIT-CNRS, 
  LIP-INRIA-Rhone-Alpes, and PARALLAB.

  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

  Up-to-date copies of the MUMPS package can be obtained
  from the Web pages http://www.enseeiht.fr/apo/MUMPS/
  or http://graal.ens-lyon.fr/MUMPS


   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.


  User documentation of any code that uses this software can
  include this complete notice. You can acknowledge (using
  references [1], [2], and [3] the contribution of this package
  in any scientific publication dependent upon the use of the
  package. You shall use reasonable endeavours to notify
  the authors of the package of this publication.

   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
   Multifrontal parallel distributed symmetric and unsymmetric solvers,
   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).

   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
   A fully asynchronous multifrontal solver using distributed dynamic
   scheduling, SIAM Journal of Matrix Analysis and Applications,
   Vol 23, No 1, pp 15-41 (2001).

   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
   S. Pralet (2005), Hybrid scheduling for the parallel solution
   of linear systems. Accepted to Parallel Computing.

*/
/* $Id: smumps_orderings.h,v 1.2 2003/04/25 12:25:46 jylexcel Exp $ */
int esmumps (const int n, const int iwlen, int * const pe, const int pfree, int * const len, int * const iw, int * const nv, int * const elen, int * const last);

int   smumps_pord(int, int, int *, int *, int *);

#define FORTRAN(nu,nl,pl,pc)                     \
void nu ();                                      \
void nl pl                                       \
{ nu pc; }                                       \
void nl##_ pl                                    \
{ nu pc; }                                       \
void nl##__ pl                                   \
{ nu pc; }                                       \
void nu pl
