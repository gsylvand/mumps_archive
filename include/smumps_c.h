/*

   THIS FILE IS PART OF MUMPS VERSION 4.3.2
   This Version was built on Wed Nov 12 16:57:09 2003

 COPYRIGHT (c) 1996-2003 P. R. Amestoy, I. S. Duff, J. Koster,
                       J.-Y. L'Excellent

  CERFACS      , Toulouse    (France)  (http://www.cerfacs.fr)
  ENSEEIHT-IRIT, Toulouse    (France)  (http://www.enseeiht.fr)
  INRIA                      (France)  (http://www.inria.fr)
  PARALLAB     , Bergen      (Norway)  (http://www.parallab.uib.no)

 All rights reserved.

  Your use or distribution of the package implies that you agree
  with this License. Up-to-date copies of the MUMPS package can be
  obtained from the Web page http://www.enseeiht.fr/apo/MUMPS/

  This package is provided to you free of charge. It was
  initially based on public domain software developed during
  the European Esprit IV project PARASOL (1996-1999).

  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

  Permission is hereby granted to use or copy this
  package provided that the Copyright and this License is
  retained on all copies and that the package is used
  under the same terms and conditions. User documentation
  of any code that uses this software should include this
  complete Copyright notice and this License.

  You can modify this code but, at no time shall the right
  or title to all or any part of this package pass to you.
  All information relating to any alteration or addition
  made to this package for the purposes of extending the
  capabilities or enhancing the performance of this package
  shall be made available free of charge to the authors for
  any purpose.

  You shall acknowledge (using references [1] and [2])
  the contribution of this package in any publication
  of material dependent upon the use of the package.
  You shall use reasonable endeavours to notify
  the authors of the package of this publication.


  [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
  Multifrontal parallel distributed symmetric and unsymmetric solvers,
  in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
  An early version appeared as a Technical Report ENSEEIHT-IRIT (1998)
  and is available at http://www.enseeiht.fr/apo/MUMPS/.

  [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
  A fully asynchronous multifrontal solver using distributed dynamic
  scheduling, SIAM Journal of Matrix Analysis and Applications,
  Vol 23, No 1, pp 15-41 (2001).
  An  early version appeared as a Technical Report ENSEEIHT-IRIT,
  RT/APO/99/2 (1999)  and is available at http://www.enseeiht.fr/apo/MUMPS/.

  None of the text from the Copyright notice up to and
  including this line shall be removed or altered in any way.

*/
/* $Id: smumps_c.h,v 1.14 2003/04/25 12:25:46 jylexcel Exp $ */
/* Mostly written in march 2002 (JYL) */

#if ( ! defined SMUMPS_C_H )
#define SMUMPS_C_H

/* Complex datatypes */
typedef struct {float r,i;} smumps_complex;
typedef struct {double r,i;} smumps_double_complex;

/* Next line defines F_INT, F_DOUBLE and F_DOUBLE2 */
#include "smumps_prec.h"
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    F_INT sym, par, job;
    F_INT comm_fortran;    /* Fortran communicator */
    F_INT icntl[40];
    F_DOUBLE2 cntl[5];
    F_INT n;
    /* Assembled entry */
    F_INT nz; F_INT *irn; F_INT *jcn; F_DOUBLE *a;
    /* Distributed entry */
    F_INT nz_loc; F_INT *irn_loc; F_INT *jcn_loc; F_DOUBLE *a_loc;
    /* Element entry */
    F_INT nelt; F_INT *eltptr; F_INT *eltvar; F_DOUBLE *a_elt;

    /* Ordering, if given by user */
    F_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    F_INT *sym_perm;
    /* column permutation */
    F_INT *uns_perm;

    /* Scaling (input only in this version) */
    F_DOUBLE *colsca; F_DOUBLE *rowsca;

    /* Ouptput data and statistics */
    F_DOUBLE *rhs;
    F_INT info[40],infog[40];
    F_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    F_INT deficiency; F_DOUBLE * nullspace; F_INT * mapping;
    /* Schur */
    F_INT size_schur; F_INT *listvar_schur; F_DOUBLE *schur;
    /* Internal parameters */
    F_INT instance_number;
  } SMUMPS_STRUC_C;


#if defined(UPPER)
  #define smumps_f77_ SMUMPS_F77
  #define smumps_affect_mapping_ SMUMPS_AFFECT_MAPPING
  #define smumps_affect_nullspace_ SMUMPS_AFFECT_NULLSPACE
  #define smumps_affect_colsca_ SMUMPS_AFFECT_COLSCA
  #define smumps_affect_rowsca_ SMUMPS_AFFECT_ROWSCA 
  #define smumps_affect_uns_perm_     SMUMPS_AFFECT_UNS_PERM
  #define smumps_affect_sym_perm_     SMUMPS_AFFECT_SYM_PERM
  #define smumps_nullify_c_mapping_   SMUMPS_NULLIFY_C_MAPPING
  #define smumps_nullify_c_nullspace_ SMUMPS_NULLIFY_C_NULLSPACE
  #define smumps_nullify_c_sym_perm_  SMUMPS_NULLIFY_C_SYM_PERM
  #define smumps_nullify_c_uns_perm_  SMUMPS_NULLIFY_C_UNS_PERM
  #define smumps_nullify_c_colsca_    SMUMPS_NULLIFY_C_COLSCA
  #define smumps_nullify_c_rowsca_    SMUMPS_NULLIFY_C_ROWSCA
#elif defined(Add__)
  #define smumps_f77_ smumps_f77__
  #define smumps_affect_mapping_ smumps_affect_mapping__
  #define smumps_affect_nullspace_ smumps_affect_nullspace__
  #define smumps_affect_colsca_ smumps_affect_colsca__
  #define smumps_affect_rowsca_ smumps_affect_rowsca__
  #define smumps_affect_uns_perm_     smumps_affect_uns_perm__     
  #define smumps_affect_sym_perm_     smumps_affect_sym_perm__     
  #define smumps_nullify_c_mapping_   smumps_nullify_c_mapping__    
  #define smumps_nullify_c_nullspace_ smumps_nullify_c_nullspace__  
  #define smumps_nullify_c_sym_perm_  smumps_nullify_c_sym_perm__   
  #define smumps_nullify_c_uns_perm_  smumps_nullify_c_uns_perm__   
  #define smumps_nullify_c_colsca_    smumps_nullify_c_colsca__     
  #define smumps_nullify_c_rowsca_    smumps_nullify_c_rowsca__    
#elif defined(Add_)
/* Default. Nothing to do */
#else
/* Name without underscore is used */
  #define smumps_f77_ smumps_f77
  #define smumps_affect_mapping_ smumps_affect_mapping
  #define smumps_affect_nullspace_ smumps_affect_nullspace
  #define smumps_affect_colsca_ smumps_affect_colsca
  #define smumps_affect_rowsca_ smumps_affect_rowsca
  #define smumps_affect_uns_perm_     smumps_affect_uns_perm     
  #define smumps_affect_sym_perm_     smumps_affect_sym_perm     
  #define smumps_nullify_c_mapping_   smumps_nullify_c_mapping    
  #define smumps_nullify_c_nullspace_ smumps_nullify_c_nullspace  
  #define smumps_nullify_c_sym_perm_  smumps_nullify_c_sym_perm   
  #define smumps_nullify_c_uns_perm_  smumps_nullify_c_uns_perm   
  #define smumps_nullify_c_colsca_    smumps_nullify_c_colsca     
  #define smumps_nullify_c_rowsca_    smumps_nullify_c_rowsca
#endif

void smumps_c(SMUMPS_STRUC_C * smumps_par);
void smumps_affect_mapping_(F_INT * f77mapping);
void smumps_affect_nullspace_(F_DOUBLE * f77nullspace);
void smumps_affect_uns_perm_(F_INT * f77sym_perm);
void smumps_affect_sym_perm_(F_INT * f77uns_perm);
void smumps_nullify_c_mapping_();
void smumps_nullify_c_nullspace_();
void smumps_nullify_c_sym_perm_();
void smumps_nullify_c_uns_perm_();
#ifdef return_scaling
void smumps_affect_colsca_(F_DOUBLE * f77colsca);
void smumps_affect_rowsca_(F_DOUBLE * f77rowsca);
void smumps_nullify_c_colsca_();
void smumps_nullify_c_rowsca_();
#endif

void smumps_f77_(F_INT *job, F_INT *sym, F_INT *par, F_INT * comm_fortran,
F_INT *n, F_INT *icntl,
F_DOUBLE2 *cntl, F_INT *nz, F_INT *irn, F_INT *irn_avail,
F_INT *jcn, F_INT *jcn_avail, F_DOUBLE *a, F_INT *a_avail,
  F_INT *nz_loc, F_INT * irn_loc, F_INT * irn_loc_avail,
  F_INT * jcn_loc, F_INT * jcn_loc_avail, F_DOUBLE *a_loc,
  F_INT * a_loc_avail, F_INT * nelt, F_INT * eltptr,
  F_INT * eltptr_avail, F_INT * eltvar, F_INT * eltvar_avail,
  F_DOUBLE * a_elt, F_INT * a_elt_avail, F_INT * perm_in,
  F_INT * perm_in_avail, F_DOUBLE * rhs, F_INT * rhs_avail,
  F_INT * info, F_DOUBLE2 * rinfo, F_INT * infog, F_DOUBLE2 * rinfog,
  F_INT * deficiency, F_INT * size_schur, F_INT * listvar_schur,
  F_INT * listvar_schur_avail, F_DOUBLE * schur,
  F_INT * schur_avail, F_DOUBLE * colsca, F_INT * colsca_avail,
  F_DOUBLE * rowsca, F_INT * rowsca_avail, F_INT * instance_number);
#endif
