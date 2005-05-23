/*

   THIS FILE IS PART OF MUMPS VERSION 4.4.2
   This Version was built on Mon May 23 09:41:55 2005

  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  european project PARASOL (1996-1999). It has been partially supported
  by the European Community, and by CERFACS, ENSEEIHT-IRIT, INRIA
  Rhone-Alpes, and LBNL.

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
   of linear systems. Submitted to Parallel Computing.

*/
/* $Id: dmumps_c.h,v 1.21 2005/04/18 10:58:44 jylexcel Exp $ */
/* Mostly written in march 2002 (JYL) */

#if ( ! defined DMUMPS_C_H )
#define DMUMPS_C_H

/* Complex datatypes */
typedef struct {float r,i;} dmumps_complex;
typedef struct {double r,i;} dmumps_double_complex;

/* Next line defines F_INT, F_DOUBLE and F_DOUBLE2 */
#include "dmumps_prec.h"
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
    /* RHS, solution, ouptput data and statistics */
    F_DOUBLE *rhs, *rhs_sparse, *sol_loc;
    F_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    F_INT nrhs, lrhs, nz_rhs, lsol_loc;
  F_INT schur_mloc, schur_nloc, schur_lld;
  F_INT mblock, nblock, nprow, npcol;
    F_INT info[40],infog[40];
    F_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    F_INT deficiency; F_DOUBLE * nullspace; F_INT * mapping;
    /* Schur */
    F_INT size_schur; F_INT *listvar_schur; F_DOUBLE *schur;
    /* Internal parameters */
    F_INT instance_number;
  } DMUMPS_STRUC_C;


#if defined(UPPER)
  #define dmumps_f77_ DMUMPS_F77
  #define dmumps_affect_mapping_ DMUMPS_AFFECT_MAPPING
  #define dmumps_affect_nullspace_ DMUMPS_AFFECT_NULLSPACE
  #define dmumps_affect_colsca_ DMUMPS_AFFECT_COLSCA
  #define dmumps_affect_rowsca_ DMUMPS_AFFECT_ROWSCA 
  #define dmumps_affect_uns_perm_     DMUMPS_AFFECT_UNS_PERM
  #define dmumps_affect_sym_perm_     DMUMPS_AFFECT_SYM_PERM
  #define dmumps_nullify_c_mapping_   DMUMPS_NULLIFY_C_MAPPING
  #define dmumps_nullify_c_nullspace_ DMUMPS_NULLIFY_C_NULLSPACE
  #define dmumps_nullify_c_sym_perm_  DMUMPS_NULLIFY_C_SYM_PERM
  #define dmumps_nullify_c_uns_perm_  DMUMPS_NULLIFY_C_UNS_PERM
  #define dmumps_nullify_c_colsca_    DMUMPS_NULLIFY_C_COLSCA
  #define dmumps_nullify_c_rowsca_    DMUMPS_NULLIFY_C_ROWSCA
#elif defined(Add__)
  #define dmumps_f77_ dmumps_f77__
  #define dmumps_affect_mapping_ dmumps_affect_mapping__
  #define dmumps_affect_nullspace_ dmumps_affect_nullspace__
  #define dmumps_affect_colsca_ dmumps_affect_colsca__
  #define dmumps_affect_rowsca_ dmumps_affect_rowsca__
  #define dmumps_affect_uns_perm_     dmumps_affect_uns_perm__     
  #define dmumps_affect_sym_perm_     dmumps_affect_sym_perm__     
  #define dmumps_nullify_c_mapping_   dmumps_nullify_c_mapping__    
  #define dmumps_nullify_c_nullspace_ dmumps_nullify_c_nullspace__  
  #define dmumps_nullify_c_sym_perm_  dmumps_nullify_c_sym_perm__   
  #define dmumps_nullify_c_uns_perm_  dmumps_nullify_c_uns_perm__   
  #define dmumps_nullify_c_colsca_    dmumps_nullify_c_colsca__     
  #define dmumps_nullify_c_rowsca_    dmumps_nullify_c_rowsca__    
#elif defined(Add_)
/* Default. Nothing to do */
#else
/* Name without underscore is used */
  #define dmumps_f77_ dmumps_f77
  #define dmumps_affect_mapping_ dmumps_affect_mapping
  #define dmumps_affect_nullspace_ dmumps_affect_nullspace
  #define dmumps_affect_colsca_ dmumps_affect_colsca
  #define dmumps_affect_rowsca_ dmumps_affect_rowsca
  #define dmumps_affect_uns_perm_     dmumps_affect_uns_perm     
  #define dmumps_affect_sym_perm_     dmumps_affect_sym_perm     
  #define dmumps_nullify_c_mapping_   dmumps_nullify_c_mapping    
  #define dmumps_nullify_c_nullspace_ dmumps_nullify_c_nullspace  
  #define dmumps_nullify_c_sym_perm_  dmumps_nullify_c_sym_perm   
  #define dmumps_nullify_c_uns_perm_  dmumps_nullify_c_uns_perm   
  #define dmumps_nullify_c_colsca_    dmumps_nullify_c_colsca     
  #define dmumps_nullify_c_rowsca_    dmumps_nullify_c_rowsca
#endif

void dmumps_c(DMUMPS_STRUC_C * dmumps_par);
void dmumps_affect_mapping_(F_INT * f77mapping);
void dmumps_affect_nullspace_(F_DOUBLE * f77nullspace);
void dmumps_affect_uns_perm_(F_INT * f77sym_perm);
void dmumps_affect_sym_perm_(F_INT * f77uns_perm);
void dmumps_nullify_c_mapping_();
void dmumps_nullify_c_nullspace_();
void dmumps_nullify_c_sym_perm_();
void dmumps_nullify_c_uns_perm_();
#ifdef return_scaling
void dmumps_affect_colsca_(F_DOUBLE * f77colsca);
void dmumps_affect_rowsca_(F_DOUBLE * f77rowsca);
void dmumps_nullify_c_colsca_();
void dmumps_nullify_c_rowsca_();
#endif

void dmumps_f77_(F_INT *job, F_INT *sym, F_INT *par, F_INT * comm_fortran,
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
  F_DOUBLE * rowsca, F_INT * rowsca_avail, F_INT * instance_number,
  F_INT * nrhs, F_INT * lrhs, F_DOUBLE * rhs_sparse, F_INT * rhs_sparse_avail,
  F_DOUBLE * sol_loc, F_INT * sol_loc_avail, F_INT * irhs_sparse,
  F_INT * irhs_sparse_avail, F_INT * irhs_ptr, F_INT * irhs_ptr_avail,
  F_INT * isol_loc, F_INT * isol_loc_avail, F_INT * nz_rhs, F_INT * lsol_loc
  ,F_INT * schur_mloc, F_INT *schur_nloc, F_INT * schur_lld,
  F_INT * schur_mblock, F_INT * schur_nblock, F_INT * schur_nprow,
  F_INT * schur_npcol
  );
#endif

