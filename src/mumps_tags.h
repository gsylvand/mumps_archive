C
C   THIS FILE IS PART OF MUMPS VERSION 4.3.2
C   This Version was built on Wed Nov 12 16:57:09 2003
C
C COPYRIGHT (c) 1996-2003 P. R. Amestoy, I. S. Duff, J. Koster,
C                       J.-Y. L'Excellent
C
C  CERFACS      , Toulouse    (France)  (http://www.cerfacs.fr)
C  ENSEEIHT-IRIT, Toulouse    (France)  (http://www.enseeiht.fr)
C  INRIA                      (France)  (http://www.inria.fr)
C  PARALLAB     , Bergen      (Norway)  (http://www.parallab.uib.no)
C
C All rights reserved.
C
C  Your use or distribution of the package implies that you agree
C  with this License. Up-to-date copies of the MUMPS package can be
C  obtained from the Web page http://www.enseeiht.fr/apo/MUMPS/
C
C  This package is provided to you free of charge. It was
C  initially based on public domain software developed during
C  the European Esprit IV project PARASOL (1996-1999).
C
C  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
C
C  Permission is hereby granted to use or copy this
C  package provided that the Copyright and this License is
C  retained on all copies and that the package is used
C  under the same terms and conditions. User documentation
C  of any code that uses this software should include this
C  complete Copyright notice and this License.
C
C  You can modify this code but, at no time shall the right
C  or title to all or any part of this package pass to you.
C  All information relating to any alteration or addition
C  made to this package for the purposes of extending the
C  capabilities or enhancing the performance of this package
C  shall be made available free of charge to the authors for
C  any purpose.
C
C  You shall acknowledge (using references [1] and [2])
C  the contribution of this package in any publication
C  of material dependent upon the use of the package.
C  You shall use reasonable endeavours to notify
C  the authors of the package of this publication.
C
C
C  [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
C  Multifrontal parallel distributed symmetric and unsymmetric solvers,
C  in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
C  An early version appeared as a Technical Report ENSEEIHT-IRIT (1998)
C  and is available at http://www.enseeiht.fr/apo/MUMPS/.
C
C  [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
C  A fully asynchronous multifrontal solver using distributed dynamic
C  scheduling, SIAM Journal of Matrix Analysis and Applications,
C  Vol 23, No 1, pp 15-41 (2001).
C  An  early version appeared as a Technical Report ENSEEIHT-IRIT,
C  RT/APO/99/2 (1999)  and is available at http://www.enseeiht.fr/apo/MUMPS/.
C
C  None of the text from the Copyright notice up to and
C  including this line shall be removed or altered in any way.
C
C     $Id: mumps_tags.h,v 1.3 2003/04/25 12:25:46 jylexcel Exp $
C -----------------------------------------
C  This file contains the definition
C  of all tags.
C -----------------------------------------
C
C -----------------------------------------
C  Tag for arrowheads distribution
C -----------------------------------------
      INTEGER ARROWHEAD, ARR_INT, ARR_REAL, ELT_INT, ELT_REAL
      PARAMETER ( ARROWHEAD = 20,
     *            ARR_INT = 29,
     *            ARR_REAL = 30,
     *            ELT_INT = 31,
     *            ELT_REAL = 32 )

C ----------------------------------------------------
C   Tags for collecting distributed integer info
C   for analysis in case of initial distributed matrix
C ----------------------------------------------------
      INTEGER COLLECT_NZ, COLLECT_IRN, COLLECT_JCN
      PARAMETER( COLLECT_NZ  = 35,
     *           COLLECT_IRN = 36,
     *           COLLECT_JCN = 37 )

C -----------------------------------------
C   Tags for factorization
C -----------------------------------------
      INTEGER RACINE,
     *        NOEUD,
     *        TERREUR,
     *        MAITRE_DESC_BANDE,
     *        MAITRE2,
     *        BLOC_FACTO,
     *        CONTSIPERENIV2,
     *        MAPLIG,
     *        FACTOR
      PARAMETER ( RACINE            = 2,
     *            NOEUD             = 3,
     *            MAITRE_DESC_BANDE = 4,
     *            MAITRE2           = 5,
     *            BLOC_FACTO        = 6,
     *            CONTSIPERENIV2    = 7,
     *            MAPLIG            = 8,
     *            FACTOR            = 9,
     *            TERREUR           = 99 )

C -----------------------------------------
C   Tags for assembly of root (in facto)
C -----------------------------------------
      INTEGER ROOT_NELIM_INDICES,
     *        ROOT_CONT_STATIC,
     *        ROOT_NON_ELIM_CB,
     *        ROOT_2SLAVE,
     *        ROOT_2SON
       PARAMETER( ROOT_NELIM_INDICES = 15,
     *        ROOT_CONT_STATIC       = 16,
     *        ROOT_NON_ELIM_CB       = 17,
     *        ROOT_2SLAVE            = 18,
     *        ROOT_2SON              = 19 )

C -----------------------------------------
C   Tags for solve
C -----------------------------------------
      INTEGER RACINE_SOLVE,
     *        ContVec,
     *        Master2Slave
      PARAMETER( RACINE_SOLVE = 10,
     *           ContVec      = 11,
     *           Master2Slave = 12 )

C -----------------------------------------
C   Tags for backsolve
C -----------------------------------------
      INTEGER FEUILLE,
     *        BACKSLV_UPDATERHS,
     *        BACKSLV_MASTER2SLAVE
      PARAMETER( FEUILLE = 21,
     *           BACKSLV_UPDATERHS = 22,
     *           BACKSLV_MASTER2SLAVE = 23 )
C ------------------------
C   Tag for symmetrization
C ------------------------
      INTEGER SYMMETRIZE
      PARAMETER ( SYMMETRIZE = 24 )

C ----------------------------
C   Tags specific to symmetric
C ----------------------------
      INTEGER BLOC_FACTO_SYM,
     *        BLOC_FACTO_SYM_SLAVE, END_NIV2_LDLT

      PARAMETER ( BLOC_FACTO_SYM = 25,
     *            BLOC_FACTO_SYM_SLAVE = 26 , 
     *            END_NIV2_LDLT = 33 )
C -------------------------------------
C   Tags specific to dynamic scheduling
C -------------------------------------
      INTEGER UPDATE_LOAD
      PARAMETER ( UPDATE_LOAD = 27 )

C   To send deficientcy
      INTEGER DEFIC_TAG
      PARAMETER(  DEFIC_TAG = 28 )
C   To send Schur
      INTEGER TAG_SCHUR
      PARAMETER( TAG_SCHUR = 38 )
C   To clean up IRECV
      INTEGER TAG_DUMMY
      PARAMETER( TAG_DUMMY = 39 )
