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
      RECURSIVE SUBROUTINE SMUMPS_210( COMM_LOAD, ASS_IRECV,
     *  BUFR, LBUFR, LBUFR_BYTES,
     *
     *  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     *  NFRONT_PERE, NASS_PERE, LMAP, TROW,
     *  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *  LRLUS, N, IW,
     *  LIW, A, LA, NIRBDU,
     *  PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     *  ICNTL, KEEP,
     *  root, OPASSW, OPELIW,
     *  ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR, ND, FRERE,
     *  LPTRAR, NELT, FRTPTR, FRTELT, 
     *
     *  ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *  )
      USE SMUMPS_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC ) :: ROOT
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ICNTL( 40 ), KEEP(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER SLAVEF, NBFIN, POSFAC
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER COMP
      INTEGER NSTK( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER INODE_PERE, ISON
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE( * )
      INTEGER LMAP 
      INTEGER TROW( LMAP )
      REAL OPASSW, OPELIW
      REAL DBLARR(MAX(1,KEEP(13)))
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER NOSLA, I, ISTCHK
      INTEGER INDICE_PERE, LREQI, LREQA
      INTEGER PDEST, PDEST_MASTER
      INTEGER NFRONT
      LOGICAL FLAG
      INTEGER NSLSON, VALSON, ROWLIST, NBCOLS, NROW,
     *        NPIV, NROWS_TO_STACK, II, COLLIST, POSROW
      INTEGER PDEST_MASTER_ISON, IPOS_IN_SLAVE
      LOGICAL DESCLU, FREEMASTER, SLAVE_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED 
      INTEGER MSGSOU, MSGTAG
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        write(*,*) myid, ' : PB allocation SLAVES_PERE in SMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      IF (NSLAVES_PERE.GT.0) 
     &SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = SMUMPS_275( STEP(INODE_PERE),
     *                 PROCNODE_STEPS, SLAVEF )
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        write(*,*) myid, ' : PB allocation NBROW in SMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        write(*,*) myid, ' : PB allocation LMAP in SMUMPS_210'
        IFLAG  =-13
        IERROR = LMAP
        GOTO 700
      endif
      MAP( 1 : LMAP ) = TROW( 1 : LMAP )
      PDEST_MASTER_ISON = SMUMPS_275(STEP(ISON),
     *                    PROCNODE_STEPS,SLAVEF)
      SLAVE_ISON = PDEST_MASTER_ISON .NE. MYID
      IF (SLAVE_ISON) THEN
        DO WHILE ( PTRIST(STEP( ISON )) .EQ. 0 )
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL SMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    PDEST_MASTER_ISON, MAITRE_DESC_BANDE,
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE, LPTRAR,
     *    NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
          IF ( IFLAG .LT. 0 ) GOTO 600
        END DO
        DO WHILE (
     *     ( IW( PTRIST(STEP(ISON)) + 1 ) .NE.
     *       IW( PTRIST(STEP(ISON)) + 3 ) ) .OR.
     *     ( KEEP(50) .NE. 0 .AND.
     *       IW( PTRIST(STEP(ISON)) + 6 ) .NE. 0 ) )
          IF ( KEEP(50).eq.0) THEN
            MSGSOU = PDEST_MASTER_ISON
            MSGTAG = BLOC_FACTO
          ELSE
            IF ( IW( PTRIST(STEP(ISON)) + 1 ) .NE.
     *           IW( PTRIST(STEP(ISON)) + 3 ) ) THEN
              MSGSOU = PDEST_MASTER_ISON
              MSGTAG = BLOC_FACTO_SYM
            ELSE
              MSGSOU = MPI_ANY_SOURCE
              MSGTAG = BLOC_FACTO_SYM_SLAVE
            END IF
          END IF
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL SMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, MSGTAG,
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE, LPTRAR,
     *    NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 600
        END DO
      ENDIF
      IF ( NSLAVES_PERE .EQ. 0 ) THEN
        NBROW( 0 ) = LMAP
      ELSE
        DO I = 0, NSLAVES_PERE
          NBROW( I ) = 0
        END DO
        DO I = 1, LMAP
          INDICE_PERE = MAP( I )
          CALL SMUMPS_47(
     *         KEEP, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *           NASS_PERE,
     *           NFRONT_PERE - NASS_PERE,
     *           NSLAVES_PERE,
     *           INDICE_PERE,
     *           NOSLA,
     *           IPOS_IN_SLAVE )
          NBROW( NOSLA ) = NBROW( NOSLA ) + 1
        END DO
      END IF
      PDEST_MASTER = SLAVES_PERE(0)
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .EQ. MYID ) THEN 
            NBPROCFILS(STEP(INODE_PERE)) =
     *                 NBPROCFILS(STEP(INODE_PERE)) - 1
            IF ( PDEST .EQ. PDEST_MASTER ) THEN
              NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
            ENDIF
            ISTCHK     = PTRIST(STEP(ISON))
            NBCOLS     = IW(ISTCHK)
            NROW       = IW(ISTCHK+2)
            NPIV       = IW(ISTCHK+3)
            NSLSON     = IW(ISTCHK+5)
            NFRONT     = NPIV + NBCOLS
            LREQI = NBROW(I)
            LREQA = NBROW(I) * NBCOLS 
            IF ( LRLU .LT. LREQA .OR. IWPOS+LREQI-1 .GT. IWPOSCB ) THEN
              IF ( LRLUS .LT. LREQA ) THEN
                IFLAG = -9
                IERROR = LREQA - LRLUS
                GOTO 700
              END IF
              CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *            NIRBDU, LRLU, IPTRLU,
     *            IWPOS, IWPOSCB, PTRIST, PTRAST,
     *            STEP, PIMASTER, PAMASTER, ITLOC)
              COMP = COMP+1
              IF ( LRLU .NE. LRLUS ) THEN
                   WRITE(*,*) 'PB compress ass..maplig:LRLU,LRLUS=',
     * LRLU, LRLUS
                   IFLAG = -9
                   IERROR = LREQA - LRLUS
                   GOTO 700
              END IF
              IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
                IFLAG  = -8
                IERROR = IWPOS + LREQI - 1 - IWPOSCB
                GOTO 700
              END IF
            END IF
            LRLU   = LRLU - LREQA
            LRLUS  = LRLUS - LREQA
            VALSON = POSFAC
            POSFAC = POSFAC + LREQA
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,LREQA,KEEP)
            ROWLIST = IWPOS
            IWPOS   = IWPOS + LREQI
            NROWS_TO_STACK = 0
            DO II = 1, LMAP
              INDICE_PERE = MAP( II )
              CALL SMUMPS_47(
     *         KEEP, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *         NASS_PERE,
     *         NFRONT_PERE - NASS_PERE,
     *         NSLAVES_PERE,
     *         INDICE_PERE,
     *         NOSLA,
     *         IPOS_IN_SLAVE )
              INDICE_PERE = IPOS_IN_SLAVE
              IF ( NOSLA .EQ. I .OR. NBROW(I) .EQ. LMAP ) THEN
                NROWS_TO_STACK = NROWS_TO_STACK + 1
                IW( ROWLIST + NROWS_TO_STACK - 1) = INDICE_PERE
                POSROW=PTRAST(STEP(ISON)) + (II-1) * NFRONT + NPIV
                A( VALSON+(NROWS_TO_STACK-1) * NBCOLS:
     *               VALSON + NROWS_TO_STACK * NBCOLS - 1) =
     *               A( POSROW: POSROW + NBCOLS - 1 )
              END IF
            ENDDO
            IF ( NROWS_TO_STACK .NE. NBROW(I) ) THEN
              WRITE(*,*) 'Internal error 1 in psl_maplig'
              stop
            endif
            IF (PDEST.EQ.PDEST_MASTER) THEN 
              FREEMASTER = .TRUE.
#if defined(ROW_BY_ROW)
            write(*,*) MYID, ':row by row NBROW(I)=',NBROW(I)
            IF (NBROW(I)==0) THEN
              CALL SMUMPS_39(N, INODE_PERE, IW, LIW, 
     *          A, LA, ISON, NBROW(I), NBCOLS, IW(ROWLIST),
     *          A(VALSON), PTRIST, PTLUST_S, PTRAST,
     *          STEP, PIMASTER, PAMASTER, NIRBDU, OPASSW, OPELIW,
     *          POSFAC, LRLU, IPTRLU, IWPOSCB, LRLUS,
     *          NBPROCFILS, IPOOL, LPOOL, LEAF, FREEMASTER,
     *          PROCNODE_STEPS, SLAVEF, MYID, KEEP, FILS, ND,
     *          COMM_LOAD)
            ELSE
              FREEMASTER = .false.
              DO II = 1, NBROW(I)
               if (II==NBROW(i)) FREEMASTER=.TRUE.
               CALL SMUMPS_39(N, INODE_PERE, IW, LIW, 
     *          A, LA, ISON, 1, NBCOLS, IW(ROWLIST+II-1),
     *          A(VALSON+(II-1)*NBCOLS), PTRIST, PTLUST_S, PTRAST,
     *          STEP, PIMASTER, PAMASTER, NIRBDU, OPASSW, OPELIW,
     *          POSFAC, LRLU, IPTRLU, IWPOSCB, LRLUS,
     *          NBPROCFILS, IPOOL, LPOOL, LEAF, FREEMASTER,
     *          PROCNODE_STEPS, SLAVEF, MYID, KEEP, FILS, ND,
     *          COMM_LOAD)
              ENDDO
            ENDIF
#else
              CALL SMUMPS_39(N, INODE_PERE, IW, LIW, 
     *          A, LA, ISON, NBROW(I), NBCOLS, IW(ROWLIST),
     *          A(VALSON), PTRIST, PTLUST_S, PTRAST,
     *          STEP, PIMASTER, PAMASTER, NIRBDU, OPASSW, OPELIW,
     *          POSFAC, LRLU, IPTRLU, IWPOSCB, LRLUS,
     *          NBPROCFILS, IPOOL, LPOOL, LEAF, FREEMASTER,
     *          PROCNODE_STEPS, SLAVEF, MYID, KEEP, FILS, ND,
     *          COMM_LOAD)
#endif
            ELSE
               ISTCHK  = PTRIST(STEP(ISON))
               COLLIST = ISTCHK + 6 + IW( ISTCHK + 5 ) + NROW + NPIV
               IF ( KEEP(55) .eq. 0 ) THEN
                 CALL SMUMPS_40(N, INODE_PERE, IW, LIW,
     *           A, LA, NBROW(I), NBCOLS, IW( ROWLIST ),
     *           IW( COLLIST ), A( VALSON ),
     *           OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *           FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL, KEEP,
     *           MYID )
               ELSE
                 CALL SMUMPS_123(NELT, FRTPTR, FRTELT,
     *           N, INODE_PERE, IW, LIW,
     *           A, LA, NBROW(I), NBCOLS, IW( ROWLIST ),
     *           IW( COLLIST ), A( VALSON ),
     *           OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *           FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL, KEEP,
     *           MYID )
               END IF
            ENDIF
            LRLU = LRLU + LREQA
            LRLUS = LRLUS + LREQA
            POSFAC = POSFAC - LREQA
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,-LREQA,KEEP)
            IWPOS = IWPOS - LREQI
        END IF
      END DO
      DO I = NSLAVES_PERE, 0, -1
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
 95         DESCLU = .FALSE.
            IF ( PTRIST(STEP(ISON)) .lt.0 .or.
     *         IW ( PTRIST(STEP(ISON) ) ) .GT. N ) THEN
              WRITE(*,*) MYID,': Internal error in Maplig'
              WRITE(*,*) MYID,': PTRIST(STEP(ISON))=',
     *                           PTRIST(STEP(ISON))
              WRITE(*,*) MYID,': I, NBROW(I)=',I, NBROW(I)
              WRITE(*,*) MYID,': NSLAVES_PERE=',NSLAVES_PERE
              WRITE(*,*) MYID,': ISON, INODE_PERE=',ISON,INODE_PERE
              WRITE(*,*) MYID,': Son header=',
     *                   IW(PTRIST(STEP(ISON)): PTRIST(STEP(ISON))+5)
              CALL SMUMPS_ABORT()
            END IF
            CALL SMUMPS_67( DESCLU, INODE_PERE,
     *      NFRONT_PERE, NASS_PERE, NSLAVES_PERE, ISON,
     *      NBROW(I), LMAP, MAP, IW( PTRIST(STEP(ISON))),
     *      A(PTRAST(STEP(ISON))), I, PDEST, PDEST_MASTER, 
     *      COMM, IERR, 
     *
     *      KEEP, STEP, N, SLAVEF,
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE)
            IF ( IERR .EQ. -2 ) THEN
              IFLAG  = -17
              IERROR =  (NBROW(I) + 3 )* KEEP( 34 ) +
     *        NBROW(I) * IW(PTRIST(STEP(ISON))) * KEEP( 35 )
              GO TO 700
            END IF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL SMUMPS_329( COMM_LOAD,
     *          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *          MPI_ANY_SOURCE, MPI_ANY_TAG,
     *          STATUS,
     *          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *          IWPOS, IWPOSCB, IPTRLU,
     *          LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, 
     *          PTRIST, PTLUST_S,
     *          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *          IFLAG, IERROR, COMM,
     *          NBPROCFILS,
     *          IPOOL, LPOOL, LEAF,
     *          NBFIN, MYID, SLAVEF,
     *
     *          root, OPASSW, OPELIW, ITLOC, FILS, 
     *          PTRARW, PTRAIW,
     *          INTARR, DBLARR, ICNTL, KEEP, ND, FRERE, LPTRAR,
     *          NELT, FRTPTR, FRTELT, 
     *          ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 600
              GO TO 95
            END IF
        END IF
      END DO
      CALL SMUMPS_314( N, ISON,
     *    PTRIST, PTRAST, PTLUST_S, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP)
 600  CONTINUE
      DEALLOCATE(SLAVES_PERE)
      DEALLOCATE(NBROW)
      DEALLOCATE(MAP)
      RETURN
 700  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_210
      SUBROUTINE SMUMPS_211( COMM_LOAD, ASS_IRECV, 
     *  BUFR, LBUFR, LBUFR_BYTES,
     *
     *  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     *  NFRONT_PERE, NASS_PERE, LMAP, TROW,
     *  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *  LRLUS, N, IW,
     *  LIW, A, LA, NIRBDU,
     *  PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     *  ICNTL, KEEP, root,
     *  OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *  ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     *
     *  ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *  )
      USE SMUMPS_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER ICNTL( 40 ), KEEP(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER SLAVEF, NBFIN, POSFAC
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      REAL A( LA )
      INTEGER COMP
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER INODE_PERE, ISON
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE(NSLAVES_PERE)
      INTEGER NELIM, LMAP, TROW( LMAP )
      REAL OPASSW, OPELIW
      REAL DBLARR(MAX(1,KEEP(13)))
      INTEGER LPTRAR, NELT
      INTEGER IW( LIW )
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK( KEEP(28) ), ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER PTRIST(KEEP(28)), PTLUST_S( KEEP(28) ), PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER APOS
      INTEGER NOSLA, I, ISTCHK
      INTEGER INDICE_PERE, LREQI, LREQA
      INTEGER PDEST, PDEST_MASTER, NFRONT
      LOGICAL FLAG, SAME_PROC, FREE, DESCLU, FREEMASTER
      INTEGER SIZFI, SIZFR
      INTEGER NSLSON, VALSON, ROWLIST, NBCOLS, NROW,
     *        NPIV, NROWS_TO_STACK, II, HS,
     *        POSROW, NCOL, IPOS_IN_SLAVE
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP
      if (NSLAVES_PERE.le.0) then
       write(6,*) ' error 2 in maplig_fils_niv1 ', NSLAVES_PERE
       stop
      endif
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        write(*,*) myid, ' : PB allocation NBROW in SMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat =allocok)
      IF ( allocok .GT. 0 ) THEN
        write(*,*) myid, ' : PB allocation SLAVES_PERE in SMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = SMUMPS_275( STEP(INODE_PERE),
     *      PROCNODE_STEPS, SLAVEF )
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        write(*,*) myid, ' : PB allocation LMAP in SMUMPS_210'
        IFLAG  =-13
        IERROR = LMAP
        GOTO 700
      endif
      MAP( 1 : LMAP ) = TROW( 1 : LMAP )
      DO I = 0, NSLAVES_PERE
        NBROW( I ) = 0
      END DO
      DO I = 1, LMAP
        INDICE_PERE = MAP( I )
        CALL SMUMPS_47(
     *         KEEP, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *         NASS_PERE,
     *         NFRONT_PERE - NASS_PERE,
     *         NSLAVES_PERE,
     *         INDICE_PERE,
     *         NOSLA,
     *         IPOS_IN_SLAVE )
        NBROW( NOSLA ) = NBROW( NOSLA ) + 1
      END DO
      PDEST_MASTER = MYID
      IF ( SLAVES_PERE(0) .NE. MYID ) THEN
        WRITE(*,*) 'Error 1 in MAPLIG_FILS_NIV1:',MYID, SLAVES_PERE
        stop
      END IF
      PDEST        = PDEST_MASTER
        I = 0
        NBPROCFILS(STEP(INODE_PERE))=NBPROCFILS(STEP(INODE_PERE))-1
        NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK)
        NELIM      = IW(ISTCHK+1)
        NROW       = IW(ISTCHK+2)
        NPIV       = IW(ISTCHK+3)
        IF (NPIV.LT.0) THEN
         write(6,*) ' Error 2 in SMUMPS_211 ', NPIV
         CALL SMUMPS_ABORT()
        ENDIF
        NSLSON     = IW(ISTCHK+5)
        NFRONT     = NPIV + NBCOLS
        LREQI = NBROW(I)
        LREQA = NBROW(I) * NBCOLS 
        IF ( LRLU .LT. LREQA .OR. IWPOS+LREQI-1 .GT. IWPOSCB ) THEN
          IF ( LRLUS .LT. LREQA ) THEN
            IFLAG = -9
            IERROR = LREQA - LRLUS
            GO TO 700
          END IF
          CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
          COMP = COMP+1
          IF ( LRLU .NE. LRLUS ) THEN
               WRITE(*,*) 'PB compress ass..maplig_fils_niv1'
               IFLAG = -9
               IERROR = LREQA - LRLUS
               GO TO 700
          END IF
          IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
            IFLAG  = -8
            IERROR = IWPOS + LREQI - 1 - IWPOSCB
            GO TO 700
          END IF
        END IF
        LRLU   = LRLU - LREQA
        LRLUS  = LRLUS - LREQA
        VALSON = POSFAC
        POSFAC = POSFAC + LREQA
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,LREQA,KEEP)
        ROWLIST = IWPOS
        IWPOS   = IWPOS + LREQI
        NROWS_TO_STACK = 0
        DO II = 1, LMAP
          INDICE_PERE = MAP(II)
          CALL SMUMPS_47(
     *         KEEP, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *         NASS_PERE,
     *         NFRONT_PERE - NASS_PERE,
     *         NSLAVES_PERE,
     *         INDICE_PERE,
     *         NOSLA,
     *         IPOS_IN_SLAVE )
          INDICE_PERE = IPOS_IN_SLAVE
          IF ( NOSLA .EQ. I .OR. NBROW(I) .EQ. LMAP ) THEN
            NROWS_TO_STACK = NROWS_TO_STACK + 1
            IW( ROWLIST + NROWS_TO_STACK - 1) = INDICE_PERE
            POSROW = PAMASTER(STEP(ISON)) + (NELIM+II-1)*NBCOLS
            A( VALSON+(NROWS_TO_STACK-1) * NBCOLS:
     *           VALSON + NROWS_TO_STACK * NBCOLS - 1) =
     *       A( POSROW: POSROW + NBCOLS - 1 )
          END IF
        ENDDO
        IF ( NROWS_TO_STACK .NE. NBROW(I) ) THEN
          WRITE(*,*) 'error 1 in psl_maplig_fils_niv1'
          stop
        endif
        FREEMASTER = .FALSE.
        CALL SMUMPS_39(N, INODE_PERE, IW, LIW, 
     *    A, LA, ISON, NBROW(I), NBCOLS, IW(ROWLIST),
     *    A(VALSON), PTRIST, PTLUST_S, PTRAST,
     *    STEP, PIMASTER, PAMASTER, NIRBDU, OPASSW, OPELIW,
     *    POSFAC, LRLU, IPTRLU, IWPOSCB, LRLUS,
     *    NBPROCFILS, IPOOL, LPOOL, LEAF, FREEMASTER,
     *    PROCNODE_STEPS, SLAVEF, MYID, KEEP, FILS, ND,
     *    COMM_LOAD)
        LRLU = LRLU + LREQA
        LRLUS = LRLUS + LREQA
        POSFAC = POSFAC - LREQA
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,-LREQA,KEEP)
        IWPOS = IWPOS - LREQI
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
 95        NFRONT = IW(PIMASTER(STEP(ISON)))
           NELIM  = IW(PIMASTER(STEP(ISON))+1)
           APOS = PAMASTER(STEP(ISON)) + NFRONT * NELIM
           DESCLU = .TRUE.
           CALL SMUMPS_67( DESCLU, INODE_PERE,
     *       NFRONT_PERE, NASS_PERE, NSLAVES_PERE,
     *      ISON, NBROW(I), LMAP, MAP, IW(PIMASTER(STEP(ISON))),
     *      A(APOS), I, PDEST, PDEST_MASTER, COMM, IERR,
     *
     *      KEEP, STEP, N, SLAVEF,
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE)
            IF ( IERR .EQ. -2 ) THEN
              IFLAG  = -17
              IERROR =  (NBROW(I) + 3 )* KEEP( 34 ) +
     *        NBROW(I) *  KEEP( 35 )
              GO TO 700
            END IF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .FALSE.
              MESSAGE_RECEIVED = .FALSE.
              CALL SMUMPS_329( COMM_LOAD,
     *          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *          MPI_ANY_SOURCE, MPI_ANY_TAG,
     *          STATUS, 
     *          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *          IWPOS, IWPOSCB, IPTRLU,
     *          LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *          PTLUST_S,
     *          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *          IFLAG, IERROR, COMM,
     *          NBPROCFILS,
     *          IPOOL, LPOOL, LEAF,
     *          NBFIN, MYID, SLAVEF,
     *          root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *          INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *          LPTRAR, NELT, FRTPTR, FRTELT, 
     *          ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 600
              GO TO 95
            END IF
        END IF
      END DO
      ISTCHK = PTRIST(STEP(ISON))
      PTRIST(STEP( ISON )) = -77777777
            IF ( IW(ISTCHK) .GE. 0 ) THEN
              WRITE(*,*) 'error 3 in SMUMPS_211'
              CALL SMUMPS_ABORT()
            ENDIF
      CALL SMUMPS_152( MYID, N, ISTCHK,
     *     PAMASTER(STEP(ISON)),
     *     IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *     IWPOSCB, LA, KEEP
     *     )
 600  CONTINUE
      DEALLOCATE(NBROW)
      DEALLOCATE(MAP)
      DEALLOCATE(SLAVES_PERE)
      RETURN
 700  CONTINUE
      CALL SMUMPS_44(MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_211
      SUBROUTINE SMUMPS_209( N, FRERE, FILS, NFSIZ )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: N
      INTEGER, intent( in    )  :: NFSIZ( N )
      INTEGER, intent( inout )  :: FRERE( N ), FILS( N )
      INTEGER INODE, IROOT, IFILS, IN, IROOTLAST, SIZE
      IROOT = -9999
      SIZE  = 0
      DO INODE = 1, N
        IF ( FRERE( INODE ) .EQ. 0 )  THEN
          IF ( NFSIZ( INODE ) .GT. SIZE ) THEN
            SIZE  = NFSIZ( INODE )
            IROOT = INODE
          END IF
        ENDIF
      END DO
      IN = IROOT
      DO WHILE ( FILS( IN ) .GT. 0 )
        IN = FILS( IN )
      END DO
      IROOTLAST = IN
      IFILS     = - FILS ( IN )
      DO INODE = 1, N
        IF ( FRERE( INODE ) .eq. 0 .and. INODE .ne. IROOT ) THEN
          IF ( IFILS .eq. 0 ) THEN
            FILS( IROOTLAST ) = - INODE
            FRERE( INODE )    = -IROOT
            IFILS             = INODE
          ELSE
            FRERE( INODE ) = -FILS( IROOTLAST )
            FILS( IROOTLAST ) = - INODE
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_209
      SUBROUTINE SMUMPS_22( PROCESS_BANDE, MYID,N, KEEP, 
     *   IW, LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP,PIMASTER,PAMASTER, ITLOC,
     *   LREQ, LREQCB,
     *   COMP, LRLUS, IFLAG, IERROR )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N,LIW,LA,NIRBDU,LRLU, KEEP(150),
     &        IPTRLU,IWPOS,IWPOSCB
      INTEGER IW(LIW),PTRIST(KEEP(28)),PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     * PAMASTER(KEEP(28)), ITLOC(N)
      INTEGER MYID
      REAL A(LA)
      LOGICAL PROCESS_BANDE
      INTEGER COMP, LREQ, LREQCB, LRLUS, IFLAG, IERROR
      LOGICAL DONE
      DONE =.FALSE.
      IF ((IPTRLU.LT.LREQCB).OR.(LRLU.LT.LREQCB)) THEN
        IF (LRLUS.LT.LREQCB) THEN
          GOTO 620
        ELSE
          CALL SMUMPS_94(N,KEEP(28),IW,LIW,A,LA,NIRBDU,
     *                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     *                    PTRIST,PTRAST,
     *                    STEP, PIMASTER,PAMASTER,ITLOC)
          IF ( LRLU .NE. LRLUS ) THEN
            WRITE(*,*) 'PB compress... alloc_cb',
     *      'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 620
          END IF
          DONE = .TRUE.
          COMP = COMP + 1
        ENDIF
      ENDIF
      IF (IWPOSCB-IWPOS+1 .LT. LREQ) THEN
       IF (DONE) GOTO 600
                 CALL SMUMPS_94(N,KEEP(28),IW,LIW,A,LA,NIRBDU,
     *                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     *                    PTRIST,PTRAST,
     *                    STEP, PIMASTER,PAMASTER,ITLOC)
                 IF ( LRLU .NE. LRLUS ) THEN
                   WRITE(*,*) 'PB compress... alloc_cb',
     *             'LRLU,LRLUS=',LRLU,LRLUS
                   GOTO 620
                 END IF
          COMP = COMP + 1
          IF (IWPOSCB-IWPOS+1 .LT. LREQ) GOTO 600
      ENDIF
      IWPOSCB = IWPOSCB - LREQ
      IPTRLU = IPTRLU - LREQCB
      LRLU   = LRLU - LREQCB
      LRLUS  = LRLUS - LREQCB
#if ! defined(OLD_LOAD_MECHANISM)
      CALL SMUMPS_471(PROCESS_BANDE,LA-LRLUS,0,LREQCB,KEEP)
#else
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,LREQCB,KEEP)
#endif
      RETURN
 600  IFLAG  = -8
      IERROR = LREQ
      RETURN
 620  IFLAG  = -9
      IERROR = LREQCB - LRLUS
      RETURN
      END SUBROUTINE SMUMPS_22
      SUBROUTINE SMUMPS_93(MYID,N,IOLDPS,TYPE,IW, LIW,
     *           A, LA, POSFAC, LRLU, LRLUS, 
     *           IWPOS, PTRAST, STEP, KEEP)
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INTEGER MYID
      INTEGER IOLDPS, TYPE, LIW, LA, N, KEEP(150)
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER POSFAC, LRLU, LRLUS, IWPOS, LDLT
      INTEGER PTRAST(KEEP(28)), STEP( N )
      INTEGER LCONT, NELIM, NROW, NFRONT, NPIV,
     *        IAPOS, NSLAVES, INTSIZ
      INTEGER IPS, IPSIZE, SIZELU, SIZECB, I
      LOGICAL MOVEPTRAST
      LDLT = KEEP(50)
      IF ( IW( IOLDPS ) < 0 ) THEN
        write(*,*) ' ERROR 1 compressLU:Should not point to a band.'
        stop
      ELSE IF ( IW( IOLDPS + 2 ) < 0 ) THEN
        write(*,*) ' ERROR 2 compressLU:Stack not performed yet'
        stop
      ENDIF
      LCONT  = IW( IOLDPS )
      NELIM  = IW( IOLDPS + 1 )
      NROW   = IW( IOLDPS + 2 )
      NPIV   = IW( IOLDPS + 3 )
      IAPOS  = IW( IOLDPS + 4 )
      NSLAVES= IW( IOLDPS + 5 )
      NFRONT = LCONT + NPIV
      INTSIZ = 6 + NSLAVES + 2 * NFRONT
      IF ( (NSLAVES > 0  .AND. TYPE .NE. 2) .OR. 
     *   (NSLAVES .eq. 0 .AND. TYPE .EQ. 2 ) ) THEN
          WRITE(*,*) ' ERROR 3 compressLU: problem with level of inode'
          stop
      END IF
      IF (LDLT.EQ.0) THEN
        SIZELU = ( LCONT + NROW ) * NPIV
      ELSE
        SIZELU =   NROW * NPIV
      ENDIF
      IF ( TYPE .EQ. 2 ) THEN
        IF (LDLT.EQ.0) THEN
         SIZECB = NELIM * LCONT
        ELSE
         SIZECB = NELIM * ( NELIM + NPIV )
        ENDIF
      ELSE
        IF (LDLT.EQ.0) THEN
         SIZECB = LCONT * LCONT
        ELSE
         SIZECB = NROW * LCONT
        ENDIF
      END IF
      IF (SIZECB.EQ.0) GOTO 500
      IF ( IOLDPS + INTSIZ .NE. IWPOS ) THEN
         IPS = IOLDPS + INTSIZ
         MOVEPTRAST = .FALSE.
         DO WHILE ( IPS .NE. IWPOS )
           IF ( IW( IPS + 2 ) < 0 ) THEN
             NFRONT = IW( IPS )
             IPSIZE = 6 + 2 * NFRONT + IW( IPS + 5 )
             MOVEPTRAST = .TRUE.
           ELSE IF ( IW( IPS ) < 0 ) THEN
             IPSIZE = - IW( IPS ) + IW ( IPS + 1 ) + 4
             IW( IPS + 3 ) = IW( IPS + 3 ) - SIZECB
           ELSE
             NFRONT = IW( IPS ) + IW( IPS + 3 )
             IPSIZE = 6 + 2 * NFRONT + IW( IPS + 5 )
             IW( IPS + 4 ) = IW( IPS + 4 ) - SIZECB
           END IF
           IPS = IPS + IPSIZE
         END DO
         IF (SIZECB .NE. 0) THEN
           DO I=IAPOS + SIZELU, POSFAC - SIZECB - 1
             A( I ) = A( I + SIZECB )
           END DO
         END IF
         IF ( MOVEPTRAST ) THEN
           DO I = 1, KEEP(28)
             IF ( PTRAST( I ) < POSFAC .AND.
     *            PTRAST( I ) > IAPOS ) THEN
               PTRAST( I ) = PTRAST( I ) - SIZECB
             END IF
           END DO
         END IF
      END IF
      POSFAC = POSFAC - SIZECB
      LRLU   = LRLU   + SIZECB
      LRLUS  = LRLUS  + SIZECB
 500  CONTINUE
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,SIZELU,-SIZECB,KEEP)
      RETURN
      END SUBROUTINE SMUMPS_93
      SUBROUTINE SMUMPS_314( N, ISON, 
     *    PTRIST, PTRAST, PTLUST_S, IW, LIW, A, LA, 
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP)
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N, ISON, LIW, LA, LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC,
     *        COMP, NIRBDU, IPTRLU, IFLAG, IERROR, SLAVEF, MYID, COMM
      INTEGER KEEP(150)
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), STEP(N), 
     * PIMASTER(KEEP(28)), PAMASTER(KEEP(28)), IW(LIW)
      INTEGER PTLUST_S(KEEP(28)), ITLOC(N)
      REAL OPELIW
      REAL FLOP1
      REAL A( LA )
      INTEGER  NFRONT, NCOL_L, NROW_L, LREQI, LREQA,NSLAVES_L,
     *         POSA, POSI, POSALOC, ISTCHK, SIZFR, SIZFI, I,
     *         NSLSON, HS, IROW_L, ICOL_L, OLDPOS, NROW, NCOL,
     *         LDA_BAND, NASS
      LOGICAL SAME_PROC, FREE
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      FLOP1 = ZERO
      NCOL_L = IW( PTRIST(STEP( ISON )) + 3 )
      NROW_L = IW( PTRIST(STEP( ISON )) + 2 )
      NSLAVES_L = IW( PTRIST(STEP( ISON )) + 5 )
      LDA_BAND = NCOL_L + IW( PTRIST(STEP( ISON )) )
      IF  ( KEEP(50) .eq. 0 ) THEN
        NFRONT = LDA_BAND
      ELSE
        NFRONT = IW( PTRIST(STEP( ISON )) + 7 )
      END IF
      IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) GOTO 80
      LREQI   = 4 + NCOL_L + NROW_L
      LREQA   = NCOL_L * NROW_L
      IF ( LRLU .LT. LREQA .OR.
     *  IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LREQA ) THEN
          IFLAG  = -9
          IERROR = lREQA - LRLUS
          GO TO 700
        END IF
        CALL SMUMPS_94( N,KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB compress Stack_band:LRLU,LRLUS=',
     *            LRLU, LRLUS
                  IFLAG = -9
                  IERROR = LREQA - LRLU
                  GOTO 700
        END IF
        IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
          IFLAG  = -8
          IERROR = IWPOS + LREQI - 1 - IWPOSCB
          GOTO 700
        END IF
      END IF
      POSA = POSFAC
      POSFAC = POSFAC + LREQA
      LRLU = LRLU - LREQA
      LRLUS = LRLUS - LREQA
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,LREQA,LREQA,KEEP)
      POSI = IWPOS
      IWPOS = IWPOS + LREQI
      PTLUST_S(STEP( ISON )) = POSI
      IW( POSI     ) = - NCOL_L
      IW( POSI + 1 ) =   NROW_L
      IW( POSI + 2 ) =   NFRONT - NCOL_L
      IW( POSI + 3 ) =   POSA
      IROW_L = PTRIST(STEP( ISON )) + 6 + NSLAVES_L
      ICOL_L = PTRIST(STEP( ISON )) + 6 + NROW_L + NSLAVES_L
      IW( POSI + 4 : POSI + 3 + NROW_L ) =
     *                       IW( IROW_L: IROW_L + NROW_L - 1 )
      IW( POSI + NROW_L + 4: POSI + NROW_L + NCOL_L + 3 ) =
     *                       IW( ICOL_L: ICOL_L + NCOL_L - 1 )
      POSALOC = POSA
      DO I = 1, NROW_L
        OLDPOS =  PTRAST( STEP(ISON)) + ( I - 1 ) * LDA_BAND
        A( POSALOC : POSALOC + NCOL_L - 1 ) = 
     *                   A( OLDPOS: OLDPOS + NCOL_L - 1 )
        POSALOC = POSALOC + NCOL_L
      END DO
      IF ( KEEP(50) .eq. 0 ) THEN
         FLOP1 = REAL( NCOL_L * NROW_L) +
     *     REAL(NROW_L*NCOL_L)*REAL(2*NFRONT-NCOL_L-1)
      ELSE
         FLOP1 = REAL( NCOL_L ) * REAL( NROW_L )
     *            * ( 2 * LDA_BAND - NROW_L - NCOL_L + 1)
      END IF
      OPELIW = OPELIW + FLOP1
      NASS = IW( PTRIST(STEP( ISON )) + 1 )
      IF ( NCOL_L .NE. NASS ) THEN
        IF ( KEEP(50).eq.0 ) THEN
           FLOP1 = REAL( NASS * NROW_L) +
     *     REAL(NROW_L*NASS)*REAL(2*NFRONT-NASS-1)
        ELSE
           FLOP1 = REAL( NASS ) * REAL( NROW_L ) *
     *     REAL( 2 * LDA_BAND - NROW_L - NASS + 1)
        END IF
      END IF
      CALL SMUMPS_190(.FALSE.,-FLOP1,KEEP)
  80  CONTINUE
      ISTCHK = PTRIST(STEP(ISON))
      PTRIST(STEP(ISON)) = -8888888
      CALL SMUMPS_152( MYID, N, ISTCHK,
     *     PTRAST(STEP(ISON)),
     *     IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *     IWPOSCB, LA, KEEP
     *     )
      PTRIST(STEP( ISON )) = -9999888
      PTRAST(STEP( ISON )) = -9999888
      RETURN
 700  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_314
      INTEGER FUNCTION SMUMPS_330(INODE, PROCNODE, SLAVEF)
      IMPLICIT NONE
      INTEGER INODE, SLAVEF 
      INTEGER PROCNODE(*)
      INTEGER MYPROC, TPN
      MYPROC = PROCNODE( INODE )
      TPN = (MYPROC-1+2*SLAVEF)/SLAVEF - 1
      IF ( TPN .LT. 1 ) TPN = 1
      SMUMPS_330 = TPN
      RETURN 
      END
      INTEGER FUNCTION SMUMPS_275(INODE, PROCNODE, SLAVEF)
      IMPLICIT NONE
      INTEGER INODE, SLAVEF 
      INTEGER PROCNODE(*)
      INTEGER MYPROC
      MYPROC = PROCNODE( INODE )
      SMUMPS_275=MOD(2*SLAVEF+MYPROC-1,SLAVEF)
      RETURN
      END
      LOGICAL FUNCTION SMUMPS_283( INODE, PROCNODE, SLAVEF )
      IMPLICIT NONE
      INTEGER INODE, SLAVEF
      INTEGER  PROCNODE(*)
      INTEGER TPN, MYPROC
      MYPROC = PROCNODE( INODE )
      TPN = (MYPROC-1+2*SLAVEF)/SLAVEF - 1
      SMUMPS_283 = ( TPN .eq. 0 )
      RETURN
      END
      LOGICAL FUNCTION SMUMPS_167( INODE, PROCNODE, SLAVEF )
      IMPLICIT NONE
      INTEGER INODE, SLAVEF
      INTEGER  PROCNODE(*)
      INTEGER TPN, MYPROC
      MYPROC = PROCNODE( INODE )
      TPN = (MYPROC-1+SLAVEF+SLAVEF)/SLAVEF - 1
      SMUMPS_167 = ( TPN .eq. -1 )
      RETURN 
      END
      LOGICAL FUNCTION SMUMPS_170
     *        ( INODE, PROCNODE, SLAVEF )
      IMPLICIT NONE
      INTEGER INODE, SLAVEF
      INTEGER  PROCNODE(*)
      INTEGER TPN, MYPROC
      MYPROC = PROCNODE( INODE )
      TPN = (MYPROC-1+SLAVEF+SLAVEF)/SLAVEF - 1
      SMUMPS_170 =
     *           ( TPN .eq. -1 .OR. TPN .eq. 0 )
      RETURN
      END
      LOGICAL FUNCTION SMUMPS_358( MYID, SLAVEF, INODE,
     *                 NMB_PAR2, ISTEP_TO_INIV2 , K71, STEP, N, 
     *                 CANDIDATES, KEEP24 )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, INODE, NMB_PAR2, KEEP24, I
      INTEGER K71, N
      INTEGER ISTEP_TO_INIV2 ( K71 ), STEP ( N )
      INTEGER CANDIDATES(SLAVEF+1, MAX(NMB_PAR2,1))
      INTEGER NCAND, POSINODE
      SMUMPS_358 = .FALSE.
      IF (KEEP24 .eq. 0) RETURN
      POSINODE = ISTEP_TO_INIV2 ( STEP (INODE) )
      NCAND = CANDIDATES( SLAVEF+1, POSINODE )
      DO I = 1, NCAND
        IF (MYID .EQ. CANDIDATES( I, POSINODE ))
     *     SMUMPS_358 = .TRUE.
      END DO
      RETURN
      END
      SUBROUTINE SMUMPS_291(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()
      RETURN
      END SUBROUTINE SMUMPS_291
      SUBROUTINE SMUMPS_292(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()-T
      RETURN
      END SUBROUTINE SMUMPS_292
      SUBROUTINE SMUMPS_312( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      REAL VAL( N )
      INTEGER I, ISWAP
      REAL SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_312
      SUBROUTINE SMUMPS_158
     *( N, MYID, STEP, PROCNODE_STEPS, K28, SLAVEF, NA, LNA,  MYLEAF )
      IMPLICIT NONE
      INTEGER, intent(in) :: N, SLAVEF, MYID, K28, LNA
      INTEGER, intent(in) :: PROCNODE_STEPS( K28 ), STEP(N), NA( LNA )
      INTEGER, intent(out) :: MYLEAF
      INTEGER INODE, NBLEAF, I
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      MYLEAF = 0
      NBLEAF = NA(1)
      DO I = 1, NBLEAF
        INODE = NA(2+I)
        IF ( SMUMPS_275(STEP(INODE),PROCNODE_STEPS,SLAVEF)
     *       .eq. MYID ) THEN
          MYLEAF = MYLEAF + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_158
#if defined (PESSL)
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     $                     LLD, INFO )
      INTEGER            ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
#if defined(DESC8)
      PARAMETER          ( DLEN_ = 8, DTYPE_ = 1,
     $                     CTXT_ = 7, M_ = 1, N_ = 2, MB_ = 3, NB_ = 4,
     $                     RSRC_ = 5, CSRC_ = 6, LLD_ = 8 )
#else
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
#endif
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           BLACS_GRIDINFO, PXERBLA
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      INTRINSIC          MAX, MIN
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( MB.LT.1 ) THEN
         INFO = -4
      ELSE IF( NB.LT.1 ) THEN
         INFO = -5
      ELSE IF( IRSRC.LT.0 .OR. IRSRC.GE.NPROW ) THEN
         INFO = -6
      ELSE IF( ICSRC.LT.0 .OR. ICSRC.GE.NPCOL ) THEN
         INFO = -7
      ELSE IF( NPROW.EQ.-1 ) THEN
         INFO = -8
      ELSE IF( LLD.LT.MAX( 1, NUMROC( M, MB, MYROW, IRSRC,
     $                                NPROW ) ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 )
     $   CALL PXERBLA( ICTXT, 'DESCINIT', -INFO )
#ifndef DESC8
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D
#endif
      DESC( M_ )  = MAX( 0, M )
      DESC( N_ )  = MAX( 0, N )
      DESC( MB_ ) = MAX( 1, MB )
      DESC( NB_ ) = MAX( 1, NB )
      DESC( RSRC_ ) = MAX( 0, MIN( IRSRC, NPROW-1 ) )
      DESC( CSRC_ ) = MAX( 0, MIN( ICSRC, NPCOL-1 ) )
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ )  = MAX( LLD, MAX( 1, NUMROC( DESC( M_ ), DESC( MB_ ),
     $                              MYROW, DESC( RSRC_ ), NPROW ) ) )
      RETURN
      END SUBROUTINE DESCINIT
      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
      INTEGER            ICTXT, INFO
      CHARACTER*(*)      SRNAME
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           BLACS_GRIDINFO
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     $        ' parameter number', I4, ' had an illegal value' )
      END SUBROUTINE PXERBLA
#endif
      SUBROUTINE SMUMPS_256( N, NZ, IRN, ICN, ASPK, X, Y,
     *                         LDLT, MTYPE, MAXTRANS, PERM )
      INTEGER N, NZ, LDLT, MTYPE, MAXTRANS
      INTEGER IRN( NZ ), ICN( NZ ) 
      INTEGER PERM( N )
      REAL ASPK( NZ ), X( N ), Y( N )
      INTEGER K, I, J
      REAL PX( N )
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      Y = REAL(ZERO)
      IF ( MAXTRANS .eq. 1 .and. MTYPE .eq. 1) THEN
        DO I = 1, N
          PX(I) = X( PERM( I ) )
        END DO
      ELSE
        PX = X
      END IF
      IF ( LDLT .eq. 0 ) THEN
       IF (MTYPE .EQ. 1) THEN
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
        ENDDO
       ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y(J) = Y(J) + ASPK(K) * PX(I)
        ENDDO
       ENDIF
      ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
          IF (J.NE.I) THEN
            Y(J) = Y(J) + ASPK(K) * PX(I)
          ENDIF
        ENDDO
      END IF
      IF ( MAXTRANS .EQ. 1 .AND. MTYPE .eq. 0 ) THEN
      PX = Y
      DO I = 1, N
        Y( PERM( I ) ) = PX( I )
      END DO
      END IF
      RETURN
      END SUBROUTINE SMUMPS_256
      SUBROUTINE SMUMPS_214( KEEP, MYID, N, NELT, LNA, NSLAVES,
     *           MEMORY_BYTES, EFF )
      IMPLICIT NONE
      INCLUDE 'smumps_struc.h'
      INTEGER MEMORY_BYTES
      INTEGER KEEP(150)
      INTEGER MYID, N, NELT, NSLAVES, LNA
      LOGICAL EFF
      LOGICAL I_AM_SLAVE, I_AM_MASTER
      INTEGER NB_INT, NB_REAL, PERLU, NBRECORDS
      REAL DMEMORY_BYTES, TEMP, NB_BYTES
      INTEGER  SMUMPS_LBUF_INT, SMUMPS_LBUFR_BYTES, SMUMPS_LBUF
      INTEGER TEMPI, TEMPR, NBUFS
      INTEGER PERLU_SAVE
      PERLU     = KEEP(12)
      PERLU_SAVE= PERLU
      NBRECORDS = KEEP(39)
      I_AM_MASTER = ( MYID .eq. 0 )
      I_AM_SLAVE  = ( KEEP(46).eq. 1 .or. MYID .ne. 0 )
      TEMP    = 0
      NB_REAL = 0
      NB_BYTES = 0.0
      NB_INT  = 0
      NB_INT = NB_INT + 4*KEEP(28)
      NB_INT = NB_INT +KEEP(28)+KEEP(56)*(NSLAVES+2)
      NB_INT = NB_INT + 3*N
      IF (KEEP(23).ne.0 .and. I_AM_MASTER) NB_INT=NB_INT + N
      IF (KEEP(55).eq.0) THEN
        NB_INT = NB_INT + 2 * N
      ELSE
        NB_INT = NB_INT + 2 * ( NELT+1 )
      ENDIF
      IF (KEEP(55) .ne. 0 ) THEN
        NB_INT = NB_INT + N+1 + NELT
      END IF
      NB_INT = NB_INT + LNA
      IF ( .NOT. EFF ) THEN
        IF (KEEP(50).EQ.1) THEN
         NB_REAL = NB_REAL + KEEP(16)+PERLU*(KEEP(16) / 100 + 1 )
     *                  + KEEP(44) + PERLU * ( KEEP(44) / 100 )
        ELSE
         NB_REAL = NB_REAL + KEEP(16)+PERLU*(KEEP(16) / 100 + 1 )
     *                  + KEEP(44) + 4*PERLU * ( KEEP(44) / 100 )
        ENDIF
      ELSE
        NB_REAL = NB_REAL + KEEP(67) + KEEP(44)
      ENDIF
      NB_REAL = NB_REAL + KEEP(13)
      IF ( .not. ( I_AM_SLAVE .and. I_AM_MASTER .and. KEEP(52) .eq. 0
     *         .and. KEEP(55) .ne. 0 ) ) THEN
        NB_INT  = NB_INT  + KEEP(14)
      END IF
      IF ( I_AM_SLAVE .and. KEEP(38) .ne. 0 ) THEN
        NB_INT = NB_INT + 2 * N
      END IF
      TEMPI= 0
      TEMPR = 0
      IF ( KEEP(54) .eq. 0 ) THEN
        IF ( I_AM_MASTER ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            NBUFS = NSLAVES 
          ELSE
            NBUFS = NSLAVES - 1
            IF (KEEP(55) .eq. 0 )
     *      TEMPI = TEMPI + 2 * N
          END IF
          TEMPI = TEMPI + 2 * NBRECORDS * NBUFS
          TEMPR = TEMPR + NBRECORDS * NBUFS
        ELSE
          IF ( KEEP(55) .eq. 0 )THEN
            TEMPI = TEMPI + 2 * NBRECORDS
            TEMPR = TEMPR + NBRECORDS
          END IF
        END IF
      ELSE
        IF ( I_AM_SLAVE ) THEN
          TEMPI = TEMPI + ( 1 + 4 * NSLAVES ) * NBRECORDS
          TEMPR = TEMPR + ( 1 + 2 * NSLAVES ) * NBRECORDS
        END IF
      END IF
      TEMP = MAX( NB_BYTES + REAL(NB_INT +TEMPI) * KEEP(34)
     *           + REAL(NB_REAL+TEMPR) * KEEP(35) 
     *            , TEMP )
      IF ( I_AM_SLAVE ) THEN
        SMUMPS_LBUFR_BYTES = KEEP( 44 ) * KEEP( 35 )
        SMUMPS_LBUFR_BYTES = MAX( SMUMPS_LBUFR_BYTES,
     *                      100000 )
        IF ( KEEP( 50 ) .EQ. 1 .OR. EFF ) PERLU = 0
        SMUMPS_LBUFR_BYTES = SMUMPS_LBUFR_BYTES
     *                 + INT( 4.0E0 * REAL(MAX(PERLU,0))*
     *                   REAL(SMUMPS_LBUFR_BYTES)/100E0)
        NB_BYTES = NB_BYTES + SMUMPS_LBUFR_BYTES
        IF (NSLAVES.EQ.2) THEN
         SMUMPS_LBUF = INT( 1.01E0 * KEEP( 43 ) *KEEP( 35 )  )
        ELSE
         SMUMPS_LBUF = INT( 2.01E0 * KEEP( 43 ) *KEEP( 35 )  )
        ENDIF
        SMUMPS_LBUF = MAX( SMUMPS_LBUF, 100000 )
        SMUMPS_LBUF = SMUMPS_LBUF
     *                 + INT( 2.0E0 * REAL(MAX(PERLU,0))*
     *                   REAL(SMUMPS_LBUF)/100E0)
        NB_BYTES = NB_BYTES + SMUMPS_LBUF
        IF ( KEEP( 50 ) .EQ. 1 .OR. EFF) PERLU = PERLU_SAVE
        SMUMPS_LBUF_INT = 
     *         NSLAVES * NSLAVES * 4
     *               * KEEP(34)
        NB_BYTES = NB_BYTES + SMUMPS_LBUF_INT
        IF ( EFF ) THEN
          NB_INT = NB_INT + KEEP(15)
        ELSE
          NB_INT = NB_INT + 
     *           KEEP(15) + 2 * MAX(PERLU,10) *
     *           ( KEEP(15) / 100 + 1 )
        ENDIF
        NB_INT = NB_INT + KEEP(28)
        NB_INT = NB_INT + N + 7 * KEEP(28) + 1
      END IF
      DMEMORY_BYTES = NB_BYTES + REAL(NB_INT) * KEEP(34) +
     *               REAL(NB_REAL) * KEEP(35)
      DMEMORY_BYTES = MAX( DMEMORY_BYTES, TEMP )
      MEMORY_BYTES = INT( DMEMORY_BYTES / 1 000 000 )
      RETURN
      END SUBROUTINE SMUMPS_214
      SUBROUTINE SMUMPS_243(MYID, COMM, INFO, INFOG, IRANK)
      IMPLICIT NONE
      INTEGER MYID, COMM, IRANK, INFO(1), INFOG(2)
      INCLUDE 'mpif.h'
      INTEGER IERR_MPI, MASTER
      INTEGER TEMP1(2), TEMP2(2)
      PARAMETER( MASTER = 0 )
      CALL MPI_REDUCE( INFO(1), INFOG(1), 1, MPI_INTEGER,
     *                 MPI_MAX, MASTER, COMM, IERR_MPI )
      CALL MPI_REDUCE( INFO(1), INFOG(2), 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, COMM, IERR_MPI )
      TEMP1(1) = INFO(1)
      TEMP1(2) = MYID
      CALL MPI_REDUCE( TEMP1, TEMP2, 1, MPI_2INTEGER,
     *                 MPI_MAXLOC, MASTER, COMM, IERR_MPI )
      IF ( MYID.eq. MASTER ) THEN
        IF ( INFOG(1) .ne. TEMP2(1) ) THEN
          write(*,*) 'Error in SMUMPS_243'
          stop
        END IF
        IRANK    = TEMP2(2)
      ELSE
        IRANK    = -1
      END IF
      RETURN
      END SUBROUTINE SMUMPS_243
      SUBROUTINE SMUMPS_362(N, LEAF, NBROOT, NROOT_LOC,
     *           MYID_NODES,
     *           SLAVEF, NA, LNA, KEEP, STEP,
     *           PROCNODE_STEPS, IPOOL, LPOOL)
      IMPLICIT NONE
      INTEGER N, LEAF, NROOT_LOC, NBROOT, MYID_NODES,
     *        SLAVEF, LPOOL, LNA
      INTEGER KEEP(150)
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28)), NA(LNA),
     *        IPOOL(LPOOL)
      INTEGER NBLEAF, INODE, I
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      NBLEAF = NA(1)
      NBROOT = NA(2)
      LEAF = 1
      DO I = 1, NBLEAF
        INODE = NA(I+2)
        IF (SMUMPS_275(STEP(INODE),PROCNODE_STEPS,SLAVEF)
     *   .EQ.MYID_NODES) THEN
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
          ENDIF
      ENDDO
      NROOT_LOC = 0
      DO I = 1, NBROOT
        INODE = NA(I+2+NBLEAF)
        IF (SMUMPS_275(STEP(INODE),PROCNODE_STEPS,
     *    SLAVEF).EQ.MYID_NODES) THEN
            NROOT_LOC = NROOT_LOC + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_362
      LOGICAL FUNCTION SMUMPS_438(TAB1,TAB2,LEN1,LEN2)
      IMPLICIT NONE
      INTEGER LEN1 , LEN2 ,I
      INTEGER TAB1(LEN1)
      INTEGER TAB2(LEN2)
      SMUMPS_438=.FALSE.
      IF(LEN1 .NE. LEN2) THEN
         RETURN
      ENDIF
      DO I=1 , LEN1
         IF(TAB1(I) .NE. TAB2(I)) THEN
            RETURN
         ENDIF
      ENDDO
      SMUMPS_438=.TRUE.
      RETURN
      END FUNCTION SMUMPS_438
      SUBROUTINE SMUMPS_463( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_463
      SUBROUTINE SMUMPS_466( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .LT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_466
      SUBROUTINE SMUMPS_ABORT()
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, IERRCODE
      IERRCODE = -99
      CALL MPI_ABORT(MPI_COMM_WORLD, IERRCODE, IERR)
      END SUBROUTINE SMUMPS_ABORT
      SUBROUTINE SMUMPS_244(N, NZ, NSTEPS,
     * A, LA, IW, LIW, SYM_PERM, NA, LNA,
     * NE_STEPS, NFSIZ, FILS,
     * STEP, FRERE, CAND, 
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     * PTRAR, LDPTRAR,
     * PTRIST, PTLUST_S, IW1, IW2, POOL, LPOOL,  
     * CNTL, ICNTL, INFO, RINFO, KEEP,PROCNODE_STEPS,
     * SLAVEF,
     * COMM_NODES, MYID, MYID_NODES,
     * BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,
     * root, NELT, FRTPTR, FRTELT, COMM_LOAD,
     * ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     * MEM_DISTRIB)
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER N,NZ,NSTEPS,LA,LIW,LPOOL,SLAVEF,COMM_NODES
      INTEGER MYID, MYID_NODES,LNA
      REAL A(LA)
      REAL RINFO(20)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER BUFR( LBUFR )
      INTEGER NELT, LDPTRAR
      INTEGER FRTPTR(*), FRTELT(*)
      REAL CNTL(5)
      INTEGER   ICNTL(40)
      INTEGER   INFO(40), KEEP(150)
      INTEGER   IW(LIW), SYM_PERM(N), NA(LNA),
     *          NE_STEPS(KEEP(28)), FILS(N),
     *          FRERE(KEEP(28)), NFSIZ(KEEP(28))
      INTEGER   CAND(SLAVEF+1, MAX(1,KEEP(56)))
      INTEGER   STEP(N)
      INTEGER   PTRAR(LDPTRAR,2)
      INTEGER   PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER   IW1(5*KEEP(28)), IW2(N), POOL(LPOOL)
      INTEGER   PROCNODE_STEPS(KEEP(28))
      INTEGER   COMM_LOAD, ASS_IRECV
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER   INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
      REAL SEUIL, SEUIL_LDLT_NIV2
      REAL UULOC 
      INTEGER LP, MP
      INTEGER NSTK,PTRAST, NBPROCFILS
      INTEGER PIMASTER, PAMASTER
      INTEGER NIRBDU
      LOGICAL PROK
      REAL ZERO, ONE
      DATA ZERO /0.0E0/
      DATA ONE /1.0E0/
      KEEP(41)=0
      KEEP(42)=0
      INFO(11)= 0
      INFO(12)= 0
      INFO(13)= 0
      INFO(14)  = 0
      RINFO(2) = ZERO
      RINFO(3) = ZERO
      NSTEPS   = 0
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      PROK   = (MP.GT.0)
      UULOC = CNTL(1)
      IF (UULOC.GT.ONE)   UULOC=ONE
      IF (UULOC.LT.ZERO)  UULOC=ZERO
      PTRAST     = 1
      PIMASTER   = PTRAST + KEEP(28)
      PAMASTER   = PIMASTER + KEEP(28) 
      NSTK       = PAMASTER + KEEP(28)
      NBPROCFILS = NSTK + KEEP(28)
      NIRBDU = LIW
      CALL SMUMPS_250(N, NZ, NSTEPS, A, LA, IW,
     * LIW, SYM_PERM, NA, LNA, NE_STEPS, IW1(NSTK),
     * IW1( NBPROCFILS ),
     * INFO(11),NFSIZ,FILS,STEP, FRERE, CAND, 
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     * PTRAR(1,2),PTRAR(1,1),IW1(PTRAST),PTRIST,PTLUST_S,
     * IW1(PIMASTER), IW1(PAMASTER),
     * IW2,NIRBDU,
     * POOL(1), LPOOL ,
     * UULOC,ICNTL, INFO, RINFO, KEEP,PROCNODE_STEPS,
     * SLAVEF,
     * COMM_NODES,MYID,MYID_NODES,
     * BUFR,LBUFR,LBUFR_BYTES,
     * INTARR,DBLARR,root,NELT, FRTPTR, 
     * FRTELT, LDPTRAR, COMM_LOAD, ASS_IRECV, SEUIL,
     * SEUIL_LDLT_NIV2, MEM_DISTRIB)
      IF (PROK) THEN
        WRITE (MP,99980) INFO(1), INFO(2),
     *       KEEP(28), INFO(9), INFO(10), INFO(11), INFO(12),
     *       INFO(13), INFO(14), RINFO(2), RINFO(3)
      ENDIF
      RETURN
99980 FORMAT (/' LEAVING FACTORIZATION PHASE WITH ...'/
     1      ' INFO (1)                                      =',I12/
     2      '  --- (2)                                      =',I12/
     5      '           NUMBER OF NODES IN THE TREE         =',I12/
     6      ' INFO (9)  REAL SPACE FOR FACTORS              =',I12/
     7      '  --- (10) INTEGER SPACE FOR FACTORS           =',I12/
     8      '  --- (11) MAXIMUM SIZE OF FRONTAL MATRICES    =',I12/
     9      '  --- (12) NUMBER OF OFF DIAGONAL PIVOTS       =',I12/
     1      '  --- (13) NUMBER OF DELAYED PIVOTS            =',I12/
     2      '  --- (14) NUMBER OF MEMORY COMPRESSES         =',I12/
     3  ' RINFO(2)  OPERATIONS DURING NODE ASSEMBLY     =',1PD10.3/
     4  ' -----(3)  OPERATIONS DURING NODE ELIMINATION  =',1PD10.3)
      END SUBROUTINE SMUMPS_244
      SUBROUTINE SMUMPS_269( MYID,KEEP, 
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS,
     *           N, IW, LIW, A, LA, NIRBDU,
     *           PTRIST, PTRAST, STEP, PIMASTER, PAMASTER,
     *           NSTK_S, COMP,
     *           FPERE, FLAG, IFLAG, IERROR, COMM, ITLOC )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER MYID
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER KEEP(150), BUFR( LBUFR )
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST( KEEP(28) ), PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP, FPERE
      LOGICAL FLAG
      INTEGER NSTK_S( KEEP(28) ), ITLOC( N )
      INTEGER IFLAG, IERROR, COMM
      INTEGER POSITION, FINODE, FLCONT, LREQ, LREQCB
      FLAG = .FALSE.
      POSITION = 0
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *                FINODE, 1, MPI_INTEGER,
     *                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *                FPERE, 1, MPI_INTEGER, 
     *                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *                FLCONT, 1, MPI_INTEGER,
     *                COMM, IERR)
      LREQ    = 2 * FLCONT + 6
      LREQCB  = FLCONT * FLCONT
      IF (IPTRLU.LT.0) WRITE(*,*) 'before alloc_cb:IPTRLU = ',IPTRLU
      CALL SMUMPS_22( .FALSE.,
     *   MYID,N, KEEP, IW, LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER,PAMASTER, ITLOC,
     *   LREQ, LREQCB,
     *   COMP, LRLUS, IFLAG, IERROR )
      IF (IPTRLU.LT.0) WRITE(*,*) 'after alloc_cb:IPTRLU = ',IPTRLU
      IF ( IFLAG .LT. 0 ) RETURN
      PIMASTER(STEP( FINODE )) = IWPOSCB + 1
      PAMASTER(STEP( FINODE )) = IPTRLU  + 1
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        IW(IWPOSCB + 1), LREQ,
     *        MPI_INTEGER, COMM, IERR)
      IF ( LREQCB .ne. 0 ) THEN
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION, A(IPTRLU + 1),
     *        LREQCB, MPI_REAL, COMM, IERR)
      END IF
      NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
      IF ( NSTK_S(STEP(FPERE)).EQ.0 ) THEN
          FLAG = . TRUE.
      END IF
      RETURN
      END SUBROUTINE SMUMPS_269
      RECURSIVE SUBROUTINE SMUMPS_322(
     *    COMM_LOAD, ASS_IRECV,
     *    MSGSOU, MSGTAG,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER MSGSOU, MSGTAG
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER KEEP(150), ICNTL( 40 )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28))
      INTEGER STEP(N),
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER MYID, SLAVEF, NBFIN
      REAL OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1), FRTELT( NELT )
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      REAL DBLARR( MAX(1,KEEP(13)) )
      INTEGER INIV2, ISHIFT, IBEG
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      LOGICAL FLAG
      INTEGER MP, LP
      INTEGER TMP( 2 )
      INTEGER NBRECU, POSITION, INODE, ISON, IROOT
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE,
     *LASTBL_PERE, LMAP, FPERE, NELIM,
     *TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      REAL FLOP1
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      MP = ICNTL(2)
      LP = ICNTL(1)
      IF ( MSGTAG .EQ. RACINE ) THEN
          NBRECU = BUFR( 1 )
          NBFIN =  NBFIN - NBRECU
      ELSEIF ( MSGTAG .EQ. NOEUD ) THEN
          CALL SMUMPS_269( MYID,KEEP,
     *    BUFR, LBUFR, LBUFR_BYTES,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTRAST,
     *    STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, FPERE, FLAG, IFLAG, IERROR, COMM, ITLOC )
          IF ( IFLAG .LT. 0 ) GO TO 500
          IF ( FLAG ) THEN
            CALL SMUMPS_507(N, IPOOL, LPOOL,
     *           PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     *           STEP, FPERE )
            IF (KEEP(47) .GE. 3) THEN
               CALL SMUMPS_500(
     $              IPOOL, LPOOL,
     *              PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *              MYID, STEP, N, ND, FILS )
            ENDIF
            CALL SMUMPS_137( FPERE, N,
     *           PROCNODE_STEPS,SLAVEF,
     &           ND, FILS, FRERE, STEP, PIMASTER,
     *           KEEP(28), KEEP(50), FLOP1,
     &           IW, LIW )
            CALL SMUMPS_190(.FALSE.,FLOP1,KEEP)
          ENDIF
      ELSEIF ( MSGTAG .EQ. END_NIV2_LDLT ) THEN
          INODE = BUFR( 1 )
          CALL SMUMPS_507(N, IPOOL, LPOOL,
     *         PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     *         STEP, -INODE )
          IF (KEEP(47) .GE. 3) THEN
             CALL SMUMPS_500(
     $            IPOOL, LPOOL,
     *            PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *            MYID, STEP, N, ND, FILS )
          ENDIF
      ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          IFLAG  = -001
          IERROR = MSGSOU
          GOTO 100
      ELSEIF ( MSGTAG .EQ. MAITRE_DESC_BANDE ) THEN
        CALL SMUMPS_266( MYID,BUFR, LBUFR,
     *    LBUFR_BYTES, IWPOS,
     *    IWPOSCB,
     *    IPTRLU, LRLU, LRLUS, NBPROCFILS,
     *    N, IW, LIW, A, LA, NIRBDU,
     *    PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, COMP,
     *    KEEP, ITLOC,
     *    IFLAG, IERROR )
        IF ( IFLAG .LT. 0 ) GO to 500
      ELSEIF ( MSGTAG .EQ. MAITRE2           ) THEN
        CALL SMUMPS_268(MYID,BUFR, LBUFR, LBUFR_BYTES,
     *    PROCNODE_STEPS, SLAVEF, IWPOS, IWPOSCB,
     *    IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *    PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM, COMM_LOAD, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    KEEP, ND, FILS, FRERE, ITLOC,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF ( IFLAG .LT. 0 ) GO to 500
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO        ) THEN
        CALL SMUMPS_264( COMM_LOAD, ASS_IRECV,
     *   BUFR,  LBUFR, LBUFR_BYTES,
     *   PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM , IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, root, OPASSW, OPELIW, ITLOC, FILS,
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP, IPOOL, LPOOL, LEAF, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO_SYM_SLAVE    ) THEN
        CALL SMUMPS_263( COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, root, OPASSW, OPELIW, ITLOC, FILS,
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP, IPOOL, LPOOL, LEAF, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      ELSEIF ( MSGTAG .EQ. BLOC_FACTO_SYM    ) THEN
        CALL SMUMPS_274( COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, root, OPASSW, OPELIW, ITLOC, FILS,
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP, IPOOL, LPOOL, LEAF, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      ELSEIF ( MSGTAG .EQ. CONTSIPERENIV2    ) THEN
        CALL SMUMPS_265( COMM_LOAD, ASS_IRECV,
     *       BUFR, LBUFR,
     *       LBUFR_BYTES, PROCNODE_STEPS,
     *       SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, POSFAC,
     *       N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S, PTRAST,
     *       STEP, PIMASTER, PAMASTER, NBPROCFILS, COMP, root,
     *       OPASSW, OPELIW, ITLOC, NSTK_S,
     *       FILS, PTRARW, PTRAIW, INTARR, DBLARR, NBFIN, MYID, COMM,
     *       ICNTL, KEEP, IFLAG, IERROR, IPOOL, LPOOL, LEAF,
     *       ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF ( IFLAG .LT. 0 ) GO TO 100
      ELSEIF ( MSGTAG .EQ. MAPLIG            ) THEN
         INODE        = BUFR( 1 )
         ISON         = BUFR( 2 )
         NSLAVES_PERE = BUFR( 3 )
         NFRONT_PERE  = BUFR( 4 )
         NASS_PERE    = BUFR( 5 )
         LMAP         = BUFR( 6 )
         IF ( (NSLAVES_PERE.NE.0).AND.(KEEP(48).NE.0) ) THEN
             INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
             ISHIFT = NSLAVES_PERE+1
             TAB_POS_IN_PERE(1:NSLAVES_PERE+1, INIV2) =
     *         BUFR(7:7+NSLAVES_PERE)
             TAB_POS_IN_PERE(SLAVEF+2, INIV2) = NSLAVES_PERE
         ELSE
             ISHIFT = 0
         ENDIF
         IBEG = 7+ISHIFT
         CALL SMUMPS_210( COMM_LOAD, ASS_IRECV,
     *    BUFR, LBUFR, LBUFR_BYTES,
     *    INODE, ISON, NSLAVES_PERE,
     *    BUFR(IBEG),
     *    NFRONT_PERE, NASS_PERE, LMAP,
     *    BUFR(IBEG+NSLAVES_PERE),
     *    PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB,
     *    IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *    PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP,
     *    IFLAG, IERROR, MYID, COMM, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, ICNTL, KEEP, root,
     *    OPASSW, OPELIW,
     *    ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *    ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
         IF ( IFLAG .LT. 0 ) GO TO 100
      ELSE IF ( MSGTAG .EQ. FACTOR ) THEN
        CALL SMUMPS_267( BUFR, LBUFR, LBUFR_BYTES,
     *       N, MSGSOU, MYID,
     *       STEP, PTLUST_S, KEEP(10), KEEP(28),
     *       SLAVEF, IW, LIW, A, LA, COMM  )
      ELSE IF ( MSGTAG .EQ. ROOT_CONT_STATIC ) THEN
        CALL SMUMPS_272(
     *        BUFR, LBUFR, LBUFR_BYTES,
     *        root, N, IW, LIW, A, LA, NIRBDU, NBPROCFILS,
     *        LRLU, IPTRLU, IWPOS, IWPOSCB,
     *        PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *        COMP, LRLUS, IPOOL, LPOOL, LEAF,
     *        FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *        KEEP, IFLAG, IERROR, COMM, COMM_LOAD, ITLOC,
     *        ND, PROCNODE_STEPS, SLAVEF)
        IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. ROOT_NON_ELIM_CB ) THEN
        IROOT  = KEEP( 38 )
        MSGSOU = SMUMPS_275( STEP(IROOT), PROCNODE_STEPS,
     *           SLAVEF )
        IF ( PTLUST_S( STEP(IROOT)) .EQ. 0 ) THEN
          CALL MPI_RECV( TMP, 2 * KEEP(34), MPI_PACKED,
     *                   MSGSOU, ROOT_2SLAVE,
     *                   COMM, STATUS, IERR )
          CALL SMUMPS_270( TMP( 1 ), TMP( 2 ),
     *    ROOT,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM, COMM_LOAD,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND )
          IF ( IFLAG .LT. 0 ) GOTO 500
        END IF
        CALL SMUMPS_272(
     *       BUFR, LBUFR, LBUFR_BYTES,
     *       root, N, IW, LIW, A, LA, NIRBDU, NBPROCFILS,
     *       LRLU, IPTRLU, IWPOS, IWPOSCB,
     *       PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *       COMP, LRLUS, IPOOL, LPOOL, LEAF,
     *       FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *       KEEP, IFLAG, IERROR, COMM, COMM_LOAD, ITLOC,
     *       ND, PROCNODE_STEPS, SLAVEF )
        IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. ROOT_2SON ) THEN
         ISON = BUFR( 1 )
         NELIM = BUFR( 2 )
         CALL SMUMPS_271( COMM_LOAD, ASS_IRECV,
     *    ISON, NELIM, ROOT,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
         IF ( IFLAG .LT. 0 ) GO TO 100
      ELSE IF ( MSGTAG .EQ. ROOT_2SLAVE ) THEN
          TOT_ROOT_SIZE    = BUFR( 1 )
          TOT_CONT_TO_RECV = BUFR( 2 )
          CALL SMUMPS_270( TOT_ROOT_SIZE,
     *    TOT_CONT_TO_RECV, ROOT,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM, COMM_LOAD,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND )
          IF ( IFLAG .LT. 0 ) GO TO 100
      ELSE IF ( MSGTAG .EQ. ROOT_NELIM_INDICES ) THEN
         ISON         = BUFR( 1 )
         NELIM        = BUFR( 2 )
         NSLAVES_PERE = BUFR( 3 )
         CALL SMUMPS_273( ROOT,
     *    ISON, NELIM, NSLAVES_PERE, BUFR(4), BUFR(4+BUFR(2)),
     *    BUFR(4+2*BUFR(2)),
     *
     *    PROCNODE_STEPS,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, ITLOC, COMP,
     *    IFLAG, IERROR,
     *    IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,
     *    COMM, COMM_LOAD, FILS, ND)
         IF ( IFLAG .LT. 0 ) GO TO 500
      ELSE IF ( MSGTAG .EQ. UPDATE_LOAD ) THEN
         WRITE(*,*) "Internal error 3 in SMUMPS_322"
         CALL SMUMPS_ABORT()
      ELSE IF ( MSGTAG .EQ. TAG_DUMMY   ) THEN
      ELSE
         IF ( LP > 0 )
     *     WRITE(LP,*) MYID,': Internal error, routine TRAITER_MESSAGE.'
         IFLAG = -100
         IERROR= MSGTAG
         GOTO 500
      ENDIF
 100  CONTINUE
      RETURN
 500  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_322
      RECURSIVE SUBROUTINE SMUMPS_280(
     *    COMM_LOAD, ASS_IRECV,
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT ,
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER KEEP(150), ICNTL(40)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST( KEEP(28) ),
     &        PTLUST_S(KEEP(28)), PTRAST(KEEP(28))
      INTEGER STEP(N),
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      REAL OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      REAL DBLARR( MAX(1,KEEP(13)) )
      INTEGER MSGSOU, MSGTAG, MSGLEN, IERR
      INTEGER allocok, OLDSIZE
      MSGSOU = STATUS( MPI_SOURCE )
      MSGTAG = STATUS( MPI_TAG )
      CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
      IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
        IFLAG  = -20
        IERROR = MSGLEN
         WRITE(*,*) ' RECEPTION BUF TOO SMALL, Msgtag/len=',
     *                MSGTAG,MSGLEN
        CALL SMUMPS_44( MYID, SLAVEF, COMM )
        RETURN
       ENDIF
       CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     *                 MSGTAG,
     *                 COMM, STATUS, IERR )
       CALL SMUMPS_322(
     *      COMM_LOAD, ASS_IRECV,
     *      MSGSOU, MSGTAG, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      PROCNODE_STEPS, POSFAC,
     *      IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *      IERROR, COMM,
     *      NBPROCFILS,
     *      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *
     *      root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *      INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *      LPTRAR, NELT, FRTPTR, FRTELT,
     *
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *      )
      RETURN
      END SUBROUTINE SMUMPS_280
      RECURSIVE SUBROUTINE SMUMPS_329(
     *    COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV,
     *    MESSAGE_RECEIVED, MSGSOU, MSGTAG,
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL, INTENT (IN)  :: BLOCKING
      LOGICAL, INTENT (IN)  :: SET_IRECV
      LOGICAL, INTENT (INOUT) :: MESSAGE_RECEIVED
      INTEGER, INTENT (IN) :: MSGSOU, MSGTAG
      INTEGER KEEP(150), ICNTL(40)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST( KEEP(28) ),
     &        PTLUST_S(KEEP(28)), PTRAST(KEEP(28))
      INTEGER STEP(N),
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      REAL OPASSW, OPELIW
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      REAL DBLARR( MAX(1,KEEP(13)) )
       LOGICAL FLAG, RIGHT_MESS, FLAGbis
       INTEGER LP, MSGSOU_LOC, MSGTAG_LOC
       INTEGER IERR
       INTEGER STATUS_BIS( MPI_STATUS_SIZE )
       INTEGER, SAVE :: RECURS = 0
      CALL SMUMPS_467(COMM_LOAD, KEEP)
      RECURS = RECURS + 1
      LP = ICNTL(1)
      IF ( MESSAGE_RECEIVED ) THEN
        MSGSOU_LOC = MPI_ANY_SOURCE
        MSGTAG_LOC = MPI_ANY_TAG
        GOTO 250
      ENDIF
      IF ( ASS_IRECV .NE. MPI_REQUEST_NULL) THEN
      RIGHT_MESS = .TRUE.
       IF (BLOCKING) THEN
         CALL MPI_WAIT(ASS_IRECV,
     *                STATUS, IERR)
         FLAG = .TRUE.
         IF ( ( (MSGSOU.NE.MPI_ANY_SOURCE) .OR.
     *      (MSGTAG.NE.MPI_ANY_TAG) )  ) THEN
           IF ( MSGSOU.NE.MPI_ANY_SOURCE) THEN
             RIGHT_MESS = MSGSOU.EQ.STATUS(MPI_SOURCE)
           ENDIF
           IF ( MSGTAG.NE.MPI_ANY_TAG) THEN
             RIGHT_MESS =
     &       ( (MSGTAG.EQ.STATUS(MPI_TAG)).AND.RIGHT_MESS )
           ENDIF
           IF (.NOT.RIGHT_MESS) THEN
             CALL MPI_PROBE(MSGSOU,MSGTAG,
     *           COMM, STATUS_BIS, IERR)
           ENDIF
         ENDIF
       ELSE
        CALL MPI_TEST(ASS_IRECV,
     *             FLAG, STATUS, IERR)
       ENDIF
       IF (IERR.LT.0) THEN
        IFLAG = -20
        IF (LP.GT.0)
     &  write(LP,*) ' Error return from MPI_TEST ',
     &     IFLAG, ' in SMUMPS_329'
        CALL SMUMPS_44( MYID, SLAVEF, COMM )
        RETURN
       ENDIF
       IF ( FLAG ) THEN
         MESSAGE_RECEIVED = .TRUE.
         MSGSOU_LOC = STATUS( MPI_SOURCE )
         MSGTAG_LOC = STATUS( MPI_TAG )
           IF (.NOT.RIGHT_MESS) RECURS = RECURS + 10
         CALL SMUMPS_322( COMM_LOAD, ASS_IRECV,
     *      MSGSOU_LOC, MSGTAG_LOC, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      PROCNODE_STEPS, POSFAC,
     *      IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA,
     *      NIRBDU, PTRIST, PTLUST_S,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *      IERROR, COMM,
     *      NBPROCFILS,
     *      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *
     *      root, OPASSW, OPELIW, ITLOC, FILS,
     *      PTRARW, PTRAIW,
     *      INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *      LPTRAR, NELT, FRTPTR, FRTELT,
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
           IF (.NOT.RIGHT_MESS) RECURS = RECURS - 10
          IF ( IFLAG .LT. 0 ) RETURN
           IF (.NOT.RIGHT_MESS) THEN
              IF (ASS_IRECV .NE. MPI_REQUEST_NULL) THEN
                stop
              ENDIF
             CALL MPI_IPROBE(MSGSOU,MSGTAG,
     *           COMM, FLAGbis, STATUS, IERR)
             IF (FLAGbis) THEN
               MSGSOU_LOC = STATUS( MPI_SOURCE )
               MSGTAG_LOC = STATUS( MPI_TAG )
               CALL SMUMPS_280( COMM_LOAD, ASS_IRECV,
     *            STATUS, BUFR, LBUFR,
     *            LBUFR_BYTES,
     *            PROCNODE_STEPS, POSFAC,
     *            IWPOS, IWPOSCB, IPTRLU,
     *            LRLU, LRLUS, N, IW, LIW, A, LA,
     *            NIRBDU, PTRIST, PTLUST_S,
     *            PTRAST, STEP, PIMASTER, PAMASTER,
     *            NSTK_S, COMP, IFLAG,
     *            IERROR, COMM,
     *            NBPROCFILS,
     *            IPOOL, LPOOL, LEAF,
     *            NBFIN, MYID, SLAVEF,
     *
     *            root, OPASSW, OPELIW, ITLOC,
     *            FILS, PTRARW, PTRAIW,
     *            INTARR, DBLARR, ICNTL,
     *            KEEP, ND, FRERE,
     *            LPTRAR, NELT, FRTPTR, FRTELT,
     *            ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
                  IF ( IFLAG .LT. 0 ) RETURN
             ENDIF
           ENDIF
       ENDIF
      ELSE
         IF (BLOCKING) THEN
           CALL MPI_PROBE(MSGSOU,MSGTAG,
     *           COMM, STATUS, IERR)
           FLAG = .TRUE.
         ELSE
           CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,
     *           COMM, FLAG, STATUS, IERR)
         ENDIF
         IF (FLAG) THEN
          MSGSOU_LOC = STATUS( MPI_SOURCE )
          MSGTAG_LOC = STATUS( MPI_TAG )
          MESSAGE_RECEIVED = .TRUE.
          CALL SMUMPS_280( COMM_LOAD, ASS_IRECV,
     *      STATUS, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      PROCNODE_STEPS, POSFAC,
     *      IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA,
     *      NIRBDU, PTRIST, PTLUST_S,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *      IERROR, COMM,
     *      NBPROCFILS,
     *      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *
     *      root, OPASSW, OPELIW, ITLOC,
     *      FILS, PTRARW, PTRAIW,
     *      INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *      LPTRAR, NELT, FRTPTR, FRTELT,
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
          IF ( IFLAG .LT. 0 ) RETURN
         ENDIF
      ENDIF
 250  CONTINUE
      RECURS  = RECURS - 1
      IF ( NBFIN .EQ. 0 ) RETURN
      IF ( RECURS .GT. 5 ) RETURN
      IF ( KEEP(36).EQ.1 .AND. SET_IRECV  .AND.
     *      (ASS_IRECV.EQ.MPI_REQUEST_NULL) .AND.
     *    MESSAGE_RECEIVED ) THEN
       CALL MPI_IRECV ( BUFR(1),
     *      LBUFR_BYTES, MPI_PACKED, MPI_ANY_SOURCE,
     *      MPI_ANY_TAG, COMM,
     *      ASS_IRECV, IERR )
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_329
      SUBROUTINE SMUMPS_255( INFO1,
     *    ASS_IRECV, STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES,
     *    COMM,
     *    MYID, SLAVEF)
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER COMM
      INTEGER MYID, SLAVEF, INFO1
      LOGICAL FLAG
      INTEGER LP, MSGSOU_LOC, MSGTAG_LOC
      INTEGER IERR, DUMMY
      IF (ASS_IRECV.EQ.MPI_REQUEST_NULL) RETURN
      CALL MPI_TEST(ASS_IRECV,
     *             FLAG, STATUS, IERR)
      IF (.NOT.FLAG) THEN
         DUMMY = 1
         CALL SMUMPS_62
     *    (DUMMY, MYID, TAG_DUMMY, COMM, IERR)
         CALL MPI_WAIT(ASS_IRECV,
     *                STATUS, IERR)
      ENDIF
      MSGSOU_LOC = STATUS( MPI_SOURCE )
      MSGTAG_LOC = STATUS( MPI_TAG )
      RETURN
      END SUBROUTINE SMUMPS_255
      SUBROUTINE SMUMPS_180(
     *    INFO1, STATUS, BUFR, LBUFR, LBUFR_BYTES,
     *    COMM_NODES, COMM_LOAD, SLAVEF )
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER COMM_NODES, COMM_LOAD, SLAVEF, INFO1
      LOGICAL FLAG, BUFFERS_EMPTY, BUFFERS_EMPTY_ON_ALL_PROCS
      INTEGER MSGSOU_LOC, MSGTAG_LOC, COMM_EFF
      INTEGER IERR, DUMMY
      IF (SLAVEF.EQ.1) RETURN
 10   CONTINUE
      FLAG = .TRUE.
      DO WHILE ( FLAG )
        COMM_EFF = COMM_NODES
        CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,
     *       COMM_NODES, FLAG, STATUS, IERR)
        IF ( .NOT. FLAG ) THEN
          COMM_EFF = COMM_LOAD
          CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     *         COMM_LOAD, FLAG, STATUS, IERR)
        END IF
        IF (FLAG) THEN
            MSGSOU_LOC = STATUS( MPI_SOURCE )
            MSGTAG_LOC = STATUS( MPI_TAG )
               CALL MPI_RECV( BUFR, LBUFR_BYTES,
     *             MPI_PACKED, MSGSOU_LOC,
     *             MSGTAG_LOC, COMM_EFF, STATUS, IERR )
           ENDIF
         END DO
        CALL SMUMPS_469(BUFFERS_EMPTY)
        CALL MPI_ALLREDUCE(BUFFERS_EMPTY,
     *                     BUFFERS_EMPTY_ON_ALL_PROCS,
     *                     1, MPI_LOGICAL, MPI_LAND,
     *                     COMM_NODES, IERR)
        IF ( .NOT. BUFFERS_EMPTY_ON_ALL_PROCS ) THEN
          write(*,*) "treating pending requests..."
          GOTO 10
        END IF
      RETURN
      END SUBROUTINE SMUMPS_180
      SUBROUTINE SMUMPS_245(root, N, A, LA, IW, LIW, W, LWC, 
     * IWCB,LIWW,RHS,RHSSOL,NA,LNA,NE_STEPS, W2, MTYPE, ICNTL,
     * STEP, FRERE,FILS, PTRIST,IW1,LIW1,
     * PROCNODE_STEPS, SLAVEF,
     * INFO, KEEP, COMM, COMM_NODES, MYID,
     * MYID_NODES,
     * BUFR, LBUFR, LBUFR_BYTES,
     * 
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE
     * )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE ( SMUMPS_ROOT_STRUC ) :: root
      INTEGER LA,LWC,N,LIW,MTYPE,LIW1,LIWW,LNA
      INTEGER ICNTL(40),INFO(40), KEEP(150)
      INTEGER IW(LIW),IW1(LIW1),NA(LNA),NE_STEPS(KEEP(28)),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)), FILS(N), PTRIST(KEEP(28))
      REAL    A(LA), W(LWC), RHS(N),
     *        RHSSOL(N), W2(KEEP(133))
      INTEGER SLAVEF, COMM, COMM_NODES, MYID, MYID_NODES
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER MP, LP, LDIAG, LWC_LOC
      INTEGER K,I
      INTEGER LPOOL,MYLEAF
      INTEGER NSTK_S,IPOOL,PTRICB,PTRACB
      INTEGER POOLSS, MTYPE_LOC
      INTEGER SIZE_ROOT, IPT_RHS_ROOT
      INTEGER IERR
      INTEGER       IOLDPS, IAPOS,
     *              LOCAL_M     ,
     *              LOCAL_N     ,
     *              MASTER_ROOT
           INTEGER INODE, IPOS, LIELL, NPIV,J1,JJ
      INTEGER IZERO
      REAL ZERO
      PARAMETER (IZERO = 0 )
      PARAMETER (ZERO = 0.0E0)
      EXTERNAL SMUMPS_248, SMUMPS_249
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      LP      = ICNTL(1)
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
      NSTK_S   = 1
      POOLSS = NSTK_S + KEEP(28)
      PTRICB = POOLSS + KEEP(28) + 1
      PTRACB = PTRICB + KEEP(28)
      IPOOL  = PTRACB + KEEP(28)
      LPOOL  = LIW1 - IPOOL + 1
      IF ( LPOOL .ne. KEEP(28) + 1 ) THEN
        WRITE(MP,*) ' warning solve: lpool.NE.n+1 :',lpool,n+1
      ENDIF
      LWC_LOC = LWC
      SIZE_ROOT  = 0
      IPT_RHS_ROOT  = LWC
      IF ( KEEP( 38 ) .NE. 0 .OR. KEEP( 20 ) .NE. 0 ) THEN
          IF ( KEEP( 38 ) .ne. 0 ) THEN
            MASTER_ROOT = SMUMPS_275(STEP( KEEP(38)),
     *                    PROCNODE_STEPS, SLAVEF )
          ELSE
            MASTER_ROOT = SMUMPS_275(STEP(KEEP(20)),
     *                    PROCNODE_STEPS, SLAVEF )
          END IF
          IF ( MASTER_ROOT .eq. MYID_NODES ) THEN
            IF ( KEEP( 38 ) .NE. 0 ) THEN
              SIZE_ROOT    = root%TOT_ROOT_SIZE
            ELSE
              SIZE_ROOT    = IW(PTRIST(STEP(KEEP(20))))
            END IF
            IPT_RHS_ROOT = LWC - SIZE_ROOT + 1
            LWC_LOC = LWC_LOC - SIZE_ROOT
          END IF
      END IF
      IF ( KEEP( 50 ) .eq. 0 ) THEN
        MTYPE_LOC = MTYPE
      ELSE
        MTYPE_LOC = 1
      ENDIF
      CALL SMUMPS_248(N, A(1), LA, IW(1), LIW, W(1),
     *           LWC_LOC, RHS, 
     *           IW1(PTRICB), IWCB, LIWW,
     *           RHSSOL,NE_STEPS, NA, LNA, STEP, FRERE,FILS,
     *           IW1(NSTK_S),IW1(IPOOL),LPOOL,PTRIST,MYLEAF,INFO,
     *           KEEP, PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           W( IPT_RHS_ROOT ), MTYPE_LOC, 
     * 
     *           ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *           )
      IF ( INFO(1) .LT. 0 ) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP,*) MYID,
     *    ': ** ERROR RETURN FROM SMUMPS_248,INFO(1:2)=',
     *    INFO(1:2)
        END IF
        RETURN
      END IF
      IF ( KEEP( 38 ) .NE. 0 .and. root%yes ) THEN
        IOLDPS = PTRIST(STEP(KEEP(38)))
        LOCAL_M = IW( IOLDPS + 2 )
        LOCAL_N = IW( IOLDPS + 1 )
        IAPOS   = IW( IOLDPS + 4 )
          IF ( .NOT. ASSOCIATED( root%QR_TAU ) )
     *          ALLOCATE( root%QR_TAU(1))
          IF ( .NOT. ASSOCIATED( root%GROW ) )
     *          ALLOCATE( root%GROW( 1 ) )
          IF ( .NOT. ASSOCIATED( root%GSIN ) )
     *          ALLOCATE( root%GSIN( 1 ) )
          IF ( .NOT. ASSOCIATED( root%GCOS ) )
     *          ALLOCATE( root%GCOS( 1 ) )
          CALL MPI_BARRIER(COMM_NODES, IERR)
#if defined(try_null_space_not_tested)
          CALL SMUMPS_352( root%DESCRIPTOR, root%DESCB,
     *       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     *       root%MBLOCK, root%NBLOCK,
     *       root%IPIV, root%LPIV, MASTER_ROOT, MYID_NODES,
     *       COMM_NODES,
     *       W( IPT_RHS_ROOT ),
     *       root%TOT_ROOT_SIZE, A( IAPOS ),
     *       INFO(1), MTYPE, KEEP(50), KEEP(19),
     *       root%QR_TAU, W(1), LWC_LOC, KEEP(17),
     *       root%MAXG, root%GIND, root%GROW, root%GCOS, root%GSIN )
#else
          CALL SMUMPS_286( root%DESCRIPTOR, root%DESCB,
     *       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     *       root%MBLOCK, root%NBLOCK,
     *       root%IPIV, root%LPIV, MASTER_ROOT, MYID_NODES,
     *       COMM_NODES,
     *       W( IPT_RHS_ROOT ),
     *       root%TOT_ROOT_SIZE, A( IAPOS ),
     *       INFO(1), MTYPE, KEEP(50))
#endif
      ELSE IF ( KEEP(20) .NE. 0 ) THEN
        IF ( MYID_NODES .eq.  SMUMPS_275( STEP(KEEP(20)),
     *        PROCNODE_STEPS, SLAVEF ) ) THEN
#if defined(try_null_space)
         IF ( KEEP(60) .eq. 0 ) THEN
           CALL SMUMPS_296(
     *     SIZE_ROOT,A( IW( PTRIST(STEP(KEEP(20))) + 4 ) ),ROOT,
     *     W( IPT_RHS_ROOT ), W(1),LWC_LOC,KEEP,MTYPE,INFO)
         ELSE
#endif
           W( IPT_RHS_ROOT : IPT_RHS_ROOT + SIZE_ROOT - 1 ) = 0.0E0
#if defined(try_null_space)
         END IF
#endif
        END IF
        CALL MPI_BARRIER( COMM_NODES, IERR )
      ELSE
        CALL MPI_BARRIER( COMM_NODES, IERR )
      END IF
      IF ( INFO(1) .LT. 0 ) RETURN
      IF ( KEEP( 50 ) .eq. 0 ) THEN
        MTYPE_LOC = MTYPE
      ELSE
        MTYPE_LOC = IZERO
      ENDIF
      CALL SMUMPS_249( N, A, LA, IW, LIW, W(1), LWC_LOC, RHS,
     *              RHSSOL, IW1(PTRICB),IW1(PTRACB),IWCB,LIWW,
     *              IW1(POOLSS),W2,
     *              NE_STEPS, NA, LNA, STEP, FRERE,FILS,
     *              IW1(IPOOL),LPOOL,PTRIST,MYLEAF,INFO,
     *              PROCNODE_STEPS, SLAVEF, COMM_NODES,MYID_NODES,
     *              BUFR, LBUFR, LBUFR_BYTES, KEEP, W( IPT_RHS_ROOT ),
     *              MTYPE_LOC, 
     *              ISTEP_TO_INIV2, TAB_POS_IN_PERE )
       IWCB(1:N)=0
       DO INODE = 1, N
         IF ( STEP(INODE) > 0 ) THEN
           IF ( MYID_NODES .eq. SMUMPS_275(
     *          STEP(INODE), PROCNODE_STEPS, SLAVEF) )
     *          THEN
                 IF ( KEEP(38).eq.INODE.or.KEEP(20).eq.INODE) THEN
                    IPOS = PTRIST(STEP(INODE)) 
                    LIELL = IW(IPOS+3)
                    NPIV = LIELL
                    IPOS= PTRIST(STEP(INODE))+5
                 ELSE
                  IPOS = PTRIST(STEP(INODE)) + 2
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(STEP(INODE)) + 5 )
                 END IF
                 IF (MTYPE.eq.1) THEN
                   J1=IPOS+1+LIELL
                 ELSE
                   J1=IPOS+1
                 END IF
                 DO JJ=J1,J1+NPIV-1
                     IWCB(IW(JJ))=1
                 ENDDO
           END IF
         END IF
       END DO
       DO I =1, N
         IF (IWCB(I) .NE.1) THEN
           RHSSOL(I)=ZERO
         END IF
       ENDDO
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        K = MIN0(10,N)
        IF (LDIAG.EQ.4) K = N
        WRITE (MP,99992)
        IF (N.GT.0) WRITE (MP,99993) (RHSSOL(I),I=1,K)
      ENDIF
      RETURN
99993 FORMAT (' RHSSOL'/(1X,1P,5E14.6))
99992 FORMAT (//' LEAVING SOLVE (MPI41C) WITH')
      END SUBROUTINE SMUMPS_245
      SUBROUTINE SMUMPS_270( TOT_ROOT_SIZE,
     *    TOT_CONT_TO_RECV, ROOT,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM, COMM_LOAD,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: ROOT
      INTEGER KEEP(150), ICNTL(40)
      INTEGER TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) ), ND( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, COMM_LOAD
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      REAL OPASSW, OPELIW
      INTEGER ITLOC(N), FILS(N), PTRARW(N), PTRAIW(N)
      INTEGER INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
      INTEGER NEW_LOCAL_M, NEW_LOCAL_N
      INTEGER OLD_LOCAL_M, OLD_LOCAL_N
      INTEGER LREQI, LREQA, IROOT
      INTEGER POSHEAD, IPOS_SON
      LOGICAL MASTER_OF_ROOT
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER NUMROC, SMUMPS_275
      EXTERNAL NUMROC, SMUMPS_275
      IROOT = KEEP( 38 )
      MASTER_OF_ROOT = ( MYID .EQ. SMUMPS_275( STEP(IROOT),
     *                   PROCNODE_STEPS, SLAVEF ) )
      NEW_LOCAL_M  = NUMROC( TOT_ROOT_SIZE, root%MBLOCK,
     *               root%MYROW, 0, root%NPROW )
      NEW_LOCAL_M  = MAX( 1, NEW_LOCAL_M )
      NEW_LOCAL_N  = NUMROC( TOT_ROOT_SIZE, root%NBLOCK,
     *               root%MYCOL, 0, root%NPCOL )
      IF ( MASTER_OF_ROOT ) THEN
        LREQI = 6 + 2 * TOT_ROOT_SIZE
      ELSE
        LREQI = 6
      END IF
      LREQA = NEW_LOCAL_M * NEW_LOCAL_N
      IF ( LRLU . LT. LREQA .OR.
     *     IWPOS + LREQI - 1. GT. IWPOSCB )THEN
           IF ( LRLUS .LT. LREQA ) THEN
             IFLAG  = -9
             IERROR = LREQA - LRLUS
             GOTO 700
           END IF
           CALL SMUMPS_94( N, KEEP(28), IW, LIW, A, LA,
     *           NIRBDU, LRLU, IPTRLU,
     *           IWPOS, IWPOSCB, PTRIST, PTRAST,
     *           STEP, PIMASTER, PAMASTER, ITLOC )
           COMP = COMP + 1
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB compress root2slave:LRLU,LRLUS=',
     *            LRLU, LRLUS
                  IFLAG = -9
                  IERROR = LREQA - LRLUS
                  GOTO 700
           END IF
           IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
              IFLAG  = -8
              IERROR = IWPOS + LREQI - 1 - IWPOSCB
              GOTO 700
           END IF
      END IF
      PTLUST_S(STEP( IROOT )) = IWPOS
      IWPOS           = IWPOS + LREQI
      PTRAST (STEP(IROOT)) = POSFAC
      POSFAC           = POSFAC + LREQA
      LRLU   = LRLU  - LREQA
      LRLUS  = LRLUS - LREQA
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,LREQA,KEEP)
      POSHEAD = PTLUST_S( STEP(IROOT))
      IW( POSHEAD ) = 0
      IW( POSHEAD + 1 ) = NEW_LOCAL_N
      IW( POSHEAD + 2 ) = NEW_LOCAL_M
      IW( POSHEAD + 4 ) = PTRAST( STEP(IROOT))
      IW( POSHEAD + 5 ) = 0
      IF ( MASTER_OF_ROOT ) THEN
        IW( POSHEAD + 3 ) = TOT_ROOT_SIZE
      ELSE
        IW( POSHEAD + 3 ) = 0
      ENDIF
      root%TOT_ROOT_SIZE = TOT_ROOT_SIZE
      IF ( KEEP( 50 ) .eq. 0 .OR. KEEP(50) .eq. 2 ) THEN
      OPELIW = OPELIW + ( REAL(2*TOT_ROOT_SIZE) *
     *         REAL(TOT_ROOT_SIZE) * REAL(TOT_ROOT_SIZE ) / REAL(3)
     *       - 0.5E0 * REAL( TOT_ROOT_SIZE ) * REAL( TOT_ROOT_SIZE )
     *       - REAL( TOT_ROOT_SIZE ) / REAL( 6 ) )
     *       / REAL( root%NPROW * root%NPCOL )
      ELSE
      OPELIW = OPELIW + ( REAL(TOT_ROOT_SIZE) *
     *         REAL( TOT_ROOT_SIZE) * 
     *         REAL( TOT_ROOT_SIZE + 1 ) )
     *          / REAL( 3 * root%NPROW * root%NPCOL )
      END IF
      IF ( PTRIST(STEP( IROOT )) .LT. 0 ) THEN
        A(PTRAST(STEP(IROOT)):PTRAST(STEP(IROOT))+LREQA-1)=REAL(ZERO)
        PTRIST(STEP( IROOT )) = 0
        PAMASTER(STEP( IROOT )) = 0
      ELSE IF ( PTRIST( STEP(IROOT)) .eq. 0 ) THEN
        A(PTRAST(STEP(IROOT)):PTRAST(STEP(IROOT))+LREQA-1)=REAL(ZERO)
        CALL SMUMPS_35( N, root, IROOT,
     *       A(PTRAST(STEP(IROOT))), NEW_LOCAL_M, NEW_LOCAL_N,
     *       FILS, PTRAIW, PTRARW,
     *       INTARR, DBLARR, KEEP,
     *       MYID )
      ELSE
        OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) )
        OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1 )
        IF ( TOT_ROOT_SIZE .eq. root%ROOT_SIZE ) THEN
          IF ( LREQA .NE. OLD_LOCAL_M * OLD_LOCAL_N ) THEN
            write(*,*) 'error 1 in PROCESS_ROOT2SLAVE',
     *     OLD_LOCAL_M, OLD_LOCAL_N
            stop
          END IF
          CALL SCOPY(LREQA, A(PAMASTER(STEP(IROOT))), 1,
     *                      A(PTRAST(STEP(IROOT))), 1)
        ELSE
          CALL SMUMPS_96( A( PTRAST(STEP(IROOT))), 
     *        NEW_LOCAL_M,
     *        NEW_LOCAL_N, A( PAMASTER( STEP(IROOT)) ), OLD_LOCAL_M,
     *        OLD_LOCAL_N )
        END IF
        IF ( PTRIST( STEP( IROOT ) ) .NE. 0 ) THEN
           IPOS_SON= PTRIST( STEP(IROOT))
           CALL SMUMPS_152( MYID,N, IPOS_SON,
     *          PTRAST(STEP(IROOT)),
     *          IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *          IWPOSCB, LA, KEEP 
     $         )
           PTRIST(STEP( IROOT ))   = 0
           PAMASTER(STEP( IROOT )) = 0
        END IF
      END IF
      NBPROCFILS(STEP(IROOT))=NBPROCFILS(STEP(IROOT)) + TOT_CONT_TO_RECV
      IF ( NBPROCFILS(STEP(IROOT)) .eq. 0 ) THEN
        CALL SMUMPS_507( N, IPOOL, LPOOL, PROCNODE_STEPS,
     *       SLAVEF, KEEP(28), KEEP(76), STEP, IROOT + N )
        IF (KEEP(47) .GE. 3) THEN
           CALL SMUMPS_500(
     $          IPOOL, LPOOL, 
     *          PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *          MYID, STEP, N, ND, FILS )
        ENDIF
      END IF
      RETURN
 700  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_270
      SUBROUTINE SMUMPS_96
     *( NEW, M_NEW, N_NEW,OLD, M_OLD, N_OLD )
      INTEGER M_NEW, N_NEW, M_OLD, N_OLD
      REAL NEW( M_NEW, N_NEW ), OLD( M_OLD, N_OLD )
      INTEGER J
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      DO J = 1, N_OLD
        NEW( 1: M_OLD, J ) = OLD( 1: M_OLD, J )
        NEW( M_OLD + 1: M_NEW, J ) = REAL(ZERO)
      END DO
      NEW( 1: M_NEW,N_OLD + 1: N_NEW ) = REAL(ZERO)
      RETURN
      END SUBROUTINE SMUMPS_96
      INTEGER FUNCTION SMUMPS_505(KEEP)
      IMPLICIT NONE
      INTEGER KEEP(150)
      SMUMPS_505 = KEEP(28) + 1 + 3
      RETURN
      END FUNCTION SMUMPS_505
      SUBROUTINE SMUMPS_506(IPOOL, LPOOL, LEAF)
      IMPLICIT NONE
      INTEGER LPOOL, LEAF
      INTEGER IPOOL(LPOOL)
      IPOOL(LPOOL-2) = 0
      IPOOL(LPOOL-1) = 0
      IPOOL(LPOOL)   = LEAF-1
      RETURN
      END SUBROUTINE SMUMPS_506
      SUBROUTINE SMUMPS_507
     *           (N, POOL, LPOOL, PROCNODE, SLAVEF,
     *           K28, K76, STEP, INODE)
      IMPLICIT NONE
      INTEGER N, INODE, LPOOL, K28, SLAVEF, K76
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(K28)
      EXTERNAL SMUMPS_170
      LOGICAL SMUMPS_170, ATM_CURRENT_NODE
      INTEGER NBINSUBTREE, NBTOP, INODE_EFF
      INTEGER IPOS1, IPOS2, ISWAP
      ATM_CURRENT_NODE = ( K76 == 2 .OR. K76 ==3)
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      IF (INODE > N ) THEN
        INODE_EFF = INODE - N
      ELSE IF (INODE < 0) THEN
        INODE_EFF = - INODE
      ELSE
        INODE_EFF = INODE
      ENDIF
      IF ( SMUMPS_170(STEP(INODE_EFF),
     *                             PROCNODE, SLAVEF) ) THEN
        POOL(NBINSUBTREE + 1 ) = INODE
        NBINSUBTREE = NBINSUBTREE + 1
      ELSE
        POOL( LPOOL - 2 - ( NBTOP + 1 ) ) = INODE
        NBTOP = NBTOP + 1
        IPOS1 = LPOOL - 2 - NBTOP
        IPOS2 = LPOOL - 2 - NBTOP + 1
 10     CONTINUE
        IF ( IPOS2 == LPOOL - 2 ) GOTO 20
        IF ( POOL(IPOS1) < 0 ) GOTO 20
        IF ( POOL(IPOS2) < 0 ) GOTO 30
        IF ( ATM_CURRENT_NODE ) THEN
          IF ( POOL(IPOS1) > N ) GOTO 20
          IF ( POOL(IPOS2) > N ) GOTO 30
        END IF
        GOTO 20
 30     CONTINUE
        ISWAP = POOL(IPOS1)
        POOL(IPOS1) = POOL(IPOS2)
        POOL(IPOS2) = ISWAP
        IPOS1 = IPOS1 + 1
        IPOS2 = IPOS2 + 1
        GOTO 10
 20     CONTINUE
      ENDIF
      POOL(LPOOL) = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      RETURN
      END SUBROUTINE SMUMPS_507
      LOGICAL FUNCTION SMUMPS_508(POOL, LPOOL)
      IMPLICIT NONE
      INTEGER LPOOL
      INTEGER POOL(LPOOL)
      INTEGER NBINSUBTREE, NBTOP
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      SMUMPS_508 = (NBINSUBTREE + NBTOP == 0)
      RETURN
      END FUNCTION SMUMPS_508
      SUBROUTINE SMUMPS_509( N, POOL, LPOOL, PROCNODE, SLAVEF,
     *           STEP, INODE, KEEP, MYID, ND )
      IMPLICIT NONE
      INTEGER INODE, LPOOL, SLAVEF, N
      INTEGER KEEP(150)
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(KEEP(28)),
     *        ND(KEEP(28))
      EXTERNAL SMUMPS_167, SMUMPS_283, SMUMPS_508
      LOGICAL SMUMPS_167, SMUMPS_283, SMUMPS_508
      INTEGER NBINSUBTREE, NBTOP, INSUBTREE, INODE_EFF, MYID
      LOGICAL LEFT, ATOMIC_SUBTREE
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF ( KEEP(76) > 3 .OR. KEEP(76) < 0 ) THEN
        WRITE(*,*) "Error 2 in SMUMPS_509: unknown strategy"
        CALL SMUMPS_ABORT()
      ENDIF
      ATOMIC_SUBTREE =  ( KEEP(76) == 1 .OR. KEEP(76) == 3)
      IF ( SMUMPS_508(POOL, LPOOL) ) THEN
        WRITE(*,*) "Error 1 in SMUMPS_509"
        CALL SMUMPS_ABORT()
      ENDIF
      IF ( .NOT. ATOMIC_SUBTREE ) THEN
        LEFT = (NBTOP == 0)
      ELSE
        IF ( INSUBTREE == 1 ) THEN
          IF (NBINSUBTREE == 0) THEN
            WRITE(*,*) "Error 3 in SMUMPS_509"
            CALL SMUMPS_ABORT()
          ENDIF
          LEFT = .TRUE.
        ELSE
          LEFT = ( NBTOP == 0)
        ENDIF
      ENDIF
      IF ( LEFT ) THEN
        INODE = POOL( NBINSUBTREE )
        NBINSUBTREE = NBINSUBTREE - 1
        IF ( INODE < 0 ) THEN
          INODE_EFF = -INODE
        ELSE IF ( INODE > N ) THEN
          INODE_EFF = INODE - N
        ELSE
          INODE_EFF = INODE
        ENDIF
        IF ( SMUMPS_167( STEP(INODE_EFF), PROCNODE, SLAVEF) ) THEN
          INSUBTREE = 1
        ELSE IF ( SMUMPS_283( STEP(INODE_EFF), PROCNODE,
     *              SLAVEF)) THEN
          INSUBTREE = 0
        END IF
      ELSE
        IF (NBTOP < 1 ) THEN
          WRITE(*,*) "Error 5 in SMUMPS_509", NBTOP
          CALL SMUMPS_ABORT()
        ENDIF
        INODE = POOL( LPOOL - 2 - NBTOP )
        NBTOP = NBTOP - 1
        IF ( INODE < 0 ) THEN
          INODE_EFF = -INODE
        ELSE IF ( INODE > N ) THEN
          INODE_EFF = INODE - N
        ELSE
          INODE_EFF = INODE
        ENDIF
      END IF
      POOL(LPOOL)     = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      POOL(LPOOL - 2) = INSUBTREE
      RETURN
      END SUBROUTINE SMUMPS_509
