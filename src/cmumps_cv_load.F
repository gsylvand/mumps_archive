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
      MODULE CMUMPS_LOAD
      implicit none
      PRIVATE
      PUBLIC :: CMUMPS_188, CMUMPS_185,
     *       CMUMPS_189,
     *       CMUMPS_184, CMUMPS_190,
     *       CMUMPS_183, CMUMPS_187, 
     *       CMUMPS_186,
     *       CMUMPS_409,
     *       CMUMPS_384, CMUMPS_461,
     *       CMUMPS_467, CMUMPS_471,
     *       CMUMPS_472,
     *       CMUMPS_500,
     $       CMUMPS_501
      REAL, DIMENSION(:),
     *       ALLOCATABLE, SAVE :: LOAD_FLOPS
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: BUF_LOAD_RECV
      INTEGER, SAVE :: LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES
      INTEGER, SAVE :: K50, K69, K35
      LOGICAL, SAVE :: BDC_MEM, BDC_POOL, BDC_SBTR
      REAL, DIMENSION(:),
     *       ALLOCATABLE, TARGET, SAVE :: WLOAD
#if ! defined(OLD_LOAD_MECHANISM)
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: FUTURE_NIV2
      REAL, SAVE :: DELTA_LOAD, DELTA_MEM
#else
      REAL, SAVE :: LAST_LOAD_SENT, DM_LAST_MEM_SENT
#endif
      INTEGER, SAVE :: CHECK_MEM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: IDWLOAD
      REAL, SAVE :: COST_SUBTREE
      REAL, SAVE :: ALPHA
      REAL, SAVE :: BETA
      INTEGER, SAVE :: MYID, NPROCS, COMM_LD
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: POOL_MEM
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: SBTR_MEM
      REAL, SAVE :: POOL_LAST_COST_SENT
      REAL, SAVE :: MIN_DIFF
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::MEM_SUBTREE
      INTEGER NB_SUBTREES,NIV1_FLAG,INDICE
      INTEGER INSIDE_SUBTREE
      PUBLIC :: NB_SUBTREES,MEM_SUBTREE,INSIDE_SUBTREE,NIV1_FLAG
      REAL, SAVE :: DM_SUMLU,
     *                   DM_THRES_MEM
      REAL, DIMENSION(:),
     *   ALLOCATABLE, SAVE :: DM_MEM
      CONTAINS
      SUBROUTINE CMUMPS_188( COST_SUBTREE_ARG, K64, K66 )
      IMPLICIT NONE
      REAL COST_SUBTREE_ARG
      INTEGER K64, K66
      REAL T64, T66
      T64 = MAX ( REAL(K64), REAL(1) )
      T64 = MIN ( T64, REAL(1000)  )
      T66 = MAX (REAL(K66), REAL(100))
      MIN_DIFF     =  ( T64 / REAL(1000)  )* 
     *                  T66 * REAL(1000000)
      DM_THRES_MEM = 30000.0E0
      COST_SUBTREE = COST_SUBTREE_ARG
      END SUBROUTINE CMUMPS_188
      SUBROUTINE CMUMPS_472(
     *  SLAVEF,
     *  KEEP,
     *  CAND_OF_NODE,
     *  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
       IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(150),SIZE_SLAVES_LIST
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER I
      IF ( KEEP(48) == 0 .OR. KEEP(48) .EQ. 3 ) THEN
         CALL CMUMPS_499(
     *        SLAVEF,
     *        KEEP,
     *        CAND_OF_NODE,
     *        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      ELSE IF ( KEEP(48) == 4 ) THEN
         CALL CMUMPS_504(
     *        SLAVEF,
     *        KEEP,
     *        CAND_OF_NODE,
     *        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID)
         DO I=1,NSLAVES_NODE
            IF(TAB_POS(I+1)-TAB_POS(I).LE.0)THEN
               WRITE(*,*)'probleme de partition dans 
     *CMUMPS_LOAD_SET_PARTI_ACTV_MEM_K10'
               CALL CMUMPS_ABORT()
            ENDIF
         ENDDO
      ELSE
        WRITE(*,*) "Strategy 5 not implemented"
        CALL CMUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_472
      SUBROUTINE CMUMPS_499(
     *  SLAVEF,
     *  KEEP,
     *  CAND_OF_NODE,
     *  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(150),SIZE_SLAVES_LIST
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER ITEMP, NMB_OF_CAND, NSLAVES_LESS
      REAL MSG_SIZE
      LOGICAL force_cand
      INTEGER  CMUMPS_12
      EXTERNAL CMUMPS_12
      IF ( KEEP(48) == 0 .AND. KEEP(50) .NE. 0) THEN
      write(*,*) "Internal error 2 in CMUMPS_499."
      CALL CMUMPS_ABORT()
      END IF
      IF ( KEEP(48) == 3 .AND. KEEP(50) .EQ. 0) THEN
      write(*,*) "Internal error 3 in CMUMPS_499."
      CALL CMUMPS_ABORT()
      END IF
      MSG_SIZE = REAL( NFRONT - NCB ) * REAL(NCB)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        force_cand = .FALSE.
      ELSE
        force_cand = (mod(KEEP(24),2).eq.0)
      END IF
      IF (force_cand) THEN
        ITEMP=CMUMPS_409
     *       (MEM_DISTRIB(0:SLAVEF-1),
     *        CAND_OF_NODE,
     *
     *        KEEP(69), SLAVEF, MSG_SIZE,
     *        NMB_OF_CAND )
      ELSE
        ITEMP=CMUMPS_186(KEEP(69),MEM_DISTRIB,MSG_SIZE)
        NMB_OF_CAND = SLAVEF - 1
      END IF
      NSLAVES_LESS = MAX(ITEMP,1)
      NSLAVES_NODE = CMUMPS_12(KEEP(10), KEEP(48),
     *          KEEP(50),SLAVEF,
     *          NCB, NFRONT, NSLAVES_less, NMB_OF_CAND)
      CALL CMUMPS_441(
     *            KEEP, SLAVEF,
     *            TAB_POS,
     *            NSLAVES_NODE, NFRONT, NCB
     *             )
      IF (force_cand) THEN
        CALL CMUMPS_384(MEM_DISTRIB(0),
     *       CAND_OF_NODE, SLAVEF, NSLAVES_NODE,
     *       SLAVES_LIST)
      ELSE
        CALL CMUMPS_189(MEM_DISTRIB(0),
     *       MSG_SIZE, SLAVES_LIST, NSLAVES_NODE)
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_499
      SUBROUTINE CMUMPS_185( COMM_LD_ARG,
     *                            SYM_ARG, KEEP, 
     *                            MEM_DISTRIB,SIZE_COMM_LD,
     *                            INFO 
     *                            ,MEM_SUB,NB_SUB
#if ! defined(OLD_LOAD_MECHANISM)
     *                            ,FUTURE_NIV2_ARG
#endif
     * )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER COMM_LD_ARG, KEEP(150),
     *        INFO(40), SYM_ARG
      INTEGER SIZE_COMM_LD
      INTEGER, DIMENSION(0:SIZE_COMM_LD - 1) :: MEM_DISTRIB
      INTEGER NB_SUB,MEM_SUB(MAX(1,NB_SUB))
#if ! defined(OLD_LOAD_MECHANISM)
      INTEGER  FUTURE_NIV2_ARG(SIZE_COMM_LD)
#endif
      INTEGER allocok, IERR, I, BUF_LOAD_SIZE
      REAL ZERO
      PARAMETER( ZERO=0.0E0 )
      K35  = KEEP(35)
      K69  = KEEP(69)
      IF ( KEEP(47) .le. 0 .OR. KEEP(47) .gt. 4 ) THEN
        write(*,*) "Internal error in CMUMPS_185"
        CALL CMUMPS_ABORT()
      END IF
      BDC_MEM  = ( KEEP(47) >= 2 )
      BDC_POOL = ( KEEP(47) >= 3 )
      BDC_SBTR = ( KEEP(47) >= 4 )
      K50      = SYM_ARG
      COMM_LD = COMM_LD_ARG
      CALL MPI_COMM_RANK( COMM_LD, MYID, IERR )
      NPROCS = SIZE_COMM_LD
      DM_SUMLU=ZERO
#if ! defined(OLD_LOAD_MECHANISM)
      ALLOCATE(FUTURE_NIV2(NPROCS), stat=allocok)
      IF (allocok > 0 ) THEN
         WRITE(*,*) 'PB allocation in CMUMPS_185'
         INFO(1) = -13
         INFO(2) = NPROCS
         RETURN
      ENDIF
      DO I = 1, NPROCS
        FUTURE_NIV2(I) = FUTURE_NIV2_ARG(I)
      ENDDO
      DELTA_MEM=ZERO
      DELTA_LOAD=ZERO
#endif
      CHECK_MEM=0
      IF (KEEP(47).EQ.4)THEN
         NB_SUBTREES=NB_SUB
         ALLOCATE(MEM_SUBTREE(NB_SUB),stat=allocok)
         DO I=NB_SUB,1,-1
            MEM_SUBTREE(NB_SUB-I+1)=MEM_SUB(I)
         ENDDO
         INSIDE_SUBTREE=0
         NIV1_FLAG=0
         INDICE=NB_SUB
      ENDIF
      IF ( ALLOCATED(LOAD_FLOPS) ) DEALLOCATE( LOAD_FLOPS )
      ALLOCATE( LOAD_FLOPS( 0: NPROCS - 1 ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in CMUMPS_185'
         INFO(1) = -13
         INFO(2) = NPROCS
         RETURN
      END IF
      IF ( ALLOCATED(WLOAD) ) DEALLOCATE( WLOAD )
      ALLOCATE( WLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in CMUMPS_185'
         INFO(1) = -13
         INFO(2) = NPROCS
         RETURN
      END IF
      IF ( ALLOCATED(IDWLOAD) ) DEALLOCATE( IDWLOAD )
      ALLOCATE( IDWLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in CMUMPS_185'
         INFO(1) = -13
         INFO(2) = NPROCS
         RETURN
      END IF
      IF ( BDC_MEM ) THEN
        IF ( ALLOCATED(DM_MEM) ) DEALLOCATE( DM_MEM )
        ALLOCATE( DM_MEM( 0:NPROCS-1 ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in CMUMPS_185'
           INFO(1) = -13
           INFO(2) = NPROCS
           RETURN
        END IF
      END IF
      IF ( BDC_POOL ) THEN
        IF ( ALLOCATED(POOL_MEM) ) DEALLOCATE(POOL_MEM)
        ALLOCATE( POOL_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in CMUMPS_185'
           INFO(1) = -13
           INFO(2) = NPROCS
           RETURN
        END IF
        POOL_MEM = 0.0E0 ! Array instruction
        POOL_LAST_COST_SENT = 0.0E0
      END IF
      IF ( BDC_SBTR ) THEN
        IF ( ALLOCATED(SBTR_MEM) ) DEALLOCATE(SBTR_MEM)
        ALLOCATE( SBTR_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in CMUMPS_185'
           INFO(1) = -13
           INFO(2) = NPROCS
           RETURN
        END IF
        SBTR_MEM = 0.0E0 ! Array instruction
      END IF
      BUF_LOAD_SIZE = KEEP(34) * 2 * ( NPROCS - 1 ) +
     *                NPROCS * ( KEEP(35) + KEEP(34) )
      IF (BDC_MEM) THEN
        BUF_LOAD_SIZE = BUF_LOAD_SIZE + NPROCS * KEEP(35)
      END IF
      LBUF_LOAD_RECV = (BUF_LOAD_SIZE+KEEP(34))/KEEP(34)
      LBUF_LOAD_RECV_BYTES = LBUF_LOAD_RECV * KEEP(34)
      IF ( ALLOCATED(BUF_LOAD_RECV) ) DEALLOCATE(BUF_LOAD_RECV)
      ALLOCATE( BUF_LOAD_RECV( LBUF_LOAD_RECV), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in CMUMPS_185'
        INFO(1) = -13
        INFO(2) = LBUF_LOAD_RECV
        RETURN
      ENDIF 
      BUF_LOAD_SIZE = BUF_LOAD_SIZE * 20
      CALL CMUMPS_54( BUF_LOAD_SIZE, IERR )
      IF ( IERR .LT. 0 ) THEN
         INFO(1) = -13
         INFO(2) = BUF_LOAD_SIZE
         RETURN
      END IF
      DO I = 0, NPROCS - 1
         LOAD_FLOPS( I ) = ZERO
      END DO
#if defined(OLD_LOAD_MECHANISM)
      LOAD_FLOPS( MYID ) = COST_SUBTREE
      LAST_LOAD_SENT = ZERO
#endif
      IF ( BDC_MEM ) THEN
        DO I = 0, NPROCS - 1
          DM_MEM( I )=ZERO
        END DO
#if defined(OLD_LOAD_MECHANISM)
        DM_LAST_MEM_SENT=ZERO
#endif
      ENDIF
      CALL CMUMPS_425(KEEP(69))
      END SUBROUTINE CMUMPS_185
      SUBROUTINE CMUMPS_190( PROCESS_BANDE, INC_LOAD, KEEP )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      REAL INC_LOAD
      INTEGER KEEP(150)
      LOGICAL PROCESS_BANDE
      INTEGER IERR
      REAL ZERO, SEND_MEM, SEND_LOAD
      PARAMETER( ZERO=0.0E0 )
      INTRINSIC MAX, ABS
      IF (INC_LOAD == 0) RETURN
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( PROCESS_BANDE ) RETURN! Info sent by master of node already
#endif
      LOAD_FLOPS( MYID ) = MAX( LOAD_FLOPS( MYID ) + INC_LOAD, ZERO)
#if ! defined(OLD_LOAD_MECHANISM)
      DELTA_LOAD = DELTA_LOAD + INC_LOAD
      IF ( DELTA_LOAD > MIN_DIFF .OR. DELTA_LOAD < -MIN_DIFF ) THEN
         SEND_LOAD = DELTA_LOAD
         IF (BDC_MEM) THEN
           SEND_MEM = DELTA_MEM
         ELSE
           SEND_MEM = ZERO
         END IF
#else
      IF ( ABS( LOAD_FLOPS ( MYID ) -
     *     LAST_LOAD_SENT ) .GT. MIN_DIFF )  THEN
         IERR = 0
         SEND_LOAD  = LOAD_FLOPS( MYID )
         IF ( BDC_MEM ) THEN
           SEND_MEM = DM_MEM(MYID)
         ELSE
           SEND_MEM = ZERO
         END IF
#endif
 111     CONTINUE
           CALL CMUMPS_77( BDC_MEM,
     *        COMM_LD, NPROCS,
     *        SEND_LOAD,
     *        SEND_MEM,
#if ! defined(OLD_LOAD_MECHANISM)
     *        FUTURE_NIV2,
#endif
     *        MYID, IERR )
           IF ( IERR == -1 )THEN
             CALL CMUMPS_467(COMM_LD, KEEP)
             GOTO 111
           ELSE IF ( IERR == -2 ) THEN
             WRITE(*,*) "Internal Error in CMUMPS_190"
             CALL CMUMPS_ABORT()
           ENDIF
         IF ( IERR .EQ. 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
           DELTA_LOAD = ZERO
           IF (BDC_MEM) DELTA_MEM  = ZERO
#else          
           LAST_LOAD_SENT = LOAD_FLOPS( MYID )
           IF ( BDC_MEM ) DM_LAST_MEM_SENT = DM_MEM( MYID )
#endif
         END IF
      ENDIF
      END SUBROUTINE CMUMPS_190
      SUBROUTINE CMUMPS_471(
     *           PROCESS_BANDE, MEM_VALUE, NEW_LU, INC_MEM_ARG,
     *           KEEP)
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MEM_VALUE, INC_MEM_ARG, NEW_LU
      LOGICAL, INTENT(IN) :: PROCESS_BANDE
      INTEGER IERR, KEEP(150)
      REAL ZERO, SEND_MEM, TMP_MEM
      PARAMETER( ZERO=0.0E0 )
      INTRINSIC MAX, ABS
      INTEGER INC_MEM
      INC_MEM = INC_MEM_ARG
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( PROCESS_BANDE .AND. NEW_LU .NE. 0) THEN
        WRITE(*,*) " Internal Error in CMUMPS_471."
        WRITE(*,*) " NEW_LU must be zero if called from PROCESS_BANDE"
        CALL CMUMPS_ABORT()
      ENDIF
#endif
      DM_SUMLU = DM_SUMLU + REAL(NEW_LU)
      CHECK_MEM = CHECK_MEM + INC_MEM
      IF ( MEM_VALUE .NE. CHECK_MEM ) THEN
         WRITE(*,*)
     *   'Problem with increments in CMUMPS_471',
     *   CHECK_MEM, MEM_VALUE, INC_MEM
         CALL CMUMPS_ABORT()
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (PROCESS_BANDE) RETURN
#endif
      IF ( .NOT. BDC_MEM ) RETURN
#if defined(OLD_LOAD_MECHANISM)
      DM_MEM( MYID ) = REAL(MEM_VALUE) - DM_SUMLU
      TMP_MEM = DM_MEM(MYID)
#endif
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( NEW_LU > 0 ) THEN
        INC_MEM = INC_MEM - NEW_LU
      ENDIF
      DM_MEM( MYID ) = DM_MEM(MYID) + INC_MEM
      DELTA_MEM = DELTA_MEM + REAL(INC_MEM)
      IF ( ABS(DELTA_MEM) > DM_THRES_MEM ) THEN
        SEND_MEM = DELTA_MEM
#else
      IF ( ABS( TMP_MEM-DM_LAST_MEM_SENT) >
     *     DM_THRES_MEM )  THEN
         IERR = 0
         SEND_MEM = TMP_MEM
#endif
 111     CONTINUE
         CALL CMUMPS_77( BDC_MEM, COMM_LD,
     *        NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     *        DELTA_LOAD,
#else
     *        LOAD_FLOPS( MYID ), 
#endif
     *        SEND_MEM,
#if ! defined(OLD_LOAD_MECHANISM)
     *        FUTURE_NIV2,
#endif
     $        MYID,IERR )
         IF ( IERR == -1 )THEN
             CALL CMUMPS_467(COMM_LD, KEEP)
             GOTO 111
         ELSE IF ( IERR == -2 ) THEN
             WRITE(*,*) "Internal Error in CMUMPS_471"
             CALL CMUMPS_ABORT()
         ENDIF
         IF ( IERR .EQ. 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
           DELTA_LOAD = ZERO
           DELTA_MEM  = ZERO
#else
           LAST_LOAD_SENT   = LOAD_FLOPS  ( MYID )
           DM_LAST_MEM_SENT = TMP_MEM
#endif
         END IF
      ENDIF
      END SUBROUTINE CMUMPS_471
      INTEGER FUNCTION CMUMPS_186( K69, MEM_DISTRIB,MSG_SIZE )
      IMPLICIT NONE
      INTEGER I, NLESS, K69 
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      REAL LREF
      REAL MSG_SIZE
      NLESS = 0
      DO I=1,NPROCS
         IDWLOAD(I) = I - 1
      ENDDO
      IF(K69 .gt. 1) THEN
         CALL CMUMPS_426(MEM_DISTRIB,MSG_SIZE,IDWLOAD,NPROCS)
      ELSE
         WLOAD(1:NPROCS) = LOAD_FLOPS(0:NPROCS-1)
      ENDIF
      LREF = LOAD_FLOPS(MYID)
      DO I=1, NPROCS
         IF (WLOAD(I).LT.LREF) NLESS=NLESS+1
      ENDDO
      CMUMPS_186 = NLESS
      RETURN
      END FUNCTION CMUMPS_186
      SUBROUTINE CMUMPS_189(MEM_DISTRIB,MSG_SIZE,DEST,
     *     NSLAVES)
      IMPLICIT NONE
      INTEGER NSLAVES
      INTEGER DEST(NSLAVES)
      INTEGER, DIMENSION(0:NPROCS - 1) :: MEM_DISTRIB
      INTEGER I,J,NBDEST
      REAL MSG_SIZE
      IF ( NSLAVES.eq.NPROCS-1 ) THEN
        J = MYID+1
        DO I=1,NSLAVES
           J=J+1
           IF (J.GT.NPROCS) J=1
           DEST(I) = J - 1
        ENDDO
      ELSE
        DO I=1,NPROCS
           IDWLOAD(I) = I - 1
        ENDDO
        CALL CMUMPS_312(NPROCS, WLOAD, IDWLOAD)
         NBDEST = 0
         DO I=1, NSLAVES
            J = IDWLOAD(I)
            IF (J.NE.MYID) THEN
               NBDEST = NBDEST+1
               DEST(NBDEST) = J
            ENDIF
         ENDDO
         IF (NBDEST.NE.NSLAVES) THEN
            DEST(NSLAVES) = IDWLOAD(NSLAVES+1)
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_189
      SUBROUTINE CMUMPS_184( IPROC, LOAD_IPROC )
      IMPLICIT NONE
      INTEGER, intent(in) :: IPROC
      REAL, intent(out) :: LOAD_IPROC
      LOAD_IPROC = LOAD_FLOPS( IPROC )
      RETURN
      END SUBROUTINE CMUMPS_184
      SUBROUTINE CMUMPS_183( INFO1, IERR )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER INFO1, IERR
      DEALLOCATE( LOAD_FLOPS )
      DEALLOCATE( WLOAD )
      DEALLOCATE( IDWLOAD )
#if ! defined(OLD_LOAD_MECHANISM)
      DEALLOCATE(FUTURE_NIV2)
#endif
      IF ( BDC_MEM ) DEALLOCATE( DM_MEM )
      IF ( BDC_POOL) DEALLOCATE( POOL_MEM )
      IF ( BDC_SBTR) DEALLOCATE( SBTR_MEM )
      CALL CMUMPS_58( IERR )
      DEALLOCATE(BUF_LOAD_RECV)
      END SUBROUTINE CMUMPS_183
      SUBROUTINE CMUMPS_467(COMM, KEEP)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MSGTAG, MSGLEN, MSGSOU,COMM
      INTEGER KEEP(150)
      INTEGER STATUS(MPI_STATUS_SIZE)
      LOGICAL FLAG
 10   CONTINUE
      CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *     FLAG, STATUS, IERR )
      IF (FLAG) THEN
        KEEP(65)=KEEP(65)+1
        MSGTAG = STATUS( MPI_TAG )
        MSGSOU = STATUS( MPI_SOURCE )
        IF ( MSGTAG .NE. UPDATE_LOAD) THEN
          write(*,*) "Internal error 1 in CMUMPS_467",
     *    MSGTAG
          CALL CMUMPS_ABORT()
        ENDIF
        CALL MPI_GET_COUNT(STATUS, MPI_PACKED, MSGLEN, IERR)
        IF ( MSGLEN > LBUF_LOAD_RECV_BYTES ) THEN
          write(*,*) "Internal error 2 in CMUMPS_467",
     *    MSGLEN, LBUF_LOAD_RECV_BYTES
          CALL CMUMPS_ABORT()
        ENDIF
        CALL MPI_RECV( BUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES,
     *    MPI_PACKED, MSGSOU, MSGTAG, COMM_LD, STATUS, IERR)
        CALL CMUMPS_187( MSGSOU, BUF_LOAD_RECV,
     *  LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES )
        GOTO 10
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_467
      SUBROUTINE CMUMPS_187( MSGSOU,
     *           BUFR, LBUFR, LBUFR_BYTES )
      IMPLICIT NONE
      INTEGER MSGSOU, LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INCLUDE 'mpif.h'
      INTEGER POSITION, IERR, WHAT, NSLAVES, I
      REAL LOAD_RECEIVED
      INTEGER, POINTER, DIMENSION (:) :: LIST_SLAVES
      REAL, POINTER, DIMENSION (:) :: LOAD_INCR
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     WHAT, 1, MPI_INTEGER,
     *     COMM_LD, IERR )
      IF ( WHAT == 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
#else
#endif
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_REAL,
     *     COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
      LOAD_FLOPS( MSGSOU ) = LOAD_FLOPS(MSGSOU) + LOAD_RECEIVED
#else
      LOAD_FLOPS( MSGSOU ) = LOAD_RECEIVED
#endif
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *       LOAD_RECEIVED, 1, MPI_REAL,
     *       COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
          DM_MEM(MSGSOU)  = DM_MEM(MSGSOU) + LOAD_RECEIVED
#else
          DM_MEM(MSGSOU)  = LOAD_RECEIVED
#endif
        END IF
      ELSEIF ( WHAT == 1) THEN
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     NSLAVES, 1, MPI_INTEGER,
     *     COMM_LD, IERR )
        LIST_SLAVES => IDWLOAD
        LOAD_INCR => WLOAD
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LIST_SLAVES, NSLAVES, MPI_INTEGER,
     *     COMM_LD, IERR)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_INCR, NSLAVES, MPI_REAL,
     *     COMM_LD, IERR)
        DO I = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
          IF ( LIST_SLAVES(I) /= MYID ) THEN
#endif
            LOAD_FLOPS(LIST_SLAVES(I)) =
     *      LOAD_FLOPS(LIST_SLAVES(I)) + 
     *      LOAD_INCR(I)
#if defined(OLD_LOAD_MECHANISM)
          END IF
#endif
        END DO
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_INCR, NSLAVES, MPI_REAL,
     *     COMM_LD, IERR)
          DO I = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
            IF ( LIST_SLAVES(I) /= MYID ) THEN
#endif
              DM_MEM(LIST_SLAVES(I)) = DM_MEM(LIST_SLAVES(I)) + 
     *        LOAD_INCR(I)
#if defined(OLD_LOAD_MECHANISM)
            END IF
#endif
          END DO
        END IF
        NULLIFY( LIST_SLAVES )
        NULLIFY( LOAD_INCR )
      ELSE IF (WHAT == 2 ) THEN
        IF ( .not. BDC_POOL ) THEN
          WRITE(*,*) "Internal error 2 in CMUMPS_187"
          CALL CMUMPS_ABORT()
        END IF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_REAL,
     *     COMM_LD, IERR )
        POOL_MEM(MSGSOU)=LOAD_RECEIVED
      ELSE IF ( WHAT == 3 ) THEN
        IF ( .NOT. BDC_SBTR) THEN
          WRITE(*,*) "Internal error 3 in CMUMPS_187"
          CALL CMUMPS_ABORT()
        ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_REAL,
     *     COMM_LD, IERR )
        SBTR_MEM(MSGSOU)=LOAD_RECEIVED
#if ! defined(OLD_LOAD_MECHANISM)
      ELSE IF (WHAT == 4) THEN
        FUTURE_NIV2(MSGSOU+1)=0
#endif
      ELSE
          WRITE(*,*) "Internal error 1 in CMUMPS_187"
          CALL CMUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE CMUMPS_187
      integer function CMUMPS_409
     *                 (MEM_DISTRIB,CAND,
     *                  K69,
     *                  SLAVEF,MSG_SIZE,
     *                  NMB_OF_CAND )
      implicit none      
      integer, intent(in) :: K69, SLAVEF
      INTEGER, intent(in) :: CAND(SLAVEF+1)
      INTEGER, DIMENSION(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      INTEGER, intent(out) :: NMB_OF_CAND
      integer i,nless
      REAL lref
      REAL MSG_SIZE
      nless = 0 
      NMB_OF_CAND=CAND(SLAVEF+1)
      IF(K69 .gt. 1) THEN
         CALL CMUMPS_426(MEM_DISTRIB,MSG_SIZE,
     *        cand,NMB_OF_CAND)
      ELSE
         do i=1,NMB_OF_CAND
            WLOAD(I)=LOAD_FLOPS(CAND(I))
         end do
      ENDIF
      lref = LOAD_FLOPS(myid)
      do i=1, NMB_OF_CAND
         if (wload(i).lt.lref) nless=nless+1
      end do 
      CMUMPS_409 = nless
      return
      end function CMUMPS_409
      subroutine CMUMPS_384
     *           (MEM_DISTRIB,CAND, 
     *
     *            SLAVEF,
     *            nslaves_inode, dest)
      implicit none
      integer, intent(in) :: nslaves_inode, SLAVEF
      integer, intent(in) :: CAND(SLAVEF+1)
      integer, dimension(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      integer, intent(out) :: dest(nslaves_inode)
      integer i,j,NMB_OF_CAND,deficit,nextproc
      REAL lref
      external CMUMPS_312
      NMB_OF_CAND = CAND(SLAVEF+1)
      if(nslaves_inode.ge.NPROCS .or.
     *   nslaves_inode.gt.NMB_OF_CAND) then
         write(*,*)'Internal error in CMUMPS_384',
     *   nslaves_inode, nprocs, NMB_OF_CAND
         CALL CMUMPS_ABORT()
      end if
      if (nslaves_inode.eq.NPROCS-1) then
         j=myid+1
         do i=1,nslaves_inode
            if(j.ge.nprocs) j=0
            dest(i)=j
            j=j+1
         end do
      else
        do i=1,NMB_OF_CAND
               IDWLOAD(i)=i
        end do
        call CMUMPS_312(NMB_OF_CAND,
     *       wload(1:NMB_OF_CAND),IDWLOAD(1:NMB_OF_CAND) )
        do i=1,nslaves_inode
           dest(i)= cand(IDWLOAD(i))
        end do
      end if   
      return
      end subroutine CMUMPS_384
      SUBROUTINE CMUMPS_425(K69)
      IMPLICIT NONE
      INTEGER K69
      IF (K69 .LE. 4) THEN
         ALPHA = 0.0
         BETA = 0.0
         RETURN
      ENDIF
      IF (K69 .EQ. 5) THEN
         ALPHA = 0.5
         BETA = 50000.0
         RETURN
      ENDIF
      IF (K69 .EQ. 6) THEN
         ALPHA = 0.5
         BETA = 100000.0
         RETURN
      ENDIF
      IF (K69 .EQ. 7) THEN
         ALPHA = 0.5
         BETA = 150000.0
         RETURN
      ENDIF
      IF (K69 .EQ. 8) THEN
         ALPHA = 1.0
         BETA = 50000.0
         RETURN
      ENDIF
      IF (K69 .EQ. 9) THEN
         ALPHA = 1.0
         BETA = 100000.0
         RETURN
      ENDIF
      IF (K69 .EQ. 10) THEN
         ALPHA = 1.0
         BETA = 150000.0
         RETURN
      ENDIF
      IF (K69 .EQ. 11) THEN
         ALPHA = 1.5
         BETA = 50000.0
         RETURN
      ENDIF
      IF (K69 .EQ. 12) THEN
         ALPHA = 1.5
         BETA = 100000.0
         RETURN
      ENDIF
      ALPHA = 1.5
      BETA = 150000.0
      RETURN
      END SUBROUTINE CMUMPS_425
      SUBROUTINE CMUMPS_426(MEM_DISTRIB,MSG_SIZE,ARRAY_ADM,LEN)
      IMPLICIT NONE
      INTEGER I,J,LEN,COUNTMORE,COUNTLOCLESS
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      REAL MSG_SIZE,FORBIGMSG
      INTEGER ARRAY_ADM(LEN)
      FORBIGMSG = 1.0E0
      IF (K69 .lt.2) THEN
         DO I = 1,LEN
            WLOAD(I) = LOAD_FLOPS(ARRAY_ADM(I))
         ENDDO
         RETURN
      ENDIF
      IF((MSG_SIZE * K35 ) .gt. 3200000.0) THEN
         FORBIGMSG = 2.0E0
      ENDIF
      IF (K69 .le. 4) THEN
         DO I = 1,LEN
            IF ((MEM_DISTRIB(ARRAY_ADM(I)) .EQ. 1) .AND.
     *      LOAD_FLOPS(ARRAY_ADM(I)) .LT. LOAD_FLOPS(MYID) ) THEN
              WLOAD(I) = LOAD_FLOPS(ARRAY_ADM(I))/LOAD_FLOPS(MYID)
            ELSE
              IF ( MEM_DISTRIB(ARRAY_ADM(I)) == 1 ) THEN
                WLOAD(I) = LOAD_FLOPS(ARRAY_ADM(I))
              ELSE ! MEM_DISTRIB equals 4 and FORBIGMSG = 4
                WLOAD(I) = LOAD_FLOPS(ARRAY_ADM(I)) *
     *              MEM_DISTRIB(ARRAY_ADM(I))
     *              * FORBIGMSG
     *              + REAL(2)
              ENDIF
            ENDIF
         ENDDO
         RETURN
      ENDIF
      DO I = 1,LEN
         IF ((MEM_DISTRIB(ARRAY_ADM(I)) .EQ. 1) .AND.
     *        LOAD_FLOPS(ARRAY_ADM(I)) .LT. LOAD_FLOPS(MYID) ) THEN
            WLOAD(I) = LOAD_FLOPS(ARRAY_ADM(I)) /  LOAD_FLOPS(MYID)
         ELSE
            IF(MEM_DISTRIB(ARRAY_ADM(I)) .EQ. 1) THEN
               WLOAD(I) = LOAD_FLOPS(ARRAY_ADM(I))
            ELSE     
               WLOAD(I) = (LOAD_FLOPS(ARRAY_ADM(I)) +
     *              ALPHA * MSG_SIZE * K35  +
     *              BETA) * FORBIGMSG
            ENDIF
         ENDIF
      ENDDO  
      RETURN
      END SUBROUTINE CMUMPS_426
      SUBROUTINE CMUMPS_461(MYID, SLAVEF, COMM,
     *     TAB_POS, NASS, KEEP, LIST_SLAVES, NSLAVES)
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: MYID, SLAVEF, COMM, NASS, NSLAVES
      INTEGER, INTENT (IN) :: TAB_POS(SLAVEF+2)
      INTEGER, INTENT (IN) :: LIST_SLAVES( NSLAVES )
      INTEGER KEEP(150)
      INTEGER NCB, NFRONT, NBROWS_SLAVE
      INTEGER I, IERR
      REAL MEM_INCREMENT ( NSLAVES )
      REAL FLOPS_INCREMENT( NSLAVES )
#if ! defined(OLD_LOAD_MECHANISM)
      FUTURE_NIV2(MYID+1) = FUTURE_NIV2(MYID+1) - 1
      IF ( FUTURE_NIV2(MYID+1) < 0 ) THEN
        WRITE(*,*) "Internal error in CMUMPS_461"
        CALL CMUMPS_ABORT()
      ENDIF
      IF ( FUTURE_NIV2(MYID + 1) == 0 ) THEN
 112    CONTINUE
        CALL CMUMPS_502(COMM,MYID,SLAVEF,IERR)
        IF (IERR == -1 ) THEN
          CALL CMUMPS_467(COMM_LD, KEEP)
          GOTO 112
        ELSE IF ( IERR == -2 ) THEN
          WRITE(*,*) "Internal Error in CMUMPS_461"
          CALL CMUMPS_ABORT()
        ENDIF
      ENDIF
#endif
      IF ( NSLAVES /= TAB_POS(SLAVEF + 2) ) THEN
        write(*,*) "Error 1 in CMUMPS_461",
     +             NSLAVES, TAB_POS(SLAVEF+2)
        stop
      ENDIF
      NCB = TAB_POS(NSLAVES+1) - 1
      NFRONT = NCB + NASS
      DO I = 1, NSLAVES
        NBROWS_SLAVE = TAB_POS(I+1) - TAB_POS(I)
        IF ( BDC_MEM ) THEN
          IF ( KEEP(50) == 0 ) THEN
            MEM_INCREMENT( I ) = REAL(NBROWS_SLAVE) * REAL(NFRONT)
          ELSE
            MEM_INCREMENT( I ) = REAL(NBROWS_SLAVE) *
     *                          REAL( NASS + TAB_POS(I+1) - 1 )
          END IF
        END IF
        IF ( KEEP(50) == 0 ) THEN
          FLOPS_INCREMENT( I ) = REAL(NBROWS_SLAVE * NASS ) +
     *    REAL(NBROWS_SLAVE * NASS ) * REAL(2*NFRONT-NASS-1)
        ELSE
          FLOPS_INCREMENT( I ) = REAL(NBROWS_SLAVE * NASS ) *
     *      REAL( 2 * ( NASS + TAB_POS(I+1) - 1 ) 
     *      - NBROWS_SLAVE - NASS + 1 )
        ENDIF
      END DO
 111  CONTINUE
      CALL CMUMPS_460(BDC_MEM, COMM, MYID, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     *     FUTURE_NIV2,
#endif
     *     NSLAVES, LIST_SLAVES,
     *     MEM_INCREMENT,
     *     FLOPS_INCREMENT, IERR)
        IF ( IERR == -1 ) THEN
          CALL CMUMPS_467(COMM_LD, KEEP)
          GOTO 111
        ELSE IF ( IERR == -2 ) THEN
          WRITE(*,*) "Internal Error in CMUMPS_461"
          CALL CMUMPS_ABORT()
        ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (FUTURE_NIV2(MYID+1) .NE. 0) THEN
#endif
        DO I = 1, NSLAVES
          LOAD_FLOPS(LIST_SLAVES(I)) = LOAD_FLOPS(LIST_SLAVES(I))
     *       +  FLOPS_INCREMENT(I)
          IF ( BDC_MEM ) THEN
            DM_MEM(LIST_SLAVES(I)) = DM_MEM(LIST_SLAVES(I))
     *       +  MEM_INCREMENT(I)
          END IF
        ENDDO
#if ! defined(OLD_LOAD_MECHANISM)
      ENDIF
#endif
      RETURN
      END SUBROUTINE CMUMPS_461
      SUBROUTINE CMUMPS_500(
     $     POOL, LPOOL,
     *     PROCNODE, KEEP, SLAVEF, COMM, MYID, STEP, N,
     *     ND, FILS )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER LPOOL, SLAVEF, COMM, MYID 
      INTEGER N, KEEP(150)
      INTEGER POOL( LPOOL ), PROCNODE( KEEP(28) ), STEP( N )
      INTEGER ND( KEEP(28) ), FILS( N )
      INTEGER I, INODE, NELIM, NFR, LEVEL, IERR, WHAT
      REAL COST
      INTEGER NBINSUBTREE,NBTOP,INSUBTREE
      LOGICAL CMUMPS_167
      INTEGER CMUMPS_330
      EXTERNAL CMUMPS_167, CMUMPS_330
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF(KEEP(76).EQ.0)THEN
         IF(NBTOP.NE.0)THEN
            DO I = LPOOL-NBTOP-2, MIN(LPOOL-3,LPOOL-NBTOP-2+3)
               INODE = POOL( I )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=0.0E0 
            GOTO 30
         ELSE
            DO I = NBINSUBTREE, MAX(1,NBINSUBTREE-3), -1
               INODE = POOL( I )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=0.0E0 
            GOTO 30
         ENDIF
      ELSE
         IF(KEEP(76).EQ.1)THEN
            IF(INSUBTREE.EQ.1)THEN
               DO I = NBINSUBTREE, MAX(1,NBINSUBTREE-3), -1
                  INODE = POOL( I )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=0.0E0 
               GOTO 30
            ELSE
               DO I = LPOOL-NBTOP-2, MIN(LPOOL-3,LPOOL-NBTOP-2+3)
                  INODE = POOL( I )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=0.0E0 
               GOTO 30
            ENDIF
         ELSE
            WRITE(*,*)
     *      'Internal error: Unknown pool management strategy'
            CALL CMUMPS_ABORT()
         ENDIF
      ENDIF
 20   CONTINUE
        I = INODE
        NELIM = 0
 10     CONTINUE
        IF ( I > 0 ) THEN
          NELIM = NELIM + 1
          I = FILS(I)
          GOTO 10
        ENDIF
        NFR = ND( STEP(INODE) )
        LEVEL = CMUMPS_330( STEP(INODE), PROCNODE, SLAVEF )
        IF (LEVEL .EQ. 1) THEN
          COST = REAL( NFR * NFR )
        ELSE
          IF ( KEEP(50) == 0 ) THEN
            COST = REAL( NFR * NELIM )
          ELSE
            COST = REAL( NELIM * NELIM )
          ENDIF
        ENDIF
 30   CONTINUE
      IF ( POOL_LAST_COST_SENT /= COST ) THEN
        WHAT = 2
 111    CONTINUE
        CALL CMUMPS_498( WHAT,
     *         COMM, SLAVEF,
     *         COST, MYID, iERR  )
        POOL_LAST_COST_SENT = COST
        IF ( IERR == -1 )THEN
          CALL CMUMPS_467(COMM_LD, KEEP)
          GOTO 111
        ELSE IF ( IERR == -2 ) THEN
          WRITE(*,*) "Internal Error in CMUMPS_500"
          CALL CMUMPS_ABORT()
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_500
      SUBROUTINE CMUMPS_501(
     $     POOL,LPOOL,MYID,slavef,COMM,KEEP)
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER LPOOL,MYID,SLAVEF,COMM
      INTEGER POOL(LPOOL),KEEP(150)
      INTEGER WHAT,IERR
      REAL ZERO
      PARAMETER(ZERO = 0.0E0)
      IF(POOL(LPOOL-2).EQ.1)THEN
         IF(INSIDE_SUBTREE.EQ.0)THEN
            WHAT = 3
            IF(MEM_SUBTREE(INDICE).GE.DM_THRES_MEM)THEN
 111          CONTINUE
              CALL CMUMPS_498(
     *           WHAT, COMM, SLAVEF,
     *           REAL(MEM_SUBTREE(INDICE)), MYID, iERR  )
              IF ( IERR == -1 )THEN
               CALL CMUMPS_467(COMM_LD, KEEP)
               GOTO 111
              ELSE IF ( IERR == -2 ) THEN
               WRITE(*,*)
     *     "Internal Error 1 in CMUMPS_501"
               CALL CMUMPS_ABORT()
              ENDIF
            ENDIF
            INDICE=INDICE-1
            INSIDE_SUBTREE=1
         ENDIF
      ELSE
         IF(POOL(LPOOL-2).EQ.0)THEN
            IF(INSIDE_SUBTREE.EQ.1)THEN
               WHAT = 3
 112           CONTINUE
               CALL CMUMPS_498(
     *              WHAT, COMM, SLAVEF,
     *              ZERO, MYID, iERR  )
               IF ( IERR == -1 )THEN
               CALL CMUMPS_467(COMM_LD, KEEP)
               GOTO 112
               ELSE IF ( IERR == -2 ) THEN
               WRITE(*,*)
     *     "Internal Error 2 in CMUMPS_501"
               CALL CMUMPS_ABORT()
               ENDIF
               INSIDE_SUBTREE=0
            ENDIF
         ELSE
            WRITE(*,*)
     *     "Internal Error 3 in CMUMPS_501"
            CALL CMUMPS_ABORT()
         ENDIF 
      ENDIF
      END SUBROUTINE CMUMPS_501
      SUBROUTINE CMUMPS_504
     &      (SLAVEF,KEEP,PROCS,MEM_DISTRIB,NCB,NFRONT,
     $       NSLAVES_NODE,TAB_POS,
     $       SLAVES_LIST,SIZE_SLAVES_LIST,MYID)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(150),SIZE_SLAVES_LIST
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER NUMBER_OF_PROCS,K47,K10, K48, K50
      INTEGER J, NSLAVES, NSLAVES_REAL, ISLAVE
      INTEGER KMIN, KMAX, BLSIZE, NPIV,
     *        NSLAVES_ref, KMINSURF, NSLAVES_max
      REAL WK_MASTER, WK_SLAVE
      REAL  CMUMPS_45
      EXTERNAL CMUMPS_45
      INTEGER OTHERS,CHOSEN,SMALL_SET,ACC
      REAL SOMME,AFFECTED
      INTEGER ADDITIONNAL_ROWS,BLOCK_SIZE,I,X,TOTAL_MEM,REF,POS
      LOGICAL force_cand
      REAL TEMP(SLAVEF),SUB(SLAVEF),PEAK
      INTEGER TEMP_ID(SLAVEF),NB_ROWS(SLAVEF),SUB_ID(SLAVEF)
      EXTERNAL MPI_WTIME
      REAL MPI_WTIME
      IF (KEEP(10) .GT. 0) THEN
      write(*,*)
     * "Internal Error 1 in CMUMPS_504"
      CALL CMUMPS_ABORT()
      ENDIF
      K10=NFRONT*NCB
      K50=KEEP(50)
      K48=KEEP(48)
      K47=KEEP(47)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        force_cand = .FALSE.
      ELSE
        force_cand = (mod(KEEP(24),2).eq.0)
      END IF
      IF(K48.NE.4)THEN
         WRITE(*,*)'CMUMPS_COMPUTE_PARTI_ACTV_MEM_K10
     $      should be called with KEEP(48) different from 4'
         CALL CMUMPS_ABORT()
      ENDIF
         KMIN=1
         KMAX=K10/NCB
         IF(FORCE_CAND)THEN
            DO I=1,PROCS(SLAVEF+1)
               WLOAD(I)=DM_MEM(PROCS(I))
               IDWLOAD(I)=PROCS(I)
            ENDDO
            NUMBER_OF_PROCS=PROCS(SLAVEF+1)
            OTHERS=NUMBER_OF_PROCS
         ELSE
            NUMBER_OF_PROCS=SLAVEF
            WLOAD(1:SLAVEF) = DM_MEM(0:NUMBER_OF_PROCS-1)
            DO I=1,NUMBER_OF_PROCS
               IDWLOAD(I) = I - 1
            ENDDO
            OTHERS=NUMBER_OF_PROCS-1
         ENDIF
         NB_ROWS=0
         CALL CMUMPS_312(NUMBER_OF_PROCS, WLOAD, IDWLOAD)
         TOTAL_MEM=NCB*NFRONT
         SOMME=0
         J=1
         DO I=1,NUMBER_OF_PROCS
            IF((IDWLOAD(I).NE.MYID))THEN
               TEMP_ID(J)=IDWLOAD(I)
               TEMP(J)=WLOAD(I)
                IF(KEEP(47).GE.3) THEN
                   TEMP(J)=TEMP(J)+POOL_MEM(TEMP_ID(J))
                ENDIF
                IF(KEEP(47).EQ.4) THEN
                   SUB(J)=SBTR_MEM(IDWLOAD(I))
                   SUB_ID(J)=IDWLOAD(I)
                ENDIF
               J=J+1
            ENDIF
         ENDDO               
         NUMBER_OF_PROCS=J-1
         CALL CMUMPS_312(NUMBER_OF_PROCS, TEMP, TEMP_ID)
         PEAK=TEMP(NUMBER_OF_PROCS)
         SOMME=0
         DO I=1,NUMBER_OF_PROCS
           SOMME=SOMME+TEMP(OTHERS)-TEMP(J)
         ENDDO               
         IF(SOMME.LE.REAL(TOTAL_MEM)) THEN
            GOTO 096
         ENDIF
         IF(KEEP(47).EQ.4)THEN
             CALL CMUMPS_312(NUMBER_OF_PROCS, SUB, SUB_ID)
             DO J=1,NUMBER_OF_PROCS
             ENDDO
             OTHERS=0
             J=1
             DO I=1,NUMBER_OF_PROCS
                IDWLOAD(I)=SUB_ID(I)
                WLOAD(I)=DM_MEM(SUB_ID(I))
                IF(SUB(I).EQ.0)THEN
                   TEMP_ID(J)=SUB_ID(I)
                   TEMP(J)=WLOAD(I)
                   J=J+1
                   OTHERS=OTHERS+1
                ENDIF
             ENDDO
             ACC=OTHERS
          ENDIF
  097     CONTINUE
          IF(KEEP(47).EQ.4)THEN
             CALL CMUMPS_312(OTHERS,TEMP,TEMP_ID)
             SOMME=0
             X=0
             DO I=1,ACC
                SOMME=SOMME+TEMP(OTHERS)-TEMP(I)
                IF((TEMP(OTHERS)-TEMP(I)).GT.K10) X=1
             ENDDO
             IF((SOMME.LT.TOTAL_MEM).OR.(OTHERS.EQ.0).OR.(X.EQ.1))THEN
                IF(OTHERS .EQ. NUMBER_OF_PROCS) GOTO 096
                ACC=OTHERS+1
                DO I=1,ACC
                   TEMP_ID(I)=SUB_ID(I)
                   TEMP(I)=WLOAD(I)
                ENDDO
                IF(OTHERS.LT.NUMBER_OF_PROCS) THEN
                   OTHERS=ACC
                   GOTO 097
                ELSE
                   CALL CMUMPS_312(OTHERS,TEMP,TEMP_ID)
                   GOTO 096
                ENDIF
             ENDIF
             OTHERS=ACC
          ENDIF
 096     CONTINUE
         SOMME=0
         DO I=1,OTHERS
            SOMME=SOMME+TEMP(OTHERS)-TEMP(I)
         ENDDO
         IF(REAL(TOTAL_MEM).GE.SOMME) THEN
            AFFECTED=0
            CHOSEN=0
            ACC=0
            DO I=1,OTHERS
               IF(K50.EQ.0)THEN
                  IF((TEMP(OTHERS)-TEMP(I)).GT.K10)THEN
                     AFFECTED=K10
                  ELSE
                     AFFECTED=TEMP(OTHERS)-TEMP(I)
                  ENDIF
                  x=INT(AFFECTED/REAL(NFRONT))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
               ENDIF
               IF(K50.EQ.1)THEN
                  IF((TEMP(OTHERS)-TEMP(I)).GT.K10)THEN
                     AFFECTED=K10
                  ELSE
                     AFFECTED=TEMP(OTHERS)-TEMP(I)
                  ENDIF
                  x=INT((-(NFRONT-NCB+ACC)
     $                 +SQRT(REAL(((NFRONT-NCB+ACC)*
     $                 (NFRONT-NCB+ACC))+REAL(4)*
     *                 (AFFECTED))))/
     $                 REAL(2))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
                  IF(X.LE.0) THEN
                    WRITE(*,*)'probleme de nb lignes'
                    CALL CMUMPS_ABORT()
                  ENDIF
               ENDIF
               NB_ROWS(I)=X
               CHOSEN=CHOSEN+1
               ACC=ACC+X
               IF(NCB-ACC.LT.KMIN) GOTO 111
               IF(NCB.EQ.ACC) GOTO 111
               ENDDO
 111           CONTINUE
               IF((ACC.GT.NCB))THEN
                  X=0
                  DO I=1,OTHERS
                     X=X+NB_ROWS(I)
                  ENDDO
                  WRITE(*,*)'NCB=',NCB,',SOMME=',X
                  WRITE(*,*)'probleme dans get_nslaves_mod'
                  CALL CMUMPS_ABORT()
               ENDIF
               WRITE(*,*)'CHOSEN=',CHOSEN
               IF((NCB.NE.ACC))THEN
                  IF(K50.EQ.1)THEN
                     IF(CHOSEN.NE.0)THEN
                        WRITE(*,*)'ici!!!!!'
                        ADDITIONNAL_ROWS=NCB-ACC
                        NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+ADDITIONNAL_ROWS
                     ELSE
                        TOTAL_MEM=TOTAL_MEM/REAL(NUMBER_OF_PROCS)
                        AFFECTED=0
                        CHOSEN=0
                        ACC=0
                        DO I=1,OTHERS
                           x=INT((-(NFRONT-NCB+ACC)
     $                          +SQRT(REAL(((NFRONT-NCB+ACC)*
     $                          (NFRONT-NCB+ACC))+REAL(4)*
     *                          (AFFECTED))))/
     $                          REAL(2))
                           IF((ACC+X).GT.NCB) X=NCB-ACC
                           NB_ROWS(I)=X
                           CHOSEN=CHOSEN+1
                           ACC=ACC+X
                           IF(NCB-ACC.LT.KMIN) GOTO 002
                           IF(NCB.EQ.ACC) GOTO 002
                        ENDDO
 002                    CONTINUE
                        IF(ACC.LT.NCB)THEN
                           NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+(NCB-ACC)
                        ENDIF
                     ENDIF
                     GOTO 333
                  ENDIF
                  ADDITIONNAL_ROWS=NCB-ACC
                  DO I=CHOSEN,1,-1
                     IF(INT(REAL(ADDITIONNAL_ROWS)/REAL(I)).NE.0)THEN
                        GOTO 222 
                     ENDIF
                  ENDDO
 222              CONTINUE
                  X=INT(REAL(ADDITIONNAL_ROWS)/REAL(I))
                  DO J=1,I
                     NB_ROWS(J)=NB_ROWS(J)+X
                     ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                  ENDDO 
                  IF(ADDITIONNAL_ROWS.NE.0) THEN
                     NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                  ENDIF    
               ENDIF
 333           CONTINUE
            ELSE
               DO I=OTHERS,1,-1
                  SOMME=0
                  DO J=1,I
                     SOMME=SOMME+TEMP(J)
                  ENDDO
                  SOMME=(REAL(I)*TEMP(I))-SOMME
                  IF(TOTAL_MEM.GE.INT(SOMME)) GOTO 444
               ENDDO
 444           CONTINUE
               REF=I
               DO J=1,I
                  IF(TEMP(J).EQ.TEMP(I)) THEN
                     SMALL_SET=J
                     GOTO 123
                  ENDIF
               ENDDO
 123           CONTINUE
               IF(I.EQ.1)THEN 
                  NB_ROWS(I)=NCB
                  CHOSEN=1
                  GOTO 666
               ENDIF
 323           CONTINUE
               AFFECTED=0
               CHOSEN=0
               ACC=0
               DO I=1,SMALL_SET
                  IF(K50.EQ.0)THEN
                     IF((TEMP(SMALL_SET)-TEMP(I)).GT.K10)THEN
                        AFFECTED=K10
                     ELSE
                        AFFECTED=TEMP(SMALL_SET)-TEMP(I)
                     ENDIF
                     x=INT(AFFECTED/REAL(NFRONT))
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  IF(K50.EQ.1)THEN
                     IF((TEMP(SMALL_SET)-TEMP(I)).GT.K10)THEN
                        AFFECTED=K10
                     ELSE
                        AFFECTED=TEMP(SMALL_SET)-TEMP(I)
                     ENDIF
                      x=INT((-(NFRONT-NCB+ACC)
     $                  +SQRT(REAL(((NFRONT-NCB+ACC)*
     $                  (NFRONT-NCB+ACC))+REAL(4)*
     *                  (AFFECTED))))/
     $                  REAL(2))
                     IF(X.LT.0)THEN
                        WRITE(*,*)
     *              'Internal error 1 in CMUMPS_504'
                        CALL CMUMPS_ABORT()
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  NB_ROWS(I)=X
                  ACC=ACC+X
                  CHOSEN=CHOSEN+1
                  IF(NCB-ACC.LT.KMIN) GOTO 888
                  IF(NCB.EQ.ACC) GOTO 888
                  IF(ACC.GT.NCB) THEN
                    WRITE(*,*)
     *              'Internal error 2 in CMUMPS_504'
                    CALL CMUMPS_ABORT()
                  ENDIF
               ENDDO
 888           CONTINUE
               SOMME=0
               X=NFRONT-NCB
               IF((ACC.GT.NCB))THEN
                  WRITE(*,*)
     *            'Internal error 3 in CMUMPS_504'
                  CALL CMUMPS_ABORT()
               ENDIF
               IF((ACC.LT.NCB))THEN
                  IF(K50.EQ.1)THEN
                     IF(SMALL_SET.LT.OTHERS)THEN
                       SMALL_SET=REF+1
                       REF=SMALL_SET
                       GOTO 323
                     ELSE
                       NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+NCB-ACC
                       GOTO 666
                     ENDIF
                 ENDIF
                 ADDITIONNAL_ROWS=NCB-ACC
                 DO I=CHOSEN,1,-1
                    IF(INT(REAL(ADDITIONNAL_ROWS)/REAL(I)).NE.0)THEN
                       GOTO 555
                    ENDIF
                 ENDDO
 555             CONTINUE
                 X=INT(REAL(ADDITIONNAL_ROWS)/REAL(I))
                 DO J=1,I
                    IF(NB_ROWS(J)+X.GT.K10/NCB)THEN
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                      ((K10/NCB)-NB_ROWS(J))
                       NB_ROWS(J)=(K10/NCB)
                    ELSE
                       NB_ROWS(J)=NB_ROWS(J)+X
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                    ENDIF
                    IF(((TEMP(J)+REAL(NB_ROWS(J)*NFRONT))
     *               .GT. TEMP(OTHERS))
     *               .AND.(SMALL_SET.LT.OTHERS))THEN
                       SMALL_SET=SMALL_SET+1
                       GOTO 323
                    ENDIF
                 ENDDO
                 SOMME=0
                 DO J=1,CHOSEN
                   SOMME=SOMME+NB_ROWS(J)
                 ENDDO
                 IF(ADDITIONNAL_ROWS.NE.0) THEN
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                        .AND.(I.LE.NUMBER_OF_PROCS))
                          IF(NB_ROWS(I)+ADDITIONNAL_ROWS.LT.K10/NCB)THEN
                             NB_ROWS(I)=NB_ROWS(I)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ELSE
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-(K10/NCB
     $                            -NB_ROWS(I))
                             NB_ROWS(I)=K10/NCB
                          ENDIF
                          I=I+1
                       ENDDO
                       CHOSEN=I-1
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                         DO I=CHOSEN,1,-1
                          IF(INT(REAL(ADDITIONNAL_ROWS)/REAL(I))
     *                          .NE.0)THEN
                             GOTO 372 
                          ENDIF
                         ENDDO
 372                     CONTINUE
                         X=INT(REAL(ADDITIONNAL_ROWS)/REAL(I))
                         DO J=1,I
                           NB_ROWS(J)=NB_ROWS(J)+X
                           ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                         ENDDO 
                         IF(ADDITIONNAL_ROWS.NE.0) THEN
                           NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                         ENDIF    
                    ENDIF
                    IF(((TEMP(I)+REAL(NB_ROWS(I)*NFRONT))
     *               .GT.TEMP(OTHERS)).AND.
     *               (SMALL_SET.LT.OTHERS)) THEN
                       SMALL_SET=SMALL_SET+1
                       GOTO 323
                    ENDIF 
                 ENDIF
              ENDIF
 666          CONTINUE
              SOMME=0
              X=0
              DO I=1,CHOSEN
                 IF(K50.EQ.1) THEN
                    IF((TEMP(I)+REAL(NB_ROWS(I)
     *                   *(X+NB_ROWS(I)+NFRONT-NCB)))
     *                   .GT.PEAK)THEN
                       IF((K47.EQ.4).AND.
     *                 (CHOSEN.LT.NUMBER_OF_PROCS))THEN
                         ACC=OTHERS+1
                         OTHERS=OTHERS+1
                         IF(OTHERS.GT.NUMBER_OF_PROCS)THEN
                         ENDIF
                         GOTO 097
                       ENDIF
                       SMALL_SET=SMALL_SET+1
                       CALL CMUMPS_ABORT()
                    ENDIF
                 ENDIF
                 IF(K50.EQ.0) THEN
                    IF((TEMP(I)+REAL(NB_ROWS(I)*NFRONT))
     *                   .GT.PEAK)THEN
                       IF((K47.EQ.4).AND.
     *                    (CHOSEN.LT.NUMBER_OF_PROCS))THEN
                         ACC=OTHERS+1
                         OTHERS=OTHERS+1
                         CHOSEN=0
                         NB_ROWS=0
                         IF(OTHERS.GT.NUMBER_OF_PROCS)THEN
                         ENDIF
                         GOTO 097
                       ENDIF
                       SMALL_SET=SMALL_SET+1
                       CALL CMUMPS_ABORT()
                    ENDIF
                 ENDIF
                 X=X+NB_ROWS(I)
                 SOMME=SOMME+ REAL(NB_ROWS(I))
              ENDDO
           ENDIF
           J=CHOSEN
           X=0
           DO I=J,1,-1
             IF(NB_ROWS(I).EQ.0)THEN
                IF(X.EQ.1)THEN
                  CALL CMUMPS_ABORT()
                ENDIF
                CHOSEN=CHOSEN-1
             ELSE 
                  IF(NB_ROWS(I).GT.0)THEN
                    X=1
                  ELSE
                    WRITE(*,*)
     *            'Internal error 6 in CMUMPS_504'
                    CALL CMUMPS_ABORT()
                  ENDIF
             ENDIF
          ENDDO
           NSLAVES_NODE=CHOSEN
           TAB_POS(NSLAVES_NODE+1)= NCB+1
           TAB_POS(SLAVEF+2) = CHOSEN
           POS=1
           DO I=1,CHOSEN
              SLAVES_LIST(I)=TEMP_ID(I)
              TAB_POS(I)=POS
              POS=POS+NB_ROWS(I) 
              IF(NB_ROWS(I).LE.0)THEN
                WRITE(*,*)
     *          'Internal error 7 in CMUMPS_504'
                 CALL CMUMPS_ABORT()
              ENDIF
           ENDDO
           IF(POS.NE.(NCB+1))THEN
              WRITE(*,*)
     *        'Internal error 8 in CMUMPS_504'
             CALL CMUMPS_ABORT()
           ENDIF
      END SUBROUTINE CMUMPS_504
      END MODULE CMUMPS_LOAD
