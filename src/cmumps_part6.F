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
      RECURSIVE SUBROUTINE CMUMPS_41(
     *     BLOQ, FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, III, IISSAR,
     *     NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *     IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *     PLEFTW, KEEP,
     *     PTRIST, IW, LIW, A, LA, W2, MYLEAFE, RHSSOL, RHS )
      IMPLICIT NONE
      LOGICAL BLOQ, FLAG
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      COMPLEX W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, III, IISSAR, LPOOL
      INTEGER IPOOL( LPOOL ), POOLSS( N )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 150 )
      INTEGER PROCNODE_STEPS( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER LIW, LA
      INTEGER PTRIST( KEEP(28) ), IW( LIW )
      COMPLEX A( LA ), W2( N ), RHSSOL( N ), RHS(N)
      INTEGER MYLEAFE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGSOU, MSGTAG, MSGLEN
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     *                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
      END IF
      IF (FLAG) THEN
         MSGSOU=STATUS(MPI_SOURCE)
         MSGTAG=STATUS(MPI_TAG)
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL CMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV(BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     *                   MSGTAG, COMM, STATUS, IERR)
           CALL CMUMPS_42( MSGTAG, MSGSOU,
     *                BUFR, LBUFR, LBUFR_BYTES,
     *                MYID, SLAVEF, COMM,
     *                N, IWCB, LIWW, POSIWCB,
     *                W, LWC, POSWCB,
     *                IIPOOL, III, IISSAR,
     *                NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *                IPOOL, LPOOL, STEP,
     *                FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,
     *                PTRIST, IW, LIW, A, LA, W2, MYLEAFE, 
     *                RHSSOL, RHS )
         END IF
      END IF
      RETURN
      END SUBROUTINE CMUMPS_41
      RECURSIVE SUBROUTINE CMUMPS_42(
     *                MSGTAG, MSGSOU,
     *                BUFR, LBUFR, LBUFR_BYTES,
     *                MYID, SLAVEF, COMM,
     *                N, IWCB, LIWW, POSIWCB,
     *                W, LWC, POSWCB,
     *                IIPOOL, III, IISSAR,
     *                NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *                IPOOL, LPOOL, STEP,
     *                FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,
     *                PTRIST, IW, LIW, A, LA, W2, MYLEAFE, 
     *                RHSSOL, RHS )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MSGTAG, MSGSOU
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      COMPLEX W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, III, IISSAR, LPOOL
      INTEGER IPOOL( LPOOL ), POOLSS( N )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 150 )
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER FRERE(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LIW, LA
      INTEGER IW( LIW ), PTRIST( KEEP(28) )
      COMPLEX A( LA ), W2( N ), RHSSOL( N ), RHS(N)
      INTEGER MYLEAFE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER POSITION, IF, INODE, IERR, LONG, DUMMY(1)
      INTEGER P_UPDATE, P_SOL_MAS, LIELL
      INTEGER NPIV, NPIV_CHECK, NROW_L, APOS, IPOS, NROW_RECU
      INTEGER I, JJ, IN, PROCDEST
      LOGICAL FLAG
      COMPLEX ZERO, ALPHA
      PARAMETER( ZERO = 0.0E0, ALPHA = -1.0E0 )
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      IF (MSGTAG .EQ. FEUILLE) THEN
          NBFINF = NBFINF - 1
      ELSE IF (MSGTAG .EQ. NOEUD) THEN
          POSITION = 0
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        INODE, 1, MPI_INTEGER,
     *        COMM, IERR)
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        LONG, 1, MPI_INTEGER,
     *        COMM, IERR)
          IF (   POSIWCB - LONG - 2 .LT. 0
     *      .OR. POSWCB - PLEFTW + 1 .LT. LONG ) THEN
            CALL CMUMPS_95(N, KEEP(28), IWCB, LIWW, W, LWC,
     *      POSWCB, POSIWCB, PTRICB, PTRACB)
            IF ((POSIWCB - LONG - 2 ) .LT. 0) THEN
              INFO(1)=-14
              INFO(2)=-POSIWCB + LONG + 2
              GOTO 260
            END IF
            IF ( POSWCB - PLEFTW + 1 .LT. LONG ) THEN
              INFO(1) = -11
              INFO(2) = LONG + PLEFTW - POSWCB - 1
              GOTO 260
            END IF
          ENDIF
          POSIWCB = POSIWCB - LONG - 2
          PTRICB(STEP(INODE)) = POSIWCB + 1
          IWCB(POSIWCB + 1) = LONG
          IWCB(POSIWCB + 2) = 1
          POSWCB = POSWCB - LONG
          PTRACB(STEP(INODE)) = POSWCB + 1
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          IWCB(POSIWCB + 3), 
     *          LONG, MPI_INTEGER, COMM, IERR)
            CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          W(POSWCB + 1), LONG, 
     *          MPI_COMPLEX, COMM, IERR)
           DO JJ=0, LONG-1
            RHSSOL(IWCB(POSIWCB+3+JJ)) = W(POSWCB+1+JJ)
           ENDDO
           IWCB(POSIWCB + 2) = 0
           CALL CMUMPS_151(N, KEEP(28), IWCB, LIWW, W, LWC,
     *        POSWCB,POSIWCB,PTRICB,PTRACB)
          ENDIF
          IF = FRERE( STEP(INODE) )
          DO WHILE ( IF .GT. 0 )
            IF ( CMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *           SLAVEF) .eq. MYID ) THEN
              IPOOL( IIPOOL ) = IF
              IIPOOL = IIPOOL + 1
            END IF
            IF = FRERE( STEP( IF ) )
          END DO
          IPOOL( IIPOOL ) = INODE
          IIPOOL = IIPOOL + 1
      ELSE IF ( MSGTAG .EQ. BACKSLV_MASTER2SLAVE ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   INODE, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NROW_RECU, 1, MPI_INTEGER, COMM, IERR )
        IPOS   = PTRIST( STEP(INODE) )
        NPIV   = - IW( IPOS     )
        NROW_L =   IW( IPOS + 1 )
        APOS   =   IW( IPOS + 3 )
        IF ( NROW_L .NE. NROW_RECU ) THEN
          WRITE(*,*) 'Error1 in 41S : NROW L/RECU=',NROW_L, NROW_RECU
          CALL CMUMPS_ABORT()
        END IF
        LONG = NROW_L + NPIV
        IF ( POSWCB - LONG .LT. PLEFTW - 1 ) THEN
           CALL CMUMPS_95(N, KEEP(28), IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
           IF ( POSWCB - LONG .LT. PLEFTW - 1 ) THEN
             INFO(1) = -11
             INFO(2) = LONG - POSWCB
             GOTO 260
           END IF
        END IF
        P_UPDATE  = PLEFTW
        P_SOL_MAS = PLEFTW + NPIV
        PLEFTW    = P_SOL_MAS + NROW_L
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   W( P_SOL_MAS),NROW_L, MPI_COMPLEX,
     *                   COMM, IERR )
        CALL CGEMV( 'N', NPIV, NROW_L, ALPHA, A( APOS ), NPIV,
     *              W( P_SOL_MAS ), 1, ZERO, W( P_UPDATE ), 1 )
        PLEFTW = PLEFTW - NROW_L
 100    CALL CMUMPS_63( INODE, W(P_UPDATE), NPIV,
     *                                MSGSOU, 
     *                                BACKSLV_UPDATERHS,
     *                                COMM, IERR )
        IF ( IERR .EQ. -1 ) THEN
          CALL CMUMPS_41(
     *     .FALSE., FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, III, IISSAR,
     *     NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *     IPOOL, LPOOL, STEP,
     *     FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,
     *     PTRIST, IW, LIW, A, LA, W2, MYLEAFE, RHSSOL, RHS )
          IF ( INFO( 1 ) .LT. 0 ) GOTO 270
          GOTO 100
        ELSE IF ( IERR .EQ. -2 ) THEN
          INFO( 1 ) = -17
          INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
          GOTO 260
        END IF
        PLEFTW = PLEFTW - NPIV
      ELSE IF ( MSGTAG .EQ. BACKSLV_UPDATERHS ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   INODE, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NPIV, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   W2, NPIV, MPI_COMPLEX,
     *                   COMM, IERR )
        W( PTRACB(STEP(INODE)): PTRACB(STEP(INODE)) + NPIV - 1 )
     *    =   W( PTRACB(STEP(INODE)): PTRACB(STEP(INODE)) + NPIV - 1 )
     *    +   W2( 1: NPIV )
        IWCB( PTRICB(STEP(INODE)) + 1 )
     *    = IWCB( PTRICB(STEP(INODE)) + 1 ) - 1
        IF ( IWCB( PTRICB(STEP(INODE)) + 1 ) .EQ. 1 ) THEN
          IPOS  = PTRIST(STEP(INODE)) + 2
          LIELL = IW(IPOS-2)+IW(IPOS+1)
          IPOS  = IPOS + 1
          NPIV_CHECK = IW( IPOS )
          IF ( NPIV .NE. NPIV_CHECK ) THEN
            WRITE(*,*) 'Error 2 in usedbympima41s.F'
            CALL CMUMPS_ABORT()
          END IF
          IPOS  = IPOS + 1
          APOS  = IW( IPOS )
          IF ( KEEP(50) .eq. 0 ) THEN
            CALL CTRSV( 'U', 'N', 'U', NPIV, A(APOS), LIELL,
     *                  W(PTRACB(STEP(INODE))),1)
          ELSE
            CALL CTRSV( 'U', 'N', 'U', NPIV, A(APOS), NPIV,
     *                  W(PTRACB(STEP(INODE))),1)
          END IF
          IPOS = PTRIST(STEP(INODE)) + 5
          IPOS = IW( IPOS ) + IPOS + 1
          DO I = 1, NPIV
            JJ = IW( IPOS + I - 1 )
            RHSSOL( JJ ) = W( PTRACB(STEP(INODE)) + I - 1 )
          END DO
          IN = INODE
  200     IN = FILS(IN)
          IF (IN .GT. 0) GOTO 200
          IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
              CALL CMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                       FEUILLE, SLAVEF )
              NBFINF = NBFINF - 1
            ENDIF
            IWCB( PTRICB(STEP(INODE)) + 1 )
     *        = IWCB( PTRICB(STEP(INODE)) + 1 ) - 1
            CALL CMUMPS_151(N, KEEP(28), IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
            GOTO 270
          ENDIF
          DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
          END DO
          IN = -IN
 300      CONTINUE
            IF (CMUMPS_275(STEP(IN),PROCNODE_STEPS,
     *          SLAVEF) .EQ. MYID) THEN
              IPOOL(IIPOOL ) = IN
              IIPOOL = IIPOOL + 1
            ELSE
              PROCDEST = CMUMPS_275( STEP(IN), PROCNODE_STEPS,
     *                   SLAVEF )
              IF ( .NOT. DEJA_SEND( PROCDEST ) ) THEN
 110            CALL CMUMPS_78( IN, 0, 0, LIELL,
     *          IWCB( PTRICB(STEP(INODE))+2 ),
     *          W( PTRACB(STEP(INODE))),
     *          PROCDEST, NOEUD, COMM, IERR )
                IF ( IERR .EQ. -1 ) THEN
                  CALL CMUMPS_41(
     *            .FALSE., FLAG,
     *            BUFR, LBUFR, LBUFR_BYTES,
     *            MYID, SLAVEF, COMM,
     *            N, IWCB, LIWW, POSIWCB,
     *            W, LWC, POSWCB,
     *            IIPOOL, III, IISSAR,
     *            NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *            IPOOL, LPOOL, STEP,
     *            FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,
     *            PTRIST, IW, LIW, A, LA, W2, MYLEAFE, RHSSOL, 
     *            RHS)
                  IF ( INFO( 1 ) .LT. 0 ) GOTO 270
                  GOTO 110
                ELSE IF ( IERR .eq. -2 ) THEN
                  INFO(1) = -17
                  INFO(2) = LIELL * KEEP(35) +
     *                    ( LIELL + 2 ) * KEEP(34)
                  GOTO 260
                END IF
                DEJA_SEND( PROCDEST ) = .TRUE.
              END IF
            END IF
          IN = FRERE( STEP( IN ) )
          IF ( IN .GT. 0 ) GOTO 300
          IWCB( PTRICB(STEP( INODE )) + 1 )
     *      = IWCB( PTRICB(STEP( INODE )) + 1 ) - 1
          CALL CMUMPS_151(N, KEEP(28), IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
        END IF
      ELSE IF (MSGTAG.EQ.TERREUR) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GO TO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     *      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
      ELSE
          INFO(1) = -100
          INFO(2) = MSGTAG
          GOTO 260
      ENDIF
      GO TO 270
 260  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
 270  CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_42
      SUBROUTINE CMUMPS_324(A, NFRONT, NPIV, NBROW, K59)
      IMPLICIT NONE
      INTEGER NFRONT, NPIV, NBROW, K59
      COMPLEX A(NFRONT*(NBROW+NPIV))
      INTEGER IOLD, INEW, I , J
      IF ((NPIV.EQ.0).OR.(NFRONT.EQ.NPIV)) GOTO 500
      IOLD = NFRONT + 1
      INEW = NPIV + 1
      IF ( NPIV < K59 ) THEN
        DO I = 1, NBROW+NPIV-1
          DO J = 0, NPIV-1
            A(INEW+J) = A(IOLD+J)
          ENDDO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        ENDDO
      ELSE
        DO I = 1, NPIV - 1
          DO J = 0, I
            A( INEW + J ) = A( IOLD + J )
          END DO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        END DO
        DO I = NPIV, NBROW+NPIV-1
          DO J = 0, NPIV - 1
            A( INEW + J ) = A( IOLD + J )
          END DO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        ENDDO
      END IF
 500  RETURN
      END SUBROUTINE CMUMPS_324
      SUBROUTINE CMUMPS_325(A, NPIV, NBROW)
      IMPLICIT NONE
      INTEGER NPIV, NBROW
      COMPLEX A((NBROW+NPIV)*(NBROW+NPIV))
      INTEGER IOLD, INEW, I , J, LDA
      IF ((NPIV.EQ.0).OR.(NBROW.EQ.0)) GOTO 500
      LDA  = NBROW+NPIV
      IOLD = LDA + 1
      INEW = NPIV + 1
      DO I = 1, NBROW+NPIV-1
        DO J = 0, NPIV-1
         A(INEW+J) = A(IOLD+J)
        ENDDO
        INEW = INEW + NPIV
        IOLD = IOLD + LDA
      ENDDO
 500  RETURN
      END SUBROUTINE CMUMPS_325
      SUBROUTINE CMUMPS_308( N, PERM, IW, LIW )
      IMPLICIT NONE
      INTEGER N, LIW
      INTEGER PERM( N ), IW( LIW )
      INTEGER I, SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, LIW - 1
          IF ( PERM( IW( I ) ) .GT. PERM( IW( I + 1 ) ) ) THEN
            DONE = .FALSE.
            SWAP  = IW( I + 1 )
            IW( I + 1 ) = IW( I )
            IW( I ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_308
      SUBROUTINE CMUMPS_309( N, NASS1, PERM, ITLOC,
     *                             SMALL, LSMALL,
     *                             LARGE, LLARGE,
     *                             MERGE, LMERGE )
      IMPLICIT NONE
      INTEGER N, NASS1, LSMALL, LLARGE, LMERGE
      INTEGER PERM( N ), ITLOC( N )
      INTEGER SMALL(LSMALL), LARGE(LLARGE), MERGE(LMERGE)
      INTEGER PSMALL, PLARGE, PMERGE, VSMALL, VLARGE, VMERGE
      PSMALL = 1
      PLARGE = 1
      PMERGE = 1
      DO WHILE ( PSMALL .LE. LSMALL .or. PLARGE.LE. LLARGE )
        IF ( PSMALL .GT. LSMALL ) THEN
          VMERGE = LARGE( PLARGE )
          PLARGE = PLARGE + 1
        ELSE IF ( PLARGE .GT. LLARGE ) THEN
          VMERGE = SMALL( PSMALL )
          PSMALL = PSMALL + 1
        ELSE
          VSMALL = SMALL( PSMALL )
          VLARGE = LARGE( PLARGE )
          IF ( PERM( VSMALL ) .LT. PERM( VLARGE ) ) THEN
            VMERGE = VSMALL
            PSMALL   = PSMALL + 1
          ELSE
            VMERGE = VLARGE
            PLARGE   = PLARGE + 1
          END IF
        END IF
        MERGE( PMERGE ) = VMERGE
        ITLOC( VMERGE ) = PMERGE + NASS1
        PMERGE = PMERGE + 1
      END DO
      PMERGE = PMERGE - 1
      RETURN
      END SUBROUTINE CMUMPS_309
      SUBROUTINE CMUMPS_88( MYID, N, PERM,
     *                             ITLOC,
     *                             IW, LIW, NASS1 )
      IMPLICIT NONE
      INTEGER MYID, N, LIW, NASS1
      INTEGER PERM( N ), ITLOC( N ), IW( LIW )
      INTEGER I, J
      IF ( LIW .eq. 0 ) RETURN
      IF ( NASS1 .eq. 0 ) THEN
      END IF
      DO I = 1, LIW - 1
        IF ( PERM( IW( I ) ) .GE. PERM( IW( I + 1 ) ) ) THEN
          WRITE(*,*) MYID,':Error in CHECK_ORDER_LIST: not sorted'
          WRITE(*,*) MYID,':List = ',IW(1:LIW)
          WRITE(*,*) MYID,':Permuted list =', (PERM(IW(J)),J=1,LIW)
          CALL CMUMPS_ABORT()
        END IF
        IF (ITLOC( IW( I + 1 ) ) .NE. ITLOC( IW( I ) ) + 1 ) THEN
          WRITE(*,*) MYID,':Error in MUMP_CHECK_ORDER_LIST: bad ITLOC'
          WRITE(*,*) MYID,':ITLOC( IW ) =',(ITLOC(IW(J)),J=1,LIW)
          CALL CMUMPS_ABORT() 
        END IF
      END DO
      IF ( ITLOC( IW( 1 ) )  .NE. NASS1 + 1 ) THEN
        WRITE(*,*) MYID,': Bad First ITLOC ITLOC(IW(1)),NASS1+1=',
     *              ITLOC( IW( 1 ) ), NASS1 + 1
        CALL CMUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE CMUMPS_88
      SUBROUTINE CMUMPS_181( N, NBROOT, NBLEAF, POOL, 
     &          LPOOL, PERM, FILS, FRERE, NSTK )
      IMPLICIT NONE
      INTEGER N, NBROOT, NBLEAF, LPOOL
      INTEGER POOL( LPOOL ), PERM( N ), FILS( N ),
     *        FRERE( N ), NSTK( N )
      INTEGER IPERM, INODE, IN
      IPERM = 1
      DO WHILE ( NBLEAF .NE. 0 )
        INODE = POOL( NBLEAF )
        NBLEAF = NBLEAF - 1
        IN = INODE
        DO WHILE ( IN .GT. 0 )
          PERM ( IN ) = IPERM
          IPERM = IPERM + 1
          IN = FILS( IN )
        END DO
        IN = FRERE( INODE )
        DO WHILE( IN .GT. 0 )
          IN = FRERE( IN )
        END DO
        IF ( IN .eq. 0 ) THEN
          NBROOT = NBROOT - 1
        ELSE
          IN = - IN
          NSTK( IN ) = NSTK( IN ) - 1
          IF ( NSTK( IN ) .eq. 0 ) THEN
            NBLEAF = NBLEAF + 1
            POOL( NBLEAF ) = IN
          END IF
        END IF
      END DO
      END SUBROUTINE CMUMPS_181
      SUBROUTINE CMUMPS_81(MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF,
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     *           INTARR, ITLOC, FILS, FRERE, 
     *           SON_LEVEL2, NIV1, NBPROCFILS, KEEP, IFLAG)
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG
      INTEGER KEEP(150), NFRONT_EFF
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PTRAIW(N), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE(KEEP(28)),
     *        NBPROCFILS(KEEP(28))
      INTEGER INTARR(MAX(1,KEEP(14)))
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, INEW, IOLDP2, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K
      LOGICAL LEVEL1
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1
      NEWEL = IOLDP2 + NASS1
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW) = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      IF (NUMSTK .NE. 0) THEN
        NTOTFS = NUMORG
        ISON = IFSON
        ICT11 = IOLDP2 + NFRONT
        DO 100 IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2)
          NELIM = IW(J2 + 1)
          NPIVS = IW(J2 + 3)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON = IW(J2 + 5)
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     *                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2)
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          IF (NELIM .EQ. 0) GOTO 70
          DO 60 JJ = J1, J3
            NTOTFS = NTOTFS + 1
            JT1 = IW(JJ)
            IW(ICT11 + NTOTFS) = JT1
            ITLOC(JT1) = NTOTFS
            IW(JJ) = NTOTFS
            IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
   60     CONTINUE
   70     J1 = J3 + 1
          IF (NASS1 .NE. NFRONT) THEN
            DO 80 JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN
                NEWEL = NEWEL + 1
                NFRONT_EFF = NFRONT_EFF + 1
                IW(NEWEL) = J
                IW(JJ) = NFRONT_EFF
                ITLOC(J) = NFRONT_EFF
              ELSE
                IW(JJ) = ITLOC(J)
              ENDIF
   80       CONTINUE
          ELSE
            DO 90 JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
   90       CONTINUE
          ENDIF
          ISON = FRERE(STEP(ISON))
  100   CONTINUE
      ENDIF
      IBROT = INODE
      DO 120 IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO 110 JJ = J1, J2
            J = INTARR(JJ)
            IF (ITLOC(J) .EQ. 0) THEN
              NEWEL = NEWEL + 1
              NFRONT_EFF = NFRONT_EFF + 1
              IW(NEWEL) = J
              INTARR(JJ) = NFRONT_EFF
              ITLOC(J) = NFRONT_EFF
            ELSE
              INTARR(JJ) = ITLOC(J)
            ENDIF
  110     CONTINUE
        ENDIF
  120 CONTINUE
      IF (NFRONT.NE.NFRONT_EFF) THEN
       IF (NUMORG.EQ.NASS1) THEN
        IP1 = IOLDPS + HF
        IP2 = IOLDPS + HF + NFRONT_EFF - 1
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ELSE
        IP1 = IOLDPS + NFRONT + HF + NUMORG
        IP2 = IOLDPS + HF + NFRONT_EFF + NUMORG
        IW(IP2:IP2+NASS-1) = IW(IP1:IP1+NASS-1)
        IP1 = IOLDPS + NASS1 + HF
        IP2 = IOLDPS + HF + NFRONT - 1
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
        IP1 = IOLDPS + HF
        IP2 = IOLDPS + HF + NUMORG - 1
        DO I = IP1, IP2
          IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ENDIF
      ELSE
       IP1 = IOLDPS + NASS1 + HF
       IP2 = IOLDPS + HF + NFRONT - 1
       DO I = IP1, IP2
        IW(I + NFRONT) = IW(I)
       ENDDO
       IP1 = IOLDPS + HF
       IP2 = IOLDPS + HF + NUMORG - 1
       DO I = IP1, IP2
         IW(I + NFRONT) = IW(I)
       ENDDO
      ENDIF
      K1 = IOLDPS + HF + NUMORG
      K2 = K1 + NFRONT_EFF - 1 + NASS
      DO 150 K = K1, K2
        I = IW(K)
        ITLOC(I) = 0
  150 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_81
      SUBROUTINE CMUMPS_86(MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF, PERM, 
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     *           INTARR, ITLOC, FILS, FRERE_STEPS, 
     *           SON_LEVEL2, NIV1, NBPROCFILS, KEEP, IFLAG)
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG, 
     *        NFRONT_EFF
      INTEGER KEEP(150)
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PTRAIW(N), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE_STEPS(KEEP(28)),
     *        NBPROCFILS(KEEP(28)), PERM(N)
      INTEGER INTARR(MAX(1,KEEP(14)))
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K, JDEBROW, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     *        LAST_J_ASS, JMIN, MIN_PERM
      LOGICAL LEVEL1
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      ALLOCATE(PTTRI(NUMSTK+1), stat=allocok)
      IF (allocok .GT. 0) THEN
       IFLAG = -13
       GOTO 800
      ENDIF
      ALLOCATE(PTLAST(NUMSTK+1), stat=allocok)
      IF (allocok .GT. 0) THEN
       IFLAG = -13
       GOTO 800
      ENDIF
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      JDEBROW = PTRAIW(INODE)+3
      PTTRI(NUMSTK+1)  = JDEBROW
      PTLAST(NUMSTK+1) = JDEBROW + INTARR(JDEBROW-3) - 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW)         = JT1
      IW (INEW+NFRONT) = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      NTOTFS = NUMORG
      IF (NUMSTK .NE. 0) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK   = IW(J2    )
          NELIM  = IW(J2 + 1)
          NPIVS  = IW(J2 + 3)
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5)
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS  = NPIVS + LSTK
          NROWS  = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     *                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2)
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          IF (NELIM .NE. 0) THEN
            DO JJ = J1, J3
              NTOTFS = NTOTFS + 1
              JT1 = IW(JJ)
              IW(ICT11 + NTOTFS) = JT1
              ITLOC(JT1) = NTOTFS
              IW(JJ) = NTOTFS
              IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
            ENDDO
          ENDIF
          PTTRI(IELL)  = J2+1
          PTLAST(IELL) = J2
          J1 = J3 + 1
          IF (NASS1 .NE. NFRONT) THEN
            DO JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN 
                PTTRI(IELL) = JJ
                EXIT
              ENDIF
            ENDDO
          ELSE
            DO JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
            ENDDO
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IF (NFRONT.EQ.NASS1) GOTO 600
 199  IF ( PTTRI( NUMSTK + 1 ) .LE. PTLAST( NUMSTK + 1 ) ) THEN
      IF ( ITLOC( INTARR( PTTRI( NUMSTK + 1 ) ) ) .NE. 0 ) THEN
       PTTRI( NUMSTK + 1 ) = PTTRI( NUMSTK + 1 ) + 1
       GOTO 199
      END IF
      END IF
      MIN_PERM = N + 1
      DO IELL = 1, NUMSTK 
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
      IELL = NUMSTK + 1
      ILOC =  PTTRI( IELL )
      IF ( ILOC .LE. PTLAST( IELL ) ) THEN
        IF ( PERM( INTARR( ILOC ) ) .LT. MIN_PERM ) THEN
         JMIN        = INTARR( ILOC )
         MIN_PERM = PERM( JMIN )
        END IF
      END IF
      NEWEL = IOLDP2 + NASS1 + NFRONT
      DO WHILE ( MIN_PERM .NE. N + 1 )
          NEWEL  = NEWEL + 1
          NFRONT_EFF = NFRONT_EFF + 1
          IW( NEWEL ) = JMIN
          ITLOC( JMIN ) = NFRONT_EFF
          LAST_J_ASS = JMIN
          MIN_PERM = N + 1
          DO IELL = 1,  NUMSTK
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
              IF ( IW( PTTRI( IELL ) ) .eq. LAST_J_ASS )
     *        PTTRI( IELL ) = PTTRI( IELL ) + 1
            ENDIF
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
             IF ( PERM(IW( PTTRI( IELL )) ) .LT. MIN_PERM ) THEN
                JMIN        = IW( PTTRI( IELL ) )
                MIN_PERM = PERM( JMIN )
             END IF
            END IF
          END DO
          IELL = NUMSTK + 1
 145      CONTINUE
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
            IF ( INTARR( PTTRI( IELL ) ) .eq. LAST_J_ASS ) THEN
              PTTRI( IELL ) = PTTRI( IELL ) + 1 
              GOTO 145
            END IF
          END IF
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
            IF (PERM(INTARR( PTTRI(IELL) )) .LT. MIN_PERM) THEN
              JMIN        = INTARR( PTTRI(IELL) )
              MIN_PERM = PERM( JMIN )
            END IF
          END IF
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT) THEN 
       IBROT = FILS( INODE )
       DO IORG = 2, NUMORG
         J1    = PTRAIW(IBROT) + 2
         IBROT = FILS( IBROT )
         J2    = J1 + INTARR(J1 - 2)
         J1    = J1 + 1
         DO JJ = J1, J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
       ENDDO
      ENDIF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
       IW( IOLDP2+NASS1+1 : IOLDP2+NFRONT_EFF ) = 
     *  IW( ICT11+NASS1+1 : ICT11+NFRONT_EFF ) 
      ELSE
        CALL CMUMPS_308( N, PERM, 
     *           IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
          CALL CMUMPS_309( N, NASS1, PERM, ITLOC,
     *               IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE,
     *               IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1,
     *               IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 )
      END IF
      IF ( NFRONT_EFF .NE. NEWEL1_SAVE ) THEN
        IP1 = IOLDPS + NASS1 + HF
        IP2 = IOLDPS + HF + NFRONT - 1
        DO I = IP1, IP2
          IW(I + NFRONT) = IW(I)
        ENDDO
      END IF
      IF (NFRONT.NE.NFRONT_EFF) THEN
        IP1 = IOLDPS + NFRONT + HF 
        IP2 = IOLDPS + NFRONT_EFF + HF
        IW (IP2:IP2+NFRONT_EFF) = IW (IP1:IP1+NFRONT_EFF)
      ENDIF
  600 CONTINUE
      IF ((NUMSTK .NE. 0).AND.(NFRONT.GT.NASS1)) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2)
          NELIM = IW(J2 + 1)
          NPIVS = IW(J2 + 3)
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5)
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2)
          ENDIF
          HS = NSLSON + 6
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          J1 = J3 + 1
          DO JJ = J1, J2
              J = IW(JJ)
                IW(JJ) = ITLOC(J)
          ENDDO
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IBROT = INODE
      DO IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO JJ = J1, J2
            J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
          ENDDO
        ENDIF
      ENDDO
        K1 = IOLDPS + HF + NUMORG
        K2 = K1 + NFRONT_EFF - 1 + NASS
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
  800 CONTINUE
      IF (ALLOCATED(PTTRI)) DEALLOCATE(PTTRI)
      IF (ALLOCATED(PTLAST)) DEALLOCATE(PTLAST)
      RETURN
      END SUBROUTINE CMUMPS_86
      SUBROUTINE CMUMPS_140( N, INODE, IW, LIW, A, LA,
     *                           IOLDPS, POSELT, IFLAG,
     *                           UU, NNEG, NPVW,
     *                           KEEP, MYID, SEUIL )
      IMPLICIT NONE
      INTEGER N, INODE, LIW, LA, IFLAG, NNEG, NPVW
      INTEGER MYID, SLAVEF, IOLDPS, POSELT
      INTEGER KEEP( 150 )
      REAL UU, UUTEMP, SEUIL
      INTEGER POSFAC, NIRBDU, COMP
      COMPLEX A( LA )
      INTEGER IW( LIW )
      INTEGER INOPV, IFINB, NFRONT, NPIV, LAELL, NBTLKJ
      INTEGER NASS, NEL1, LDA
      EXTERNAL CMUMPS_330, CMUMPS_222, CMUMPS_234, 
     *         CMUMPS_230, CMUMPS_226, 
     *         CMUMPS_237
      INTEGER  CMUMPS_330
      NFRONT = IW(IOLDPS)
      LDA    = NFRONT
      NASS   = IABS(IW(IOLDPS+2))
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5), NASS )
      ENDIF
      UUTEMP = UU
 50   CONTINUE
      CALL CMUMPS_222(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *                INOPV, NNEG, IFLAG,IOLDPS,POSELT,UUTEMP,
     *                SEUIL)
      IF (IFLAG.LT.0) GOTO 500
      IF (INOPV.EQ.1) THEN
         CALL CMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT, KEEP, UUTEMP )
         GOTO 500
      END IF
      IF (INOPV.EQ.2) THEN
         CALL CMUMPS_234(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT,NBTLKJ,KEEP(4),UUTEMP,
     *            KEEP)
         GOTO 50
      ENDIF
      NPVW = NPVW + 1
      IF (NASS.LE.1) THEN
       CALL CMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     *                 IOLDPS,POSELT)
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       GO TO 500
      ENDIF
       CALL CMUMPS_226(NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     *             LDA, UUTEMP, IOLDPS,POSELT,IFINB,
     *             NBTLKJ,KEEP(4))
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       IF (IFINB.EQ.0) GOTO 50
       NPIV   = IW(IOLDPS+1)
       NEL1   = NASS - NPIV
       CALL CMUMPS_234(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT,NBTLKJ,KEEP(4),UUTEMP,
     *            KEEP)
       IF (IFINB.EQ.-1) THEN 
         CALL CMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT, KEEP, UUTEMP)
         GOTO 500
       ENDIF
      GO TO 50
 500  CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_140
      SUBROUTINE CMUMPS_222 
     *   (NFRONT,NASS,N,INODE,IW,LIW,
     *    A,LA, INOPV,
     *    NNEG, 
     *    IFLAG,IOLDPS,POSELT,UU, SEUIL)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV,
     *        IOLDPS, POSELT, NNEG
      COMPLEX A(LA) 
      REAL UU, UULOC, SEUIL
      INTEGER IW(LIW)
       include 'mpif.h'
      INTEGER APOS, LDA
      REAL ZERO, RMAX
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,J1,JJ,J2,IDIAG
      INTRINSIC MAX
      DATA ZERO /0.0E0/
        LDA   = NFRONT
        UULOC = UU
        INOPV   = 0
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
        NASSW   = IABS(IW(IOLDPS+3))
        DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + LDA*(IPIV-1) + NPIV
            IDIAG = APOS + IPIV - NPIVP1
            IF (UULOC.EQ.ZERO) THEN 
              IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
              GO TO 420
            ENDIF
            RMAX = ABS(A(IDIAG))
            J1 = APOS
            J2 = IDIAG - 1
            IF ( J1.LE.J2 ) THEN
             DO JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
             ENDDO
            END IF
            J1 = IDIAG + LDA
            J2 = NFRONT - IPIV 
            DO JJ=1,J2
              RMAX = MAX(ABS(A(J1)),RMAX)
              J1 = J1 + LDA
            ENDDO
            IF (ABS(A(IDIAG)).LE. MAX(UULOC*RMAX, SEUIL)) GOTO 460
            IF (IPIV.EQ.NPIVP1) GO TO 420
            CALL CMUMPS_319( A, LA, IW, LIW,
     *                     IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     *                     LDA, NFRONT, 1 )
            GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
  420 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_222
      SUBROUTINE CMUMPS_226(NFRONT,NASS,N,INODE,IW,LIW,
     *     A,LA,LDA, UU,
     *     IOLDPS,POSELT,IFINB,LKJIB,LKJIT)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFINB,LKJIB,LDA,
     &        NPBEG
      COMPLEX    A(LA)
      REAL   UU
      INTEGER IW(LIW)
      COMPLEX    VALPIV
      INTEGER APOS, IOLDPS, POSELT, K1POS, NCB1
      INTEGER LKJIT
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,LPOS,NEL
      COMPLEX ONE, ALPHA
      REAL  ZERO
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0, ZERO=0.0E0)
      NPIV   = IW(IOLDPS+1)
      NPIVP1 = NPIV + 1
      NEL    = NFRONT - NPIVP1
      IFINB  = 0
      IF (IW(IOLDPS+3).LE.0) THEN
         IW(IOLDPS+3) = MIN0(NASS,LKJIB)
      ENDIF
      JROW2 = IW(IOLDPS+3)
      NPBEG = IW(IOLDPS+4)
      NEL2   = JROW2 - NPIVP1
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      APOS   = POSELT + NPIV*(NFRONT + 1)
      VALPIV = ONE/A(APOS)
      LPOS   = APOS + LDA
      IF ( UU.NE.ZERO) THEN
        CALL CCOPY(NFRONT-NPIVP1, A(LPOS), LDA, A(APOS+1), 1)
      ELSE
        CALL CCOPY(NASS-NPIVP1, A(LPOS), LDA, A(APOS+1), 1)
      END IF
      CALL CMUMPS_CSYR('U', NEL2, -VALPIV, A(LPOS), LDA, A(LPOS+1), 
     &             LDA)
      IF (UU.NE.ZERO) THEN
         CALL CSCAL(NFRONT-NPIVP1, VALPIV, A(LPOS), LDA)
      ELSE
         CALL CSCAL(NASS-NPIVP1, VALPIV, A(LPOS), LDA)
      ENDIF
      IF (NEL2.GT.0) THEN
           K1POS = LPOS + NEL2*LDA
           IF ( UU.NE.ZERO ) THEN
             NCB1  = NFRONT - JROW2
           ELSE
             NCB1  = NASS   - JROW2
           END IF
           CALL CGERU(NEL2, NCB1 , ALPHA, A(APOS+1), 1, 
     &          A(K1POS), LDA, A(K1POS+1), LDA)
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_226
      SUBROUTINE CMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     &           IOLDPS,POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LA,LIW
      COMPLEX    A(LA)
      INTEGER IW(LIW)
      COMPLEX    VALPIV
      INTEGER APOS, POSELT
      INTEGER IOLDPS,NEL
      INTEGER LPOS,JROW
      INTEGER LDA
      COMPLEX ONE
      PARAMETER(ONE=1.0E0)
        NEL    = NFRONT - 1
        APOS   = POSELT 
        IF (NEL.EQ.0) GO TO 500
        VALPIV = ONE/A(APOS)
        LDA    = NFRONT
        LPOS   = APOS + LDA
        CALL CMUMPS_CSYR('U',NEL, -VALPIV, 
     &             A(LPOS), NFRONT, A(LPOS+1), NFRONT)
          DO JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
          END DO
  500   RETURN
        END SUBROUTINE CMUMPS_230
      SUBROUTINE CMUMPS_234(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *    LDA,
     *    IOLDPS,POSELT,LKJIB,LKJIT,UUTEMP,KEEP)
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,LA,LIW
      COMPLEX    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE, KEEP(150)
      INTEGER POSELT, LDA
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1, NEL11
      INTEGER LBP, HF
      INTEGER LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER I, Block
      INTEGER BLSIZE
      REAL UUTEMP
      COMPLEX ONE, ALPHA
      REAL ZERO
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0, ZERO = 0.0E0)
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1)
      JROW2  = IABS(IW(IOLDPS+3))
      NPBEG  = IW(IOLDPS+4)
      HF     = 6 + IW(IOLDPS+5)
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      NEL11  = NFRONT - NPIV
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB + NONEL
          IW(IOLDPS+3)= MIN0(NPIV+LKJIB,NASS)
          LKJIB       = MIN0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3) = NASS
        ENDIF
        IW(IOLDPS+4) = NPIV + 1
      ELSEIF (JROW2.LT.NASS) THEN
          IW(IOLDPS+4) = NPIV + 1
          IW(IOLDPS+3) = MIN0(JROW2+LKJIB,NASS)
          LKJIB  = MIN0(LKJIB,NASS-NPIV)
      ENDIF
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = MIN( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + ( IROW - 1 ) * LDA + NPBEG - 1
           UPOS = POSELT + ( NPBEG - 1 ) * LDA + IROW - 1
           APOS = POSELT + ( IROW - 1 ) * LDA + IROW - 1
           CALL CGEMM( 'N','N', Block, NASS - IROW + 1, LKJIW,
     *                ALPHA, A( UPOS ), LDA,
     *                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         END DO
        END IF
       LPOS = POSELT + NASS*LDA + NPBEG - 1
       UPOS = POSELT + (NPBEG-1) * LDA + JROW2
       APOS = POSELT + NASS*LDA + JROW2
       IF ( UUTEMP. NE. ZERO ) THEN
         CALL CGEMM('N', 'N', NEL1, NFRONT-NASS, LKJIW, ALPHA, 
     *              A(UPOS), LDA, A(LPOS), LDA, ONE, 
     *              A(APOS), LDA)
       END IF
      ENDIF
      LBP = JROW2 - NPIV
      IF ( LBP  .GT. 0     .AND.
     *     NEL1 .GE. LKJIB .AND.
     *     NEL1 .GE. LBP + LBP ) THEN
      DO I = 1, LBP
        CALL CMUMPS_319( A, LA, IW, LIW,
     *                 IOLDPS, NPIV+I, NASS-LBP+I, POSELT, NASS,
     *                 LDA, NFRONT, 1 )
      END DO
      END IF
  500 RETURN
      END SUBROUTINE CMUMPS_234
      SUBROUTINE CMUMPS_321(NFRONT,N,INODE,IW,LIW,A,LA,
     *                          IOLDPS,POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LA,LIW
      COMPLEX    A(LA)
      INTEGER IW(LIW)
      COMPLEX    VALPIV
      INTEGER APOS, POSELT
      INTEGER IOLDPS,NPIV,NEL
      INTEGER LPOS,JROW
      COMPLEX ONE
      DATA ONE /1.0E0/
        NPIV   = IW(IOLDPS+1)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + (NPIV)*NFRONT + NPIV
        IF (NEL.EQ.0) GO TO 500
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        CALL CMUMPS_CSYR('U', NEL, -VALPIV, A(LPOS), NFRONT, A(LPOS+1), 
     &             NFRONT)
        DO JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
        ENDDO
  500   RETURN
        END SUBROUTINE CMUMPS_321
        SUBROUTINE CMUMPS_319( A, LA, IW, LIW,
     *                       IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     *                       LDA, NFRONT, LEVEL )
        IMPLICIT NONE
      INTEGER LA, LIW, IOLDPS, NPIVP1, POSELT, IPIV
      INTEGER LDA, NFRONT, NASS, LEVEL
      COMPLEX A( LA )
      INTEGER IW( LIW )
      INTEGER ISW, ISWPS1, ISWPS2, APOS, IDIAG, HF
      COMPLEX SWOP
            APOS = POSELT + LDA*(IPIV-1) + NPIVP1-1
            IDIAG = APOS + IPIV - NPIVP1
            HF = 6 + IW( IOLDPS + 5 )
            ISWPS1 = IOLDPS + HF + NPIVP1 - 1
            ISWPS2 = IOLDPS + HF + IPIV - 1
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            ISW = IW(ISWPS1+NFRONT)
            IW(ISWPS1+NFRONT) = IW(ISWPS2+NFRONT)
            IW(ISWPS2+NFRONT) = ISW
            IF ( LEVEL .eq. 2 ) THEN
              CALL CSWAP( NPIVP1 - 1,
     *            A( POSELT + NPIVP1 - 1 ), LDA,
     *            A( POSELT + IPIV - 1 ), LDA )
            END IF
            CALL CSWAP( NPIVP1-1, A( POSELT+(NPIVP1-1) * LDA ), 1,
     *                          A( POSELT + (IPIV-1) * LDA ), 1 )
             CALL CSWAP( IPIV - NPIVP1 - 1,
     *                  A( POSELT + NPIVP1 * LDA + NPIVP1 - 1),
     *                  LDA, A( APOS + 1 ), 1 )
            SWOP = A( IDIAG )
            A( IDIAG  ) = A( POSELT+(NPIVP1-1)*LDA+NPIVP1-1 )
            A( POSELT+(NPIVP1-1)*LDA+NPIVP1-1) = SWOP
            CALL CSWAP( NASS - IPIV, A( APOS + LDA ), LDA,
     *                  A( IDIAG + LDA ), LDA )
            IF ( LEVEL .eq. 1 ) THEN
              CALL CSWAP( NFRONT - NASS,
     *        A( APOS  + ( NASS - IPIV + 1 ) * LDA ), LDA,
     *        A( IDIAG + ( NASS - IPIV + 1 ) * LDA ), LDA )
            END IF
        RETURN
        END SUBROUTINE CMUMPS_319
      SUBROUTINE CMUMPS_237(NFRONT,NASS,N,INODE,
     *    IW,LIW,A,LA,
     *    LDA,
     *    IOLDPS,POSELT,KEEP,UU)
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,INODE,LA,LIW
      COMPLEX    A(LA)
      REAL   UU
      INTEGER IW(LIW) 
      INTEGER KEEP(150)
      INTEGER POSELT, LDA
      INTEGER IOLDPS
      INTEGER BLSIZE, Block, IROW, NPIV, I
      INTEGER UPOS, APOS, LPOS
      COMPLEX  ONE, ALPHA
      REAL ZERO
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0,ZERO=0.0E0)
      IF ( NFRONT - NASS > KEEP(57) ) THEN
        BLSIZE = KEEP(58)
      ELSE
        BLSIZE = NFRONT - NASS
      END IF
      NPIV = IW( IOLDPS + 1 )
      IF ( NFRONT - NASS .GT. 0 ) THEN
       IF ( UU.eq.ZERO ) THEN
         CALL CTRSM( 'L', 'L', 'N', 'N',
     *               NPIV, NFRONT-NPIV, ONE,
     *               A( POSELT ), LDA,
     *               A( POSELT + LDA * NPIV ), LDA )
       ENDIF
       DO IROW = 1, NFRONT - NASS, BLSIZE
        Block = MIN( BLSIZE, NFRONT - NASS - IROW + 1 )
        LPOS = POSELT + NASS*LDA + (IROW-1) * LDA
        APOS = POSELT + NASS*LDA + (IROW-1) * LDA + NASS + IROW - 1
        UPOS = POSELT + NASS
        DO I = 1, NPIV
          CALL CCOPY( Block, A(LPOS+I-1), LDA, A(UPOS+(I-1)*LDA),1)
          CALL CSCAL( Block, A(POSELT+(LDA+1)*(I-1)),
     *                A( UPOS + ( I - 1 ) * LDA ), 1 )
        ENDDO
        CALL CGEMM( 'N', 'N', Block, NFRONT-NASS-IROW+1, NPIV,
     *              ALPHA,  A( UPOS ), LDA,
     *              A( LPOS ), LDA,
     *              ONE,
     *              A( APOS ), LDA )
       END DO
      END IF
      END SUBROUTINE CMUMPS_237
      SUBROUTINE CMUMPS_320( BUF, BLOCK_SIZE,
     *                           MYROW, MYCOL, NPROW, NPCOL,
     *                           A, LOCAL_M, LOCAL_N, N, MYID, COMM )
      IMPLICIT NONE
      INTEGER BLOCK_SIZE, NPROW, NPCOL, LOCAL_M, LOCAL_N, N, COMM
      INTEGER MYROW, MYCOL, MYID
      COMPLEX BUF( BLOCK_SIZE * BLOCK_SIZE )
      COMPLEX A( LOCAL_M, LOCAL_N )
      INTEGER NBLOCK, IBLOCK, JBLOCK, IBLOCK_SIZE, JBLOCK_SIZE
      INTEGER ROW_SOURCE, ROW_DEST, COL_SOURCE, COL_DEST
      INTEGER IGLOB, JGLOB
      INTEGER IROW_LOC_SOURCE, JCOL_LOC_SOURCE
      INTEGER IROW_LOC_DEST, JCOL_LOC_DEST
      INTEGER PROC_SOURCE, PROC_DEST
      NBLOCK = ( N - 1 ) / BLOCK_SIZE + 1
      DO IBLOCK = 1, NBLOCK
        IF ( IBLOCK .NE. NBLOCK
     *    ) THEN
          IBLOCK_SIZE = BLOCK_SIZE
        ELSE
          IBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
        END IF
        ROW_SOURCE = MOD( IBLOCK - 1, NPROW ) 
        COL_DEST   = MOD( IBLOCK - 1, NPCOL )
        IGLOB = ( IBLOCK - 1 ) * BLOCK_SIZE + 1
        IROW_LOC_SOURCE = BLOCK_SIZE *
     *                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( IGLOB - 1, BLOCK_SIZE ) + 1
        JCOL_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( IGLOB - 1, BLOCK_SIZE ) + 1
        DO JBLOCK = 1, IBLOCK
          IF ( JBLOCK .NE. NBLOCK
     *      ) THEN
            JBLOCK_SIZE = BLOCK_SIZE
          ELSE
            JBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
          END IF
          COL_SOURCE = MOD( JBLOCK - 1, NPCOL )
          ROW_DEST   = MOD( JBLOCK - 1, NPROW )
          PROC_SOURCE = ROW_SOURCE * NPCOL + COL_SOURCE
          PROC_DEST   = ROW_DEST   * NPCOL + COL_DEST
          IF ( PROC_SOURCE .eq. PROC_DEST ) THEN
           IF ( MYID .eq. PROC_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     *                  ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            IF ( IBLOCK .eq. JBLOCK ) THEN
              IF ( IBLOCK_SIZE .ne. JBLOCK_SIZE ) THEN
                WRITE(*,*) MYID,': Error in calling transdiag:unsym'
                CALL CMUMPS_ABORT()
              END IF
              CALL CMUMPS_327( A( IROW_LOC_SOURCE,
     *                 JCOL_LOC_SOURCE),
     *                 IBLOCK_SIZE, LOCAL_M )
            ELSE
              CALL CMUMPS_326(
     *           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ),
     *           A( IROW_LOC_DEST, JCOL_LOC_DEST ),
     *           IBLOCK_SIZE, JBLOCK_SIZE, LOCAL_M )
            END IF
           END IF
          ELSE IF (  MYROW .eq. ROW_SOURCE 
     *    .AND. MYCOL .eq. COL_SOURCE ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL CMUMPS_293( BUF,
     *           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ), LOCAL_M,
     *           IBLOCK_SIZE, JBLOCK_SIZE, COMM, PROC_DEST )
          ELSE IF ( MYROW .eq. ROW_DEST 
     *    .AND.     MYCOL .eq. COL_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL CMUMPS_281( BUF,
     *           A( IROW_LOC_DEST, JCOL_LOC_DEST ), LOCAL_M,
     *           JBLOCK_SIZE, IBLOCK_SIZE, COMM, PROC_SOURCE )
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_320
      SUBROUTINE CMUMPS_293( BUF, A, LDA, M, N, COMM, DEST )
      IMPLICIT NONE
      INTEGER M, N, LDA, DEST, COMM
      COMPLEX BUF(*), A(LDA,*)
      INTEGER I, IBUF, IERR
      INTEGER J
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IBUF = 1
      DO J = 1, N
        BUF( IBUF: IBUF + M - 1 ) = A( 1 : M, J )
        DO I = 1, M
        END DO
        IBUF = IBUF + M
      END DO
      CALL MPI_SEND( BUF, M * N, MPI_COMPLEX,
     *     DEST, SYMMETRIZE, COMM, IERR )
      RETURN
      END SUBROUTINE CMUMPS_293
      SUBROUTINE CMUMPS_281( BUF, A, LDA, M, N, COMM, SOURCE )
      IMPLICIT NONE
      INTEGER LDA, M, N, COMM, SOURCE
      COMPLEX BUF(*), A( LDA, *)
      INTEGER I, IBUF, IERR
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      CALL MPI_RECV( BUF(1), M * N, MPI_COMPLEX, SOURCE,
     *               SYMMETRIZE, COMM, STATUS, IERR )
      IBUF = 1
      DO I = 1, M
        CALL CCOPY( N, BUF(IBUF), 1, A(I,1), LDA )
        IBUF = IBUF + N
      END DO
      RETURN
      END
      SUBROUTINE CMUMPS_327( A, N, LDA )
      IMPLICIT NONE
      INTEGER N,LDA
      COMPLEX A( LDA, * )
      INTEGER I, J
      DO I = 2, N
        DO J = 1, I - 1
          A( J, I ) = A( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_327
      SUBROUTINE CMUMPS_326( A1, A2, M, N, LD )
      IMPLICIT NONE
      INTEGER M,N,LD
      COMPLEX A1( LD,* ), A2( LD, * )
      INTEGER I, J
      DO J = 1, N
        DO I = 1, M
          A2( J, I ) = A1( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_326
      RECURSIVE SUBROUTINE CMUMPS_274( 
     *   COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, root, OPASSW, OPELIW, ITLOC, FILS,  
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 150 )
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, LIW, LA
      INTEGER COMP
      INTEGER NIRBDU, IFLAG, IERROR, POSFAC, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)), ITLOC(N), FILS(N), ND(KEEP(28))
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      REAL OPASSW, OPELIW
      COMPLEX DBLARR( MAX(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      LOGICAL FLAG
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER INODE, POSITION, NPIV, IERR, IERR_MPI
      INTEGER LAELL, NCOL, POSBLOCFACTO, NROW
      INTEGER IOLDPS, POSELT, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11, LPOS, LPOS2, DPOS, UPOS
      INTEGER I, IPOS, KPOS, IPIV, FPERE, NSLAVES_TOT,
     *        NSLAVES_FOLLOW,  NB_BLOC_FAC
      INTEGER LCONT,NELIM,NASS, LDA, NCOL_TO_SEND,
     *        SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON, SHIFT_VAL_SON
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok
      COMPLEX, DIMENSION(:),ALLOCATABLE :: UIP21K
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      COMPLEX ONE,ALPHA
      COMPLEX VALPIV
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      FPERE    = -1
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     *                 MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     *                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LAELL = NPIV * NCOL
      IF ( NPIV.GT.0 ) THEN
       IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          IERROR = LAELL - LRLU
          GOTO 700
        END IF
        CALL CMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *       NIRBDU, LRLU, IPTRLU,
     *       IWPOS, IWPOSCB, PTRIST, PTRAST,
     *       STEP, PIMASTER, PAMASTER, ITLOC)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     *       ,LRLU,LRLUS
             IFLAG = -9
             IERROR = LAELL -LRLU
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
       END IF
       LRLU  = LRLU - LAELL
       LRLUS = LRLUS - LAELL
      ENDIF
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL CMUMPS_471(.FALSE.,LA-LRLUS,0,LAELL,KEEP)
      IF ( NPIV.GT.0 ) THEN
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 IW( IPIV ), NPIV,
     *                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *              A(POSBLOCFACTO), NPIV*NCOL, MPI_COMPLEX,
     *              COMM, IERR )
      ENDIF
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
        WRITE(*,*) MYID,
     *   ': Internal ERROR 1  in CMUMPS_274',
     *   ' INODE =', INODE,
     *   ' MAITRE_DESC_BANDE not yet received from ', MSGSOU
        CALL CMUMPS_ABORT()
      ENDIF
      DO WHILE ( NBPROCFILS(STEP(INODE)) .NE. 0 ) 
        BLOCKING = .TRUE.
        SET_IRECV=.FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL CMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, CONTSIPERENIV2,
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
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END  DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL CMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, MPI_ANY_TAG, 
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
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS )
      NASS1  = IW( IOLDPS + 1 )
      NROW1  = IW( IOLDPS + 2 )
      NPIV1  = IW( IOLDPS + 3 )
      NSLAV1 = IW( IOLDPS + 5 )
      NSLAVES_FOLLOW = NSLAV1 - 2
      HS     = 6 + NSLAV1
      NCOL1  = LCONT1 + NPIV1
      IF ( LASTBL ) THEN
        IW( IOLDPS + 6 ) = IW( IOLDPS + 6 ) -
     *  ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     *  NB_BLOC_FAC
      END IF
      IF (NPIV.EQ.0) GOTO 200
      IF ( NPIV1 + NCOL .NE. NASS1 ) THEN
        WRITE(*,*) 'SymBLFC Error: NPIV1 + NCOL .NE. NASS1 :',
     *             NPIV1,NCOL,NASS1
        CALL CMUMPS_ABORT()
      END IF
      ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
      DO 100 I = 1, NPIV
        IF (IW(IPIV+I-1).EQ.I) GOTO 100
        ISW = IW(ICT11+I)
        IW(ICT11+I) = IW(ICT11+IW(IPIV+I-1))
        IW(ICT11+IW(IPIV+I-1)) = ISW
        IPOS = POSELT + NPIV1 + I - 1
        KPOS = POSELT + NPIV1 + IW(IPIV+I-1) - 1
        CALL CSWAP(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
 100  CONTINUE
      ALLOCATE( UIP21K( NPIV * NROW1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IFLAG = -13
        IERROR = NPIV * NROW1
        GOTO 700
      END IF
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
        ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     *            stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          IFLAG = -13
          IERROR = NSLAVES_FOLLOW
          GOTO 700
        END IF
        LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     *  IW(IOLDPS+8:IOLDPS+7+NSLAVES_FOLLOW)
      END IF
      CALL CTRSM( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     *             A( POSBLOCFACTO ), NCOL, A(POSELT+NPIV1), NCOL1 )
      LPOS = POSELT + NPIV1
      UPOS = 1
      DO I = 1, NROW1
        UIP21K( UPOS: UPOS + NPIV-1 ) = A( LPOS: LPOS+NPIV-1)
        LPOS = LPOS + NCOL1
        UPOS = UPOS + NPIV
      END DO
      LPOS = POSELT + NPIV1
      DPOS = POSBLOCFACTO
      DO I = 1, NPIV
        VALPIV = ONE / A( DPOS )
        CALL CSCAL( NROW1, VALPIV, A( LPOS ), NCOL1 )
        LPOS = LPOS + 1
        DPOS = DPOS + NCOL + 1
      END DO
      LPOS2 = POSELT + NPIV1
      UPOS = POSBLOCFACTO+NPIV
      LPOS  = LPOS2 + NPIV
      CALL CGEMM('N','N', NCOL-NPIV,NROW1,NPIV,ALPHA,A(UPOS),NCOL,
     *           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
      DPOS = POSELT + NCOL1 - NROW1
      IF ( NROW1 .GT. KEEP(7) ) THEN
        BLSIZE = KEEP(8)
      ELSE
        BLSIZE = NROW1
      ENDIF
      IF ( NROW1 .GT. 0 ) THEN
      DO IROW = 1, NROW1, BLSIZE
        Block = MIN( BLSIZE, NROW1 - IROW + 1 )
        DPOS  = POSELT + NCOL1 - NROW1 + ( IROW - 1 ) * ( NCOL1 + 1 )
        LPOS2 = POSELT + NPIV1         + ( IROW - 1 ) * NCOL1
        UPOS  = ( IROW - 1 ) * NPIV + 1
        DO I = 1, Block
          CALL CGEMV( 'T', NPIV, Block-I+1, ALPHA,
     *                A( LPOS2 + (I - 1) * NCOL1 ), NCOL1,
     *                UIP21K( UPOS + NPIV * ( I - 1 ) ), 1, ONE,
     *                A(DPOS+(NCOL1+1)*(I-1)),NCOL1 )
        END DO
       IF ( NROW1-IROW+1-Block .ne. 0 )
     * CALL CGEMM( 'T', 'N', Block, NROW1-IROW+1-Block, NPIV, ALPHA,
     *             UIP21K( UPOS ), NPIV,
     *             A( LPOS2 + Block * NCOL1 ), NCOL1, ONE,
     *             A( DPOS + Block * NCOL1 ), NCOL1 )
      ENDDO
      ENDIF
 200  CONTINUE
      IW(IOLDPS) = IW(IOLDPS) - NPIV
      IW(IOLDPS + 3) = IW(IOLDPS+3) + NPIV
      IF (LASTBL) IW(IOLDPS+1) = IW(IOLDPS + 3)
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      POSFAC = POSFAC - LAELL
      IWPOS = IWPOS - NPIV
      CALL CMUMPS_471(.FALSE.,LA-LRLUS,0,-LAELL,KEEP)
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
          IERR = -1
           DO WHILE ( IERR .eq. -1 )
            CALL CMUMPS_64(
     *                    INODE, NPIVSENT, FPERE,
     *                    IPOSK, JPOSK,
     *                    UIP21K, NROW1,
     *                    NSLAVES_FOLLOW,
     *                    LIST_SLAVES_FOLLOW(1),
     *                    COMM, IERR )
            IF (IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV= .FALSE.
              MESSAGE_RECEIVED = .FALSE.
              CALL CMUMPS_329( COMM_LOAD, ASS_IRECV,
     *         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *         MPI_ANY_SOURCE, MPI_ANY_TAG,
     *         STATUS, 
     *         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *         IWPOS, IWPOSCB, IPTRLU,
     *         LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *         IFLAG, IERROR, COMM,
     *         NBPROCFILS,
     *         IPOOL, LPOOL, LEAF,
     *         NBFIN, MYID, SLAVEF,
     *         root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *         INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *         LPTRAR, NELT, FRTPTR, FRTELT, 
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
             IF ( IFLAG .LT. 0 ) GOTO 600
            END IF
            IF ( IERR .eq. -2 ) THEN
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
            END IF
           END DO
           DEALLOCATE( LIST_SLAVES_FOLLOW )
      END IF
      IF ( NPIV .NE. 0 ) DEALLOCATE( UIP21K )
      IF ( LASTBL .AND. IW(PTRIST(STEP(INODE)) + 6 ) .eq. 0
     *   .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 1 )
     *   THEN
         DEST = CMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
         CALL CMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     *                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
           IFLAG = -99
           GOTO 700
         END IF
      END IF
      IF ( LASTBL  .AND. IW( PTRIST(STEP(INODE)) + 6 ) .eq. 0
     *     .AND. (KEEP(38).EQ.FPERE) ) THEN
       IOLDPS = PTRIST(STEP( INODE ))
       LCONT  = IW(IOLDPS)
       NROW   = IW(IOLDPS+2)
       NPIV   = IW(IOLDPS+3)
       NASS   = IW(IOLDPS+4)
       NELIM  = NASS-NPIV
       NCOL_TO_SEND =  LCONT-NELIM
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5) 
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NASS
       SHIFT_VAL_SON      = NASS
       LDA                = LCONT + NPIV
       CALL CMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE, 
     *    PTRIST, PTRAST, 
     *    root, NROW, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA, 
     *    ROOT_CONT_STATIC, MYID, COMM,
     *    
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, .FALSE., ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
       IF ( IFLAG < 0 ) GOTO 600
       IF (NELIM.EQ.0) THEN
        CALL CMUMPS_314( N, INODE,  
     *    PTRIST, PTRAST, PTLUST_S, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP)
       ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_274
      SUBROUTINE CMUMPS_141( COMM_LOAD, ASS_IRECV, 
     *           N, INODE, FPERE, IW, LIW, A, LA,
     *           UU, NOFFW,
     *           NPVW,
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, LA, NOFFW, NPVW
      COMPLEX A( LA )
      REAL UU, SEUIL
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER NBFIN, IFLAG, IERROR, LEAF, LPOOL, SLAVEF, 
     *        POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU,
     *        COMP 
      INTEGER NB_BLOC_FAC
      INTEGER ICNTL(40), KEEP(150)
      INTEGER IW( LIW )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL), ITLOC(N)
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND( KEEP(28) )
      INTEGER FRERE(KEEP(28)), FILS(N) 
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), PTLUST_S(KEEP(28)),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     *        PROCNODE_STEPS(KEEP(28)), STEP(N)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      REAL OPASSW, OPELIW
      COMPLEX  DBLARR(MAX(1,KEEP(13)))
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEGKJI, NBTLKJ
      INTEGER NASS, IEND, IOLDPS, POSELT, LDAFS,allocok
      LOGICAL LASTBL 
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      EXTERNAL CMUMPS_223, CMUMPS_235,
     *         CMUMPS_227, CMUMPS_294,
     *         CMUMPS_44
      NB_BLOC_FAC = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST( STEP( INODE ))
      NFRONT = IW(IOLDPS)
      NASS   = IABS(IW(IOLDPS+2))
      LDAFS  = NASS
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5), NASS )
      ENDIF
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID,' : FACTO_NIV2 :failed to allocate ',NASS,
     * ' integers'
        IFLAG=-13
        IERROR=NASS
        GO TO 490
      END IF
 50   CONTINUE
      IBEGKJI = IW(IOLDPS+4)
      CALL CMUMPS_223(
     *                NFRONT,NASS,IBEGKJI, NASS, IPIV,
     *                N,INODE,IW,LIW,A,LA,NOFFW,INOPV,
     *                IFLAG,IOLDPS,POSELT,UU, SEUIL)
      IF (IFLAG.LT.0) GOTO 490
      IF (INOPV.GE.1) THEN
          LASTBL = (INOPV.EQ.1)
          IEND = IW(IOLDPS+1)
          CALL CMUMPS_294( COMM_LOAD, ASS_IRECV,
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDAFS,
     *             IBEGKJI, IEND, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     *
     *             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
      ENDIF
      IF (INOPV.EQ.1) GO TO 500
      IF (INOPV.EQ.2) THEN
         CALL CMUMPS_235(
     *            NASS,N,INODE,IW,LIW,A,LA,
     *            LDAFS, 
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4),KEEP)
         GOTO 50
      ENDIF
      NPVW = NPVW + 1
      IF (NASS.LE.1) THEN
        IFINB = -1
      ELSE
         CALL CMUMPS_227(
     *             NASS, N,INODE,IW,LIW,A,LA,
     *             LDAFS, IOLDPS,POSELT,IFINB,
     *             NBTLKJ,KEEP(4))
      ENDIF
      IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       IF (IFINB.EQ.0) GOTO 50
       IF ((IFINB.EQ.1).OR.(IFINB.EQ.-1)) THEN
          LASTBL = (IFINB.EQ.-1) 
          IEND = IW(IOLDPS+1)
          CALL CMUMPS_294(COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDAFS, 
     *             IBEGKJI, IEND, IPIV, NASS, LASTBL,NB_BLOC_FAC,
     *
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
       ENDIF
       IF (IFINB.EQ.(-1)) GOTO 500
       NPIV   = IW(IOLDPS+1)
       CALL CMUMPS_235(
     *            NASS,N,INODE,IW,LIW,A,LA,
     *            LDAFS, 
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4),KEEP)
      GO TO 50
 490  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
 500  CONTINUE
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE CMUMPS_141
      SUBROUTINE CMUMPS_223( NFRONT, NASS,
     *                   IBEGKJI, LPIV, TIPIV,
     *                   N, INODE, IW, LIW,
     *                   A, LA, NNEG, 
     *                   INOPV, IFLAG,
     *                   IOLDPS, POSELT, UU, 
     *                   SEUIL )
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV
      INTEGER LPIV, IBEGKJI, NNEG
      INTEGER TIPIV( LPIV )
      COMPLEX A(LA) 
      REAL UU, UULOC, SEUIL
      INTEGER IW(LIW) 
      INTEGER   IOLDPS, POSELT
       include 'mpif.h'
      INTEGER APOS, LDAFS
      REAL ZERO, RMAX
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,J1,JJ,J2,IDIAG,ILOC
      INTRINSIC MAX
      DATA ZERO /0.0E0/
        LDAFS = NASS
        UULOC = UU
        INOPV   = 0
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
        ILOC = NPIVP1 - IBEGKJI + 1
        TIPIV( ILOC ) = ILOC
        NASSW   = IABS(IW(IOLDPS+3))
        DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + LDAFS*(IPIV-1) + NPIV
            IDIAG = APOS + IPIV - NPIVP1
            IF (UULOC.EQ.ZERO) THEN 
              IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
              GO TO 420
            ENDIF
            RMAX = ABS(A(IDIAG))
            J1 = APOS
            J2 = IDIAG - 1
            DO JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
            ENDDO
            J1 = IDIAG + LDAFS
            J2 = NASS - IPIV 
            DO JJ=1,J2
              RMAX = MAX(ABS(A(J1)),RMAX)
              J1 = J1 + LDAFS
            ENDDO
            IF (ABS(A(IDIAG)).LE.MAX(UULOC*RMAX,SEUIL)) GOTO 460
            IF (IPIV.EQ.NPIVP1) GO TO 420
            TIPIV( ILOC ) = IPIV - IBEGKJI + 1
            CALL CMUMPS_319( A, LA, IW, LIW,
     *                     IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     *                     LDAFS, NFRONT, 2 )
            GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
  420 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_223
      SUBROUTINE CMUMPS_235(
     *                 NASS, N, INODE,
     *                 IW, LIW, A, LA,
     *                 LDAFS,
     *                 IOLDPS, POSELT,
     *                 LKJIB, LKJIT, KEEP )
      IMPLICIT NONE
      INTEGER NASS,N,LA,LIW
      COMPLEX    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE, KEEP(150)
      INTEGER POSELT, LDAFS
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1
      INTEGER HF
      INTEGER LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER J, Block
      INTEGER BLSIZE
      COMPLEX ONE, ALPHA
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1)
      JROW2  = IABS(IW(IOLDPS+3))
      NPBEG  = IW(IOLDPS+4)
      HF     = 6 + IW(IOLDPS+5)
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB + 2 * NONEL
          IW(IOLDPS+3)= MIN0(NPIV+LKJIB,NASS)
          LKJIB       = MIN0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3) = NASS
        ENDIF
      ELSEIF (JROW2.LT.NASS) THEN
          IW(IOLDPS+3) = MIN0(JROW2+LKJIB,NASS)
      ENDIF
      IW(IOLDPS+4) = NPIV + 1
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = MIN( BLSIZE, NASS - IROW + 1 )
          LPOS = POSELT + ( IROW - 1 ) * LDAFS + NPBEG - 1
          UPOS = POSELT + ( NPBEG - 1 ) * LDAFS + IROW - 1
          APOS =  POSELT + (IROW-1) * LDAFS + IROW - 1
          DO J=1, Block
            CALL CGEMV( 'T', LKJIW, Block - J + 1, ALPHA,
     &                  A( LPOS ), LDAFS, A( UPOS ), LDAFS,
     &                  ONE, A( APOS ), LDAFS )
            LPOS = LPOS + LDAFS
            APOS = APOS + LDAFS + 1
            UPOS = UPOS + 1
          END DO
          LPOS = POSELT + ( IROW - 1 + Block ) * LDAFS + NPBEG - 1
          UPOS = POSELT + ( NPBEG - 1 ) * LDAFS + IROW - 1
          APOS = POSELT + ( IROW - 1 + Block ) * LDAFS + IROW - 1
          CALL CGEMM( 'N','N', Block, NASS - IROW + 1 - Block, LKJIW,
     *                ALPHA, A( UPOS ), LDAFS,
     *                A( LPOS ), LDAFS, ONE, A( APOS ), LDAFS )
         END DO
        END IF
      END IF
  500 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_235
      SUBROUTINE CMUMPS_227
     *     ( NASS, N, INODE, IW, LIW,
     *     A, LA, LDAFS, 
     *     IOLDPS,POSELT,IFINB,LKJIB,LKJIT)
      IMPLICIT NONE
      INTEGER NASS,N,LA,LIW,INODE,IFINB,LKJIB,LDAFS,
     &        NPBEG
      COMPLEX    A(LA)
      INTEGER IW(LIW)
      COMPLEX    VALPIV
      INTEGER APOS, IOLDPS, POSELT, K1POS, NCB1
      INTEGER LKJIT
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,LPOS
      COMPLEX ONE, ALPHA
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
      NPIV   = IW(IOLDPS+1)
      NPIVP1 = NPIV + 1
      IFINB  = 0
      IF (IW(IOLDPS+3).LE.0) THEN
         IW(IOLDPS+3) = MIN0(NASS,LKJIB)
      ENDIF
      JROW2 = IW(IOLDPS+3)
      NPBEG = IW(IOLDPS+4)
      NEL2   = JROW2 - NPIVP1
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      APOS   = POSELT + NPIV*(LDAFS + 1)
      VALPIV = ONE/A(APOS)
      LPOS   = APOS + LDAFS
      CALL CCOPY(NASS-NPIVP1, A(LPOS), LDAFS, A(APOS+1), 1)
      CALL CMUMPS_CSYR('U', NEL2, -VALPIV, A(LPOS), LDAFS, A(LPOS+1), 
     &             LDAFS)
      CALL CSCAL(NASS-NPIVP1, VALPIV, A(LPOS), LDAFS)
      IF (NEL2.GT.0) THEN
           K1POS = LPOS + NEL2*LDAFS
           NCB1  = NASS - JROW2
           CALL CGERU(NEL2, NCB1 , ALPHA, A(APOS+1), 1, 
     &          A(K1POS), LDAFS, A(K1POS+1), LDAFS)
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_227
      RECURSIVE SUBROUTINE CMUMPS_263(
     *   COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, root, OPASSW, OPELIW, ITLOC, FILS,  
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *     )
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 150 )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, IPTRLU, 
     &        LRLU, LRLUS, LIW, LA
      INTEGER COMP
      INTEGER NIRBDU, IFLAG, IERROR, POSFAC, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28)),
     *        NSTK_S(KEEP(28))
      INTEGER NBPROCFILS(KEEP(28)), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)), ITLOC( N ), FILS( N )
      INTEGER ND( KEEP(28) ), FRERE_STEPS( KEEP(28) )
      INTEGER INTARR( MAX(1,KEEP(14)) )
      REAL OPASSW, OPELIW
      COMPLEX  DBLARR( MAX(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      INTEGER INODE, IPOSK, JPOSK, NCOLU, NPIV, POSITION, IERR
      INTEGER LAELL, POSBLOCFACTO
      INTEGER IOLDPS, POSELT, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NELIM1, NCOL_TO_SEND
      INTEGER LDA
      INTEGER NSLAV1, HS, DEST, NSLAVES_FOLLOW
      INTEGER FPERE
      INTEGER CPOS, LPOS
      INTEGER SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON, SHIFT_VAL_SON
      LOGICAL DYNAMIC
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER allocok
      COMPLEX, ALLOCATABLE, DIMENSION(:) :: UDYNAMIC
      COMPLEX ONE,ALPHA
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      DYNAMIC = .FALSE.
      POSITION  = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPOSK, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, JPOSK, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     *                 MPI_INTEGER, COMM, IERR )
      IF ( NPIV .LE. 0 ) THEN
      NPIV = - NPIV
        WRITE(*,*) MYID,':error, received negative NPIV in BLFAC'
        CALL CMUMPS_ABORT()
      END IF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOLU, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LAELL = NPIV * NCOLU
      IF ( LRLU .LT. LAELL ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          IERROR = LAELL - LRLU
          GOTO 700
        END IF
        CALL CMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     *       ,LRLU,LRLUS
             IFLAG = -9
             IERROR = LAELL -LRLU
             GOTO 700
        END IF
      END IF
      LRLU  = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL CMUMPS_471(.FALSE.,LA-LRLUS,0, LAELL,KEEP)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 A(POSBLOCFACTO), NPIV*NCOLU,
     *                 MPI_COMPLEX,
     *                 COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) DYNAMIC = .TRUE.
      IF ( (PTRIST(STEP( INODE )).NE.0) .AND.
     &    (IPOSK + NPIV -1 .GT. IW(PTRIST(STEP(INODE))+3)) )  THEN
        DYNAMIC = .TRUE.
      ENDIF
      IF (DYNAMIC)  THEN
        ALLOCATE(UDYNAMIC(LAELL), stat=allocok)
        if (allocok .GT. 0) THEN
          write(*,*) myid, ' : PB allocation U in blfac_slave '
     *     , LAELL
          IFLAG = -13 
          IERROR = LAELL
          GOTO 700
        endif
        UDYNAMIC(1:LAELL) = A(POSBLOCFACTO:POSBLOCFACTO+LAELL-1)
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL CMUMPS_471(.FALSE.,LA-LRLUS,0,-LAELL,KEEP)
      ENDIF
      DO WHILE ( PTRIST(STEP(INODE)) .EQ. 0 )
        MSGSOU = CMUMPS_275( STEP(INODE),
     *           PROCNODE_STEPS, SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL CMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, MAITRE_DESC_BANDE,
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
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      ENDDO
      DO WHILE ( IPOSK + NPIV -1 .GT.
     *                 IW( PTRIST(STEP( INODE )) + 3 ) )
        MSGSOU = CMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL CMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, BLOC_FACTO_SYM,
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
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL CMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, MPI_ANY_TAG, 
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
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IOLDPS  = PTRIST(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      LCONT1 = IW( IOLDPS )
      NROW1  = IW( IOLDPS + 2 )
      NPIV1  = IW( IOLDPS + 3 )
      NSLAV1 = IW( IOLDPS + 5 )
      HS     = 6 + NSLAV1
      NCOL1  = LCONT1 + NPIV1
      CPOS = POSELT + JPOSK - 1
      LPOS = POSELT + IPOSK - 1
      IF ( NPIV .GT. 0 ) THEN
       IF (DYNAMIC) THEN
        CALL CGEMM( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     *            UDYNAMIC(1), NPIV,
     *            A( LPOS ), NCOL1, ONE,
     *            A( CPOS ), NCOL1 )
       ELSE
        CALL CGEMM( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     *            A( POSBLOCFACTO ), NPIV,
     *            A( LPOS ), NCOL1, ONE,
     *            A( CPOS ), NCOL1 )
       ENDIF
      ENDIF
      IW( IOLDPS + 6 ) = IW( IOLDPS + 6 ) + 1
      IF (DYNAMIC) THEN
       DEALLOCATE(UDYNAMIC)
      ELSE
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL CMUMPS_471(.FALSE.,LA-LRLUS,0,-LAELL,KEEP)
      ENDIF
      NSLAVES_FOLLOW = IW( IOLDPS + 5 ) - 2
      IF ( IW( IOLDPS + 6 ) .eq. 0 .and.
     *  KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 1 ) THEN
         DEST = CMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
         CALL CMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     *                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_BLFAC_SLAVE.'
           IFLAG = -99
           GOTO 700
         END IF
      END IF
      IF ( IW( IOLDPS + 6 ) .eq. 0 .AND. (KEEP(38).EQ.FPERE) ) THEN
       NPIV1  = IW( IOLDPS + 3 )
       NASS1  = IW(IOLDPS+4)
       NELIM1  = NASS1-NPIV1
       NCOL_TO_SEND =  LCONT1-NELIM1
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5) 
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW1 + NASS1
       SHIFT_VAL_SON      = NASS1
       LDA                = LCONT1 + NPIV1
       CALL CMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE, 
     *    PTRIST, PTRAST, 
     *    root, NROW1, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA, 
     *    ROOT_CONT_STATIC, MYID, COMM,
     *    
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, .FALSE., ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
       IF ( IFLAG < 0 ) GOTO 600
       IF (NELIM1.EQ.0) THEN
        CALL CMUMPS_314( N, INODE,  
     *    PTRIST, PTRAST, PTLUST_S, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP)
       ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_263
      SUBROUTINE CMUMPS_179(MYID, N, STEP, FRERE, FILS,
     *     NA, LNA, NE, ND, PROCNODE, SLAVEF,
     *     NRLADU, NIRADU, NIRNEC, NRLNEC, NRLBDU, NIRBDU,
     *     MAXFR, OPSA, PERLU, 
     *     KEEP, LOCAL_M, LOCAL_N,
     *     SBUF_SEND, SBUF_REC, OPS_SUBTREE, NSTEPS,
     *     NMB_PAR2, ISTEP_TO_INIV2, CANDIDATES, 
     *     IFLAG, IERROR)
      IMPLICIT NONE
      INTEGER  MYID, N, LNA, IFLAG, IERROR
      INTEGER  NRLADU, NIRADU, NRLNEC, NIRNEC
      INTEGER  MAXFR, PERLU, NSTEPS
      INTEGER STEP(N)
      INTEGER FRERE(NSTEPS), FILS(N), NA(LNA), NE(NSTEPS),
     *        ND(NSTEPS), PROCNODE(NSTEPS)
      INTEGER  SLAVEF, KEEP(150), LOCAL_M, LOCAL_N
      INTEGER  SBUF_SEND, SBUF_REC
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      INTEGER  CANDIDATES( SLAVEF+1, MAX(NMB_PAR2,1) )
      INTEGER  NRLBDU, NIRBDU, IDUMMY
      REAL OPSA
      REAL OPS_SUBTREE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSTKR, TNSTK, IPOOL, 
     *                                      LSTKI 
      INTEGER I,ISTKR,ITOP,NELIM,NFR
      INTEGER ISTKI, ISTKIM, STKI, INDX
      INTEGER K,LSTK,NSTK, IFATH
      INTEGER INODE, LEAF, NBLEAF, NBROOT, III, IN
      INTEGER LEVEL, MAXTEMPCB, MAXITEMPCB
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF 
      INTEGER LEVELF, NCB, SIZECB, SIZECBI
      INTEGER NBROWMAX, NSLAVES, TEMP,
     *         NELIMF, NFRF, NCBF, CBMAXR, CBMAXS,
     *         NSLAVESF, NBROWMAXF, NCOLMAX, LKJIB,
     *         LKJIBT, KMAX
      INTEGER LWK_QR, LIWK_QR, SIZE_ROOT
      INTEGER IROOT
      INTEGER allocok
      LOGICAL ROOT_OWNER
      REAL OPS_NODE
       INTEGER WHAT
      INTRINSIC MIN, INT, REAL
      INTEGER CMUMPS_275, CMUMPS_330, CMUMPS_497
      INTEGER CMUMPS_52
      EXTERNAL CMUMPS_503, CMUMPS_52
      LOGICAL CMUMPS_170, CMUMPS_358
      EXTERNAL CMUMPS_275, CMUMPS_330, 
     *         CMUMPS_170, CMUMPS_358
      logical :: force_cand , concerned, UPDATES, STACKCB, MASTERSON
      integer :: nmb_of_cand,istat, IFSON, LEVELSON
      integer, allocatable :: cand_of_node(:)
      ALLOCATE( LSTKR(NSTEPS), TNSTK(NSTEPS), IPOOL(NSTEPS),
     *          LSTKI(NSTEPS) , stat=allocok)
      if (allocok .GT. 0) THEN
        IFLAG  =-13
        IERROR = 4*NSTEPS
        RETURN
      endif
      IDUMMY = 0
      LKJIB = MAX(KEEP(5),KEEP(6))
      TNSTK = NE
      LEAF = NA(1)+1
      IPOOL(1:LEAF-1) = NA(3:3+LEAF-2)
      NBROOT = NA(2)
      ISTKR  = 0
      ISTKI  = 0
      ISTKIM = 0
      OPSA   = 0.0E0
      OPS_SUBTREE = 0.0E0
      NRLADU = 0
      NIRADU = 0
      NRLNEC = 0
      NIRNEC = 0
      MAXFR  = 0
      ITOP = 0
      MAXTEMPCB  = 0
      MAXITEMPCB = 0
      SBUF_SEND  = 1
      SBUF_REC   = 1
      IF (KEEP(38) .NE. 0) THEN
        INODE  = KEEP(38)
        NRLADU = LOCAL_M*LOCAL_N
        NRLNEC = NRLADU
        IF (CMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *                                       .EQ. MYID) THEN
          NIRADU = 6+2*ND(STEP(INODE))
        ELSE
          NIRADU = 6
        ENDIF
        NIRNEC = NIRADU
      ENDIF
      IF((KEEP(24).eq.0).OR.(KEEP(24).eq.1)) THEN
         force_cand=.FALSE.           
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      END IF
      III = 1
 90   CONTINUE
        IF (LEAF .NE. 1) THEN
           LEAF = LEAF - 1
           INODE = IPOOL(LEAF)
        ELSE 
           WRITE(MYID+6,*) ' ERROR 1 in file mpima41md '
           CALL CMUMPS_ABORT()
        ENDIF
 95     CONTINUE 
        MASTER = CMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *           .EQ. MYID
        LEVEL  = CMUMPS_330(STEP(INODE),PROCNODE,SLAVEF)
      UPDATE=.FALSE.
      if(.NOT.force_cand) then
         UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
      else
         if(master.and.(level.ne.3)) then
            UPDATE = .TRUE.
         else if(level.eq.2) then
            if ( CMUMPS_358( MYID, SLAVEF, inode,
     *           NMB_PAR2, ISTEP_TO_INIV2, KEEP(71), STEP, N,
     *           CANDIDATES, KEEP(24) ) ) then
              UPDATE = .TRUE.
            end if
         end if
      end if
        NFR    = ND(STEP(INODE))
        NSTK   = NE(STEP(INODE))
        NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
        IFSON = -IN
        IF (LEVEL.EQ.3) THEN
         NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+
     *                 LOCAL_M*LOCAL_N)
         IF (MASTER) THEN 
          IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
        ENDIF
        NCB      = NFR-NELIM
        SIZECB   =  NCB*NCB
        SIZECBI  = 2* NCB  + 6 + 3
        IF (LEVEL.EQ.2) THEN
          IF (MASTER) THEN
             SBUF_SEND = MAX(SBUF_SEND, NELIM*LKJIB+NELIM+6)
          ELSE if (UPDATE) then
             SBUF_REC = MAX( SBUF_REC, NELIM*LKJIB+NELIM+6 )
             IF (KEEP(50).EQ.1) THEN
               LKJIBT  = LKJIB
             ELSE
               LKJIBT  = NELIM
             ENDIF
             WHAT = 1 ! Compute NBROWMAX only
             CALL CMUMPS_503(WHAT, KEEP,
     &       NCB, NFR, SLAVEF, NBROWMAX, IDUMMY)
             SBUF_SEND = MAX(SBUF_SEND,
     *                       LKJIBT*NBROWMAX+6)
             SBUF_REC = MAX( SBUF_REC, NBROWMAX*LKJIBT+6 )
          ENDIF
        ENDIF
        IN = INODE
 660    IN = FRERE(STEP(IN))
        IF (IN.GT.0) GO TO 660
        IFATH = -IN
        IF ( UPDATE ) THEN
          IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
            SIZECB = NCB*NCB
             NIRADU = NIRADU + 6 + 2*NFR
             NRLADU = NRLADU + NELIM*NFR
            SIZECBI    = 2* NCB  + 6 + 3
          ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
             NIRADU = NIRADU + 6 + SLAVEF - 1 + 2*NFR
             NRLADU = NRLADU + NELIM*NELIM
             SIZECB = 0
             SIZECBI    = NCB + 5 +  SLAVEF - 1
            ELSE
             WHAT = 2 ! Both NBROWMAX and SIZECB
             CALL CMUMPS_503(WHAT, KEEP,
     &            NCB, NFR, SLAVEF, NBROWMAX, SIZECB )
             NIRADU   = NIRADU+4+NELIM+NBROWMAX 
             NRLADU   = NRLADU + NELIM*NBROWMAX
             NSLAVES = CMUMPS_52( SLAVEF,
     *       KEEP(48), KEEP(10), KEEP(50), NFR, NCB)
             SIZECBI      = 4 + NBROWMAX + NCB + NSLAVES + 1
            ENDIF
          ENDIF
          NIRNEC = MAX0(NIRNEC,NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
          NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+SIZECB+
     *                  MAX0(SIZECB,MAXTEMPCB) )
          IF (NFR.GT.MAXFR) MAXFR = NFR
          IF (NSTK.GT.0) THEN
           DO 70 K=1,NSTK
            LSTK = LSTKR(ITOP)
            ISTKR = ISTKR - LSTK
            STKI = LSTKI( ITOP )
            ISTKI = ISTKI - STKI
            ITOP = ITOP - 1
            IF (ITOP.LT.0) THEN
             write(MYID+6,*) ' ERROR 2 in file mpima41md'
             CALL CMUMPS_ABORT()
            ENDIF
   70      CONTINUE
          ENDIF  
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = CMUMPS_275(STEP(IFSON),PROCNODE,SLAVEF)
     *                  .EQ.MYID
            LEVELSON  = CMUMPS_330(STEP(IFSON),PROCNODE,SLAVEF)
            if(.NOT.force_cand) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     *                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( CMUMPS_358( MYID, SLAVEF, IFSON,
     *                 NMB_PAR2, ISTEP_TO_INIV2, KEEP(71), STEP, N,
     *                 CANDIDATES, KEEP(24)) ) then
                    UPDATES = .TRUE.
                  end if
               end if
            end if
            IF (UPDATES) THEN
              LSTK = LSTKR(ITOP)
              ISTKR = ISTKR - LSTK
              STKI = LSTKI( ITOP )
              ISTKI = ISTKI - STKI
              ITOP = ITOP - 1
              IF (ITOP.LT.0) THEN
                write(*,*) MYID,': ERROR 2 in mpima41md. ITOP = ',ITOP
               CALL CMUMPS_ABORT()
              ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) ! process next son
         END DO
        ENDIF
        IF (MASTER.OR.LEVEL.EQ.3) THEN
         IF (INODE.NE.KEEP(20).OR.(KEEP(60).EQ.0)) THEN
          OPS_NODE = REAL(NELIM) * (
     *        REAL( NFR * NFR ) + REAL( NFR )
     *      - REAL( NFR*NELIM + NELIM + 1 )   )
     *      + ( REAL(NELIM)*REAL(NELIM+1)*REAL(2*NELIM+1) )
     *      / 6 
         ELSE
           OPS_NODE = 0.0E0
         ENDIF
         IF (CMUMPS_170(STEP(INODE),
     *   PROCNODE, SLAVEF) ) THEN
           OPS_SUBTREE = OPS_SUBTREE + OPS_NODE
         ENDIF
         IF ( LEVEL .eq. 3 ) THEN
            OPSA = OPSA + OPS_NODE / REAL( SLAVEF )
         ELSE IF (MASTER) THEN
            OPSA = OPSA + OPS_NODE
         ENDIF
        ENDIF
        IF (IFATH .EQ. 0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
        ELSE
         LEVELF = CMUMPS_330(STEP(IFATH),PROCNODE,SLAVEF)
         MASTERF= CMUMPS_275(STEP(IFATH),PROCNODE,SLAVEF)
     *            .EQ. MYID
         UPDATEF=.FALSE.
         if(.NOT.force_cand) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(masterf.and.(levelf.ne.3)) then
               UPDATEF = .TRUE.
            else if(levelf.eq.2) then
               if ( CMUMPS_358( MYID, SLAVEF, ifath,
     *              NMB_PAR2, ISTEP_TO_INIV2, KEEP(71), STEP, N,
     *              CANDIDATES, KEEP(24) )) THEN
                 UPDATEF = .TRUE.
               end if
            end if
         end if
         concerned  = UPDATEF .OR. UPDATE
         NCB    = NFR-NELIM
         IF (LEVELF.EQ.3) THEN
           TEMP       = MIN(SIZECB, LOCAL_M*LOCAL_N) 
           MAXTEMPCB  = MAX(MAXTEMPCB, TEMP)
           MAXITEMPCB = MAX(MAXITEMPCB,SIZECBI)
           TEMP       = MIN(SIZECB, LOCAL_M*LOCAL_N) + SIZECBI
           SBUF_REC   = MAX(SBUF_REC, TEMP)
         ENDIF
        IF (concerned) THEN
         IN = IFATH
         NELIMF = 0
         DO WHILE (IN.GT.0)
            IN = FILS(IN)
            NELIMF = NELIMF+1
         ENDDO
         NFRF = ND(STEP(IFATH))
         NCBF = NFRF - NELIMF
         IF (LEVELF.EQ.2) THEN
           WHAT = 1 ! NBROWMAX only
           CALL CMUMPS_503(WHAT,
     *     KEEP, NCBF, NFRF, SLAVEF, NBROWMAXF,IDUMMY)
           NCOLMAX  = NCB
           if (update.or.updatef) then 
              CBMAXS = MAX(MIN(NELIMF,NCB)*NCOLMAX, 
     *          MIN(NBROWMAXF,NCB)*NCOLMAX)
           else
              CBMAXS = 0
           end if
           IF (MASTERF) THEN
            CBMAXR = MIN(NELIMF,NCB)*NCOLMAX
           ELSE IF (UPDATEF) THEN
            CBMAXR = MIN(NBROWMAXF,NCB)*NCOLMAX
           else
            cbmaxr = 0
           ENDIF
         ELSEIF (LEVELF.EQ.3) THEN
           CBMAXR = LOCAL_M*LOCAL_N
           CBMAXS = LOCAL_M*LOCAL_N
         ELSE
           CBMAXR = NFRF*NFRF
           CBMAXS = NFRF*NFRF
         ENDIF
         IF (UPDATE) THEN
           SIZECB = MIN(CBMAXS, SIZECB)
           IF ( .not. ( LEVELF .eq. 1 .AND. UPDATEF ) )
     *     SBUF_SEND = MAX(SBUF_SEND, SIZECB+SIZECBI)
         ENDIF
         STACKCB = .FALSE.
         IF (UPDATEF) THEN 
          STACKCB = .TRUE.
          SIZECBI = 2 * NFR + 6 
          IF (LEVEL.EQ.1) THEN
             IF (MASTER) THEN
               SIZECB  = NCB*NCB
               SIZECBI = 2
             ELSE IF (LEVELF.EQ.1) THEN
              SIZECBI = 2 * NCB + 6 + 3
              SIZECB  = MIN(CBMAXR,NCB*NCB)
              SBUF_REC = MAX(SBUF_REC, SIZECBI+SIZECB)
             ELSE 
              SIZECB = MIN(NCB*NCB,CBMAXR)
              MAXTEMPCB  = MAX(MAXTEMPCB, SIZECB)
              MAXITEMPCB = MAX(MAXITEMPCB, 2 * NCB + 6 + 3)
              SBUF_REC = MAX(SBUF_REC, SIZECB +2 * NCB + 6 + 3)
              SIZECBI= 0
              SIZECB = 0
             ENDIF
          ELSE 
             WHAT = 2 ! Both NBROWMAX and SIZECB max are computed
             CALL CMUMPS_503( WHAT, KEEP,
     &       NCB, NFR, SLAVEF, NBROWMAX, SIZECB )
             MAXTEMPCB  = MAX(MAXTEMPCB, SIZECB)
             MAXITEMPCB = MAX(MAXITEMPCB, NBROWMAX + NCB + 6 )
             IF (MASTER) THEN
              SIZECBI =  NCB + 5 +  SLAVEF - 1
              SIZECB  = 0
             ELSE IF (UPDATE) THEN
              SIZECBI =  NFR + 6 + SLAVEF - 1 
              NSLAVES = CMUMPS_52(
     &            SLAVEF, KEEP(48), KEEP(10), KEEP(50),
     &            NFR, NCB)
              SIZECBI = SIZECBI + NBROWMAX + NFR + 6 + NSLAVES + 1
             ELSE
              SIZECB  = 0 
              SIZECBI = 0
             ENDIF
          ENDIF
         ELSE
           IF (LEVELF.NE.3) THEN
               STACKCB = .TRUE.
               SIZECB  = 0 
               SIZECBI = 0
               IF ( (LEVEL.EQ.1) .AND. (LEVELF.NE.1) ) THEN
                  SIZECB  = NCB*NCB
                  SIZECBI = 2 * NCB + 6 + 3
               ELSE IF (LEVEL.EQ.2) THEN
                IF (MASTER) THEN
                 SIZECBI =  NCB + 5 +  SLAVEF - 1
                ELSE 
                 WHAT = 2 ! Both SIZECB and NBROWMAX computed
                 CALL CMUMPS_503(WHAT, KEEP,
     &           NCB, NFR, SLAVEF, NBROWMAX, SIZECB )
                 SIZECBI =  NFR + 6 + SLAVEF - 1
                 SIZECBI = SIZECBI + NBROWMAX + NFR + 6
                ENDIF 
               ENDIF
           ENDIF
         ENDIF
         IF (STACKCB) THEN
          IF (FRERE(STEP(INODE)).EQ.0) THEN
            write(*,*) ' ERROR 3 in LDLT_mpima41md'
            CALL CMUMPS_ABORT()
          ENDIF
           ITOP = ITOP + 1
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI+LSTKI(ITOP)
           ISTKIM=MAX( ISTKI, ISTKIM )
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
           NIRNEC = MAX0(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
         ENDIF ! end of (stackcb)
        ENDIF ! end of (concerned) 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
        ENDIF
  115 CONTINUE
      ISTKIM = ISTKIM / 2 + MAXFR + N
      IF (KEEP(38) .ne. 0 ) THEN
        IF (CMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *                                       .EQ. MYID) THEN
          ISTKIM = ISTKIM + ND(STEP( KEEP(38) ))
        END IF
      END IF
      IF ( NRLNEC - NRLADU .LT. ISTKIM ) THEN
        NRLNEC = NRLADU + ISTKIM
      END IF
      IF ( NIRNEC - NIRADU .LT. ISTKIM ) THEN
        NIRNEC = NIRADU + ISTKIM
      END IF
      IF ( KEEP(53) .NE. 0 ) THEN
        IF ( KEEP(38) .ne. 0 ) THEN
          IROOT = KEEP( 38 )
        ELSE
          IROOT = KEEP( 20 )
        END IF
        ROOT_OWNER = ( MYID .eq.
     *  CMUMPS_275( STEP(IROOT), PROCNODE, SLAVEF ) )
        SIZE_ROOT = ND(STEP(IROOT))
        CALL CMUMPS_279( KEEP(51), KEEP(51), SIZE_ROOT,
     *       LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP, LIWK_QR, LWK_QR )
        IF ( NRLNEC-NRLADU .LT. LWK_QR ) THEN
          NRLNEC = NRLADU + LWK_QR
        END IF
        IF ( NIRNEC-NIRADU .LT. LIWK_QR ) THEN
          NIRNEC = NIRNEC + LIWK_QR
        END IF
      END IF
      SBUF_REC = MAX(SBUF_REC, MAXTEMPCB+MAXITEMPCB )
      SBUF_REC = SBUF_REC   + 10
      SBUF_SEND = SBUF_SEND + 10
      NRLBDU = NRLNEC + PERLU*(NRLNEC/100+1) 
      IF (SLAVEF.EQ.1) THEN 
         SBUF_REC = 1
         SBUF_SEND= 1
      ELSE
       IF (KEEP(50).EQ.1) THEN
        NRLBDU = NRLBDU + SBUF_REC+
     &      MAX(PERLU,0)*(SBUF_REC/100)
       ELSE
        NRLBDU = NRLBDU + SBUF_REC+4*
     &     MAX(PERLU,0) *(SBUF_REC/100)
       ENDIF
      ENDIF
      NIRBDU = NIRNEC + 2 * MAX(PERLU,10)*(NIRNEC/100+1)
      DEALLOCATE( LSTKR, TNSTK, IPOOL,
     *          LSTKI )
      RETURN
      END SUBROUTINE CMUMPS_179
      SUBROUTINE CMUMPS_137( INODE, N, PROCNODE_STEPS,
     *           SLAVEF,
     &           ND, FILS, FRERE_STEPS, STEP, PIMASTER,
     *           KEEP28, KEEP50,
     *           FLOP1,
     &           IW, LIW )
      IMPLICIT NONE
      INTEGER INODE, N, KEEP50, LIW, SLAVEF, KEEP28
      INTEGER PROCNODE_STEPS(KEEP28), ND(KEEP28),
     &        FILS(N), FRERE_STEPS(KEEP28),
     *        STEP(N), 
     * PIMASTER(KEEP28),
     *  IW( LIW )
      REAL FLOP1
      INTEGER NUMORG, IN, NASS, IFSON, NUMSTK, NFRONT, NPIV, NCB, 
     &        LEVEL, ISON
      LOGICAL CMUMPS_170
      INTEGER CMUMPS_330
      EXTERNAL CMUMPS_170, CMUMPS_330
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      FLOP1 = 0.0E0
      IF (CMUMPS_170(STEP(INODE),
     *    PROCNODE_STEPS, SLAVEF) ) RETURN
      IN     = INODE
      NUMORG = 0
   10 NUMORG = NUMORG + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 10
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .EQ. 0) GOTO 30
   20 NUMSTK = NUMSTK + 1
      NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1)
      ISON = FRERE_STEPS(STEP(ISON))
      IF (ISON .GT. 0) GOTO 20
   30 NFRONT = ND(STEP(INODE)) + NASS
      NPIV  = NASS + NUMORG
      NCB   = NFRONT - NPIV
      LEVEL = CMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF (KEEP50.EQ.0) THEN
        IF (LEVEL.EQ.1) THEN
          FLOP1 = REAL(2) * REAL(NFRONT) * REAL(NPIV) *
     *      REAL(NFRONT - NPIV - 1) +
     *      REAL(NPIV * (NPIV + 1)) * REAL(2 * NPIV + 1) / REAL(3)
          FLOP1 = FLOP1 + REAL(2 * NFRONT - NPIV - 1) * REAL(NPIV) /
     *      REAL(2)
        ELSEIF (LEVEL.EQ.2) THEN
          FLOP1 = REAL(2*NPIV)*REAL(NFRONT) -
     *          REAL(NPIV+NFRONT)*REAL(NPIV+1)
          FLOP1 = REAL(NPIV)*FLOP1 +
     *     REAL(((2 * NPIV - NPIV - 1) * NPIV) / 2) +
     *     REAL(NPIV) * REAL(NPIV + 1) * REAL(2 * NPIV + 1) /REAL(3)
        ENDIF
      ELSE
        IF (LEVEL.EQ.1) THEN
          FLOP1 = REAL(NPIV) * (
     *          REAL( NFRONT ) * REAL( NFRONT ) + REAL( NFRONT )
     *        - REAL( NFRONT*NPIV + NPIV + 1 )    )
     *        + ( REAL(NPIV)*REAL(NPIV+1)*REAL(2*NPIV+1)  )
     *        / REAL(6)
        ELSE
          FLOP1 = REAL(NPIV) * (
     *          REAL( NPIV ) * REAL( NPIV ) + REAL( NPIV )
     *        - REAL( NPIV*NPIV +  NPIV + 1 )    )
     *        + ( REAL(NPIV)*REAL(NPIV+1)*REAL(2*NPIV+1) )
     *        / REAL( 6 )
        ENDIF
      ENDIF
      RETURN
      END
