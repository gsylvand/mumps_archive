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
      RECURSIVE SUBROUTINE CMUMPS_323
     *         ( BUFR, LBUFR, LBUFR_BYTES,
     *           MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *           N, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *           IWCB, LIWCB,
     *           WCB, LWCB, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB,
     *           INFO, KEEP, STEP, PROCNODE_STEPS, W1 )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIW, LA
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER INFO( 40 ), KEEP( 150 )
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ),  NSTK_S( N )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      INTEGER PTRICB(KEEP(28)),PTRIST(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      COMPLEX WCB( LWCB ), A( LA ), W1( N )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER FINODE, FPERE, LONG, NCB, POSITION, NCV, NPIV
      INTEGER PTRX, PTRY, APOS, PDEST, I
      LOGICAL FLAG
      EXTERNAL CMUMPS_275
      INTEGER  CMUMPS_275
      COMPLEX ALPHA, ONE
      PARAMETER( ONE = 1.0E0, ALPHA = -1.0E0 )
      IF ( MSGTAG .EQ. RACINE_SOLVE ) THEN
          NBFIN = NBFIN - 1
          IF ( NBFIN .eq. 0 ) GOTO 270
      ELSE  IF (MSGTAG .EQ. ContVec ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NCB, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         LONG, 1, MPI_INTEGER, COMM, IERR )
          IF ( NCB .eq. 0 ) THEN
             PTRICB(STEP(FINODE)) = -1
             NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
               IPOOL( LEAF ) = FPERE
               LEAF = LEAF + 1
               IF ( LEAF > LPOOL ) THEN
                 WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                 CALL CMUMPS_ABORT()
               END IF
             END IF
          ELSE
          IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
             PTRICB(STEP(FINODE)) = NCB + 1
          END IF
          IF ( ( POSIWCB - LONG ) .LT. 0 ) THEN
            INFO( 1 ) = -14
            INFO( 2 ) = LONG
            GOTO 260
          END IF
          IF ( POSWCB - PLEFTWCB + 1 .LT. LONG ) THEN
            INFO( 1 ) = -11
            INFO( 2 ) = PLEFTWCB - POSWCB - 1 + LONG
            GOTO 260
          END IF
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          IWCB( 1 ),
     *          LONG, MPI_INTEGER, COMM, IERR )
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *               WCB( PLEFTWCB ),
     *               LONG, MPI_COMPLEX, COMM, IERR )
            DO I = 1, LONG
              W1(IWCB(I)) = W1(IWCB(I)) + WCB( PLEFTWCB + I - 1)
            ENDDO
            PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - LONG
          ENDIF
          IF ( PTRICB(STEP(FINODE)) == 1 ) THEN
               NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
          END IF
          IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
            IPOOL( LEAF ) = FPERE
            LEAF = LEAF + 1
            IF ( LEAF > LPOOL ) THEN
              WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
              CALL CMUMPS_ABORT()
            END IF
          ENDIF
        END IF
        ELSEIF ( MSGTAG .EQ. MASTER2SLAVE ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NCV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NPIV, 1, MPI_INTEGER, COMM, IERR )
          PTRY = PLEFTWCB
          PTRX = PLEFTWCB + NCV
          PLEFTWCB = PLEFTWCB + NPIV + NCV
          IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
            INFO(1) = -11
            INFO(2) = -POSWCB + PLEFTWCB -1
            GO TO 260
          END IF
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         WCB( PTRY ), NCV, MPI_COMPLEX, COMM, IERR )
          IF ( NPIV .GT. 0 ) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         WCB( PTRX ), NPIV, MPI_COMPLEX, COMM, IERR )
          END IF
          APOS = IW( PTRIST( STEP(FINODE)) + 3 )
          CALL CGEMV( 'T', NPIV, NCV, ALPHA, A( APOS ), NPIV,
     *                WCB( PTRX ), 1, ONE, WCB( PTRY ), 1 )
          PLEFTWCB = PLEFTWCB - NPIV
          PDEST = CMUMPS_275( STEP(FPERE),
     *            PROCNODE_STEPS, SLAVEF )
          IF ( PDEST .EQ. MYID ) THEN
            IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
              NCB = IW( PTRIST(STEP(FINODE)) + 2 )
              PTRICB(STEP(FINODE)) = NCB + 1
            END IF
            DO I = 1, NCV
              W1(IW(PTRIST(STEP(FINODE))+3+I))=
     *        W1(IW(PTRIST(STEP(FINODE))+3+I)) + WCB(PTRY + I - 1)
            END DO
            PTRICB(STEP(FINODE)) =
     *      PTRICB(STEP(FINODE)) - NCV
            IF ( PTRICB( STEP( FINODE ) ) == 1 ) THEN
                 NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
            END IF
            IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
              IPOOL( LEAF ) = FPERE
              LEAF = LEAF + 1
              IF ( LEAF > LPOOL ) THEN
                WRITE(*,*) 'INTERNAL Error 41r: Pool is too small.'
                CALL CMUMPS_ABORT()
              END IF
            ENDIF
          ELSE
 210       CONTINUE
           CALL CMUMPS_78( FINODE, FPERE,
     *       IW(PTRIST(STEP( FINODE )) + 2 ), NCV,
     *       IW(PTRIST(STEP(FINODE))+4),
     *       WCB( PTRY ), PDEST, ContVec, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
              CALL CMUMPS_303( .FALSE., FLAG,
     *               BUFR, LBUFR, LBUFR_BYTES,
     *               MYID, SLAVEF, COMM,
     *               N, IPOOL, LPOOL, III, LEAF,
     *               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *               IWCB, LIWCB,
     *               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *               PTRICB, INFO, KEEP, STEP,
     *               PROCNODE_STEPS, W1 )
              IF ( INFO( 1 )  .LT. 0 )  GOTO 270
              GOTO 210
            ELSE IF ( IERR .EQ. -2 ) THEN
              INFO( 1 ) = -17
              INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     *                      NCV * KEEP( 35 )
              GOTO 260
            END IF
          END IF
          PLEFTWCB = PLEFTWCB - NCV
        ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GOTO 270
        ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     *      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
        ELSE
          INFO(1)=-100
          INFO(2)=MSGTAG
          GO TO 260
        ENDIF
        GO TO 270
 260    CONTINUE
        CALL CMUMPS_44( MYID, SLAVEF, COMM )
 270    CONTINUE
        RETURN
        END SUBROUTINE CMUMPS_323
      SUBROUTINE CMUMPS_302( INODE,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *           N, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S,
     *           IWCB, LIWCB,
     *           WCB, LWCB, A, LA, IW, LIW, RHS, W1, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB, PTRIST, PROCNODE_STEPS,
     *           FILS, STEP, FRERE,
     *           MYROOT,
     *           INFO, KEEP, RHS_ROOT, MTYPE,
     *
     *           ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *           )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER INODE, LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIWCB, LWCB, LIW, LA, POSWCB, PLEFTWCB, POSIWCB
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER MYROOT
      INTEGER INFO( 40 ), KEEP( 150 )
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ), NSTK_S(KEEP(28))
      INTEGER IWCB( LIWCB ), IW( LIW )
      COMPLEX WCB( LWCB ), A( LA ), RHS( N ), W1( N ),
     *                 RHS_ROOT( * )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER FILS( N ), STEP( N ), FRERE(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      EXTERNAL CGEMV, CTRSV, CMUMPS_275
      INTEGER CMUMPS_275
      COMPLEX ALPHA,ONE,ZERO
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0, ZERO=0.0E0)
      INTEGER I, J, IPOS, NSLAVES, J1, J2, J3, FPERE, NPIV, NCB, IERR,
     *        APOS, APOS1, IF, IFR, IPOSCB, APOSCB, LIELL, IN, JJ,
     *        NELIM, PLEFT
      INTEGER Effective_CB_SIZE, NUPDATE, ISLAVE, PDEST, FirstIndex
      LOGICAL FLAG
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY( 1 )
      IF ( INODE .eq. KEEP(38 ) .OR. INODE .eq.KEEP( 20 ) ) THEN
        LIELL = IW( PTRIST( STEP(INODE)) + 3 )
        NPIV  = LIELL
        NELIM = 0
        NSLAVES = 0
        IPOS = PTRIST( STEP(INODE)) + 5
      ELSE
        IPOS = PTRIST(STEP(INODE)) + 2
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        NELIM = IW(IPOS-1)
        NSLAVES = IW( PTRIST(STEP(INODE)) +5 )
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        APOS = IW(IPOS)
        NSLAVES = IW( PTRIST(STEP(INODE)) +5 )
        IPOS = IPOS + 1 + NSLAVES
      END IF
      IF ( MTYPE .EQ. 1 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + LIELL
        J3 = IPOS + NPIV
      ELSE
        J1 = IPOS + LIELL + 1
        J2 = IPOS + 2 * LIELL
        J3 = IPOS + LIELL + NPIV
      END IF
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq. KEEP(20) ) THEN
        IFR = 0
        DO JJ = J1, J3
          J = IW( JJ )
          IFR = IFR + 1
          RHS_ROOT( IFR ) = RHS( J ) + W1( J )
          W1( J ) = ZERO
        END DO
        IF ( NPIV .LT. LIELL ) THEN
          WRITE(*,*) ' Internal error in SOLVE_NODE for Root node'
          CALL CMUMPS_ABORT()
        END IF
        MYROOT = MYROOT - 1
        IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
              DUMMY (1) = 1
              CALL CMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     *             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
        END IF
        GO TO 270
      END IF
      PLEFT    = PLEFTWCB
      PLEFTWCB = PLEFTWCB + LIELL
      IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
        INFO(1) = -11
        INFO(2) = PLEFTWCB - POSWCB - 1
        GO TO 260
      END IF
      IFR = PLEFT - 1
      DO 130 JJ = J1, J3
        J = IW(JJ)
        IFR = IFR + 1
        WCB(IFR) = RHS(J) + W1(J)
        W1(J) = ZERO
  130 CONTINUE
      IF (NPIV .LT. LIELL) THEN
        DO 140 JJ = J3 + 1, J2
          J = IW(JJ)
          IFR = IFR + 1
          WCB(IFR) = W1(J)
          W1(J) = ZERO
  140   CONTINUE
      ENDIF
      IF ( NPIV .NE. 0 ) THEN
        IF ( KEEP(50).NE.0) THEN
          CALL CTRSV( 'U', 'T', 'U', NPIV, A(APOS), NPIV,
     *             WCB(PLEFT), 1 )
        ELSE
          IF ( MTYPE .eq. 1 ) THEN
                   CALL CTRSV( 'U', 'T', 'U', NPIV, A(APOS), LIELL, 
     *             WCB(PLEFT), 1 )
          ELSE
                   CALL CTRSV( 'L', 'N', 'N', NPIV, A(APOS), LIELL,
     *             WCB(PLEFT), 1 )
          END IF
        END IF
      END IF
      NCB   = LIELL - NPIV
      IF ( MTYPE .EQ. 1 ) THEN
        IF ( KEEP(50) .eq. 0 ) THEN
          APOS1 = APOS  + NPIV * LIELL
        ELSE
          APOS1 = APOS + NPIV * NPIV
        END IF
        IF ( NSLAVES .EQ. 0 .OR. NPIV .eq. 0 ) THEN
          NUPDATE = NCB
        ELSE
          NUPDATE = NELIM
        END IF
      ELSE
        APOS1 = APOS + NPIV
        NUPDATE = NCB
      END IF
      IF ( NPIV .NE. 0 .AND. NUPDATE.NE.0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
                CALL CGEMV('T', NPIV, NUPDATE, ALPHA, A(APOS1),
     *         NPIV,  WCB(PLEFT), 1, ONE, WCB(PLEFT + NPIV), 1)
        ELSE
                CALL CGEMV('N',NUPDATE, NPIV, ALPHA, A(APOS1),
     *        LIELL, WCB(PLEFT), 1, ONE, WCB(PLEFT + NPIV ), 1 )
        END IF
      END IF
      IFR = PLEFT - 1
      IF ( KEEP(50) .eq. 0 ) THEN
        DO 150 JJ = J1, J3
          J = IW(JJ)
          IFR = IFR + 1
          RHS(J) = WCB(IFR)
  150   CONTINUE
      ELSE
        APOS1 = APOS
        DO 155 JJ = J1, J3
          J = IW( JJ )
          IFR = IFR + 1
          RHS( J ) = WCB( IFR ) / A( APOS1 )
          APOS1 = APOS1 + NPIV + 1
  155   CONTINUE
      END IF
      IN = INODE
      DO WHILE ( IN .GT. 0 )
        IN = FRERE(STEP( IN ))
      END DO
      FPERE = -IN
      IF ( FPERE .EQ. 0 ) THEN
          MYROOT = MYROOT - 1
          PLEFTWCB = PLEFTWCB - LIELL
          IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
              DUMMY (1) = 1
              CALL CMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     *             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
          END IF
          GO TO 270
      ENDIF
      IF ( NUPDATE .NE. 0 .OR. NCB.eq.0 ) THEN
        IF (CMUMPS_275(STEP(FPERE),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
         IF ( NCB .ne. 0 ) THEN
          PTRICB(STEP(INODE)) = NCB + 1
          DO 190 I = 1, NUPDATE
            W1( IW(J3 + I) ) = W1( IW(J3 + I) ) + WCB(PLEFT-1+NPIV+I)
  190     CONTINUE
          PTRICB(STEP( INODE )) = PTRICB(STEP( INODE )) - NUPDATE
          IF ( PTRICB(STEP(INODE)) == 1 ) THEN
            NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
            IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
              IPOOL( LEAF ) = FPERE
              LEAF = LEAF + 1
            ENDIF
          END IF
         ELSE
          PTRICB(STEP( INODE )) = -1
          NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
          IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
            IPOOL( LEAF ) = FPERE
            LEAF = LEAF + 1
          ENDIF
         ENDIF
        ELSE
 210      CONTINUE
          CALL CMUMPS_78( INODE, FPERE, NCB, NUPDATE,
     *                          IW( J3 + 1 ), WCB( PLEFT + NPIV ),
     *                          CMUMPS_275(STEP(FPERE),
     *                          PROCNODE_STEPS,SLAVEF),
     *                          ContVec,
     *                          COMM, IERR )
          IF ( IERR .EQ. -1 ) THEN
            CALL CMUMPS_303( .FALSE., FLAG,
     *             BUFR, LBUFR, LBUFR_BYTES,
     *             MYID, SLAVEF, COMM,
     *             N, IPOOL, LPOOL, III, LEAF,
     *             NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *             IWCB, LIWCB,
     *             WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *             PTRICB, INFO, KEEP, STEP,
     *             PROCNODE_STEPS, W1 )
            IF ( INFO( 1 )  .LT. 0 )  GOTO 270
            GOTO 210
          ELSE IF ( IERR .EQ. -2 ) THEN
            INFO( 1 ) = -17
            INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     *                ( NUPDATE + 3 ) * KEEP( 34 )
            GOTO 260
          END IF
        ENDIF
      END IF
      IF ( NSLAVES .NE. 0 .AND. MTYPE .eq. 1
     *     .and. NPIV .NE. 0 ) THEN
        DO ISLAVE = 1, NSLAVES
          PDEST = IW( PTRIST(STEP(INODE)) + 5 + ISLAVE )
          CALL CMUMPS_49( 
     &                KEEP, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB - NELIM, 
     &                NSLAVES, 
     &                Effective_CB_Size, FirstIndex )
 222      CALL CMUMPS_72( INODE, FPERE,
     *    Effective_CB_Size, NPIV,
     *    WCB( PLEFT + NPIV + NELIM + FirstIndex - 1 ),
     *    WCB( PLEFT ),
     *    PDEST, COMM, IERR )
          IF ( IERR .EQ. -1 ) THEN
            CALL CMUMPS_303( .FALSE., FLAG,
     *               BUFR, LBUFR, LBUFR_BYTES,
     *               MYID, SLAVEF, COMM,
     *               N, IPOOL, LPOOL, III, LEAF,
     *               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *               IWCB, LIWCB,
     *               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *               PTRICB, INFO, KEEP, STEP,
     *               PROCNODE_STEPS, W1 )
            IF ( INFO( 1 )  .LT. 0 )  GOTO 270
            GOTO 222
          ELSE IF ( IERR .EQ. -2 ) THEN
            INFO( 1 ) = -17
            INFO( 2 ) = Effective_CB_Size* KEEP( 35 ) +
     *      ( Effective_CB_Size + 4 ) * KEEP( 34 )
            GOTO 260
          END IF
        END DO
      END IF
      PLEFTWCB = PLEFTWCB - LIELL
  270 CONTINUE
      RETURN
  260 CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_302
      RECURSIVE SUBROUTINE CMUMPS_303( BLOQ, FLAG,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           MYID, SLAVEF, COMM,
     *           N, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *           IWCB, LIWCB,
     *           WCB, LWCB, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB, INFO, KEEP, STEP, PROCNODE_STEPS,
     *           W1 )
      LOGICAL BLOQ
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER LIW, LA
      INTEGER INFO( 40 ), KEEP( 150 )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL)
      INTEGER NSTK_S( KEEP(28) )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      COMPLEX WCB( LWCB ), A( LA ), W1( N )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      LOGICAL FLAG
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      INTEGER MSGSOU, MSGTAG, MSGLEN
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     *                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
      END IF
      IF ( FLAG ) THEN
         MSGSOU = STATUS( MPI_SOURCE )
         MSGTAG = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL CMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED,
     *                  MSGSOU, MSGTAG, COMM, STATUS, IERR )
           CALL CMUMPS_323( BUFR, LBUFR, LBUFR_BYTES,
     *          MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *          N, IPOOL, LPOOL, III, LEAF,
     *          NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *          IWCB, LIWCB,
     *          WCB, LWCB, POSWCB,
     *          PLEFTWCB, POSIWCB,
     *          PTRICB, INFO, KEEP, STEP,
     *          PROCNODE_STEPS, W1 )
         END IF
      END IF
      RETURN
      END SUBROUTINE CMUMPS_303
      SUBROUTINE CMUMPS_250(N, NZ, NSTEPS, A, LA,
     * IW, LIW, PERM, NA, LNA, NE_STEPS,
     * NSTK_STEPS, NBPROCFILS,
     * MAXFRT, ND, FILS, STEP, FRERE_STEPS, CAND,
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     * PTRARW, PTRAIW,
     * PTRAST, PTRIST, PTLUST_S, PIMASTER, PAMASTER, ITLOC, 
     * NIRBDU, IPOOL, LPOOL,
     * UULOC, ICNTL, INFO, RINFO, KEEP,PROCNODE_STEPS,
     * SLAVEF,COMM_NODES,MYID,MYID_NODES,
     * BUFR,LBUFR,LBUFR_BYTES,
     * INTARR,DBLARR,root, NELT, FRTPTR, 
     * FRTELT, LPTRAR, COMM_LOAD, ASS_IRECV, SEUIL,
     * SEUIL_LDLT_NIV2, MEM_DISTRIB)
      IMPLICIT NONE
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER N, LIW, LA, LPTRAR, LNA
      INTEGER NZ, NSTEPS, LPOOL, KEEP(150)
      COMPLEX A(LA)
      INTEGER MAXFRT, INFO(40), NIRBDU
      INTEGER SLAVEF,COMM_NODES,MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IW(LIW)
      INTEGER PERM(N), NA(LNA),
     *        NE_STEPS(KEEP(28)), NSTK_STEPS(KEEP(28)),
     *        STEP(N),
     *        ND(KEEP(28)), FILS(N), FRERE_STEPS(KEEP(28)),
     *        PTRARW(LPTRAR), PTRAIW(LPTRAR),PTLUST_S(KEEP(28)),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        PTRAST(KEEP(28)), PTRIST(KEEP(28)),
     *        ITLOC(N), NBPROCFILS( KEEP(28) )
      INTEGER   CAND(SLAVEF+1,MAX(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER IPOOL(LPOOL), PROCNODE_STEPS(KEEP(28))
      INTEGER NELT, COMM_LOAD, ASS_IRECV
      INTEGER FRTPTR(N+1), FRTELT(NELT)
      REAL    RINFO(20)
      INTEGER ICNTL(40)
      REAL UULOC
      REAL SEUIL, SEUIL_LDLT_NIV2
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER INTARR(MAX(1,KEEP(14)))
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      INTEGER MPRINT,LP
      INTEGER I, NBLEAF
      INTEGER INODE, IERR, NTOTPVTOT
      INTEGER NTOTPV
      INTEGER LRLU, IPTRLU, LRLUS,
     *    LEAF, NBROOT, NROOT
      INTRINSIC INT,REAL,LOG
      MPRINT = ICNTL(2)
      LP     = ICNTL(1)
      IF (KEEP(4).LE.0) KEEP(4)=32
      IF (KEEP(5).LE.0) KEEP(5)=16
      IF (KEEP(5).GT.KEEP(4)) KEEP(5) = KEEP(4)
      IF (KEEP(6).LE.0) KEEP(6)=24
      IF (KEEP(3).LE.KEEP(4)) KEEP(3)=KEEP(4)*2
      IF (KEEP(6).GT.KEEP(3)) KEEP(6) = KEEP(3)
      INFO(9)   = 1
      INFO(10)  = 1
      MAXFRT = 0
        LRLU     = LA
        LRLUS    = LRLU
        KEEP(67) = 0 
        IPTRLU   = LRLU
      KEEP(18) = LRLU
      NTOTPV = 0
      NSTK_STEPS(1:KEEP(28)) = NE_STEPS(1:KEEP(28))
      CALL CMUMPS_362(N, LEAF, NBROOT, NROOT,
     *                     MYID_NODES,
     *                     SLAVEF, NA, LNA,
     *                     KEEP, STEP,
     *                     PROCNODE_STEPS,
     *                     IPOOL, LPOOL)
      CALL CMUMPS_506(IPOOL, LPOOL, LEAF)
      IF ( KEEP( 38 ) .NE. 0 ) THEN
        NBROOT = NBROOT + root%NPROW * root%NPCOL - 1
      END IF
      IF ( root%yes )  THEN 
         IF ( CMUMPS_275( STEP(KEEP(38)),PROCNODE_STEPS, SLAVEF)
     *         .NE. MYID_NODES ) THEN
             NROOT = NROOT + 1
         END IF
      END IF
      CALL CMUMPS_251(N,IW,LIW,A,LA,NSTK_STEPS,NBPROCFILS,
     *         INFO(1),ND,FILS,STEP,FRERE_STEPS, CAND,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *         MAXFRT, NTOTPV,PTRIST,PTRAST,
     *         PIMASTER, PAMASTER, PTRARW,PTRAIW,
     *         ITLOC, INFO(2), IPOOL, LPOOL, 
     *         NIRBDU, RINFO,
     *         INFO(9),INFO(10),LRLU,IPTRLU, 
     *         LRLUS, LEAF, NROOT, NBROOT,
     *         UULOC,ICNTL,PTLUST_S,NSTEPS,INFO,KEEP,
     *         PROCNODE_STEPS,SLAVEF,MYID,COMM_NODES,
     *         MYID_NODES,
     *         BUFR,LBUFR, LBUFR_BYTES,
     *         INTARR, DBLARR, root, PERM,
     *         NELT, FRTPTR, FRTELT, LPTRAR, 
     *         COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     *         MEM_DISTRIB(0),NE_STEPS)
      INFO(9)  = INFO(9)  - 1
      INFO(10) = INFO(10) - 1
      KEEP(31) = INFO(9)
      KEEP(32) = INFO(10)
      CALL MPI_ALLREDUCE(NTOTPV, NTOTPVTOT, 1, MPI_INTEGER, MPI_SUM, 
     *                COMM_NODES, IERR)
      IF ( ( (INFO(1).EQ.-10) .and. (NTOTPVTOT.EQ.N)) 
     *              .OR. ( NTOTPVTOT.GT.N ) ) THEN
       write(*,*) ' Error in mc51d NTOTPVTOT=', NTOTPVTOT 
       stop
      endif
      IF ( (INFO(1) .GE. 0 ) .AND. (KEEP(19).EQ.0 )
     *      .AND. ( KEEP(60) .eq. 0 )
     *      .AND. (NTOTPVTOT.NE.N) ) THEN
         INFO(1) = -10
         INFO(2) = NTOTPVTOT
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_250
!      
      SUBROUTINE CMUMPS_26(id)
      USE CMUMPS_LOAD
      USE CMUMPS_CV
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INCLUDE 'cmumps_struc.h'
      TYPE(CMUMPS_STRUC), TARGET :: id
      INTEGER LIW, IKEEP, FILS, FRERE, PTRAR, NFSIZ
      INTEGER NE, NA
      INTEGER I, ALLOCOK
      INTEGER MAXIS1_CHECK
      INTEGER INDX, NRECV, NB_NIV2, LAST_ALLOWED, IDEST
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LOCAL_M, LOCAL_N, K44MAX
      INTEGER NUMROC
      EXTERNAL NUMROC
      INTEGER IRANK
      INTEGER MP, LP, MPG
      LOGICAL PROK, PROKG
      INTEGER SBUF_SEND, SBUF_REC, MEMORY_SIZE
      INTEGER NSLAVES_SAVE
      INTEGER LEAF, INODE, ISTEP, INN, LPTRAR
      INTEGER NBLEAF, NBROOT, MYROW_CHECK, INIV2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: REQPTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SSARBR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PAR2_NODES
      INTEGER, DIMENSION(:), POINTER :: IWtemp 
      INTEGER, DIMENSION(:), POINTER :: XNODEL, NODEL
      INTEGER, POINTER ::  NELT, LELTVAR
      INTEGER, DIMENSION(:), POINTER :: KEEP, INFO, INFOG
      REAL, DIMENSION(:), POINTER :: RINFO
      REAL, DIMENSION(:), POINTER :: RINFOG
      INTEGER, DIMENSION(:), POINTER :: ICNTL
      LOGICAL I_AM_SLAVE
      INTEGER CMUMPS_330, CMUMPS_275
      EXTERNAL CMUMPS_330, CMUMPS_275
      INTEGER K,J
      INTEGER SIZE_TEMP_MEM
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_MEM 
      KEEP   => id%KEEP
      INFO   => id%INFO
      RINFO  => id%RINFO
      INFOG  => id%INFOG
      RINFOG => id%RINFOG
      ICNTL  => id%ICNTL
      NELT    => id%NELT
      LELTVAR => id%LELTVAR
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               id%KEEP(46) .eq. 1 ) )
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IF (PROK) WRITE( MP, 220 )
      IF (PROKG.AND.(MPG .NE. MP)) WRITE( MPG, 220 )
 220  FORMAT( /' CMUMPS Version 4.3.2 -- November 2003' )
      IF ( PROK ) THEN
        IF ( KEEP(50) .eq. 0 ) THEN
          WRITE(MP, '(A)') 'L U Solver for unsymmetric matrices'
        ELSE IF ( KEEP(50) .eq. 1 ) THEN
          WRITE(MP, '(A)') 
     * 'L D L^T Solver for symmetric positive definite matrices'
        ELSE
          WRITE(MP, '(A)') 
     * 'L D L^T Solver for general symmetric matrices'
        END IF
        IF ( KEEP(46) .eq. 1 ) THEN
          WRITE(MP, '(A)') 'Type of parallelism: Working host'
        ELSE
          WRITE(MP, '(A)') 'Type of parallelism: Host not working'
        END IF
      END IF
      IF ( PROKG .AND. (MP.NE.MPG)) THEN
        IF ( KEEP(50) .eq. 0 ) THEN
          WRITE(MPG, '(A)') 'L U Solver for unsymmetric matrices'
        ELSE IF ( KEEP(50) .eq. 1 ) THEN
          WRITE(MPG, '(A)') 
     * 'L D L^T Solver for symmetric positive definite matrices'
        ELSE
          WRITE(MPG, '(A)') 
     * 'L D L^T Solver for general symmetric matrices'
        END IF
        IF ( KEEP(46) .eq. 1 ) THEN
          WRITE(MPG, '(A)') 'Type of parallelism: Working host'
        ELSE
          WRITE(MPG, '(A)') 'Type of parallelism: Host not working'
        END IF
      END IF
      IF (PROK) WRITE( MP, 110 )
      IF (PROKG .AND. (MPG.NE.MP)) WRITE( MPG, 110 )
      IF ( (id%KEEP(24).NE.0) .AND.
     *     id%NSLAVES.eq.1 ) THEN
        id%KEEP(24) = 0
        IF ( PROKG ) THEN
          WRITE(MPG, '(A)')
     *    ' Resetting candidate strategy to 0 because NSLAVES=1'
          WRITE(MPG, '(A)') ' '
        END IF
      END IF
      IF ( (id%KEEP(24).NE.0) .AND.(id%KEEP(24).NE.1).AND.
     *     (id%KEEP(24).NE.8).AND.(id%KEEP(24).NE.10)
     *  .AND.(id%KEEP(24).NE.12).AND.(id%KEEP(24).NE.14)
     *  .AND.(id%KEEP(24).NE.16).AND.(id%KEEP(24).NE.18)) THEN
       id%KEEP(24) = 8
        IF ( PROKG ) THEN
          WRITE(MPG, '(A)')
     *    ' Resetting candidate strategy to 8 '
          WRITE(MPG, '(A)') ' '
        END IF
      END IF
      IF ( id%MYID.EQ.MASTER ) THEN
        IF ( id%ICNTL(7) .EQ.1 ) THEN
          IF ( .not. ASSOCIATED( id%PERM_IN ) ) THEN
            id%INFO(1) = -22
            id%INFO(2) = 3
            GO TO 333
          ELSE IF ( SIZE( id%PERM_IN ) < id%N ) THEN
            id%INFO(1) = -22
            id%INFO(2) = 3
            GO TO 333
          END IF
        ENDIF
        IF (KEEP(9) .LE. 1 ) KEEP(9) = 500
        IF ( KEEP(10) .GT. 0 ) THEN 
         IF ((KEEP(10).LE.1) .OR. (KEEP(10).GT.KEEP(9))) 
     *      KEEP(10) = MIN(KEEP(9),100)
        ENDIF
        IF (KEEP(48). EQ. 1 ) KEEP(48) = -12345
        LAST_ALLOWED = 3
        IF ( (KEEP(48).LT.0) .OR. (KEEP(48).GT.LAST_ALLOWED) ) THEN
         IF (KEEP(50).EQ.0) THEN 
           KEEP(48) = 0  
         ELSE
           KEEP(48) = 3  
         ENDIF
        ENDIF
        KEEP(60) = ICNTL(19)
        IF ( KEEP( 60 ) .LE. 0 ) KEEP( 60 ) = 0
        IF ( KEEP(60) .NE.0 ) THEN
         IF ( .not. ASSOCIATED( id%LISTVAR_SCHUR ) ) THEN
            id%INFO(1) = -22
            id%INFO(2) = 8
            GO TO 333
         ELSE IF (SIZE(id%LISTVAR_SCHUR)<id%SIZE_SCHUR) THEN
            id%INFO(1) = -22
            id%INFO(2) = 8
            GO TO 333
         END IF
        ENDIF
        IF ( KEEP(60) .ne. 0 .and. ICNTL(7) .eq. 1 ) THEN
          DO I = 1,id%SIZE_SCHUR
            IF (id%PERM_IN(id%LISTVAR_SCHUR(I))
     *        ==id%N-id%SIZE_SCHUR+I)
     *        CYCLE
            ICNTL(7) = 0
            IF (MPG.GT.0) THEN
              WRITE(MPG,'(A)')
     *    ' ** Ignoring user-ordering, because incompatible with Schur.'
              WRITE(MPG,'(A)') ' ** ICNTL(7) reset to 0.'
            END IF
            EXIT
          ENDDO
        END IF
        KEEP(23) = ICNTL(6)
        KEEP(54) = ICNTL(18)
        IF ( KEEP(54) .LT. 0 .or. KEEP(54).GT.3 ) THEN
          IF ( PROKG ) THEN
             WRITE(MPG, *) ' Out-of-range value for ICNTL(18).'
             WRITE(MPG, *) ' Used 0 ie matrix not distributed'
          END IF
          KEEP(54) = 0
        END IF
        KEEP(55) = ICNTL(5)
        IF ( KEEP(55) .LT. 0 .OR. KEEP(55) .GT. 2 ) THEN
          IF ( PROKG ) THEN
             WRITE(MPG, *) ' Out-of-range value for ICNTL(5).'
             WRITE(MPG, *) ' Used 0 ie matrix is assembled'
          END IF
          KEEP(55) = 0
        END IF
        IF (KEEP(23).LT.0.OR.KEEP(23).GT.7) KEEP(23) = 0
        IF ( (KEEP(23) .NE. 0) .AND. (KEEP(50) .NE. 0) ) THEN
          KEEP(23) = 0
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Max-trans not compatible with LDLT factorization '
          END IF
        END IF
        IF  ( (KEEP(23) .NE. 0) .AND. (KEEP(60) .GT. 0) ) THEN
          KEEP(23) = 0
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed because of Schur'
          END IF
        END IF
#if !defined(GRIDTLSE)
        IF ( (KEEP(23) .NE. 0) .AND. (ICNTL(7) .EQ. 1) ) THEN
          KEEP(23) = 0
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed because ordering is given'
          END IF
        END IF
#endif
        IF ( KEEP(23) .NE. 0 .AND. KEEP(54) .NE. 0 ) THEN
          KEEP(23) = 0
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed because matrix is distributed'
          END IF
        END IF
        IF ( KEEP(23) .NE. 0 .AND. KEEP(55) .ne. 0 ) THEN
          KEEP(23) = 0
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed for element matrix'
          END IF
        ENDIF
        IF ( KEEP(54) .ne. 0 .AND. KEEP(55) .ne. 0 ) THEN
          KEEP(54) = 0
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Distributed entry not available for element matrix'
          END IF
        ENDIF
      END IF
 333  CONTINUE
      CALL CMUMPS_276( ICNTL, INFO,
     *                  id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( KEEP(54), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(69), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      IF ( ASSOCIATED(id%MEM_DIST) ) DEALLOCATE( id%MEM_DIST )
      ALLOCATE( id%MEM_DIST( 0:id%NSLAVES-1 ), STAT=IERR )
      IF ( IERR .GT. 0 ) THEN
         WRITE(*,*) 'PB ALLOCATION MEM_DIST'
         INFO(1) = -13
         INFO(2) = id%NSLAVES
         RETURN
      END IF
      id%MEM_DIST(0:id%NSLAVES-1) = 0
      CALL CMUMPS_427(
     *     id%COMM,id%COMM_NODES,KEEP(69),KEEP(46),
     *     id%NSLAVES,id%MEM_DIST,INFO)
      IF ( KEEP(54) .eq. 3 ) THEN
          IF ( KEEP(46) .eq. 0 .and. id%MYID .eq. MASTER ) THEN
            id%NZ_loc = 0
          END IF
          CALL MPI_REDUCE( id%NZ_loc, id%NZ, 1, MPI_INTEGER, MPI_SUM,
     *                     MASTER, id%COMM, IERR )
          IF ( id%MYID .eq. MASTER ) THEN
            ALLOCATE( REQPTR( id%NPROCS, 3 ), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = 3 * id%NPROCS
              GOTO 13
            END IF
            ALLOCATE( id%IRN( id%NZ ), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%NZ
              GOTO 13
            END IF
            ALLOCATE( id%JCN( id%NZ ), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%NZ
              GOTO 13
            END IF
          END IF
 13       CONTINUE
          CALL CMUMPS_276( ICNTL, INFO,
     *                              id%COMM, id%MYID )
          IF ( INFO(1) < 0 ) RETURN
          IF ( id%MYID .eq. MASTER ) THEN
            DO I = 1, id%NPROCS - 1
              CALL MPI_RECV( REQPTR( I+1, 1 ), 1, 
     *        MPI_INTEGER, I,
     *        COLLECT_NZ, id%COMM, STATUS, IERR )
            END DO
            IF ( KEEP(46) .eq. 0 ) THEN
              REQPTR( 1, 1 ) = 1
            ELSE
              REQPTR( 1, 1 ) = id%NZ_loc + 1
            END IF
            DO I = 2, id%NPROCS
              REQPTR( I, 1 ) = REQPTR( I, 1 ) + REQPTR( I-1, 1 )
            END DO
          ELSE
            CALL MPI_SEND( id%NZ_loc, 1, MPI_INTEGER, MASTER,
     *        COLLECT_NZ, id%COMM, IERR )
          END IF
          IF ( id%MYID .eq. MASTER ) THEN
            NRECV = 0
            DO I = 1, id%NPROCS - 1
              IF ( REQPTR( I + 1, 1 ) - REQPTR( I, 1 ) .NE. 0 ) THEN
                NRECV = NRECV + 2
                CALL MPI_IRECV( id%IRN( REQPTR( I, 1 ) ),
     *          REQPTR( I + 1, 1 ) - REQPTR( I, 1 ), 
     *          MPI_INTEGER,
     *          I, COLLECT_IRN, id%COMM, REQPTR(I, 2), IERR )
                CALL MPI_IRECV( id%JCN( REQPTR( I, 1 ) ),
     *          REQPTR( I + 1, 1 ) - REQPTR( I, 1 ),   
     *          MPI_INTEGER,
     *          I, COLLECT_JCN, id%COMM, REQPTR(I, 3), IERR )
              ELSE
                REQPTR(I, 2) = MPI_REQUEST_NULL
                REQPTR(I, 3) = MPI_REQUEST_NULL
              END IF
            END DO
          ELSE
            IF ( id%NZ_loc .NE. 0 ) THEN
              CALL MPI_SEND( id%IRN_loc(1), id%NZ_loc, 
     *          MPI_INTEGER, MASTER,
     *          COLLECT_IRN, id%COMM, IERR )
              CALL MPI_SEND( id%JCN_loc(1), id%NZ_loc, 
     *          MPI_INTEGER, MASTER,
     *          COLLECT_JCN, id%COMM, IERR )
            END IF
          END IF
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( id%NZ_loc .NE. 0 ) THEN
              id%IRN(1:id%NZ_loc) = id%IRN_loc(1:id%NZ_loc)
              id%JCN(1:id%NZ_loc) = id%JCN_loc(1:id%NZ_loc)
            END IF
            REQPTR( id%NPROCS, 2 ) = MPI_REQUEST_NULL
            REQPTR( id%NPROCS, 3 ) = MPI_REQUEST_NULL
            DO I = 1, NRECV
              CALL MPI_WAITANY
     *        ( 2 * id%NPROCS, REQPTR( 1, 2 ), INDX, STATUS, IERR )
            END DO
            DEALLOCATE( REQPTR )
          END IF
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
 1234   CONTINUE
        IF ((KEEP(23) .NE. 0).AND.(KEEP(23).NE.7) ) THEN
          CALL CMUMPS_203(id%N, id%NZ, KEEP(23), id%IS1(1), id,
     *        ICNTL, INFO)
          IF (INFO(1) .LT. 0) THEN
            GOTO 12
          END IF
        END IF
        IF (KEEP(55) .EQ. 0) THEN
#if !defined(OLDSymbolic)
          IF ( ICNTL(7) .eq. 1 ) THEN
            LIW = 2 * id%NZ + 3 * id%N + 2
#else
          IF ( ICNTL(7) .eq. 1 ) THEN
            LIW = id%NZ + 3 * id%N + 1
#endif
          ELSE
            LIW = 2 * id%NZ + 3 * id%N + 2
          ENDIF
        ELSE
         IF (KEEP(60).NE.0) THEN
           LIW = id%N + id%N 
         ELSE
           LIW =  id%N + id%N + id%N+3 + id%N+1
         ENDIF
        ENDIF
        IF ((KEEP(23) .NE. 0).AND.(KEEP(23).NE.7)) THEN
          IKEEP = id%N + 1
        ELSE
          IKEEP = 1
        END IF
        NA      = IKEEP +     id%N
        NE      = IKEEP + 2 * id%N
        FILS    = IKEEP + 3 * id%N
        FRERE   = FILS  +     id%N
        PTRAR   = FRERE +     id%N
        IF (KEEP(55) .EQ. 0) THEN
          NFSIZ  = PTRAR + 4 * id%N
          MAXIS1_CHECK = NFSIZ + id%N - 1
        ELSE
          NFSIZ  = PTRAR + 2 * (NELT + 1)
          MAXIS1_CHECK = NFSIZ + id%N -1
        ENDIF
        IF ( id%MAXIS1 .LT. MAXIS1_CHECK ) THEN
          IF (LP.GE.0) THEN
          WRITE(LP,*) '***********************************'
          WRITE(LP,*) 'MAXIS1 and MAXIS1_CHECK are different !!'
          WRITE(LP,*) 'MAXIS1, MAXIS1_CHECK=',id%MAXIS1,
     *                 MAXIS1_CHECK
          WRITE(LP,*) 'This might cause problems ...'
          WRITE(LP,*) '***********************************'
          ENDIF
        END IF
        IF ( ICNTL( 7 ) .EQ. 1 ) THEN
          DO I = 1, id%N
            id%IS1( IKEEP + I - 1 ) = id%PERM_IN( I )
          END DO
        END IF
        KEEP(12) = ICNTL(14)
        IF ( (id%KEEP(54).NE.0).AND.
     *       (id%KEEP(50).NE.1) ) id%KEEP(12)= id%KEEP(12)*2
#if defined(try_null_space)
        KEEP(53) = ICNTL(16)
#else
        KEEP(53)=0
#endif
        IF ( KEEP(53) .LT. 0 .OR.
     *       KEEP(53) .GT. 10 ) THEN
           KEEP(53) = 0
        END IF
        IF ( KEEP(53) > 2 .and. KEEP(53) < 9 ) THEN
          IF ( MPG > 0 ) THEN
            WRITE( MPG, '(A)' )
     *    ' Null space option currently disabled.'
            WRITE( MPG, '(A)' )
     *    ' Contact authors to experiment it.'
          END IF
          KEEP(53) = 0
        END IF
        INFOG(1) = 0
        INFOG(2) = 0
        INFOG(8) = -1
        IF ( .not. ASSOCIATED( id%LISTVAR_SCHUR ) ) THEN
          id%SIZE_SCHUR = 1
          ALLOCATE( id%LISTVAR_SCHUR( 1 ), stat=allocok )
          IF ( allocok .GT. 0 ) THEN
            WRITE(*,*) 'PB allocating an array of size 1 in Schur '
            CALL CMUMPS_ABORT()
          END IF
        END IF
        IF (KEEP(55) .EQ. 0) THEN
          CALL CMUMPS_195(id%N, id%NZ, id%IRN(1), id%JCN(1),
     *        LIW, id%IS1(IKEEP),
     *        id%IS1(PTRAR), ICNTL(7), id%IS1(NFSIZ),
     *        id%IS1(FILS), id%IS1(FRERE),
     *        id%LISTVAR_SCHUR(1), id%SIZE_SCHUR,
     *        ICNTL(1), INFOG(1), KEEP(1),id%NSLAVES, 
     *        id%SYM)
          IF ( (KEEP(23).LE.-1).AND.(KEEP(23).GE.-6) ) THEN
            KEEP(23) = -KEEP(23)
            IF (.not.ASSOCIATED(id%A)) KEEP(23) = 1
            GOTO 1234
          ENDIF
        ELSE
          ALLOCATE( IWtemp ( 3*id%N ), stat = ierr )
          IF ( ierr .GT. 0 ) THEN
            INFO( 1 ) = -13
            INFO( 2 ) = 3*id%N
            GOTO 12
          ENDIF
          ALLOCATE( XNODEL ( id%N+1 ), stat = ierr )
          IF ( ierr .GT. 0 ) THEN
            INFO( 1 ) = -13
            INFO( 2 ) = id%N + 1
            GOTO 12
          ENDIF
          IF (LELTVAR.ne.id%ELTPTR(NELT+1)-1)  THEN
           INFO(1) = -2002
           INFO(2) = id%ELTPTR(NELT+1)-1
           GOTO 12
          ENDIF
          ALLOCATE( NODEL ( LELTVAR ), stat = ierr )
          IF ( ierr .GT. 0 ) THEN
            INFO( 1 ) = -13
            INFO( 2 ) = LELTVAR
            GOTO 12
          ENDIF
          CALL CMUMPS_128(id%N, NELT,
     *        id%ELTPTR(1), id%ELTVAR(1), LIW,
     *        id%IS1(IKEEP),
     *        IWtemp, ICNTL(7), id%IS1(NFSIZ), id%IS1(FILS),
     *        id%IS1(FRERE), id%LISTVAR_SCHUR, id%SIZE_SCHUR, 
     *        ICNTL(1), INFOG(1), KEEP(1),
     *        id%ELTPROC(1), id%NSLAVES, 
     *        XNODEL(1), NODEL(1))
          IF ( KEEP(24).EQ.0 ) DEALLOCATE( IWtemp )
        ENDIF
         IF ( ICNTL(8).EQ.7 ) THEN
           IF (id%SYM.EQ.0) THEN
            ICNTL(8) = 4
           ELSE
            ICNTL(8) = 0
           ENDIF
         ENDIF
        IF ( KEEP(60) .eq. 0 ) THEN
          DEALLOCATE( id%LISTVAR_SCHUR )
          NULLIFY   ( id%LISTVAR_SCHUR )
          id%SIZE_SCHUR = 0
        ELSE
         I = id%SIZE_SCHUR
         IF (id%SYM.EQ.0) THEN
          INFOG(20) = INFOG(20) - I*I
         ELSE
          INFOG(20) = INFOG(20) - (I*(I+1))/2
         ENDIF
        END IF
        INFO(1)=INFOG(1)
        INFO(2)=INFOG(2)
        KEEP(28) = INFOG(6)
        IF ( INFO(1) .LT. 0 ) THEN
          GO TO 12
        ENDIF
        ALLOCATE(SSARBR(id%N), stat=allocok)
        IF ( allocok .GT. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = id%N
          GOTO 12
        END IF
        IF (KEEP(24).EQ.0) THEN
          IF ( .NOT. associated ( id%POIDS ) ) THEN
           ALLOCATE( id%POIDS( id%NSLAVES ),stat=ierr )
           IF ( ierr .GT. 0 ) THEN
            INFO(1) = -13
            INFO(2) = id%NSLAVES
            GOTO 12
           ENDIF
           id%POIDS = 1
          END IF
        NSLAVES_SAVE = id%NSLAVES
        CALL CMUMPS_DISTRI(id%N, id%NSLAVES, id%IS1( NE ),
     *    id%IS1( NA ),
     *    id%IS1( NFSIZ ),
     *    id%IS1( FRERE ), id%IS1( FILS ),
     *    id%PROCNODE(1), id%POIDS(1),
     *    SSARBR(1),
     *    id%NBSA , KEEP(50) )
        id%NSLAVES = NSLAVES_SAVE
        CALL CMUMPS_21( id%N, id%MYID, id%PROCNODE,
     *                id%IS1( NFSIZ ),
     *                SSARBR(1), id%NBSA,
     *                id%IS1( FRERE ), id%IS1( FILS ),
     *                id%NSLAVES,
     *                KEEP(1), ICNTL(1), INFO(1) )
      ELSE
        call CMUMPS_369(id%N,id%NSLAVES,ICNTL,INFOG,
     *              id%IS1(NE:NE+id%N-1),
     *              id%IS1(NA:NA+id%N-1),id%IS1(NFSIZ:NFSIZ+id%N-1),
     *              id%IS1(FRERE:FRERE+id%N-1),
     *              id%IS1(FILS:FILS+id%N-1),
     *              KEEP,id%PROCNODE,
     *              SSARBR(1),id%NBSA,ierr)
          if(ierr.eq.-999) then 
              write(6,*) ' Internal error in CMUMPS_369'
             INFO(1) = ierr
             GOTO 12
          endif
          if(ierr.ne.0) then 
             INFO(1) = -13
             INFO(2) = ierr
             GOTO 12
          endif
        IF (KEEP(55) .EQ. 0) THEN
         CALL CMUMPS_202(id%N, id%NZ, id%IS1(IKEEP),
     *      id%IS1(FILS), id%IS1(FRERE), id%IS1(NFSIZ),
     *      id%IS1(IKEEP+2*id%N), 
     *      id%IS1(IKEEP+id%N),
     *      id%IRN(1), id%JCN(1), id%IS1(PTRAR+2*id%N) , 
     *      id%IS1(PTRAR+3*id%N), id%IS1(PTRAR),
     *      ICNTL(1), INFOG(1), KEEP(1))
        ELSE
         CALL CMUMPS_353(id%N, id%IS1(IKEEP),
     *      id%IS1(FILS), 
     *      id%IS1(FRERE), id%IS1(NFSIZ),
     *      id%IS1(IKEEP+2*id%N), 
     *      id%IS1(IKEEP+id%N), IWtemp,
     *      ICNTL(1), INFOG(1), KEEP(1))
          DEALLOCATE( IWtemp )
        ENDIF
      ENDIF
       DEALLOCATE( SSARBR)
        IF ( KEEP(54) .eq. 3 ) THEN
          DEALLOCATE( id%IRN )
          DEALLOCATE( id%JCN )
        END IF
       IF (KEEP(55).NE.0) THEN
          IF (ASSOCIATED(id%FRTPTR)) DEALLOCATE(id%FRTPTR)
          IF (ASSOCIATED(id%FRTELT)) DEALLOCATE(id%FRTELT)
          ALLOCATE( id%FRTPTR(id%N+1), id%FRTELT(id%NELT),
     *              stat=allocok )
          IF ( allocok .ne. 0) THEN
            INFO(1) = -13
            INFO(2) = id%NELT+id%N+1
            GOTO 12
          ENDIF
         CALL CMUMPS_153(
     *   id%N, NELT, id%ELTPTR(NELT+1)-1, id%IS1(FRERE),
     *   id%IS1(FILS),
     *   id%IS1(IKEEP+id%N), id%IS1(IKEEP+2*id%N), XNODEL, NODEL,
     *   id%FRTPTR(1), id%FRTELT(1), id%ELTPROC(1))
         DEALLOCATE(XNODEL)
         DEALLOCATE(NODEL)
       ENDIF
        IF ( INFO( 1 ) .LT. 0 ) GOTO 12
        IF ( KEEP(55) .ne. 0 ) THEN
           CALL CMUMPS_120(id%N, NELT, id%ELTPROC(1),id%NSLAVES,
     *                        id%PROCNODE)
        END IF
        NB_NIV2 = KEEP(56)
        IF ( NB_NIV2.GT.0 ) THEN
          ALLOCATE(PAR2_NODES(NB_NIV2),
     *             stat=allocok)
          if (allocok .gt.0) then
             info(1)= -13
             info(2)= NB_NIV2
             goto 12
          end if
        ENDIF
        IF ((NB_NIV2.GT.0) .AND. (KEEP(24).EQ.0)) THEN
          INIV2 = 0
          DO 777 INODE = 1, id%N
            IF ( ( id%IS1(FRERE+INODE-1) .NE. id%N+1 ) .AND.
     *           ( CMUMPS_330(INODE,id%PROCNODE(1),id%NSLAVES)
     *             .eq. 2) ) THEN
              INIV2 = INIV2 + 1
              PAR2_NODES(INIV2) = INODE
            END IF
 777      CONTINUE
          IF ( INIV2 .NE. NB_NIV2 ) THEN
            WRITE(*,*) "Internal Error 2 in CMUMPS_26",
     *                  INIV2, NB_NIV2
            CALL CMUMPS_ABORT()
          ENDIF
        ENDIF
        IF ( (KEEP(24) .NE. 0) .AND. (NB_NIV2.GT.0) ) THEN
          IF ( ASSOCIATED(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
          ALLOCATE( id%CANDIDATES(id%nslaves+1,NB_NIV2),
     *             stat=allocok)
          if (allocok .gt.0) then
             info(1)= -13
             info(2)= NB_NIV2*(id%nslaves+2)
             goto 12
          end if
          call CMUMPS_393
     *         (PAR2_NODES,id%CANDIDATES,ierr)
          if(ierr.ne.0)  then
             INFO(1) = -2002
             GOTO 12
          endif
          call CMUMPS_494()
          if(ierr.ne.0)  then
             INFO(1) = -2002
             GOTO 12
          endif
        ENDIF
        KEEP(129)=INFO(7)
        KEEP(130)=INFO(8)
 12     CONTINUE
      END IF
      CALL CMUMPS_276( ICNTL, INFO, id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
      CALL MPI_BCAST( id%KEEP, 100, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      CALL MPI_BCAST( id%N, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      IF ( id%KEEP(55) .EQ. 0) THEN
        CALL MPI_BCAST( id%NZ, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      ENDIF
      CALL MPI_BCAST( id%NBSA, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      IF ( id%KEEP(55) .NE. 0) THEN
        CALL MPI_BCAST( id%NELT, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      ENDIF
      IF ( ASSOCIATED( id%STEP ) ) DEALLOCATE(id%STEP)
      ALLOCATE  (id%STEP(id%N),stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = id%N
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%PROCNODE_STEPS ))
     *    DEALLOCATE(id%PROCNODE_STEPS)
      ALLOCATE(id%PROCNODE_STEPS( id%KEEP(28)), stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%NE_STEPS ) ) DEALLOCATE(id%NE_STEPS)
      ALLOCATE (id%NE_STEPS(id%KEEP(28)),stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%ND_STEPS ) ) DEALLOCATE(id%ND_STEPS)
      ALLOCATE (id%ND_STEPS(id%KEEP(28)),stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%FRERE_STEPS ) ) DEALLOCATE(id%FRERE_STEPS)
      ALLOCATE (id%FRERE_STEPS(id%KEEP(28)),stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%SYM_PERM ) ) DEALLOCATE(id%SYM_PERM)
      ALLOCATE(id%SYM_PERM(id%N),stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = id%N
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%UNS_PERM ) ) DEALLOCATE(id%UNS_PERM)
      IF ( ASSOCIATED( id%FILS ) ) DEALLOCATE(id%FILS)
      ALLOCATE(id%FILS( id%N ), stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = id%N
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%PTRAR ) ) DEALLOCATE(id%PTRAR)
      IF ( id%KEEP(55).eq.0 ) THEN
        LPTRAR = id%N+id%N
      ELSE
        LPTRAR = id%NELT+id%NELT+2
      ENDIF
      ALLOCATE(id%PTRAR( LPTRAR ), stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = LPTRAR
        GOTO 94
      ENDIF
      IF ( id%KEEP(55) /= 0 ) THEN
        IF ( id%MYID .NE. MASTER ) THEN
          IF (ASSOCIATED(id%FRTPTR)) DEALLOCATE(id%FRTPTR)
          IF (ASSOCIATED(id%FRTELT)) DEALLOCATE(id%FRTELT)
          ALLOCATE( id%FRTPTR(id%N+1), id%FRTELT(id%NELT),
     *              stat=allocok )
          IF ( allocok .ne. 0) THEN
            INFO(1) = -13
            INFO(2) = id%NELT+id%N+1
            GOTO 94
          ENDIF
        END IF
      ELSE
        ALLOCATE( id%FRTPTR(1), id%FRTELT(1),
     *            stat=allocok )
        IF ( allocok .ne. 0) THEN
          INFO(1) = -13
          INFO(2) = 2
          GOTO 94
        ENDIF
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
        id%SYM_PERM(1:id%N) = id%IS1(IKEEP:IKEEP+id%N-1)
      END IF
      IF ( id%MYID == MASTER .AND. id%KEEP(23) .NE. 0 ) THEN
        ALLOCATE(id%UNS_PERM(id%N),stat=allocok)
        IF ( allocok .ne. 0) THEN
          INFO(1) = -13
          INFO(2) = id%N
          GOTO 94
        ENDIF
        id%UNS_PERM(1:id%N) = id%IS1(1:id%N)
      ENDIF
  94  CONTINUE
      CALL CMUMPS_276( ICNTL, INFO,
     *                       id%COMM, id%MYID )
      IF ( id%MYID .EQ. MASTER ) THEN
        id%FILS(1:id%N) = id%IS1(FILS:FILS+id%N-1)
        IF ( id%KEEP(55) == 0 ) THEN
          id%PTRAR(1:LPTRAR) = id%IS1(PTRAR:PTRAR+LPTRAR-1)
        ELSE
          id%PTRAR(id%NELT+2:2*id%NELT+2)=id%ELTPTR(1:id%NELT+1)
        ENDIF
      END IF
      IF (id%MYID .EQ. MASTER ) THEN
        IF (id%N.eq.1) THEN
          NBROOT = 1
          NBLEAF = 1
        ELSE IF (id%IS1(NA+id%N-1) .LT.0) THEN
          NBLEAF = id%N
          NBROOT = id%N
        ELSE IF (id%IS1(NA+id%N-2) .LT.0) THEN
          NBLEAF = id%N-1
          NBROOT = id%IS1(NA+id%N-1)
        ELSE
          NBLEAF = id%IS1(NA+id%N-2)
          NBROOT = id%IS1(NA+id%N-1)
        ENDIF
        id%LNA = 2+NBLEAF+NBROOT
      ENDIF
      CALL MPI_BCAST( id%LNA, 1, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      IF ( ASSOCIATED( id%NA ) ) DEALLOCATE(id%NA)
      ALLOCATE(id%NA(id%LNA),stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -13
        INFO(2) = LPTRAR
        GOTO 96
      ENDIF
      IF (id%MYID .EQ.MASTER ) THEN
        id%NA(1) = NBLEAF
        id%NA(2) = NBROOT
        LEAF = 3
        IF ( id%N == 1 ) THEN
          id%NA(LEAF) = 1
          LEAF = LEAF + 1
        ELSE IF (id%IS1(NA+id%N-1) < 0) THEN
          id%NA(LEAF) = - id%IS1(NA+id%N-1)-1
          LEAF = LEAF + 1
          DO I = 1, NBLEAF - 1
            id%NA(LEAF) = id%IS1(NA+I-1)
            LEAF = LEAF + 1
          ENDDO
        ELSE IF (id%IS1(NA+id%N-2) < 0 ) THEN
          INODE = - id%IS1(NA+id%N-2) - 1
          id%NA(LEAF) = INODE
          LEAF =LEAF + 1
          IF ( NBLEAF > 1 ) THEN
            DO I = 1, NBLEAF - 1
              id%NA(LEAF) = id%IS1(NA+I-1)
              LEAF = LEAF + 1
            ENDDO
          ENDIF
        ELSE
          DO I = 1, NBLEAF
            id%NA(LEAF) = id%IS1(NA+I-1)
            LEAF = LEAF + 1
          ENDDO
        END IF
      END IF
 96   CONTINUE
      CALL CMUMPS_276( ICNTL, INFO,
     *                       id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( id%MYID .EQ. MASTER ) THEN
        ISTEP = 0
        DO I = 1, id%N
          IF ( id%IS1(FRERE+I-1) .ne. id%N + 1 ) THEN
            ISTEP = ISTEP + 1
            id%STEP(I)=ISTEP
            INN = id%IS1(FILS+I-1)
            DO WHILE ( INN .GT. 0 )
              id%STEP(INN) = - ISTEP
              INN = id%IS1(FILS + INN -1)
            END DO
            IF (id%IS1(FRERE+I-1) .eq. 0) THEN
              id%NA(LEAF) = I
              LEAF = LEAF + 1
            ENDIF
          ENDIF
        END DO
        DO I = 1, id%N
          IF (id%IS1(FRERE+I-1) .NE. id%N+1) THEN
            id%PROCNODE_STEPS(id%STEP(I)) = id%PROCNODE( I )
            id%FRERE_STEPS(id%STEP(I))    = id%IS1(FRERE+I-1)
            id%NE_STEPS(id%STEP(I))    = id%IS1(NE+I-1)
            id%ND_STEPS(id%STEP(I))    = id%IS1(NFSIZ+I-1)
          ENDIF
        ENDDO
        DEALLOCATE(id%PROCNODE)
        NULLIFY(id%PROCNODE)
        DEALLOCATE(id%IS1)
        NULLIFY(id%IS1)
       IF (( id%KEEP(70) .NE. 2 ).OR.(id%KEEP(47).EQ.4))
     *      THEN
         IF (id%KEEP(47) .EQ. 4) THEN
           IF(id%NSLAVES.GT.1) THEN
              SIZE_TEMP_MEM = id%NBSA
           ELSE
              SIZE_TEMP_MEM = 1
           ENDIF
         ELSE
           SIZE_TEMP_MEM = 1
         ENDIF
         ALLOCATE(TEMP_MEM(SIZE_TEMP_MEM,id%NSLAVES),stat=allocok) 
         IF (allocok .ne.0) then
           info(1)= -13
           info(2)= SIZE_TEMP_MEM*id%NSLAVES
           GOTO 80
         end if
         CALL CMUMPS_363(id%N, id%FRERE_STEPS(1), id %STEP(1),
     *      id%FILS(1), id%NA(1), id%LNA, id%NE_STEPS(1),
     *      id%ND_STEPS(1),
     *      id%KEEP(28), id%KEEP(70),
     *      id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(47)
     *,id%PROCNODE_STEPS,TEMP_MEM,id%NSLAVES, SIZE_TEMP_MEM
     * )
       END IF
 80     CONTINUE
      ENDIF
      CALL CMUMPS_276( ICNTL, INFO,
     *                      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      CALL MPI_BCAST( id%SYM_PERM(1), id%N, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%FILS(1), id%N, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      IF ( id%KEEP(55) .eq. 0 ) THEN
        CALL MPI_BCAST( id%PTRAR(1), LPTRAR, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
      ELSE
        CALL MPI_BCAST( id%PTRAR(id%NELT+2), id%NELT+1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
      END IF
      CALL MPI_BCAST( id%NA(1), id%LNA, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%STEP, id%N, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%PROCNODE_STEPS, id%KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%FRERE_STEPS, id% KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR)
      CALL MPI_BCAST( id%NE_STEPS, id%KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%ND_STEPS, id%KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      IF (id%KEEP(47) .EQ. 4 ) THEN
       IF(id%MYID .EQ. MASTER)THEN
         DO K=1,id%NSLAVES
           DO J=1,SIZE_TEMP_MEM
             IF(TEMP_MEM(J,K).EQ.-1) GOTO 666 
           ENDDO
 666       CONTINUE
           J=J-1
           IF (id%KEEP(46) == 1) THEN
             IDEST = K - 1
           ELSE
             IDEST = K
           ENDIF
           IF (IDEST .NE. MASTER) THEN
             CALL MPI_SEND(J,1,MPI_INTEGER,IDEST,0,
     *         id%COMM,IERR)
             CALL MPI_SEND(TEMP_MEM(1,K),J,MPI_INTEGER,
     *         IDEST, 0, id%COMM,IERR)
           ELSE
             IF(ASSOCIATED(id%MEM_SUBTREE))
     *         DEALLOCATE(id%MEM_SUBTREE)
             ALLOCATE(id%MEM_SUBTREE(J),stat=allocok)
             IF (allocok .ne.0) then
               info(1)= -13
               info(2)= J
               GOTO 87
             END IF
             id%NBSA_LOCAL = J
             id%MEM_SUBTREE(1:J)=TEMP_MEM(1:J,1)
           ENDIF
         ENDDO
         DEALLOCATE(TEMP_MEM)
       ELSE
         CALL MPI_RECV(id%NBSA_LOCAL,1,MPI_INTEGER,
     *                 MASTER,0,id%COMM,STATUS, IERR)
         IF(ASSOCIATED(id%MEM_SUBTREE))
     *     DEALLOCATE(id%MEM_SUBTREE)
         ALLOCATE(id%MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
         IF (allocok .ne.0) then
               info(1)= -13
               info(2)= id%NBSA_LOCAL
               GOTO 87
         END IF
         CALL MPI_RECV(id%MEM_SUBTREE(1),id%NBSA_LOCAL,
     *       MPI_INTEGER,MASTER,0,
     *       id%COMM,STATUS,IERR)
       ENDIF
      ELSE
         IF(id%MYID.EQ.MASTER)THEN
            IF(KEEP(70).NE.2) DEALLOCATE(TEMP_MEM)
         ENDIF
         id%NBSA_LOCAL = -999999
         IF(ASSOCIATED(id%MEM_SUBTREE))
     *     DEALLOCATE(id%MEM_SUBTREE)
         ALLOCATE(id%MEM_SUBTREE(1),stat=allocok)
         IF (allocok .ne.0) then
               info(1)= -13
               info(2)= 1
               GOTO 87
         END IF
      ENDIF
87    CONTINUE
      CALL CMUMPS_276( ICNTL, INFO,
     *                      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( id%KEEP(55) /= 0 ) THEN
        CALL MPI_BCAST( id%FRTPTR(1), id%N+1, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
        CALL MPI_BCAST( id%FRTELT(1), id%NELT,  MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      END IF
      NB_NIV2 = KEEP(56) ! KEEP(1:100) was broadcast in CMUMPS_BDCAST
      IF (  NB_NIV2.GT.0  ) THEN
        if (id%myid.ne.master) then
        IF (ASSOCIATED(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
          allocate(PAR2_NODES(NB_NIV2),
     *             id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     *             STAT=allocok)
          IF (allocok .ne.0) then
             info(1)= -13
             info(2)= NB_NIV2*(id%NSLAVES+1)
          end if
        end if
        CALL CMUMPS_276( ICNTL, INFO,
     *                       id%COMM, id%MYID )
        IF ( INFO(1).LT.0 ) RETURN
        CALL MPI_BCAST(PAR2_NODES(1),NB_NIV2,
     *         MPI_INTEGER, MASTER, id%COMM, IERR )
        IF (KEEP(24) .NE.0 ) THEN
          CALL MPI_BCAST(id%CANDIDATES,(NB_NIV2*(id%NSLAVES+1)),
     *         MPI_INTEGER, MASTER, id%COMM, IERR )
        ENDIF
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        IF ( ASSOCIATED(id%ISTEP_TO_INIV2)) THEN
          DEALLOCATE(id%ISTEP_TO_INIV2)
          NULLIFY(id%ISTEP_TO_INIV2)
        ENDIF
        IF ( ASSOCIATED(id%TAB_POS_IN_PERE)) THEN
          DEALLOCATE(id%TAB_POS_IN_PERE)
          NULLIFY(id%TAB_POS_IN_PERE)
        ENDIF
        IF (NB_NIV2.EQ.0) THEN 
         id%KEEP(71) = 1
        ELSE
         id%KEEP(71) = id%KEEP(28)
        ENDIF
        ALLOCATE(id%ISTEP_TO_INIV2(id%KEEP(71)),
     *      id%TAB_POS_IN_PERE(id%NSLAVES+2,MAX(NB_NIV2,1)),
     *      stat=allocok)
        IF (allocok .gt.0) THEN
             info(1)= -13
             IF (NB_NIV2.EQ.0) THEN
              info(2)= 2
             ELSE
              info(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
             END IF
             GOTO 321
        ENDIF
        IF ( NB_NIV2 .GT.0 ) THEN
         DO INIV2 = 1, NB_NIV2
          INN = PAR2_NODES(INIV2)
          id%ISTEP_TO_INIV2(id%STEP(INN)) = INIV2
         END DO 
        ENDIF
      END IF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (ASSOCIATED(id%FUTURE_NIV2)) THEN
        DEALLOCATE(id%FUTURE_NIV2)
        NULLIFY(id%FUTURE_NIV2)
      ENDIF
      ALLOCATE(id%FUTURE_NIV2(id%NSLAVES), stat=allocok)
      IF (allocok .gt.0) THEN
           info(1)= -13
           info(2)= id%NSLAVES
           GOTO 321
      ENDIF
      id%FUTURE_NIV2=0
      DO INIV2 = 1, NB_NIV2
        IDEST = CMUMPS_275(id%STEP(PAR2_NODES(INIV2)),
     *                       id%PROCNODE_STEPS,
     *                       id%NSLAVES)
        id%FUTURE_NIV2(IDEST+1)=id%FUTURE_NIV2(IDEST+1)+1
      ENDDO
#endif
      IF (NB_NIV2.GT.0) DEALLOCATE (PAR2_NODES)
  321 CONTINUE
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( KEEP(23).NE.0 .and. id%MYID .EQ. MASTER ) THEN
        IKEEP = id%N + 1
      ELSE
        IKEEP = 1
      END IF
      FILS   = IKEEP + 3 * id%N
        NE   = IKEEP + 2 * id%N
        NA   = IKEEP +     id%N
      FRERE  = FILS  + id%N
      PTRAR  = FRERE + id%N
      IF (KEEP(55) .EQ. 0) THEN
        IF ( id%MYID.EQ.MASTER ) THEN
          NFSIZ   = PTRAR  + 4 * id%N
        ELSE
          NFSIZ   = PTRAR  + 2 * id%N
        ENDIF
      ELSE
        NFSIZ   = PTRAR  + 2 * (NELT + 1)
      END IF
      IF ( KEEP(38) .NE. 0 ) THEN
          CALL CMUMPS_164( id%MYID,
     *    id%NSLAVES, id%N, id%root,
     *    id%COMM_NODES, KEEP( 38 ), id%FILS(1),
     *    id%KEEP(50), id%KEEP(46),
     *    id%KEEP(51) )
      ELSE
          id%root%yes = .FALSE.
      END IF
      IF ( KEEP(38) .NE. 0 .and. I_AM_SLAVE ) THEN
          CALL MPI_ALLREDUCE(id%root%MYROW, MYROW_CHECK, 1,
     *         MPI_INTEGER, MPI_MAX, id%COMM_NODES, IERR)
          IF ( MYROW_CHECK .eq. -1) THEN
            INFO(1) = -25
            INFO(2) = 0
          END IF
          IF ( id%root%MYROW .LT. -1 .OR.
     *         id%root%MYCOL .LT. -1 ) THEN
            INFO(1) = -25
            INFO(2) = 0
          END IF
          IF ( LP > 0 .AND. INFO(1) == -25 ) THEN
            WRITE(LP, '(A)') 'Problem with your version of the BLACS.'
            WRITE(LP, '(A)') 'Try using a BLACS version from netlib.'
          ENDIF
      END IF
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF (INFO(1) .LT. 0) RETURN
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( I_AM_SLAVE ) THEN
        IF (KEEP(55) .EQ. 0) THEN
          CALL CMUMPS_24( id%MYID,
     *        id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     *        id%STEP(1), id%PTRAR(1),
     *        id%PTRAR(id%N +1),
     *        id%KEEP(1), ICNTL, id )
        ELSE
          CALL CMUMPS_25( id%MYID,
     *        id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     *        id%STEP(1),
     *        id%PTRAR(1),
     *        id%PTRAR(id%NELT+2 ),
     *        id%NELT, 
     *        id%FRTPTR(1), id%FRTELT(1),
     *        id%KEEP(1), ICNTL, id%SYM )
        ENDIF
        IF ( id%root%yes ) THEN
          LOCAL_M = NUMROC( id%ND_STEPS(id%STEP(KEEP(38))),
     *              id%root%MBLOCK, id%root%MYROW, 0,
     *              id%root%NPROW )
          LOCAL_M = MAX(1, LOCAL_M)
          LOCAL_N = NUMROC( id%ND_STEPS(id%STEP(KEEP(38))),
     *              id%root%NBLOCK, id%root%MYCOL, 0,
     *              id%root%NPCOL )
        ELSE
          LOCAL_M = 0
          LOCAL_N = 0
        END IF
        IF ( .NOT. ASSOCIATED(id%CANDIDATES)) THEN
          ALLOCATE(id%CANDIDATES(id%NSLAVES+1,1))
        ENDIF
        NB_NIV2 = MAX(KEEP(56), 1) 
        IF (KEEP(50).EQ.0) THEN
          CALL CMUMPS_246( id%MYID_NODES, id%N,
     *    id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     *    id%NA, id%LNA, id%NE_STEPS(1),
     *    id%ND_STEPS(1), id%PROCNODE_STEPS(1),
     *    id%NSLAVES,
     *    INFO( 3 ), INFO( 4 ), KEEP(15),
     *    KEEP(16), INFO(8),INFO(7),INFO(5), RINFO(1),
     *    KEEP(12), KEEP(1), LOCAL_M, LOCAL_N,
     *    SBUF_SEND, SBUF_REC, id%COST_SUBTREES, KEEP(28),
     *    NB_NIV2, id%ISTEP_TO_INIV2, id%CANDIDATES, 
     *    INFO(1), INFO(2) )
        ELSE
          CALL CMUMPS_179( id%MYID_NODES, id%N, 
     *    id%STEP(1), id%FRERE_STEPS, id%FILS,
     *    id%NA, id%LNA, id%NE_STEPS,
     *    id%ND_STEPS, id%PROCNODE_STEPS,
     *    id%NSLAVES,
     *    INFO( 3 ), INFO( 4 ), KEEP(15),
     *    KEEP(16), INFO(8),INFO(7),INFO(5), RINFO(1),
     *    KEEP(12), KEEP(1), LOCAL_M, LOCAL_N,
     *    SBUF_SEND, SBUF_REC, id%COST_SUBTREES, KEEP(28),
     *    NB_NIV2, id%ISTEP_TO_INIV2, id%CANDIDATES, 
     *    INFO(1), INFO(2) )
        ENDIF
        K44MAX = 0
        CALL MPI_ALLREDUCE (SBUF_REC, K44MAX, 1, 
     *    MPI_INTEGER, MPI_MAX, 
     *    id%COMM_NODES, IERR)
        SBUF_REC  = K44MAX
        KEEP(66) = INT(RINFO(1)/1000000)
        KEEP( 43 ) = SBUF_SEND
        KEEP( 44 ) = SBUF_REC
        KEEP( 25 ) = INFO( 3 )
        KEEP( 26 ) = INFO( 4 )
        KEEP( 27 ) = INFO( 5 )
        KEEP( 29 ) = INFO( 7 )
        IF (KEEP(50).eq.1) THEN
          INFO(8) = KEEP(16) + KEEP(12) * ( KEEP(16) / 100 + 1 )
     &              + KEEP(44) + KEEP(12) * (KEEP(44) / 100 )
        ELSE
          INFO(8) = KEEP(16) + KEEP(12) * ( KEEP(16) / 100 + 1 )
     &              + KEEP(44) + 4 * KEEP(12) * (KEEP(44) / 100 )
        ENDIF
        KEEP( 30 ) = INFO( 8 )
        IF ( MP .GT. 0 ) THEN
          WRITE(MP,'(A,I10) ') 
     *    ' estimated INTEGER space for factors         :',
     *                KEEP(26)
          WRITE(MP,'(A,I10) ') 
     *    ' estimated REAL    space for factors         :',
     *                KEEP(25)
          WRITE(MP,'(A,I10) ') 
     *    ' Minimum value of MAXIS for factorization    :',
     *                KEEP(15)
          WRITE(MP,'(A,I10) ') 
     *    ' Minimum value of MAXS  for factorization    :',
     *                KEEP(16)+KEEP(44)
          WRITE(MP,'(A,I10) ') 
     *    ' Current value of space relaxation parameter :',
     *                KEEP(12)
          WRITE(MP,'(A,I10) ') 
     *    ' Relaxed Value of MAXIS for factorization    :',
     *                KEEP(29)
          WRITE(MP,'(A,I10) ') 
     *    ' Relaxed Value of MAXS  for factorization    :',
     *                KEEP(30)
          WRITE(MP,'(A,I10) ') 
     *    ' Estimated maximum SENT message size (bytes) :',
     *                KEEP(43) * KEEP(35)
          WRITE(MP,'(A,I10) ') 
     *    ' Estimated maximum RECV message size (bytes) :',
     *                KEEP(44) * KEEP(35)
        END IF
      ELSE
       KEEP(30) = 0
       KEEP(29) = 0
      END IF
      CALL CMUMPS_100(id%MYID, id%COMM, KEEP, 
     *     INFO, INFOG, RINFO, RINFOG, ICNTL)
      IF ( id%MYID. eq. MASTER .AND. KEEP(54) .eq. 1 ) THEN
        IF (ASSOCIATED( id%MAPPING))
     *        DEALLOCATE( id%MAPPING)
        ALLOCATE( id%MAPPING(id%NZ), stat=allocok)
        IF ( allocok .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = id%NZ
          IF ( LP .GT. 0 ) THEN
            WRITE(LP, '(A)') 'Could not allocate mapping'
          END IF
          GOTO 92
        END IF
        ALLOCATE(IWtemp( id%N ), stat=allocok)
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=id%N
          IF ( LP .GT. 0 ) THEN
            WRITE(LP, '(A)') 'Could not allocate a workspace of size N'
          END IF
          GOTO 92
        END IF
        CALL CMUMPS_83(
     *    id%N, id%MAPPING(1),
     *    id%NZ, id%IRN(1),id%JCN(1), id%PROCNODE_STEPS(1),
     *    id%STEP(1),
     *    id%NSLAVES, id%SYM_PERM(1),
     *    id%FILS, IWtemp, id%KEEP(1),
     *    id%root%MBLOCK, id%root%NBLOCK,
     *    id%root%NPROW, id%root%NPCOL )
        DEALLOCATE( IWtemp )
 92     CONTINUE
      END IF
      CALL CMUMPS_276( ICNTL, INFO,
     *                  id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      CALL CMUMPS_214( id%KEEP(1), id%MYID, id%N, id%NELT, id%LNA,
     *                    id%NSLAVES, MEMORY_SIZE, .FALSE. )
      IF ( MP .gt. 0 ) THEN
          WRITE(MP,'(A,I10) ')
     *    ' Estimated space in MBYTES for factorization  :',
     *                MEMORY_SIZE
      END IF
      id%INFO(15) = MEMORY_SIZE
      CALL CMUMPS_243( id%MYID, id%COMM,
     *                           id%INFO(15), id%INFOG(16), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     *' ** Rank of processor needing largest memory in facto        :',
     *  IRANK
        WRITE( MPG,'(A,I10) ')
     *' ** Estimated corresponding space in MBYTES for facto        :',
     *  id%INFOG(16)
        IF ( KEEP(46) .eq. 0 ) THEN
         WRITE( MPG,'(A,I10) ')
     *' ** Estimated avg. space in MBYTES per working proc at facto :'
     *  ,(id%INFOG(17)-id%INFO(15))/id%NSLAVES
        ELSE
         WRITE( MPG,'(A,I10) ')
     *' ** Estimated avg. space in MBYTES per working proc at facto :'
     *  ,id%INFOG(17)/id%NSLAVES
        END IF
        WRITE(MPG,'(A,I10) ')
     *' ** TOTAL     space in MBYTES for factorization              :'
     * ,id%INFOG(17)
      END IF
      RETURN
 110  FORMAT(/' ****** ANALYSIS STEP ********'/)
      END SUBROUTINE CMUMPS_26
      SUBROUTINE CMUMPS_151( N, KEEP28, IWCB, LIWW, W, LWC,
     *       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER N,LIWW,LWC,POSWCB,IWPOSCB, KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28),PTRACB(KEEP28)
      COMPLEX W(LWC)
      INTEGER SIZFI, SIZFR
      IF ( IWPOSCB .eq. LIWW ) RETURN
      DO WHILE ( IWCB( IWPOSCB + 2 ) .eq. 0 )
        SIZFR = IWCB( IWPOSCB + 1 )
        SIZFI = SIZFR + 2
        IWPOSCB = IWPOSCB + SIZFI
        POSWCB  = POSWCB  + SIZFR
        IF ( IWPOSCB .eq. LIWW ) RETURN
      END DO
      RETURN
      END SUBROUTINE CMUMPS_151
      SUBROUTINE CMUMPS_95(N,KEEP28,IWCB,LIWW,W,LWC,
     *       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER N,LIWW,LWC,POSWCB,IWPOSCB,KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28),PTRACB(KEEP28)
      COMPLEX W(LWC)
      INTEGER IPTIW,IPTA,SIZFI,SIZFR,LONGI,LONGR
      INTEGER I
      IPTIW = IWPOSCB
      IPTA  = POSWCB
      LONGI = 0
      LONGR = 0
      IF ( IPTIW .EQ. LIWW ) RETURN
10    CONTINUE
      IF (IWCB(IPTIW+2).EQ.0) THEN
        SIZFR  = IWCB(IPTIW+1)
        SIZFI  = SIZFR + 2
        IF (LONGI.NE.0) THEN
          DO 20 I=0,LONGI-1
            IWCB(IPTIW + SIZFI - I) = IWCB (IPTIW - I )
 20       CONTINUE 
          DO 30 I=0,LONGR-1
            W(IPTA + SIZFR - I)   = W(IPTA - I )
 30       CONTINUE
        ENDIF
        DO 40 I=1,KEEP28
          IF ((PTRICB(I).LE.(IPTIW+1)).AND.
     &        (PTRICB(I).GT.IWPOSCB) ) THEN
            PTRICB(I) = PTRICB(I) + SIZFI
            PTRACB(I) = PTRACB(I) + SIZFR
          ENDIF 
40      CONTINUE 
        IWPOSCB = IWPOSCB + SIZFI
        IPTIW   = IPTIW + SIZFI
        POSWCB = POSWCB + SIZFR
        IPTA   = IPTA + SIZFR     
       ELSE
        SIZFR  = IWCB(IPTIW+1)
        SIZFI  = SIZFR + 2
        IPTIW = IPTIW + SIZFI
        LONGI = LONGI + SIZFI
        IPTA  = IPTA + SIZFR
        LONGR = LONGR + SIZFR
       ENDIF
       IF (IPTIW.NE.LIWW) GOTO 10
       RETURN
       END SUBROUTINE CMUMPS_95
      SUBROUTINE CMUMPS_205(MTYPE, IFLAG, N, NZ,
     *    LHS, WRHS, W, RHS, GIVSOL, SOL, ANORM, XNORM, SCLNRM,
     *    MPRINT, ICNTL, KEEP)
      INTEGER MTYPE,N,NZ,IFLAG,ICNTL(40), KEEP(150)
      COMPLEX RHS(N),LHS(N)
      COMPLEX WRHS(N),SOL(*)
      COMPLEX ERL2,ERREL,SCLNRM,W(N)
      REAL RESMAX,RESL2,XNORM, ERMAX,MAXSOL,
     *     COMAX
      REAL ANORM,DZERO,EPSI
      LOGICAL GIVSOL,PROK
      INTEGER MPRINT, MP
      INTEGER I,J,K
      INTRINSIC ABS, MAX, SQRT
      MP = ICNTL(2)
      PROK = (MPRINT .GT. 0)
      DZERO = 0.0E0
      EPSI = 0.1E-9
      ANORM = DZERO
      RESMAX = DZERO
      RESL2 = DZERO
      DO 40 K = 1, N
        RESMAX = MAX(RESMAX, ABS(RHS(K)))
        RESL2 = RESL2 + ABS(RHS(K)) * ABS(RHS(K))
        ANORM = MAX(ANORM, ABS(W(K)))
   40 CONTINUE
      XNORM = DZERO
      DO 50 K = 1, N
        XNORM = MAX(XNORM, ABS(LHS(K)))
   50 CONTINUE
      IF (XNORM .GT. EPSI) THEN
        SCLNRM = RESMAX / (ANORM * XNORM)
      ELSE
        IFLAG = IFLAG + 2
        IF ((MP .GT. 0) .AND. (ICNTL(4) .GE. 2)) WRITE( MP, * ) 
     *' MAX-NORM of computed solut. is zero'
        SCLNRM = RESMAX / ANORM
      ENDIF
      RESL2 = SQRT(RESL2)
      ERMAX = DZERO
      COMAX = DZERO
      ERL2 = DZERO
      IF (.NOT.GIVSOL) THEN
        IF (PROK) WRITE( MPRINT, 90 ) RESMAX, RESL2, ANORM, XNORM, 
     *      SCLNRM
      ELSE
        MAXSOL = DZERO
        DO 60 K = 1, N
          MAXSOL = MAX(MAXSOL, ABS(SOL(K)))
   60   CONTINUE
        DO 70 K = 1, N
          ERL2 = (LHS(K) - SOL(K)) ** 2 + ERL2
          ERMAX = MAX(ERMAX, ABS(LHS(K) - SOL(K)))
   70   CONTINUE
        DO 80 K = 1, N
          IF (ABS(SOL(K)) .GT. EPSI) THEN
            COMAX = MAX(COMAX, (ABS(LHS(K) - SOL(K)) / ABS(SOL(K))))
          ENDIF
   80   CONTINUE
        ERL2 = SQRT(ERL2)
        IF (MAXSOL .GT. EPSI) THEN
          ERREL = ERMAX / MAXSOL
        ELSE
          IFLAG = IFLAG + 2
          IF ((MP .GT. 0) .AND. (ICNTL(4) .GE. 2)) WRITE( MP, * ) 
     *' MAX-NORM of exact solution is zero'
          ERREL = ERMAX
        ENDIF
        IF (PROK) WRITE( MPRINT, 100 ) ERMAX, ERL2, ERREL, COMAX, RESMAX
     *      , RESL2, ANORM, XNORM, SCLNRM
      ENDIF
   90  FORMAT (/' RESIDUAL IS ............ (MAX-NORM)        =',1PD9.2/
     *       '                       .. (2-NORM)          =',1PD9.2/
     *       ' RINFOG(4):NORM OF input  Matrix  (MAX-NORM)=',1PD9.2/
     *       ' RINFOG(5):NORM OF Computed SOLUT (MAX-NORM)=',1PD9.2/
     *       ' RINFOG(6):SCALED RESIDUAL ...... (MAX-NORM)=',1PD9.2)
      RETURN
  100  FORMAT (/' ERROR IS     ............ (MAX-NORM)       =',1PD9.2/
     *       '              ............ (2-NORM)         =',1PD9.2/
     *       ' RELATIVE ERROR........... (MAX-NORM)       =',1PD9.2/
     *       ' Comp. Wise ERROR......... (MAX-NORM)       =',1PD9.2/
     *       ' AND RESIDUAL IS ......... (MAX-NORM)       =',1PD9.2/
     *       '                        .. (2-NORM)         =',1PD9.2/
     *       ' NORM OF input  MATRIX ... (MAX-NORM)       =',1PD9.2/
     *       ' NORM of computed SOLUT... (MAX-NORM)       =',1PD9.2/
     *       ' SCALED RESIDUAL ......... (MAX-NORM)       =',1PD9.2)
      END SUBROUTINE CMUMPS_205
      SUBROUTINE CMUMPS_206(NZ, N, RHS,
     *    X, Y, D, W, IW, KASE,
     *    OMEGA, ERX, JOB, COND, MAXIT, NOITER, LP, KEEP, ARRET )
      INTEGER NZ, N, KASE, KEEP(150), JOB
      INTEGER IW(N,2)
      COMPLEX RHS(N)
      COMPLEX X(N), Y(N), D(N), W(N,3)
      INTEGER LP, MAXIT, NOITER
      REAL COND(2),OMEGA(2)
      REAL ARRET
      REAL CGCE, CTAU
      DATA  CTAU /1.0E3/, CGCE /0.2E0/
      LOGICAL LCOND1, LCOND2
      INTEGER IFLAG, JUMP, I, IMAX
      COMPLEX ERX, DXMAX
      REAL CONVER, OM1, OM2, DXIMAX
      REAL  ZERO, ONE,TAU, DD
      COMPLEX OLDOMG(2)
      INTEGER CMUMPS_ICAMAX
      INTRINSIC     ABS, MAX
      SAVE LCOND1, LCOND2, JUMP,  DXIMAX, DXMAX, CONVER,
     *     OM1, OLDOMG, IFLAG
      DATA ZERO /0.0E0/, ONE /1.0E0/
      IF (KASE .EQ. 0) THEN
        LCOND1 = .FALSE.
        LCOND2 = .FALSE.
        COND(1) = ONE
        COND(2) = ONE
        ERX = CMPLX(ZERO)
        OM1 = ZERO
        IFLAG = 0
        NOITER = 0
        JUMP = 1
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 30
      CASE(2)
        GOTO 10
      CASE(3)
        GOTO 110
      CASE(4)
        GOTO 150
      CASE(5)
        GOTO 35
      CASE DEFAULT
      END SELECT
   10 CONTINUE
      DO 20 I = 1, N
        X(I) = X(I) + Y(I)
   20 CONTINUE
      IF (NOITER .GT. MAXIT) THEN
        IFLAG = IFLAG + 8
        GOTO 70
      ENDIF
   30 CONTINUE
      KASE = 14
      JUMP = 5
      RETURN
   35 CONTINUE
      IMAX = CMUMPS_ICAMAX(N, X, 1)
      DXMAX = ABS(X(IMAX))
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO 40 I = 1, N
        TAU = (W(I, 2) * DXMAX + ABS(RHS(I))) * N * CTAU
        DD = W(I, 1) + ABS(RHS(I))
        IF ((DD + TAU) .GT. TAU) THEN
          OMEGA(1) = MAX(OMEGA(1), ABS(Y(I) / DD))
          IW(I, 1) = 1
        ELSE
          IF (TAU .GT. ZERO) THEN
            OMEGA(2) = MAX(OMEGA(2), ABS(Y(I) / (DD + W(I, 2) * DXMAX)))
          ENDIF
          IW(I, 1) = 2
        ENDIF
   40 CONTINUE
      OM2 = OMEGA(1) + OMEGA(2)
      IF (OM2 .LT. ARRET ) GOTO 70
      IF (MAXIT .EQ. 0) GOTO 70
      IF (NOITER .GT. 1 .AND. OM2 .GT. OM1 * CGCE) THEN
        CONVER = OM2 / OM1
        IF (OM2 .GT. OM1) THEN
          OMEGA(1) = OLDOMG(1)
          OMEGA(2) = OLDOMG(2)
          DO 50 I = 1, N
            X(I) = W(I, 3)
   50     CONTINUE
        ENDIF
        GOTO 70
      ENDIF
      DO 60 I = 1, N
        W(I, 3) = X(I)
   60 CONTINUE
      OLDOMG(1) = OMEGA(1)
      OLDOMG(2) = OMEGA(2)
      OM1 = OM2
      NOITER = NOITER + 1
      KASE = 2
      JUMP = 2
      RETURN
   70 KASE = 0
      IF (JOB .LE. 0) GOTO 170
      DO 80 I = 1, N
        IF (IW(I, 1) .EQ. 1) THEN
          W(I, 1) = W(I, 1) + ABS(RHS(I))
          W(I, 2) = CMPLX(ZERO)
          LCOND1 = .TRUE.
        ELSE
          W(I, 2) = W(I, 2) * DXMAX + W(I, 1)
          W(I, 1) = CMPLX(ZERO)
          LCOND2 = .TRUE.
        ENDIF
   80 CONTINUE
      DO 90 I = 1, N
        W(I, 3) = X(I) * D(I)
   90 CONTINUE
      IMAX = CMUMPS_ICAMAX(N, W(1, 3), 1)
      DXIMAX = ABS(W(IMAX, 3))
      IF (.NOT.LCOND1) GOTO 130
  100 CALL CMUMPS_218(N, KASE, Y, COND(1), W(1, 3), IW(1, 2))
      IF (KASE .EQ. 0) GOTO 120
      IF (KASE .EQ. 1) CALL CMUMPS_204(N, Y, D)
      IF (KASE .EQ. 2) CALL CMUMPS_204(N, Y, W)
      JUMP = 3
      RETURN
  110 CONTINUE
      IF (KASE .EQ. 1) CALL CMUMPS_204(N, Y, W)
      IF (KASE .EQ. 2) CALL CMUMPS_204(N, Y, D)
      GOTO 100
  120 IF (DXIMAX .GT. ZERO) COND(1) = COND(1) / DXIMAX
      ERX = OMEGA(1) * COND(1)
  130 IF (.NOT.LCOND2) GOTO 170
      KASE = 0
  140 CALL CMUMPS_218(N, KASE, Y, COND(2), W(1, 3), IW(1, 2))
      IF (KASE .EQ. 0) GOTO 160
      IF (KASE .EQ. 1) CALL CMUMPS_204(N, Y, D)
      IF (KASE .EQ. 2) CALL CMUMPS_204(N, Y, W(1, 2))
      JUMP = 4
      RETURN
  150 CONTINUE
      IF (KASE .EQ. 1) CALL CMUMPS_204(N, Y, W(1, 2))
      IF (KASE .EQ. 2) CALL CMUMPS_204(N, Y, D)
      GOTO 140
  160 IF (DXIMAX .GT. ZERO) THEN
        COND(2) = COND(2) / DXIMAX
      ENDIF
      ERX = ERX + OMEGA(2) * COND(2)
  170 KASE = -IFLAG
      RETURN
      END SUBROUTINE CMUMPS_206
      SUBROUTINE CMUMPS_207(A, NZ, N, IRN, ICN, Z, KEEP)
      INTEGER NZ, N, I, J, K, KEEP(150)
      INTEGER IRN(NZ), ICN(NZ)
      COMPLEX A(NZ), Z(N)
      REAL ZERO
      INTRINSIC     ABS
      DATA ZERO /0.0E0/
      DO 10 I = 1, N
        Z(I) = CMPLX(ZERO)
   10 CONTINUE
      IF (KEEP(50) .EQ.0) THEN
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K))
       ENDDO
      ELSE
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K))
        IF (J.NE.I) THEN 
          Z(J) = Z(J) + ABS(A(K))
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_207
      SUBROUTINE CMUMPS_289(A, NZ, N, IRN, ICN, Z, KEEP, 
     &             COLSCA)
      INTEGER NZ, N, I, J, K, KEEP(150)
      INTEGER IRN(NZ), ICN(NZ)
      COMPLEX A(NZ), Z(N), COLSCA(N)
      REAL  ZERO
      INTRINSIC     ABS
      DATA ZERO /0.0E0/
      DO 10 I = 1, N
        Z(I) = CMPLX(ZERO)
   10 CONTINUE
      IF (KEEP(50) .EQ.0) THEN
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K)*COLSCA(J))
       ENDDO
      ELSE
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K)*COLSCA(J))
        IF (J.NE.I) THEN
          Z(J) = Z(J) + ABS(A(K)*COLSCA(I))
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_289
      SUBROUTINE CMUMPS_208(A, NZ, N, IRN, ICN, RHS, X, R, W, KEEP)
      IMPLICIT NONE
      INTEGER NZ, N, I, K, J, KEEP(150)
      INTEGER IRN(NZ), ICN(NZ)
      COMPLEX A(NZ), RHS(N), X(N), R(N), W(N),
     *     D
      REAL ZERO
      INTRINSIC    ABS
      DATA ZERO /0.0E0/
      DO I = 1, N
        R(I) = RHS(I)
        W(I) = CMPLX(ZERO)
      ENDDO
      DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .GT. N) .OR. (J .GT. N) .OR. (I .LT. 1) .OR. (J .LT. 1))
     *      CYCLE
        D = A(K) * X(J)
        R(I) = R(I) - D
        W(I) = W(I) + ABS(D)
        IF ((I.NE.J) .AND. (KEEP(50).NE.0) ) THEN
          D = A(K) * X(I)
          R(J) = R(J) - D
          W(J) = W(J) + ABS(D)
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_208
      SUBROUTINE CMUMPS_204(N, R, W)
      INTEGER N, I
      COMPLEX  R(N), W(N)
      DO 10 I = 1, N
        R(I) = R(I) * W(I)
   10 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_204
      SUBROUTINE CMUMPS_218(N, KASE, X, EST, W, IW)
      INTEGER N, IW(N), KASE
      COMPLEX W(N), X(N)
      REAL EST
      INTRINSIC ABS, NINT, REAL, SIGN
      INTEGER CMUMPS_ICAMAX
      EXTERNAL CMUMPS_ICAMAX
      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      INTEGER I, ITER, J, JLAST, JUMP
      REAL ALTSGN
      REAL ZERO, ONE, TEMP
      SAVE ITER, J, JLAST, JUMP
      DATA ZERO /0.0E0/
      DATA ONE  /1.0E0/
      IF (KASE .EQ. 0) THEN
        DO 10 I = 1, N
          X(I) = ONE / REAL(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        RETURN
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 20
      CASE(2)
        GOTO 40
      CASE(3)
        GOTO 70
      CASE(4)
        GOTO 120
      CASE(5)
        GOTO 160
      CASE DEFAULT
      END SELECT
   20 CONTINUE
      IF (N .EQ. 1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
        GOTO 190
      ENDIF
      DO 30 I = 1, N
        X(I) = SIGN(ONE,ABS( X(I)))
        IW(I) = NINT(ABS(X(I)))
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
   40 CONTINUE
      J = CMUMPS_ICAMAX(N, X, 1)
      ITER = 2
   50 CONTINUE
      DO 60 I = 1, N
        X(I) = CMPLX(ZERO)
   60 CONTINUE
      X(J) = CMPLX(ONE)
      KASE = 1
      JUMP = 3
      RETURN
   70 CONTINUE
      DO 80 I = 1, N
        W(I) = X(I)
   80 CONTINUE
      DO 90 I = 1, N
        IF (NINT(SIGN(ONE, ABS(X(I)))) .NE. IW(I)) GOTO 100
   90 CONTINUE
      GOTO 130
  100 CONTINUE
      DO 110 I = 1, N
        X(I) = SIGN(ONE, ABS(X(I)))
        IW(I) = NINT(ABS(X(I)))
  110 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
  120 CONTINUE
      JLAST = J
      J = CMUMPS_ICAMAX(N, X, 1)
      IF ((ABS(X(JLAST)) .NE. ABS(X(J))) .AND. (ITER .LT. ITMAX)) THEN
        ITER = ITER + 1
        GOTO 50
      ENDIF
  130 CONTINUE
      EST = CMPLX(ZERO)
      DO 140 I = 1, N
        EST = EST + ABS(W(I))
  140 CONTINUE
      ALTSGN = ONE
      DO 150 I = 1, N
        X(I) = ALTSGN * (ONE + REAL(I - 1) / REAL(N - 1))
        ALTSGN = -ALTSGN
  150 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
  160 CONTINUE
      TEMP = ZERO
      DO 170 I = 1, N
        TEMP = TEMP + ABS(X(I))
  170 CONTINUE
      TEMP = 2.0 * TEMP / REAL(3 * N)
      IF (TEMP .GT. EST) THEN
        DO 180 I = 1, N
          W(I) = X(I)
  180   CONTINUE
        EST = TEMP
      ENDIF
  190 KASE = 0
      RETURN
      END SUBROUTINE CMUMPS_218
      SUBROUTINE CMUMPS_278( MTYPE, N, NZ, ASPK, IRN, ICN,
     *    LHS, WRHS, W, RHS, KEEP)
      IMPLICIT NONE
      INTEGER MTYPE, N, NZ
      INTEGER IRN( NZ ), ICN( NZ )
      INTEGER KEEP(150)
      COMPLEX ASPK( NZ )
      COMPLEX LHS( N ), WRHS( N ), W( N ), RHS( N )
      INTEGER K, I, J
      REAL DZERO
      PARAMETER(DZERO = 0.0E0)
      DO 10 K = 1, N
        W(K) = DZERO
        RHS(K) = WRHS(K)
   10 CONTINUE
      IF ( KEEP(50) .EQ. 0 ) THEN
       IF (MTYPE .EQ. 1) THEN
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          RHS(I) = RHS(I) - ASPK(K) * LHS(J)
          W(I) = W(I) + ABS(ASPK(K))
        ENDDO
       ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          RHS(J) = RHS(J) - ASPK(K) * LHS(I)
          W(J) = W(J) + ABS(ASPK(K))
        ENDDO
       ENDIF
      ELSE
       DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          RHS(I) = RHS(I) - ASPK(K) * LHS(J)
          W(I) = W(I) + ABS(ASPK(K))
          IF (J.NE.I) THEN
            RHS(J) = RHS(J) - ASPK(K) * LHS(I)
            W(J) = W(J) + ABS(ASPK(K))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_278
      SUBROUTINE CMUMPS_121( MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    LHS, WRHS, W, RHS, KEEP)
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(150)
      COMPLEX A_ELT(NA_ELT)
      COMPLEX LHS( N ), WRHS( N ), W( N ), RHS( N )
      CALL CMUMPS_257(N, NELT, ELTPTR, ELTVAR, A_ELT,
     *                         LHS, RHS, KEEP(50), MTYPE )
      RHS = WRHS - RHS
      CALL CMUMPS_119( MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    W, KEEP )
      RETURN
      END SUBROUTINE CMUMPS_121
      SUBROUTINE CMUMPS_119( MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    W, KEEP )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(150)
      COMPLEX A_ELT(NA_ELT), W(N), TEMP
      INTEGER K, I, J, IEL, SIZEI, IELPTR
      REAL DZERO
      PARAMETER(DZERO = 0.0E0)
      W = CMPLX(DZERO)
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) = 
     &           W( ELTVAR( IELPTR + I) )
     &           + ABS(A_ELT( K ))
               K = K + 1
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
               TEMP = TEMP + ABS( A_ELT(K))
               K = K + 1
              END DO
              W(ELTVAR( IELPTR + J )) = 
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) = 
     &        W(ELTVAR( IELPTR + J )) + ABS(A_ELT( K ))
          K = K + 1
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) = 
     &           W(ELTVAR( IELPTR + J )) + ABS(A_ELT( K ))
              W(ELTVAR( IELPTR + I ) ) = 
     &           W(ELTVAR( IELPTR + I )) + ABS(A_ELT( K ))
              K = K + 1
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_119
      SUBROUTINE CMUMPS_135(MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    W, KEEP, COLSCA )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(150)
      COMPLEX A_ELT(NA_ELT), W(N), COLSCA(N), 
     *                 TEMP, TEMP2
      INTEGER K, I, J, IEL, SIZEI, IELPTR
      REAL DZERO
      PARAMETER(DZERO = 0.0E0)
      W = DZERO
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              TEMP2 = COLSCA(ELTVAR( IELPTR + J) )
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) =
     &           W( ELTVAR( IELPTR + I) )
     &           + ABS(A_ELT( K )* TEMP2 )
               K = K + 1
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              TEMP2= COLSCA(ELTVAR( IELPTR + J) )
              DO I = 1, SIZEI
               TEMP = TEMP + ABS( A_ELT(K)*TEMP2)
               K = K + 1
              END DO
              W(ELTVAR( IELPTR + J )) =
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) =
     &        W(ELTVAR( IELPTR + J )) + 
     &        ABS( A_ELT( K )*COLSCA(ELTVAR( IELPTR + J)) )
          K = K + 1
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) =
     &           W(ELTVAR( IELPTR + J )) + 
     &           ABS(A_ELT( K )*COLSCA(ELTVAR( IELPTR + J)))
              W(ELTVAR( IELPTR + I ) ) =
     &           W(ELTVAR( IELPTR + I )) + 
     &           ABS(A_ELT( K )*COLSCA(ELTVAR( IELPTR + I)))
              K = K + 1
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_135
      SUBROUTINE CMUMPS_122( MTYPE, N, NELT, ELTPTR, 
     *                     LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *                     SAVERHS, X, Y, W, K50 )
      IMPLICIT NONE
      INTEGER N, NELT, K50, MTYPE, LELTVAR, NA_ELT
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( LELTVAR )
      COMPLEX A_ELT( NA_ELT ), X( N ), Y( N ), 
     *                 SAVERHS(N), W(N)
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      REAL ZERO
      COMPLEX TEMP, TEMP2
      PARAMETER( ZERO = 0.0E0 )
      Y = SAVERHS
      W = CMPLX(ZERO)
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( K50 .eq. 0 ) THEN
          IF ( MTYPE .eq. 1 ) THEN
            DO J = 1, SIZEI
              TEMP = X( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                Y( ELTVAR( IELPTR + I ) ) =
     *          Y( ELTVAR( IELPTR + I ) ) -
     *             A_ELT( K ) * TEMP
                W( ELTVAR( IELPTR + I ) ) =
     *          W( ELTVAR( IELPTR + I ) ) +
     *             ABS( A_ELT( K ) * TEMP )
                K = K + 1
              END DO
            END DO
          ELSE
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              TEMP2 = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP - 
     *          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                TEMP2 = TEMP2 +  ABS(
     *          A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
              W( ELTVAR( IELPTR + J ) ) = TEMP2
            END DO
          END IF
        ELSE
          DO J = 1, SIZEI
            Y( ELTVAR( IELPTR + J ) ) =
     *      Y( ELTVAR( IELPTR + J ) ) -
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            W( ELTVAR( IELPTR + J ) ) =
     *      W( ELTVAR( IELPTR + J ) ) + ABS(
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
            K = K + 1
            DO I = J+1, SIZEI
              Y( ELTVAR( IELPTR + I ) ) =
     *        Y( ELTVAR( IELPTR + I ) ) -
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     *        Y( ELTVAR( IELPTR + J ) ) -
     *           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              W( ELTVAR( IELPTR + I ) ) =
     *        W( ELTVAR( IELPTR + I ) ) + ABS(
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
              W( ELTVAR( IELPTR + J ) ) =
     *        W( ELTVAR( IELPTR + J ) ) + ABS(
     *           A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE CMUMPS_122
      SUBROUTINE CMUMPS_251(N,IW,LIW,A,LA,
     *             NSTK_STEPS, NBPROCFILS,IFLAG,ND,FILS,STEP,
     *             FRERE,CAND, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *             MAXFRT, NTOTPV,PTRIST,PTRAST,
     *             PIMASTER, PAMASTER, PTRARW,PTRAIW,
     *             ITLOC, IERROR,IPOOL, LPOOL, NIRBDU,
     *             RINFO, POSFAC ,IWPOS, LRLU, IPTRLU, 
     *             LRLUS, LEAF, NBROOT, NBRTOT,
     *             UU, ICNTL, PTLUST_S, NSTEPS, INFO, KEEP,
     *             PROCNODE_STEPS,SLAVEF,MYID, COMM_NODES,
     *             MYID_NODES,
     *             BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,root,
     *             PERM, NELT, FRTPTR, FRTELT, LPTRAR, 
     *             COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     *             MEM_DISTRIB, NE )
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER N,IFLAG,NTOTPV,MAXFRT,LA,LIW, LPTRAR,
     *        IERROR,NIRBDU, NSTEPS, INFO(40)
      COMPLEX A(LA)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER KEEP(150), ICNTL(40)
      INTEGER LPOOL 
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER ITLOC(N)
      INTEGER IW(LIW), NSTK_STEPS(KEEP(28)), NBPROCFILS(KEEP(28))
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND(KEEP(28))
      INTEGER FILS(N),PTRIST(KEEP(28)),PTRAST(KEEP(28))
      INTEGER STEP(N), FRERE(KEEP(28))
      INTEGER PIMASTER(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PTLUST_S(KEEP(28)), PERM(N)
      INTEGER CAND(SLAVEF+1,MAX(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER IPOOL(LPOOL)
      INTEGER NE(KEEP(28))
      REAL RINFO(20)
      INTEGER POSFAC,IWPOS,LRLU,
     *    IPTRLU, LRLUS,
     *    LEAF, NBROOT 
      INTEGER COMM_LOAD, ASS_IRECV
      REAL UU, SEUIL, SEUIL_LDLT_NIV2
      INTEGER NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER        INTARR( MAX(1,KEEP(14)) )
      COMPLEX DBLARR( MAX(1,KEEP(13)) )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      LOGICAL FLAG
      REAL ZERO
      INTEGER INODE
      INTEGER IWPOSCB
      INTEGER FPERE, TYPEF
      INTEGER MP, LP, DUMMY(1)
      INTEGER NBFIN, NBRTOT, NBROOT_TRAITEES
      INTEGER NFRONT, IOLDPS, POSELT
      INTEGER LBUFRX
      COMPLEX, ALLOCATABLE, DIMENSION(:) :: BUFRX
      REAL FLOP1
      INTEGER TYPE
      LOGICAL SON_LEVEL2, SET_IRECV, BLOCKING, 
     &        MESSAGE_RECEIVED
      EXTERNAL CMUMPS_330, CMUMPS_275 
      INTEGER CMUMPS_330, CMUMPS_275
      LOGICAL CMUMPS_167
      EXTERNAL CMUMPS_167
      LOGICAL CMUMPS_508
      EXTERNAL CMUMPS_508, CMUMPS_509
      INTEGER MAXFRW, NPVW, NOFFW, NELVAW, COMP
      REAL    OPASSW, OPELIW
      DATA ZERO /0.0E0/
       ASS_IRECV = MPI_REQUEST_NULL
       ITLOC(1:N) =0
       PTRIST  (1:KEEP(28))=0
       PTLUST_S(1:KEEP(28))=0
       PTRAST(1:KEEP(28))=0
      MP = ICNTL(2)
      LP = ICNTL(1)
      MAXFRW = 0
      NPVW   = 0
      NOFFW  = 0
      NELVAW = 0
      COMP  = 0
      OPASSW = ZERO
      OPELIW = ZERO
      IWPOSCB = NIRBDU    
      NBFIN = NBRTOT
      NBROOT_TRAITEES = 0
      NBPROCFILS(1:KEEP(28)) = 0
      IF ( KEEP(38).NE.0 ) THEN
        IF (root%yes) THEN
          CALL CMUMPS_284( 
     *        root, KEEP(38), N, IW, LIW,
     *        A, LA,
     *        NIRBDU, FILS, MYID_NODES, PTRAIW, PTRARW,
     *        INTARR, DBLARR,
     *        LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC,
     *        COMP, LRLUS, IFLAG, KEEP, IERROR )
        ENDIF
        IF ( IFLAG .LT. 0 ) GOTO 635
      END IF
 20   CONTINUE
      NIV1_FLAG=0
      SET_IRECV = .TRUE.
      BLOCKING = .FALSE.
      MESSAGE_RECEIVED = .FALSE.
      CALL CMUMPS_329(
     *      COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV, 
     *      MESSAGE_RECEIVED,
     *      MPI_ANY_SOURCE, MPI_ANY_TAG,
     *      STATUS, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      PROCNODE_STEPS, POSFAC,
     *      IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA, 
     *      NIRBDU, PTRIST, PTLUST_S,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS,
     *      COMP, IFLAG,
     *      IERROR, COMM_NODES,
     *      NBPROCFILS,
     *      IPOOL, LPOOL, LEAF, NBFIN, MYID_NODES, SLAVEF,
     *
     *      root, OPASSW, OPELIW, ITLOC, FILS, 
     *      PTRARW, PTRAIW,
     *      INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *      LPTRAR, NELT, FRTPTR, FRTELT, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      CALL CMUMPS_467(COMM_LOAD, KEEP)
      IF (MESSAGE_RECEIVED) THEN
          IF ( IFLAG .LT. 0 ) GO TO 640
          IF ( NBFIN .eq. 0 ) GOTO 640
      ELSE
          IF ( .NOT. CMUMPS_508( IPOOL, LPOOL) )THEN
            CALL CMUMPS_509( N, IPOOL, LPOOL,
     *      PROCNODE_STEPS,
     *      SLAVEF, STEP, INODE, KEEP, MYID_NODES, ND )
            IF (KEEP(47) .GE. 3) THEN
              CALL CMUMPS_500(
     $              IPOOL, LPOOL, 
     *              PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *              MYID_NODES, STEP, N, ND, FILS )
            ENDIF
            IF (KEEP(47).EQ.4) THEN
              CALL CMUMPS_501(
     *        IPOOL,LPOOL,MYID_NODES,slavef,COMM_LOAD,KEEP)
            ENDIF
            GOTO 30
          ENDIF
      ENDIF
      GO TO 20
 30   CONTINUE
      IF ( INODE .LT. 0 ) THEN
        INODE = -INODE
        FPERE = INODE
        DO WHILE (FPERE.GT.0) 
          FPERE = FRERE(STEP(FPERE))
        END DO
        FPERE = -FPERE
        GOTO 130
      ELSE IF (INODE.GT.N) THEN
       INODE = INODE - N
       IF (INODE.EQ.KEEP(38)) THEN
         NBROOT_TRAITEES = NBROOT_TRAITEES + 1
         IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
            NBFIN = NBFIN - NBROOT
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL CMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID_NODES,
     *          COMM_NODES, RACINE, SLAVEF)
            END IF
         ENDIF
         IF (NBFIN.EQ.0) GOTO 640
         GOTO 20
       ENDIF
       TYPE = CMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
       IF (TYPE.EQ.1) GOTO 100
       FPERE = INODE
       DO WHILE (FPERE.GT.0) 
        FPERE = FRERE(STEP(FPERE))
       END DO
       FPERE = -FPERE
       IF ( KEEP(50) .eq. 0 ) THEN
         CALL  CMUMPS_144( COMM_LOAD, ASS_IRECV,
     *             N, INODE, FPERE, IW, LIW, A, LA,
     *             UU, NOFFW,
     *             NPVW,
     *             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     *             NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, STEP, PIMASTER, PAMASTER,
     *             NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     *             root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 640
       ELSE
         CALL  CMUMPS_141( COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, A, LA,
     *             UU, NOFFW,
     *             NPVW,
     *             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     *             NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, STEP, PIMASTER, PAMASTER,
     *             NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     *             root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL_LDLT_NIV2, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF ( IFLAG .LT. 0 ) GOTO 640
        IF ( IW( PTLUST_S(STEP(INODE)) + 5 ) .GT. 1 ) THEN
             GOTO 20
        END IF
       END IF
       GOTO 130
      ENDIF
      IF (INODE.EQ.KEEP(38)) THEN
         CALL  CMUMPS_176( COMM_LOAD, ASS_IRECV, 
     *    ROOT, FRERE,
     *    INODE,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS, COMP,
     *    IFLAG, IERROR, COMM_NODES,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID_NODES, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, ND,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        GOTO 20
      ENDIF
      TYPE = CMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF (TYPE.EQ.1) THEN
        IF (KEEP(55).NE.0) THEN
         CALL CMUMPS_36( COMM_LOAD, ASS_IRECV,
     *        NELT, FRTPTR, FRTELT,
     *        N,INODE,IW,LIW,A,LA,
     *        IFLAG,IERROR,ND,NIRBDU,
     *        FILS,FRERE,MAXFRW,root,OPASSW, OPELIW,
     *        PTRIST,PTLUST_S,PTRAST,STEP, PIMASTER,PAMASTER,
     *        PTRARW,PTRAIW,
     *        ITLOC, NSTEPS, SON_LEVEL2,
     *        COMP, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, POSFAC, LRLUS,
     *        ICNTL, KEEP, INTARR, DBLARR,
     *    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     *    COMM_NODES, MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        ELSE
         CALL CMUMPS_252(COMM_LOAD, ASS_IRECV, 
     *        N,INODE,IW,LIW,A,LA,
     *        IFLAG,IERROR,ND,NIRBDU,
     *        FILS,FRERE,MAXFRW,root,OPASSW, OPELIW,
     *        PTRIST,PTLUST_S,PTRAST,STEP, PIMASTER,PAMASTER,
     *        PTRARW,PTRAIW,
     *        ITLOC, NSTEPS, SON_LEVEL2,
     *        COMP, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, POSFAC, LRLUS,
     *        ICNTL, KEEP, INTARR, DBLARR,
     *    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     *    COMM_NODES, MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        ENDIF
       IF ( IFLAG .LT. 0 ) GOTO 640
        IF ((NBPROCFILS(STEP(INODE)).GT.0).OR.(SON_LEVEL2)) GOTO 20
      ELSE
        IF ( KEEP(55) .eq. 0 ) THEN
          CALL CMUMPS_253(COMM_LOAD, ASS_IRECV, 
     *    N, INODE, IW, LIW, A, LA,
     *    IFLAG, IERROR,
     *    ND, NIRBDU, FILS, FRERE, CAND, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    MAXFRW,
     *    root, OPASSW, OPELIW, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     *    PTRAIW, ITLOC, NSTEPS,
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     *    ICNTL, KEEP,INTARR,DBLARR,
     *    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     *    MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES,
     *    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     *    MEM_DISTRIB(0)
     *    )
        ELSE
          CALL CMUMPS_37( COMM_LOAD, ASS_IRECV, 
     *    NELT, FRTPTR, FRTELT,
     *    N, INODE, IW, LIW, A, LA, IFLAG, IERROR,
     *    ND, NIRBDU, FILS, FRERE, CAND, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    MAXFRW,
     *    root, OPASSW, OPELIW, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     *    PTRAIW, ITLOC, NSTEPS,
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     *    ICNTL, KEEP,INTARR,DBLARR,
     *    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     *    MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES,
     *    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     *    MEM_DISTRIB(0))
        END IF
        IF (IFLAG.LT.0) GOTO 640
        GOTO 20
      ENDIF
 100  CONTINUE
      FPERE = INODE
      DO WHILE (FPERE.GT.0)
        FPERE = FRERE(STEP(FPERE))
      END DO
      FPERE = -FPERE
      IF ( INODE .eq. KEEP(20) ) THEN
        POSELT = PTRAST(STEP(INODE))
        CALL CMUMPS_87( IW(PTLUST_S(STEP(INODE))),
     *                            POSELT )
        GOTO 200
      END IF
      POSELT = PTRAST(STEP(INODE))
      IOLDPS = PTLUST_S(STEP(INODE))
      IF (KEEP(50).EQ.0) THEN
         CALL CMUMPS_143( N, INODE, IW, LIW, A, LA,
     *               IOLDPS, POSELT,
     *               IFLAG, UU, NOFFW, NPVW,
     *               KEEP, STEP, PROCNODE_STEPS, MYID_NODES, SLAVEF,
     *               SEUIL)
      ELSE
         CALL CMUMPS_140( N, INODE, IW, LIW, A, LA,
     *               IOLDPS, POSELT,
     *               IFLAG, UU, NOFFW, NPVW,
     *               KEEP, MYID_NODES, SEUIL )
      ENDIF
      IF (IFLAG.LT.0) GOTO 635
 130  CONTINUE
      TYPE  = CMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF ( FPERE .NE. 0 ) THEN
        TYPEF = CMUMPS_330(STEP(FPERE),PROCNODE_STEPS,SLAVEF)
      ELSE
        TYPEF = -9999
      END IF
      CALL CMUMPS_254( COMM_LOAD, ASS_IRECV, 
     *       N,INODE,TYPE,TYPEF,LA,IW,LIW,A,
     *       IFLAG,IERROR,OPELIW,NELVAW,
     *       PTRIST,PTLUST_S,PTRAST, STEP, PIMASTER, PAMASTER,
     *       POSFAC,LRLU,
     *       LRLUS,IPTRLU,ICNTL,KEEP,NIRBDU,COMP,IWPOS,IWPOSCB,
     *       PROCNODE_STEPS,SLAVEF,FPERE,COMM_NODES,MYID_NODES,
     *       IPOOL, LPOOL, LEAF,
     *       NSTK_STEPS, NBPROCFILS, BUFR, LBUFR, LBUFR_BYTES, NBFIN,
     * root, OPASSW, ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     * ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IF (IFLAG.LT.0) GOTO 640
 200  CONTINUE
      IF ( INODE .eq. KEEP(38) ) THEN
        WRITE(*,*) 'Error .. in CMUMPS_251: ',
     &             ' INODE == KEEP(38)'
        Stop
      END IF
      IF ( FPERE.EQ.0 ) THEN
        NBROOT_TRAITEES = NBROOT_TRAITEES + 1
        IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
            NBFIN = NBFIN - NBROOT
            IF ( NBFIN .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in CMUMPS_251: ', 
     &                   ' NBFIN=', NBFIN 
              CALL CMUMPS_ABORT()
            END IF
            IF ( NBROOT .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in CMUMPS_251: ', 
     &                   ' NBROOT=', NBROOT
              CALL CMUMPS_ABORT()
            END IF
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL CMUMPS_242( DUMMY(1), 1, MPI_INTEGER,
     *          MYID_NODES, COMM_NODES, RACINE, SLAVEF)
            END IF
        ENDIF
        IF (NBFIN.EQ.0) GOTO 640
      ELSEIF ( FPERE.NE.KEEP(38) .AND.
     *         CMUMPS_275(STEP(FPERE),PROCNODE_STEPS,SLAVEF).EQ.
     *         MYID_NODES ) THEN
        NSTK_STEPS(STEP(FPERE)) = NSTK_STEPS(STEP(FPERE))-1
        IF ( NSTK_STEPS( STEP( FPERE )).EQ.0) THEN
          CALL CMUMPS_507(N, IPOOL, LPOOL,
     *         PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     *         STEP, FPERE )
          IF (KEEP(47) .GE. 3) THEN
             CALL CMUMPS_500(
     $            IPOOL, LPOOL, 
     *            PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *            MYID_NODES, STEP, N, ND, FILS )
          ENDIF
          CALL CMUMPS_137( FPERE, N, PROCNODE_STEPS,SLAVEF,
     &           ND, FILS, FRERE, STEP, PIMASTER, KEEP(28),
     *           KEEP(50), FLOP1,
     &           IW, LIW )
          CALL CMUMPS_190(.FALSE.,FLOP1,KEEP)
        ENDIF
      ENDIF
      GO TO 20
 635  CONTINUE
      CALL CMUMPS_44( MYID_NODES, SLAVEF, COMM_NODES )
 640  CONTINUE
      IF (ASS_IRECV.NE.MPI_REQUEST_NULL) THEN
        CALL CMUMPS_255( INFO(1),
     *       ASS_IRECV, STATUS, BUFR, LBUFR,
     *       LBUFR_BYTES,
     *       COMM_NODES,
     *       MYID_NODES, SLAVEF)
      ENDIF
       CALL CMUMPS_180( INFO(1),
     *      STATUS, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      COMM_NODES, COMM_LOAD, SLAVEF )
      CALL MPI_BARRIER( COMM_NODES, IERR )
       IF ( INFO(1) .GE. 0 ) THEN
         IF( KEEP(38) .NE. 0 .OR. KEEP(20).NE.0) THEN
           IF ( KEEP(38) .NE. 0 )THEN
             IF ( LRLU .GT. LBUFR ) THEN
#if defined (try_null_space_not_tested)
               IF ( KEEP(19) .ne. 0 ) THEN
                CALL CMUMPS_147( MYID_NODES, 
     *          root, N, 
     *          KEEP(38), COMM_NODES, IW, LIW, IWPOS + 1,
     *          A, LA, PTRAST, PTLUST_S, STEP,
     *          INFO(1), KEEP(50), KEEP(19), A(POSFAC),
     *          LRLU, KEEP,A(POSFAC),LRLU,
     *          A(POSFAC+3*root%tot_root_size+3),
     *          A(POSFAC+2*root%tot_root_size+2))
               ELSE
#endif
                CALL CMUMPS_146( MYID_NODES, 
     *          root, N, 
     *          KEEP(38), COMM_NODES, IW, LIW, IWPOS + 1,
     *          A, LA, PTRAST, PTLUST_S, STEP,
     *          INFO(1), KEEP(50), KEEP(19), A(POSFAC),
     *          LRLU, KEEP )
#if defined(try_null_space_not_tested)
               END IF
#endif
             ELSE
#if defined (try_null_space_not_tested)
               IF ( KEEP(19) .ne. 0 ) THEN
                 NFRONT = IW(PTLUST_S(STEP(KEEP(38))))
                 OPELIW = OPELIW + 
     *                  ( REAL(4)*REAL(NFRONT)*REAL(NFRONT)*
     *                  REAL(NFRONT)/REAL(3) ) /
     *                  REAL(SLAVEF)
                CALL CMUMPS_147( MYID_NODES,
     *          root, N, KEEP(38),
     *          COMM_NODES, IW, LIW, IWPOS + 1,
     *          A, LA, PTRAST, PTLUST_S, STEP,
     *          INFO(1), KEEP(50), KEEP(19),
     *          BUFR, MAX(LBUFR/2,1), KEEP,
     *          BUFR, MAX(LBUFR/2,1),
     *          BUFR(2*(3*root%tot_root_size+4)),
     *          BUFR(2*(2*root%tot_root_size+3)))
               ELSE
#endif
                LBUFRX = MIN(root%MBLOCK * root%NBLOCK,
     *                   root%TOT_ROOT_SIZE*root%TOT_ROOT_SIZE )
                ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                IF (IERR.gt.0) THEN
                  INFO(1) = -9
                  INFO(2) =  LBUFRX
                  write(LP,*) ' Error allocating, real array ',
     1              'of size before CMUMPS_146',  LBUFRX
                  stop
                ENDIF
                CALL CMUMPS_146( MYID_NODES,
     *          root, N, KEEP(38),
     *          COMM_NODES, IW, LIW, IWPOS + 1,
     *          A, LA, PTRAST, PTLUST_S, STEP,
     *          INFO(1), KEEP(50), KEEP(19),
     *          BUFRX, LBUFRX, KEEP )
                DEALLOCATE ( BUFRX )
#if defined(try_null_space_not_tested)
               END IF
#endif
             ENDIF
             IF ( MYID_NODES .eq. CMUMPS_275(STEP(KEEP(38)),
     *                            PROCNODE_STEPS,SLAVEF)) THEN
               IF ( INFO(1) .eq. -10 ) THEN
                 NPVW = NPVW + INFO(2)
               ELSE
                 NPVW = NPVW + root%tot_root_size
               END IF
             END IF
            IF (root%yes) THEN
              CALL CMUMPS_471(.FALSE.,
     *             LA-LRLUS
     *             ,IW( PTLUST_S(STEP(KEEP(38))) + 2 )*
     *             IW( PTLUST_S(STEP(KEEP(38))) + 1 ),0,KEEP)
            ENDIF
           ELSE
             IF (MYID_NODES.EQ. CMUMPS_275(STEP(KEEP(20)),
     *           PROCNODE_STEPS,SLAVEF)) THEN
                 NFRONT = IW(PTLUST_S(STEP(KEEP(20))))
                 CALL CMUMPS_471(.FALSE.,LA-LRLUS,
     *           NFRONT*NFRONT,0,KEEP)
             END IF
#if defined(try_null_space)
             IF ( KEEP(60) .EQ. 0 ) THEN
             IF (MYID_NODES.EQ. CMUMPS_275(STEP(KEEP(20)),
     *           PROCNODE_STEPS,SLAVEF)) THEN
               NFRONT = IW(PTLUST_S(STEP(KEEP(20))))
               OPELIW = OPELIW + 
     1                  REAL(4)*REAL(NFRONT)*REAL(NFRONT)*
     2                  REAL(NFRONT)/REAL(3)
               IF ( LRLU .GT. 4*NFRONT ) THEN
                CALL CMUMPS_295(
     1             NFRONT,A(PTRAST(STEP(KEEP(20)))),
     2             root,
     3             A(POSFAC+NFRONT), MAX(LRLU-NFRONT,1), 
     4             KEEP,A(POSFAC),INFO, LP )
               ELSE
                LBUFRX = 4*NFRONT
                ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                IF (IERR.gt.0) THEN
                  INFO(1) = -13
                  INFO(2) =  LBUFRX
                  write(LP,*) ' Error allocating, real array ',
     1              'of size ', LBUFRX
                  stop
                ENDIF
                CALL CMUMPS_295(
     1             NFRONT,A(PTRAST(STEP(KEEP(20)))),
     2             root,
     3             BUFRX(NFRONT+1), MAX(LBUFRX-NFRONT,1), 
     4             KEEP, BUFRX(1), INFO, LP )
                DEALLOCATE ( BUFRX )
               ENDIF
             ENDIF
             CALL MPI_BCAST( KEEP(17), 1, MPI_INTEGER,
     *       CMUMPS_275(STEP(KEEP(20)),PROCNODE_STEPS,SLAVEF),
     *       COMM_NODES, IERR )
             END IF
#endif
           END IF
         END IF
       END IF
       IF ( KEEP(38) .NE. 0 ) THEN
         IF (MYID_NODES.EQ. CMUMPS_275(STEP(KEEP(38)),
     *           PROCNODE_STEPS,SLAVEF)) THEN
           MAXFRW = MAX ( MAXFRW, root%TOT_ROOT_SIZE)
         END IF
       END IF
       MAXFRT       = MAXFRW
       NTOTPV       = NPVW
       INFO(12)     = NOFFW
       RINFO(2)     = OPASSW
       RINFO(3)     = OPELIW
       INFO(13)     = NELVAW
       INFO(14)     = COMP
      RETURN
      END SUBROUTINE CMUMPS_251
      SUBROUTINE CMUMPS_87( HEADER, VALPTR )
        INTEGER HEADER( 6 )
        INTEGER VALPTR
        INTEGER NFRONT, NASS
        NFRONT = HEADER(1)
        IF ( HEADER(2) .ne. 0 ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 1 :',HEADER(2)
          CALL CMUMPS_ABORT()
        END IF
        NASS   = IABS( HEADER( 3 ) )
        IF ( NASS .NE. IABS( HEADER( 4 ) ) ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 2 :',HEADER(3:4)
          CALL CMUMPS_ABORT()
        END IF
        IF ( NASS .NE. NFRONT ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 3 : not root'
          CALL CMUMPS_ABORT()
        END IF
        HEADER( 1 ) = NFRONT
        HEADER( 2 ) = 0
        HEADER( 3 ) = NFRONT
        HEADER( 4 ) = NFRONT
        HEADER( 5 ) = VALPTR
        RETURN
      END SUBROUTINE CMUMPS_87
      SUBROUTINE CMUMPS_136( id )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      include 'mpif.h'
      include 'cmumps_struc.h'
      TYPE( CMUMPS_STRUC ) :: id
      LOGICAL I_AM_SLAVE
      INTEGER IERR, MASTER, STATUS( MPI_STATUS_SIZE )
      PARAMETER ( MASTER = 0 )
      I_AM_SLAVE = ( id%MYID .ne. MASTER .OR. id%KEEP(46) .NE. 0 )
      IF (I_AM_SLAVE) THEN
        IF ( ASSOCIATED(id%S) ) THEN
           CALL CMUMPS_150
     *     (id%MYID_NODES,id%COMM_NODES,id%S(1),id%MAXS,id%KEEP(1))
           CALL CMUMPS_150
     *     (id%MYID_NODES,id%COMM_LOAD,id%S(1),id%MAXS,id%KEEP(1))
        ENDIF
        CALL MPI_BARRIER( id%COMM_NODES, IERR )
      END IF
      IF (id%root%gridinit_done) THEN
        IF ( id%KEEP(38).NE.0 ) THEN
          CALL BLACS_GRIDEXIT( id%root%CNTXT_BLACS )
          id%root%gridinit_done = .FALSE.
        END IF
      END IF
      IF ( id%MYID .NE. MASTER .OR. id%KEEP(46) .ne. 0 ) THEN
        CALL MPI_COMM_FREE( id%COMM_NODES, IERR )
        CALL MPI_COMM_FREE( id%COMM_LOAD, IERR )
      END IF
      IF (ASSOCIATED(id%MEM_DIST))  THEN
         DEALLOCATE(id%MEM_DIST)
         NULLIFY(id%MEM_DIST)
      ENDIF
      IF (id%MYID .NE. MASTER ) THEN
        IF (ASSOCIATED(id%RHS)) THEN
          WRITE(*,*) "Warning: Freeing RHS in CMUMPS_136"
          DEALLOCATE(id%RHS)
          NULLIFY(id%RHS)
        END IF
      END IF
      IF (ASSOCIATED(id%MAPPING)) THEN
        DEALLOCATE(id%MAPPING)
        NULLIFY(id%MAPPING)
      END IF
      IF (ASSOCIATED(id%NULL_SPACE)) THEN
        DEALLOCATE(id%NULL_SPACE)
        NULLIFY(id%NULL_SPACE)
      ENDIF
      IF ( id%KEEP(52) .NE. -1 .or. id%MYID .ne. MASTER ) THEN
        IF (ASSOCIATED(id%COLSCA)) THEN
          DEALLOCATE(id%COLSCA)
          NULLIFY(id%COLSCA)
        ENDIF
        IF (ASSOCIATED(id%ROWSCA)) THEN
          DEALLOCATE(id%ROWSCA)
          NULLIFY(id%ROWSCA)
        ENDIF
      END IF
      IF ( id%KEEP(74).EQ.1 .and. id%KEEP(52) .eq. -1 ) THEN 
        IF (ASSOCIATED(id%COLSCA)) THEN
          DEALLOCATE(id%COLSCA)
          NULLIFY(id%COLSCA)
        ENDIF
        IF (ASSOCIATED(id%ROWSCA)) THEN
          DEALLOCATE(id%ROWSCA)
          NULLIFY(id%ROWSCA)
        ENDIF
      ENDIF
      IF (ASSOCIATED(id%PTLUST_S)) THEN
        DEALLOCATE(id%PTLUST_S)
        NULLIFY(id%PTLUST_S)
      END IF
      IF (ASSOCIATED(id%POIDS)) THEN
        DEALLOCATE(id%POIDS)
        NULLIFY(id%POIDS)
      ENDIF
      IF (ASSOCIATED(id%IS)) THEN
        DEALLOCATE(id%IS)
        NULLIFY(id%IS)
      ENDIF
      IF (ASSOCIATED(id%IS1)) THEN
        DEALLOCATE(id%IS1)
        NULLIFY(id%IS1)
      ENDIF
      IF (ASSOCIATED(id%STEP))      THEN
        DEALLOCATE(id%STEP)
        NULLIFY(id%STEP)
      ENDIF
      IF (ASSOCIATED(id%NE_STEPS))  THEN
        DEALLOCATE(id%NE_STEPS)
        NULLIFY(id%NE_STEPS)
      ENDIF
      IF (ASSOCIATED(id%ND_STEPS))  THEN
        DEALLOCATE(id%ND_STEPS)
        NULLIFY(id%ND_STEPS)
      ENDIF
      IF (ASSOCIATED(id%FRERE_STEPS))  THEN
        DEALLOCATE(id%FRERE_STEPS)
        NULLIFY(id%FRERE_STEPS)
      ENDIF
      IF (ASSOCIATED(id%SYM_PERM))  THEN
        DEALLOCATE(id%SYM_PERM)
        NULLIFY(id%SYM_PERM)
      ENDIF
      IF (ASSOCIATED(id%UNS_PERM))  THEN
        DEALLOCATE(id%UNS_PERM)
        NULLIFY(id%UNS_PERM)
      ENDIF
      IF (ASSOCIATED(id%FILS))      THEN
        DEALLOCATE(id%FILS)
        NULLIFY(id%FILS)
      ENDIF
      IF (ASSOCIATED(id%PTRAR))     THEN
        DEALLOCATE(id%PTRAR)
        NULLIFY(id%PTRAR)
      ENDIF
      IF (ASSOCIATED(id%FRTPTR))    THEN
        DEALLOCATE(id%FRTPTR)
        NULLIFY(id%FRTPTR)
      ENDIF
      IF (ASSOCIATED(id%FRTELT))    THEN
        DEALLOCATE(id%FRTELT)
        NULLIFY(id%FRTELT)
      ENDIF
      IF (ASSOCIATED(id%NA))        THEN
        DEALLOCATE(id%NA)
        NULLIFY(id%NA)
      ENDIF
      IF (ASSOCIATED(id%PROCNODE_STEPS)) THEN
        DEALLOCATE(id%PROCNODE_STEPS)
        NULLIFY(id%PROCNODE_STEPS)
      ENDIF
      IF (ASSOCIATED(id%PROCNODE)) THEN
        DEALLOCATE(id%PROCNODE)
        NULLIFY(id%PROCNODE)
      ENDIF
      IF (id%KEEP(46).eq.1 .and.
     *    id%KEEP(55).ne.0 .and.
     *    id%MYID .eq. MASTER .and.
     *    id%KEEP(52) .eq. 0 ) THEN
        NULLIFY(id%DBLARR)
      ELSE
        IF (ASSOCIATED(id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
      END IF
      IF (ASSOCIATED(id%INTARR))       THEN
        DEALLOCATE(id%INTARR)
        NULLIFY(id%INTARR)
      ENDIF
      IF (ASSOCIATED(id%root%RG2L_ROW))THEN
        DEALLOCATE(id%root%RG2L_ROW)
        NULLIFY(id%root%RG2L_ROW)
      ENDIF
      IF (ASSOCIATED(id%root%RG2L_COL))THEN
        DEALLOCATE(id%root%RG2L_COL)
        NULLIFY(id%root%RG2L_COL)
      ENDIF
      IF (ASSOCIATED(id%root%IPIV))    THEN
        DEALLOCATE(id%root%IPIV)
        NULLIFY(id%root%IPIV)
      ENDIF
      IF (ASSOCIATED(id%root%QR_TAU))  THEN
        DEALLOCATE(id%root%QR_TAU)
        NULLIFY(id%root%QR_TAU)
      ENDIF
      IF (ASSOCIATED(id%root%GROW))    THEN
        DEALLOCATE(id%root%GROW)
        NULLIFY(id%root%GROW)
      ENDIF
      IF (ASSOCIATED(id%root%GSIN))    THEN
        DEALLOCATE(id%root%GSIN)
        NULLIFY(id%root%GSIN)
      ENDIF
      IF (ASSOCIATED(id%root%GCOS))    THEN
        DEALLOCATE(id%root%GCOS)
        NULLIFY(id%root%GCOS)
      ENDIF
      IF (ASSOCIATED(id%root%PERM_ROW))THEN
        DEALLOCATE(id%root%PERM_ROW) 
        NULLIFY(id%root%PERM_ROW) 
      ENDIF
      IF (ASSOCIATED(id%root%PERM_COL))THEN
        DEALLOCATE(id%root%PERM_COL) 
        NULLIFY(id%root%PERM_COL) 
      ENDIF
      IF (ASSOCIATED(id%root%ELROW))   THEN
        DEALLOCATE(id%root%ELROW)
        NULLIFY(id%root%ELROW)
      ENDIF
      IF (ASSOCIATED(id%root%EUROW))   THEN
        DEALLOCATE(id%root%EUROW)
        NULLIFY(id%root%EUROW)
      ENDIF
      IF (ASSOCIATED(id%root%PTREL))   THEN
        DEALLOCATE(id%root%PTREL)
        NULLIFY(id%root%PTREL)
      ENDIF
      IF (ASSOCIATED(id%root%PTREU))   THEN
        DEALLOCATE(id%root%PTREU)
        NULLIFY(id%root%PTREU)
      ENDIF
      IF (ASSOCIATED(id%root%ELELG))   THEN
        DEALLOCATE(id%root%ELELG)
        NULLIFY(id%root%ELELG)
      ENDIF
      IF (ASSOCIATED(id%root%EUELG))   THEN
        DEALLOCATE(id%root%EUELG)
        NULLIFY(id%root%EUELG)
      ENDIF
      IF (ASSOCIATED(id%root%DL))      THEN
        DEALLOCATE(id%root%DL)
        NULLIFY(id%root%DL)
      ENDIF
      IF (ASSOCIATED(id%ELTPROC))     THEN
        DEALLOCATE(id%ELTPROC)
        NULLIFY(id%ELTPROC)
      ENDIF
      IF (I_AM_SLAVE) THEN
       IF (ASSOCIATED(id%ISTEP_TO_INIV2)) THEN
        DEALLOCATE(id%ISTEP_TO_INIV2)
        NULLIFY(id%ISTEP_TO_INIV2)
       ENDIF
       IF (ASSOCIATED(id%TAB_POS_IN_PERE)) THEN
        DEALLOCATE(id%TAB_POS_IN_PERE)
        NULLIFY(id%TAB_POS_IN_PERE)
       ENDIF
       IF (ASSOCIATED(id%FUTURE_NIV2)) THEN
        DEALLOCATE(id%FUTURE_NIV2)
        NULLIFY(id%FUTURE_NIV2)
       ENDIF
      ENDIF
      IF (ASSOCIATED(id%CANDIDATES)) THEN
        DEALLOCATE(id%CANDIDATES)
        NULLIFY(id%CANDIDATES)
      ENDIF
      IF (ASSOCIATED(id%MEM_SUBTREE)) THEN
        DEALLOCATE(id%MEM_SUBTREE)
        NULLIFY(id%MEM_SUBTREE)
      ENDIF      
      IF (ASSOCIATED(id%S))        DEALLOCATE(id%S)
      NULLIFY(id%S)
      IF (I_AM_SLAVE) THEN
        CALL CMUMPS_57( IERR )
        CALL CMUMPS_59( IERR )
      END IF
      IF ( ASSOCIATED( id%BUFR ) ) DEALLOCATE( id%BUFR )
      NULLIFY( id%BUFR )
      RETURN
      END SUBROUTINE CMUMPS_136
      SUBROUTINE CMUMPS_150(MYID,COMM,S,MAXS,KEEP)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      INTEGER COMM, MYID, MAXS
      COMPLEX S( MAXS )
      INTEGER KEEP(150)
      INTEGER MSGTAG, MSGSOU, MSGLEN
      LOGICAL FLAG
      INTEGER MAXS_PACKED
      MAXS_PACKED = MAXS * KEEP( 35 )
      FLAG = .TRUE.
      DO WHILE ( FLAG )
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
        IF (FLAG) THEN
          MSGTAG=STATUS(MPI_TAG)
          MSGSOU=STATUS(MPI_SOURCE)
          CALL MPI_GET_COUNT(STATUS,MPI_PACKED,MSGLEN,IERR)
          CALL MPI_RECV(S(1),MAXS_PACKED,MPI_PACKED,
     *    MSGSOU, MSGTAG, COMM, STATUS, IERR)
        END IF
      END DO
      RETURN
      END SUBROUTINE CMUMPS_150
      SUBROUTINE CMUMPS_254(COMM_LOAD, ASS_IRECV, 
     *    N, INODE, TYPE, TYPEF, 
     *    LA, IW, LIW, A,
     *    IFLAG, IERROR, OPELIW, NELVAW, PTRIST, PTLUST_S, PTRAST,
     *    STEP, PIMASTER, PAMASTER,
     *    POSFAC, LRLU, LRLUS, IPTRLU, ICNTL, KEEP,
     *    NIRBDU, COMP, IWPOS, IWPOSCB, PROCNODE_STEPS, SLAVEF,
     *    FPERE, COMM, MYID,
     *    IPOOL, LPOOL, LEAF, NSTK_S,
     *    NBPROCFILS,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, root,
     *    OPASSW, ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *    ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER COMM, MYID, TYPE, TYPEF
      INTEGER N, LA, LIW, INODE,IFLAG,IERROR
      INTEGER ICNTL(40), KEEP(150)
      INTEGER POSFAC, LRLU,LRLUS, IPTRLU,NIRBDU,IWPOSCB,IWPOS,
     &        FPERE, SLAVEF, NELVAW
      INTEGER IW(LIW),PROCNODE_STEPS(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      COMPLEX    A(LA)
      REAL OPASSW, OPELIW
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      INTEGER INTARR(MAX(1,KEEP(14))) 
      INTEGER ITLOC( N ), FILS( N ),
     *        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER LPOOL, LEAF, COMP
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK_S( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NBFIN
      INTEGER IWRITEFAC, IREADFAC
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER APOS, POSELT, OPSFAC, NSLAVES, NCB,
     *        H_INODE, IERR, IERR_MPI, NBCOL, NBROW, NBROW_SEND,
     *        NBROW_STACK, NELIM
      INTEGER SHIFT_LIST_ROW_SON,
     *        SHIFT_LIST_COL_SON, SHIFT_VAL_SON,
     *        LIST_ROW_SON, LIST_COL_SON, LIST_SLAVES
      INTEGER IOLDPS,NFRONT,NPIV,NASS,IOLDP1,PTROWEND,
     *        LREQI, LREQCB_SAMEPROC, LREQCB_OTHERPROC, LCONT
      INTEGER I,NPOS,J,J1,LDA, INIV2
      INTEGER MSGDEST, MSGTAG
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      REAL FLOP1
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      IOLDPS = PTLUST_S(STEP(INODE))
      NFRONT = IW(IOLDPS)
      NPIV   = IW(IOLDPS + 1)
      NASS   = IABS(IW(IOLDPS + 2))
      NSLAVES=  IW(IOLDPS+5)
      H_INODE= 6 + NSLAVES
      LCONT = NFRONT - NPIV
      NBCOL = LCONT
      IF ((FPERE.EQ.0).AND.(NASS.NE.NPIV)) THEN 
        IFLAG = -10
        GOTO 600
      ENDIF
      NBROW      = LCONT
      IF (TYPE.EQ.2) NBROW = NASS - NPIV
      IF ((KEEP(50).NE.0).AND.(TYPE.EQ.2)) THEN
        LDA = NASS
      ELSE
        LDA = NFRONT
      ENDIF
      NBROW_SEND = NBROW
      IF (TYPEF.EQ.2) NBROW_SEND = NASS - NPIV
      LREQCB_SAMEPROC  = NBROW * NBCOL
      LREQCB_OTHERPROC = (NBROW-NBROW_SEND) * NBCOL
      NELIM = NASS-NPIV
      POSELT = PTRAST(STEP(INODE))
      NELVAW = NELVAW + NASS - NPIV
      IF ( KEEP(50) .eq. 0 ) THEN
       IF ( TYPE .eq. 1 ) THEN
        FLOP1 = REAL(2) * REAL(NFRONT) * REAL(NPIV) *
     *      REAL(NFRONT - NPIV - 1) + 
     *      REAL(NPIV * (NPIV + 1)) * REAL(2 * NPIV + 1) / REAL(3)
        FLOP1 = FLOP1 + REAL(2 * NFRONT - NPIV - 1) * REAL(NPIV) /
     *      REAL(2)
       ELSE ! Node is of level 2
        FLOP1 = REAL(2*NASS)*REAL(NFRONT) -
     *          REAL(NASS+NFRONT)*REAL(NPIV+1)
        FLOP1 = REAL(NPIV)*FLOP1 + 
     *     REAL(((2 * NASS - NPIV - 1) * NPIV) / 2) + 
     *     REAL(NPIV) * REAL(NPIV + 1) * REAL(2 * NPIV + 1) /REAL(3)
       END IF
      ELSE
       IF ( TYPE .eq. 1 ) THEN
        FLOP1 = REAL(NPIV) * (
     *          REAL( NFRONT ) * REAL( NFRONT ) + REAL( NFRONT )
     *        - REAL( NFRONT*NPIV + NPIV + 1 )    )  
     *        + ( REAL(NPIV)*REAL(NPIV+1)*REAL(2*NPIV+1)  )
     *        / REAL(6)
       ELSE ! Node is of level 2
        FLOP1 = REAL(NPIV) * (
     *          REAL( NASS ) * REAL( NASS ) + REAL( NASS )
     *        - REAL( NASS*NPIV +  NPIV + 1 )    )
     *        + ( REAL(NPIV)*REAL(NPIV+1)*REAL(2*NPIV+1) )
     *        / REAL( 6 )
       END IF
      END IF
      OPELIW = OPELIW + FLOP1
      IF ( NPIV .NE. NASS ) THEN
        IF ( TYPE .eq. 1 ) THEN
          FLOP1 = REAL(NASS) * (
     *          REAL( NFRONT ) * REAL( NFRONT ) + REAL( NFRONT )
     *        - REAL( NFRONT*NASS + NASS + 1 )
     *        + REAL( NASS + 1 ) * REAL( 2 * NASS + 1 )
     *        / REAL( 6 )
     *          )
        ELSE
          FLOP1 = REAL(NASS) * (
     *          REAL( NASS ) * REAL( NASS ) + REAL( NASS )
     *        - REAL( NASS*NASS + NASS + 1 )
     *        + REAL( NASS + 1 ) * REAL( 2 * NASS + 1 )
     *        / REAL( 6 )
     *          )
        END IF
      END IF
      FLOP1=-FLOP1
      CALL CMUMPS_190(.FALSE.,FLOP1,KEEP)
      IW(IOLDPS + 4) = POSELT
      IF ( FPERE .EQ. 0 ) GOTO 200
      IF ( FPERE.EQ.KEEP(38) ) THEN
       NCB   = NFRONT - NASS
       SHIFT_LIST_ROW_SON = H_INODE + NASS
       SHIFT_LIST_COL_SON = H_INODE + NFRONT + NASS
       SHIFT_VAL_SON      = NASS*NFRONT + NASS
       IF (TYPE.EQ.1) THEN
         CALL CMUMPS_80(
     *    COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE,
     *    PTLUST_S, PTRAST,
     *    ROOT, NCB, NCB, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, NFRONT,
     *    ROOT_CONT_STATIC, MYID, COMM,
     *
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, .FALSE., ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF (IFLAG < 0 ) GOTO 500
       ENDIF
       MSGDEST=  CMUMPS_275(STEP(FPERE),PROCNODE_STEPS,SLAVEF)
       IOLDPS = PTLUST_S(STEP(INODE))
       LIST_ROW_SON = IOLDPS + H_INODE + NPIV
       LIST_COL_SON = IOLDPS + H_INODE + NFRONT + NPIV
       LIST_SLAVES  = IOLDPS + 6
       IF (MSGDEST.EQ.MYID) THEN
         CALL CMUMPS_273( root, 
     *      INODE, NELIM, NSLAVES, IW(LIST_ROW_SON), 
     *      IW(LIST_COL_SON), IW(LIST_SLAVES),
     *
     *      PROCNODE_STEPS, IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, ITLOC, COMP,
     *      IFLAG, IERROR, 
     *      IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,
     *      COMM, COMM_LOAD, FILS, ND)
         IF (IFLAG.LT.0) GOTO 600
       ELSE
        IERR = -1
        DO WHILE (IERR.EQ.-1)
         CALL CMUMPS_76( INODE, NELIM, 
     *     IW(LIST_ROW_SON), IW(LIST_COL_SON), NSLAVES, 
     *     IW(LIST_SLAVES), MSGDEST, COMM, IERR)
         IF ( IERR .EQ. -1 ) THEN
            BLOCKING   =.FALSE.
            SET_IRECV  =.TRUE.
            MESSAGE_RECEIVED = .FALSE.
            CALL CMUMPS_329( COMM_LOAD, ASS_IRECV, 
     *       BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *       MPI_ANY_SOURCE, MPI_ANY_TAG,
     *       STATUS,
     *       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *       IWPOS, IWPOSCB, IPTRLU,
     *       LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *       IFLAG, IERROR, COMM,
     *       NBPROCFILS,
     *       IPOOL, LPOOL, LEAF,
     *       NBFIN, MYID, SLAVEF,
     *
     *       root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *       INTARR, DBLARR, ICNTL, KEEP, ND, FRERE, LPTRAR, NELT,
     *       FRTPTR, FRTELT, 
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
            IF ( IFLAG .LT. 0 ) GOTO 500
            IOLDPS = PTLUST_S(STEP(INODE))
            LIST_ROW_SON = IOLDPS + H_INODE + NPIV
            LIST_COL_SON = IOLDPS + H_INODE + NFRONT + NPIV
            LIST_SLAVES  = IOLDPS + 6
         ENDIF
        ENDDO
        IF ( IERR .EQ. -2 ) THEN 
            IERROR = ( 3 + NSLAVES + 2 * NELIM ) * KEEP( 34 ) 
            GOTO 220
        ENDIF
       ENDIF
       IF (NELIM.EQ.0) THEN 
          POSELT = PTRAST(STEP(INODE))
          OPSFAC = POSELT + NPIV * NFRONT + NPIV
          GOTO 190
       ELSE
          GOTO 500
       ENDIF
      ENDIF
      OPSFAC = POSELT + NPIV * LDA + NPIV
      IF ( CMUMPS_275(STEP(FPERE),PROCNODE_STEPS,
     *     SLAVEF) .EQ. MYID ) THEN
        LREQI = 2
        CALL CMUMPS_22(.FALSE.,
     *   MYID,N,KEEP,IW, LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER,PAMASTER, ITLOC,
     *   LREQI, LREQCB_SAMEPROC,
     *   COMP, LRLUS, IFLAG, IERROR )
        IF (IFLAG.LT.0) GOTO 600
        NPOS = IPTRLU + 1
        PAMASTER(STEP(INODE)) = NPOS
        PTRAST(STEP(INODE)) = -99999999
        IF ((KEEP(50).NE.0).AND.(TYPE.EQ.2)) THEN
          IF (NBROW.GT.0) THEN
            APOS = OPSFAC
            DO I = 1, NBROW
             J1 = APOS
             DO J = 1, NBROW
               A(NPOS) = A(J1)
               NPOS = NPOS + 1
               J1 = J1 + 1
             ENDDO
             APOS = APOS + LDA
            ENDDO
          ENDIF
          PIMASTER (STEP(INODE)) = PTLUST_S(STEP(INODE))
          PTRIST(STEP(INODE)) = IWPOSCB+1
          IW(IWPOSCB+1) = MIN(-NBROW,-1)
          IW(IWPOSCB+2) = NBROW
        ELSE
          IF ((KEEP(50).NE.0).AND.(NPIV.GT.0).AND.(NFRONT.NE.NPIV)) THEN
            IWRITEFAC = POSELT + NPIV
            IREADFAC  = POSELT + NFRONT
            IF (NPIV < KEEP(59) ) THEN
             DO I = 1, NPIV - 1
              DO J = 0, NPIV-1
                A(IWRITEFAC+J)=A(IREADFAC+J)
              ENDDO
              IREADFAC = IREADFAC + NFRONT
              IWRITEFAC = IWRITEFAC + NPIV
             ENDDO
            ELSE
             DO I = 1, NPIV - 1
              DO J = 0, I
                A(IWRITEFAC+J)=A(IREADFAC+J)
              ENDDO
              IREADFAC = IREADFAC + NFRONT
              IWRITEFAC = IWRITEFAC + NPIV
             ENDDO
            ENDIF
          ENDIF
          IF (NBROW.GT.0) THEN
            APOS  = OPSFAC
            DO I = 1, NBROW
              IF ( KEEP(50) .eq. 0 ) THEN
                CALL CCOPY( NBCOL, A(APOS), 1, A(NPOS), 1 )
                IF  ( I == 1 ) THEN
                   IWRITEFAC  = OPSFAC
                   IREADFAC   = POSELT + (NPIV + 1) * NFRONT
                ELSE IF (NPIV.GT.0) THEN
                   DO J=0,NPIV-1
                      A(IWRITEFAC + J) = A(IREADFAC + J)
                   ENDDO
                   IWRITEFAC=IWRITEFAC+NPIV
                   IREADFAC=IREADFAC+NPIV+NBCOL
                END IF
              ELSE
                CALL CCOPY( I    , A(APOS), 1, A(NPOS), 1 )
                IF (NPIV.GT.0) THEN
                  DO J = 0, NPIV - 1
                    A(IWRITEFAC+J) = A(IREADFAC+J)
                  ENDDO
                  IWRITEFAC = IWRITEFAC + NPIV
                  IREADFAC = IREADFAC + NFRONT
                ENDIF
#if ! defined(ALLOW_NON_INIT)
                A(NPOS + I: NPOS + NBCOL - 1) = CMPLX(ZERO)
#endif
              END IF
              NPOS = NPOS + NBCOL
              APOS = APOS + LDA
            ENDDO
          ENDIF
          PIMASTER (STEP(INODE)) = PTLUST_S(STEP(INODE))
          PTRIST(STEP(INODE)) = IWPOSCB+1
          IW(IWPOSCB+1) = MIN( -NBCOL, -1 )
          IW(IWPOSCB+2) = NBROW
          IF ( NBCOL .eq. 0 .and. NBROW .ne. 0 ) THEN
            write(*,*) 'error: more rows than columns',NBCOL,NBROW
            stop
          END IF
          GOTO 200
        ENDIF
      ELSE
        MSGTAG =NOEUD
        MSGDEST=CMUMPS_275( STEP(FPERE), PROCNODE_STEPS, SLAVEF )
        IERR = -1 
        DO WHILE (IERR.EQ.-1)
          IF ( (TYPE.EQ.1) .AND. (TYPEF.EQ.1) ) THEN
             CALL CMUMPS_66( INODE, FPERE, NFRONT, 
     *         LCONT, NASS, NPIV, IW( IOLDPS +  H_INODE + NPIV ),
     *         IW( IOLDPS +  H_INODE + NPIV + NFRONT ),
     *         A( OPSFAC ),
     *         MSGDEST, MSGTAG, COMM, IERR )
          ELSE
             IF ( (TYPE.EQ.2).AND.(KEEP(48).NE.0)) THEN
              INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
             ELSE
              INIV2 = -9999
             ENDIF
             CALL CMUMPS_70( FPERE, INODE, 
     *             NBROW_SEND, IW(IOLDPS +  H_INODE + NPIV ), 
     *             NBCOL, IW(IOLDPS +  H_INODE + NPIV + NFRONT ),
     *             A(OPSFAC), LDA, NELIM, TYPE,
     *             NSLAVES, IW(IOLDPS+6), MSGDEST, 
     *             COMM, IERR, 
     * 
     *             SLAVEF, KEEP, INIV2, TAB_POS_IN_PERE )
          END IF
          IF ( IERR .EQ. -1 ) THEN
            BLOCKING  = .FALSE.
            SET_IRECV = .TRUE.
            MESSAGE_RECEIVED = .FALSE.
            CALL CMUMPS_329( COMM_LOAD, ASS_IRECV, 
     *       BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *       MPI_ANY_SOURCE, MPI_ANY_TAG,
     *       STATUS,
     *       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *       IWPOS, IWPOSCB, IPTRLU,
     *       LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *       PTRIST, PTLUST_S,
     *       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *       IFLAG, IERROR, COMM,
     *       NBPROCFILS,
     *       IPOOL, LPOOL, LEAF,
     *       NBFIN, MYID, SLAVEF,
     *
     *       root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *       INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *       LPTRAR, NELT, FRTPTR, FRTELT,
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
             IF ( IFLAG .LT. 0 ) GOTO 500
          ENDIF
          IOLDPS = PTLUST_S(STEP( INODE ))
        END DO
        IF ( IERR .EQ. -2 ) THEN 
          IF ( (TYPE.EQ.1) .AND. (TYPEF.EQ.1) ) THEN
            IERROR = ( 2*LCONT + 9 ) * KEEP( 34 ) + 
     *          LCONT*LCONT * KEEP( 35 )
          ELSE IF (KEEP(50).ne.0 .AND. TYPE .eq. 2 ) THEN
            IERROR = ( NBROW_SEND + NBCOL+ 5 + NSLAVES) 
     *                 * KEEP( 34 ) + 
     *          NBROW_SEND*NBROW_SEND*KEEP( 35 )
          ELSE
            IERROR = ( NBROW_SEND + NBCOL+ 5 + NSLAVES) * KEEP( 34 ) + 
     *          NBROW_SEND*NBCOL*KEEP( 35 )
          ENDIF
          GOTO 220
        ENDIF
        IF ((FPERE.EQ.0).OR.(LREQCB_OTHERPROC.EQ.0)) GOTO 190
        NBROW_STACK = NBROW-NBROW_SEND
        LREQI       = 6 + NBROW_STACK + NBCOL
        CALL CMUMPS_22(.FALSE.,
     *   MYID,N,KEEP,IW,LIW,A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER,PAMASTER, ITLOC,
     *   LREQI, LREQCB_OTHERPROC,
     *   COMP, LRLUS, IFLAG, IERROR )
        IF (IFLAG.LT.0) GOTO 600
        PTRIST(STEP(INODE)) = IWPOSCB+1
        IW(IWPOSCB+1) = NBCOL
        IW(IWPOSCB+2) = 0
        IW(IWPOSCB+3) = NBROW_STACK
        IW(IWPOSCB+4) = 0
        IW(IWPOSCB+5) = 1
        IW(IWPOSCB+6) = 0
        IOLDP1   = PTLUST_S(STEP(INODE))+H_INODE
        PTROWEND = IWPOSCB+6+NBROW_STACK
        IW(IWPOSCB+7:PTROWEND) =
     *    IW(IOLDP1+NFRONT-NBROW_STACK:IOLDP1+NFRONT-1)
        IW(PTROWEND+1:PTROWEND+NBCOL) =
     *    IW(IOLDP1+NFRONT+NPIV: IOLDP1+2*NFRONT-1)
        PTRAST(STEP(INODE)) = IPTRLU+1
        NPOS = IPTRLU+1
        APOS = OPSFAC + NBROW_SEND*NFRONT
        IF ( NELIM .ne. NBCOL-NBROW_STACK ) THEN
          write(*,*) 'Internal error in CMUMPS_254',
     *  NBCOL, NBROW_STACK, NELIM
          CALL CMUMPS_ABORT()
        ENDIF
        DO I = 1, NBROW_STACK
          IF ( KEEP(50) .eq. 0 ) THEN
            CALL CCOPY( NBCOL, A(APOS), 1, A(NPOS), 1 )
          ELSE
            CALL CCOPY( I+NELIM, A(APOS), 1, A(NPOS), 1 )
#if ! defined(ALLOW_NON_INIT)
            A(NPOS+I+NELIM:NPOS+NBCOL-1)= CMPLX(ZERO)
#endif
          END IF
          NPOS = NPOS + NBCOL
          APOS = APOS + NFRONT
        ENDDO
      END IF
  190 CONTINUE
      IF (KEEP(50).EQ.0) THEN
        IF ((NPIV .NE. 0) .AND. (NBROW .GT. 1)) THEN
          IF ( IOLDPS .NE. PTLUST_S( STEP(INODE)) ) THEN
            WRITE(*,*)
     *      ' ERROR 1 in 51G : IOLDPS,PTLUST_S(STEP(INODE))=',
     *      IOLDPS,PTLUST_S(STEP(INODE))
            CALL CMUMPS_ABORT()
          END IF
          IOLDPS = PTLUST_S(STEP(INODE))
          POSELT = IW(IOLDPS + 4)
          APOS     = POSELT + (NPIV + 1) * NFRONT
          DO 180 I = 1, NBROW - 1
            DO J = 0, NPIV-1
              A(OPSFAC + J) = A(APOS + J)
            END DO
            OPSFAC = OPSFAC + NPIV
            APOS   = APOS + NPIV + NBCOL
  180     CONTINUE
        ENDIF
      ELSE
        IOLDPS = PTLUST_S(STEP(INODE))
        POSELT = IW(IOLDPS + 4)
        IF ( TYPE .eq. 1 ) THEN
          CALL CMUMPS_324(A(POSELT), NFRONT, NPIV, NBROW,
     *                               KEEP(59))
        ELSE
          CALL CMUMPS_325(A(POSELT), NPIV, NBROW)
        ENDIF
      END IF
  200 CONTINUE
      IOLDPS = PTLUST_S(STEP(INODE))
      IW(IOLDPS)     = NBCOL
      IW(IOLDPS + 1) = NASS - NPIV
      IF (TYPE.EQ.2) THEN
        IW(IOLDPS + 2) = NASS
      ELSE
        IW(IOLDPS + 2) = NFRONT
      ENDIF
      IW(IOLDPS + 3) = NPIV
      CALL CMUMPS_93(MYID,N,IOLDPS,TYPE, IW, LIW, 
     *    A, LA, POSFAC, LRLU, LRLUS, 
     *    IWPOS, PTRAST, STEP, KEEP)
 500  CONTINUE
      RETURN
  220 CONTINUE
       IFLAG  = -17
  600 CONTINUE
      IF (IFLAG .NE. -1) CALL CMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_254
      SUBROUTINE CMUMPS_142( id)
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
#ifndef SUN_
      INTERFACE
      SUBROUTINE CMUMPS_84( id )
      INCLUDE 'cmumps_struc.h'
      TYPE (CMUMPS_STRUC), TARGET :: id
      END SUBROUTINE CMUMPS_84
      SUBROUTINE CMUMPS_27(id, ANORMINF, LSCAL)
      INCLUDE 'cmumps_struc.h'
      TYPE (CMUMPS_STRUC), TARGET :: id
      REAL, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      END SUBROUTINE CMUMPS_27
      END INTERFACE
#endif
      INCLUDE 'cmumps_struc.h'
      TYPE(CMUMPS_STRUC), TARGET :: id
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER NSEND, NSEND_TOT, LDPTRAR, NELT
      INTEGER NLOCAL, NLOCAL_TOT, KEEP13_SAVE
      INTEGER  CMUMPS_275
      EXTERNAL CMUMPS_275
      INTEGER MP, LP, MPG, allocok
      LOGICAL PROK, PROKG, LSCAL
      INTEGER CMUMPS_LBUF, CMUMPS_LBUFR_BYTES, CMUMPS_LBUF_INT
      INTEGER PTRIST, PTRWB, MAXELT_SIZE,
     *        ITLOC, IPOOL, NSTEPS, LPOOL, LIW
      INTEGER IRANK, ID_SCHUR, SIZE_SCHUR
      INTEGER KKKK
      DOUBLE PRECISION TIME
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER PERLU, MEMORY_SIZE
      INTEGER ITOT
      LOGICAL I_AM_SLAVE
      REAL :: ANORMINF, SEUIL, SEUIL_LDLT_NIV2
      INTEGER, DIMENSION(:), POINTER :: IWK
      COMPLEX, DIMENSION(:), POINTER :: WK
      INTEGER, POINTER :: MAXS, JOB, NZ, MAXIS
      REAL,DIMENSION(:),POINTER::RINFO, RINFOG
      REAL,DIMENSION(:),POINTER::    CNTL
      INTEGER,DIMENSION(:),POINTER::INFO, INFOG, KEEP
      INTEGER,DIMENSION(:),POINTER::ICNTL
      EXTERNAL CMUMPS_505
      INTEGER CMUMPS_505
#if defined(LARGEMATRICES)
      INTEGER NUMROC
      EXTERNAL NUMROC
#endif
      JOB=>id%JOB
      MAXS=>id%MAXS
      NZ=>id%NZ
      MAXIS=>id%MAXIS
      RINFO=>id%RINFO
      RINFOG=>id%RINFOG
      CNTL=>id%CNTL
      INFO=>id%INFO
      INFOG=>id%INFOG
      KEEP=>id%KEEP
      ICNTL=>id%ICNTL
      KEEP13_SAVE = id%KEEP(13)
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LP  = ICNTL( 1 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IF ( PROK ) WRITE( MP, 130 )
      IF ( PROKG ) WRITE( MPG, 130 )
      IF ( PROKG .and. KEEP(53).GT.0 ) THEN
        WRITE(MPG,'(/A,I3)') ' Null space option :', KEEP(19)
        IF ( KEEP(21) .ne. 0 ) THEN 
          WRITE( MPG, '(A,I10)') ' Max deficiency    : ', KEEP(21)
        END IF
      END IF
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               KEEP(46) .eq. 1 ) )
        IF ( KEEP(50) .eq. 1 ) THEN
          IF (id%CNTL(1) .ne. 0.0E0 ) THEN
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)')
     *' ** Warning : SPD solver called, resetting CNTL(1) to 0.0E0'
            END IF
          END IF
          id%CNTL(1) = 0.0E0
        END IF
      IF ( KEEP(60) .ne. 0 ) THEN
        id%CNTL(1) = 0.0E0
      END IF
      CALL CMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF (INFO(1).LT.0) GOTO 550
      IF ( PROKG ) THEN
          WRITE( MPG, 172 ) id%NSLAVES,
     *    KEEP(125), KEEP(126), KEEP(127), KEEP(28)
      ENDIF
      IF ( id%MYID.EQ.MASTER ) THEN
        LSCAL = ((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 6))
        IF ( LSCAL ) THEN
          IF ( KEEP(52) .eq. 5 .or. 
     &         KEEP(52) .eq. 6 ) THEN
            ITOT = 5 * id%N + NZ
          ELSE
            ITOT = 5 * id%N
          END IF
          ALLOCATE( WK( ITOT ), stat = IERR )
          IF ( IERR .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=ITOT
            GOTO 137
          END IF
          CALL CMUMPS_217(id%N, NZ, KEEP(52), id%A(1),
     *             id%IRN(1), id%JCN(1),
     *             id%COLSCA(1), id%ROWSCA(1),
     *             WK, ITOT, ICNTL, INFO )
          DEALLOCATE( WK )
        END IF
        LSCAL = (LSCAL .OR. (KEEP(52) .EQ. -1))
      END IF
 137  CONTINUE
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 550
      SEUIL = id%CNTL(3)
      CALL MPI_BCAST(SEUIL, 1, MPI_REAL,
     *               MASTER, id%COMM, IERR)
      IF (SEUIL .LT. 0) THEN
        IF (KEEP(50) .EQ. 2 .OR. KEEP(19) .NE. 0) THEN
           CALL CMUMPS_27(  id , ANORMINF, LSCAL )
           SEUIL = EPSILON( 0.0E0 ) * ANORMINF
        ELSE
           SEUIL = ZERO
        ENDIF
        SEUIL_LDLT_NIV2 = SEUIL
        IF (KEEP(19) .EQ. 0) THEN
          SEUIL = ZERO
        ENDIF
      ELSE
        SEUIL_LDLT_NIV2 = SEUIL
      ENDIF
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 550
      IF ( I_AM_SLAVE ) THEN
        CALL CMUMPS_188( id%COST_SUBTREES, KEEP(64), 
     *                                KEEP(66) )
        CALL CMUMPS_185(   id%COMM_LOAD, id%SYM, KEEP,
     *                          id%MEM_DIST, id%NSLAVES,
     *                          INFO 
     *                          ,id%MEM_SUBTREE,id%NBSA_LOCAL
#if ! defined(OLD_LOAD_MECHANISM)
     *     , id%FUTURE_NIV2
#endif
     *  )
#if ! defined(OLD_LOAD_MECHANISM)
        CALL CMUMPS_190(.FALSE.,id%COST_SUBTREES,KEEP)
#endif
        IF (INFO(1).LT.0) GOTO 111
#if defined(stephinfo)
        write(*,*) 'proc ',id%MYID,' array of dist : ',
     *       id%MEM_DIST(0:id%NSLAVES - 1)
#endif
      END IF
      PERLU = KEEP(12)
      IF ( KEEP(16) .gt. 0 ) THEN
        IF ( KEEP(50).EQ.1) THEN
          id%MAXS = MAX( id%MAXS,
     *                KEEP(16) + PERLU * ( KEEP(16) / 100 + 1 )
     *                + KEEP(44) + PERLU * ( KEEP(44) / 100 ) )
        ELSE
          id%MAXS = MAX( id%MAXS,
     *                KEEP(16) + PERLU * ( KEEP(16) / 100 + 1 )
     *                + KEEP(44) + 4*PERLU * ( KEEP(44) / 100 ) )
        ENDIF
      ELSE
        id%MAXS = 1
      END IF
      IF ( ASSOCIATED (id%S) ) DEALLOCATE(id%S)
#if defined (LARGEMATRICES)
      IF ( id%MYID .ne. MASTER ) THEN
#endif
      ALLOCATE (id%S(id%MAXS),stat=ierr)
      IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = id%MAXS
      ENDIF
#if defined (LARGEMATRICES)
      END IF
#endif
 111  CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP(55) .eq. 0 ) THEN
        IF (ASSOCIATED( id%DBLARR)) DEALLOCATE(id%DBLARR)
        IF ( I_AM_SLAVE .and. KEEP(13) .ne. 0 ) THEN
          ALLOCATE( id%DBLARR( KEEP(13) ), stat = IERR )
        ELSE
          ALLOCATE( id%DBLARR( 1 ), stat =IERR )
        END IF
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     *       ':Error allocating DBLARR : IERR = ', IERR
          INFO(1)=-13
          INFO(2)=KEEP(13)
          GOTO 100
        END IF
      ELSE
         IF ( ASSOCIATED( id%INTARR ) ) THEN
           DEALLOCATE( id%INTARR )
           NULLIFY( id%INTARR )
         END IF
         IF ( I_AM_SLAVE .and. KEEP(14) .ne. 0 ) THEN
           ALLOCATE( id%INTARR( KEEP(14) ), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = KEEP(14)
             GOTO 100
           END IF
         ELSE
           ALLOCATE( id%INTARR(1),stat=allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             GOTO 100
           END IF
         END IF
         IF (ASSOCIATED( id%DBLARR)) DEALLOCATE(id%DBLARR)
         IF ( I_AM_SLAVE ) THEN
           IF (      id%MYID_NODES .eq. MASTER
     *       .AND.   KEEP(46)   .eq. 1
     *       .AND.   KEEP(52)   .eq. 0 ) THEN
             id%DBLARR => id%A_ELT
           ELSE
             IF ( KEEP(13) .ne. 0 ) THEN
               ALLOCATE( id%DBLARR( KEEP(13) ), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = KEEP(13)
                 GOTO 100
               END IF
             ELSE
               ALLOCATE( id%DBLARR(1), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = 1
                 GOTO 100
               END IF 
             END IF
           END IF
         ELSE
           ALLOCATE( id%DBLARR(1), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             GOTO 100
           END IF
         END IF
      END IF
      IF ( KEEP(38).NE.0 .AND.  I_AM_SLAVE ) THEN
         CALL CMUMPS_165( id%N,
     *   id%root, id%FILS(1), KEEP(38), id%INFO )
      END IF
 100  CONTINUE
      CALL CMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP( 55 ) .eq. 0 ) THEN
      IF ( KEEP(54) .eq. 0 ) THEN
      IF ( id%MYID .eq. MASTER ) THEN
        ALLOCATE(IWK(id%N), stat=allocok)
        IF ( allocok .NE. 0 ) THEN
          INFO(1)=-13
          INFO(2)=id%N
        END IF
      ENDIF
      CALL CMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( id%MYID .eq. MASTER ) THEN
        IF (PROKG ) THEN
          CALL CMUMPS_291(TIME)
        END IF
        IF ( .not. associated( id%INTARR ) )
     *  ALLOCATE( id%INTARR( 1 ) )
#if defined(LARGEMATRICES)
         IF (KEEP(38).NE.0) THEN
          ITOT = NUMROC( id%root%ROOT_SIZE, id%root%MBLOCK,
     *               id%root%MYROW, 0, id%root%NPROW ) 
          ITOT = MAX( 1, ITOT )
          ITOT = ITOT*
     *          NUMROC( id%root%ROOT_SIZE, id%root%NBLOCK,
     *          id%root%MYCOL, 0, id%root%NPCOL )
          ITOT = MAX( 1, ITOT )
         ELSE
          ITOT = 1
         ENDIF
         IF ( ASSOCIATED (id%S) ) DEALLOCATE(id%S)
         ALLOCATE (WK(ITOT),stat=ierr)
         IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = ITOT
          write(6,*) ' PB1 LARGEMAT'
          CALL CMUMPS_ABORT()
         ENDIF
           CALL CMUMPS_148(id%N, NZ, id%A(1),
     *      id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     *      LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     *      id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     *      KEEP(39),
     *      LP, id%COMM, id%root, KEEP,
     *      id%FILS(1), IWK(1), ! workspace of size N
     *
     *      id%INTARR(1), id%DBLARR(1),
     *      id%PTRAR(1), id%PTRAR(id%N+1),
     *      id%FRERE_STEPS(1), id%STEP(1), WK(1), ITOT ) 
        write(6,*) '!!! A,IRN,JCN are freed during facto '
        DEALLOCATE (id%A)
        DEALLOCATE (id%IRN)
        DEALLOCATE (id%JCN)
        ALLOCATE (id%S(id%MAXS),stat=ierr)
        IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = id%MAXS
          write(6,*) ' PB2 LARGEMAT'
          CALL CMUMPS_ABORT()
        ENDIF
        id%S(id%MAXS-ITOT+1:id%MAXS) = WK(1:ITOT)
        DEALLOCATE (WK)
#else
        CALL CMUMPS_148(id%N, NZ, id%A(1),
     *   id%IRN(1), id%JCN(1), id%SYM_PERM,  
     *   LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     *   id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     *   KEEP(39),
     *   LP, id%COMM, id%root, KEEP,
     *   id%FILS(1), IWK(1),
     *
     *   id%INTARR(1), id%DBLARR(1),
     *   id%PTRAR(1), id%PTRAR(id%N+1),
     *   id%FRERE_STEPS(1), id%STEP(1), id%S(1), id%MAXS ) 
#endif
      DEALLOCATE(IWK)
        IF ( PROKG ) THEN
          CALL CMUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL CMUMPS_291(TIME)
        END IF
      ELSE
        CALL CMUMPS_145( id%N,
     *       id%DBLARR( 1 ), MAX(1,KEEP( 13 )),
     *       id%INTARR( 1 ), MAX(1,KEEP( 14 )),
     *       id%PTRAR( 1 ),
     *       id%PTRAR(id%N+1),
     *       KEEP( 1 ), id%MYID, id%COMM,
     *       id%KEEP(39),
     *
     *       id%S(1), MAXS,
     *       id%root,
     *       id%PROCNODE_STEPS(1), id%NSLAVES,
     *       id%SYM_PERM, id%FRERE_STEPS(1), id%STEP(1),
     *       id%INFO(1), id%INFO(2) )
      ENDIF
      ELSE
      IF (PROKG ) THEN
        CALL CMUMPS_291(TIME)
      END IF
      IF ( I_AM_SLAVE ) THEN
        CALL CMUMPS_282( id%N,
     * id%NZ_loc,
     * id,
     * id%DBLARR(1), KEEP(13), id%INTARR(1),
     * KEEP(14), id%PTRAR(1), id%PTRAR(id%N+1),
     * KEEP(1), id%MYID_NODES,
     * id%COMM_NODES, id%KEEP(39),
     * id%S(1), MAXS, id%root, id%PROCNODE_STEPS,
     * id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     * id%ICNTL(1), id%INFO(1), NSEND, NLOCAL )
#if defined(LARGEMATRICES)
       IF (ASSOCIATED(id%IRN_loc)) THEN
            DEALLOCATE(id%IRN_loc)
         ENDIF
         IF (ASSOCIATED(id%JCN_loc)) THEN
            DEALLOCATE(id%JCN_loc)
         ENDIF
         IF (ASSOCIATED(id%A_loc)) THEN
            DEALLOCATE(id%A_loc)
         ENDIF
       write(6,*) ' Warning :', 
     &        ' id%A_loc, IRN_loc, JCN_loc deallocated !!! '
#endif
      IF (PROK) THEN
        WRITE(MP,120) NLOCAL, NSEND
      END IF
      END IF
      IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
        NSEND  = 0
        NLOCAL = 0
      END IF
      CALL MPI_REDUCE( NSEND, NSEND_TOT, 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( NLOCAL, NLOCAL_TOT, 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
        WRITE(MPG,125) NLOCAL_TOT, NSEND_TOT
      END IF
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( PROKG ) THEN
          CALL CMUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL CMUMPS_291(TIME)
      END IF
      END IF
      ELSE
        IF (PROKG ) THEN
          CALL CMUMPS_291(TIME)
        END IF
      IF ( id%MYID.eq.MASTER)
     *CALL CMUMPS_213( id%ELTPTR(1),
     *                        id%NELT,
     *                        MAXELT_SIZE )
      CALL CMUMPS_126( id%N, id%NELT,
     *     id%COMM, id%MYID,
     *     id%NSLAVES, id%PTRAR,
     *     id%PTRAR(id%NELT+2),
     *     id%INTARR(1), id%DBLARR(1),
     *     id%KEEP(1), MAXELT_SIZE,
     *     id%FRTPTR(1), id%FRTELT(1),
     *     id%S(1), id%MAXS, id%FILS,
     *     id, id%root )
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( PROKG ) THEN
          CALL CMUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL CMUMPS_291(TIME)
      END IF
      END IF
      IF ( I_AM_SLAVE )  THEN
        CMUMPS_LBUFR_BYTES = KEEP( 44 ) * KEEP( 35 )
        CMUMPS_LBUFR_BYTES = MAX( CMUMPS_LBUFR_BYTES,
     *                      100000 )
        PERLU          = KEEP( 12 )
        IF ( KEEP( 50 ) .EQ. 1 ) PERLU = MIN(0,PERLU)
        CMUMPS_LBUFR_BYTES = CMUMPS_LBUFR_BYTES
     *                 + INT( 4.0E0 * REAL(MAX(PERLU,0))*
     *                   REAL(CMUMPS_LBUFR_BYTES)/100E0)
        IF (id%NSLAVES.EQ.2) THEN
         CMUMPS_LBUF = INT( 1.01E0 * KEEP( 43 ) *KEEP( 35 )  )
        ELSE
         CMUMPS_LBUF = INT( 2.01E0 * KEEP( 43 ) *KEEP( 35 )  )
        ENDIF
        CMUMPS_LBUF = MAX( CMUMPS_LBUF, 100000 )
        CMUMPS_LBUF = CMUMPS_LBUF
     *                 + INT( 2.0E0 * REAL(MAX(PERLU,0))*
     *                   REAL(CMUMPS_LBUF)/100E0)
        CMUMPS_LBUF_INT = ( KEEP(56) + id%NSLAVES * id%NSLAVES ) * 4
     *               * KEEP(34)
        IF ( KEEP( 38 ) .NE. 0 ) THEN
          KKKK = CMUMPS_275( id%STEP(KEEP(38)),
     *         id%PROCNODE_STEPS, id%NSLAVES )
          IF ( KKKK .EQ. id%MYID ) THEN
             CMUMPS_LBUF_INT = CMUMPS_LBUF_INT + 2 * id%NSLAVES
     *                      * KEEP( 34 )
          END IF
        END IF
        IF ( MP .GT. 0 ) THEN
          WRITE( MP, 9999 ) PERLU, CMUMPS_LBUFR_BYTES,
     *                      CMUMPS_LBUF, CMUMPS_LBUF_INT
        END IF
 9999   FORMAT( /,' Allocated buffers',/,' ------------------',/,
     *  ' Percentage relaxation of buffer size ... = ', I9,/,
     *  ' Size of reception buffer in bytes ...... = ', I9,
     *  /,
     *  ' Size of async. emission buffer (bytes).. = ', I9,/,
     *  ' Small emission buffer (bytes) .......... = ', I9)
        CALL CMUMPS_55( CMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     *   ':Error allocating small Send buffer:IERR='
     *   ,IERR
          INFO(1)= -13
          INFO(2)= CMUMPS_LBUF_INT
          GO TO 110
        END IF
        CALL CMUMPS_53( CMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating Send buffer:IERR='
     *   ,IERR
          INFO(1)= -13
          INFO(2)= CMUMPS_LBUF
          GO TO 110
        END IF
        id%LBUFR_BYTES = CMUMPS_LBUFR_BYTES
        id%LBUFR = (CMUMPS_LBUFR_BYTES+KEEP(34)-1)/KEEP(34)
        IF (ASSOCIATED(id%BUFR)) DEALLOCATE(id%BUFR)
        ALLOCATE( id%BUFR( id%LBUFR ),stat=ierr )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating BUFR:IERR='
     *   ,IERR
          INFO(1)=-13
          INFO(2)=id%LBUFR
          GO TO 110
        END IF
        PERLU          = KEEP( 12 )
        id%MAXIS = MAX( id%MAXIS, 1,
     *       KEEP(15) + 2 * MAX(PERLU,10) * ( KEEP(15) / 100 + 1 )
     *  )
        IF (ASSOCIATED(id%IS)) DEALLOCATE( id%IS )
        ALLOCATE( id%IS( id%MAXIS  ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
        WRITE(*,*) id%MYID,':Error allocating IS:IERR=',IERR
        INFO(1)=-13
        INFO(2)=id%MAXIS
        GO TO 110
        END IF
        LIW = id%MAXIS
        IF (ASSOCIATED( id%PTLUST_S )) DEALLOCATE(id%PTLUST_S)
        ALLOCATE( id%PTLUST_S( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocatingPTLUST:IERR = ',
     *    IERR
          INFO(1)=-13
          INFO(2)=id%KEEP(28)
          GOTO 100
        END IF
        PTRIST = 1
        PTRWB  = PTRIST + id%KEEP(28)
        ITLOC  = PTRWB  + 5 * id%KEEP(28)
        IPOOL  = ITLOC  + id%N
        LPOOL  = CMUMPS_505(id%KEEP(1))
        ALLOCATE( IWK(  IPOOL + LPOOL - 1 ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating IWK : IERR = ',
     *    IERR
          INFO(1)=-13
          INFO(2)=IPOOL + LPOOL - 1
          GOTO 110
        END IF
      ENDIF
 110  CONTINUE
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( I_AM_SLAVE )  THEN
        CALL CMUMPS_60( id%LBUFR_BYTES )
        KEEP(29) = KEEP(15)
        KEEP(30) = KEEP(16) + KEEP(44)
        IF (MP .GT. 0) THEN
          WRITE( MP, 170 ) MAXS, MAXIS, KEEP(30), KEEP(29), KEEP(13),
     *    KEEP(14), KEEP(25), KEEP(26), KEEP(27)
        ENDIF
      END IF
      CALL CMUMPS_214( id%KEEP(1), id%MYID, id%N, id%NELT, id%LNA,
     *     id%NSLAVES, MEMORY_SIZE, .FALSE. )
      id%INFO(16) = MEMORY_SIZE
      IF ( MP .gt. 0 ) THEN
          WRITE(MP,'(A,I10) ')
     *    ' ** Space in MBYTES used during factorization  :',
     *                MEMORY_SIZE
      END IF
      CALL CMUMPS_243( id%MYID, id%COMM,
     *                           id%INFO(16), id%INFOG(18), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     *  ' ** Memory relaxation parameter ( ICNTL(14)  )            :',
     *  KEEP(12)
        WRITE( MPG,'(A,I10) ')
     *  ' ** Rank of processor needing largest memory in facto     :',
     *  IRANK
        WRITE( MPG,'(A,I10) ')
     *  ' ** Space in MBYTES used by this processor for facto      :',
     *  id%INFOG(18)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     *  ' ** Avg. Space in MBYTES per working proc during facto    :',
     *  ( id%INFOG(19)-id%INFO(16) ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     *  ' ** Avg. Space in MBYTES per working proc during facto    :',
     *  id%INFOG(19) / id%NSLAVES
        END IF
      END IF
      IF ( I_AM_SLAVE ) THEN
        id%PTLUST_S = -1
        IF ( KEEP(55) .eq. 0 ) THEN
          LDPTRAR = id%N
        ELSE
          LDPTRAR = id%NELT + 1
        END IF
        IF ( id%KEEP(55) .NE. 0 ) THEN
          NELT = id%NELT
        ELSE
          NELT = 1
        END IF
        CALL CMUMPS_244( id%N, NZ, NSTEPS, id%S(1),
     *      id%MAXS, id%IS( 1 ), LIW,
     *      id%SYM_PERM(1), id%NA(1), id%LNA, id%NE_STEPS(1),
     *      id%ND_STEPS(1), id%FILS(1), id%STEP(1),
     *      id%FRERE_STEPS(1), id%CANDIDATES(1,1), 
     *      id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     *      id%PTRAR(1), LDPTRAR, IWK( PTRIST ),
     *      id%PTLUST_S(1), IWK( PTRWB ),
     *      IWK( ITLOC ), IWK( IPOOL ), LPOOL,
     *      id%CNTL(1), ICNTL, INFO, RINFO, KEEP,
     *      id%PROCNODE_STEPS(1),
     *      id%NSLAVES, id%COMM_NODES,
     *      id%MYID, id%MYID_NODES,
     *      id%BUFR(1),id%LBUFR,id%LBUFR_BYTES,
     *      id%INTARR(1), id%DBLARR(1), id%root,
     *      NELT, id%FRTPTR(1), 
     *      id%FRTELT(1), id%COMM_LOAD, id%ASS_IRECV, SEUIL,
     *      SEUIL_LDLT_NIV2, id%MEM_DIST(0))
        IF ( MP . GT. 0 .and. KEEP(38) .ne. 0 ) THEN
          WRITE( MP, 175 ) KEEP(49)
        END IF
        DEALLOCATE( IWK )
      ENDIF
        IF ( KEEP(55) .eq. 0 ) THEN
          IF (ASSOCIATED( id%DBLARR)) DEALLOCATE(id%DBLARR)
        ELSE
          DEALLOCATE( id%INTARR)
          IF (      id%MYID_NODES .eq. MASTER
     *      .AND.   KEEP(46)   .eq. 1
     *      .AND.   KEEP(52)   .eq. 0 ) THEN
            NULLIFY( id%DBLARR )
          ELSE
            IF (ASSOCIATED( id%DBLARR)) DEALLOCATE(id%DBLARR)
          END IF
        END IF
      CALL CMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( KEEP(19) .NE. 0 ) THEN
        IF ( KEEP(46) .NE. 1 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( KEEP(17), 1, MPI_INTEGER, 1, DEFIC_TAG,
     *                   id%COMM, STATUS, IERR )
          ELSE IF ( id%MYID .EQ. 1 ) THEN
            CALL MPI_SEND( KEEP(17), 1, MPI_INTEGER, 0, DEFIC_TAG,
     *                   id%COMM, IERR )
          END IF
        END IF
        id%Deficiency = id%KEEP(17)
      END IF
      IF ( KEEP(19) .EQ. 2 .OR. KEEP(19) .EQ. 4
     *.OR. KEEP(19) .EQ. 6 .OR. KEEP(19) .EQ. 8
     *.OR. KEEP(19) .EQ. 10
     * ) THEN
        CALL CMUMPS_84( id )
      END IF
      CALL CMUMPS_57( IERR )
      CALL CMUMPS_59( IERR )
      IF ( KEEP(60) .ne. 0 ) THEN
        ID_SCHUR =CMUMPS_275(id%STEP(id%KEEP(20)),
     *            id%PROCNODE_STEPS(1),id%NSLAVES)
        IF ( id%KEEP( 46 )  .NE. 1 ) THEN
          ID_SCHUR = ID_SCHUR + 1
        END IF
        CALL CMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
        IF ( id%INFO(1) < 0 ) GOTO 500
        IF ( id%MYID .eq. MASTER .and. ID_SCHUR .eq. MASTER ) THEN
          CALL CCOPY( id%SIZE_SCHUR * id%SIZE_SCHUR,
     *       id%S(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+4)), 1,
     *       id%SCHUR(1), 1 )
        ELSE
          IF ( id%MYID .eq. ID_SCHUR ) THEN
            SIZE_SCHUR = id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+0)
            CALL MPI_SEND( id%S(
     *              id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+4)),
     *              SIZE_SCHUR*SIZE_SCHUR,
     *              MPI_COMPLEX,
     *              MASTER, TAG_SCHUR,
     *              id%COMM, IERR )
          ELSE IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( id%SCHUR(1), id%SIZE_SCHUR*id%SIZE_SCHUR,
     *                     MPI_COMPLEX, ID_SCHUR, TAG_SCHUR,
     *                     id%COMM, STATUS, IERR )
          END IF
        END IF
      END IF
      IF ( PROKG ) THEN
        CALL CMUMPS_292(TIME)
        WRITE(MPG,180) TIME
      END IF
      IF (INFO(1).LT.0) GOTO 500
      IF (id%MYID.EQ.MASTER) THEN
        IF ( id%KEEP(46) .eq. 0 ) THEN
          RINFO(2)=ZERO
          RINFO(3)=ZERO
          INFO(9:14) = 0
        END IF
      ENDIF
      KEEP(33) = INFO(11)
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     *                 MPI_REAL,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     *                 MPI_REAL,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( INFO(9), INFOG(9), 2, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_ALLREDUCE( INFO(11), INFOG(11), 1, MPI_INTEGER,
     *                 MPI_MAX, id%COMM, IERR)
      KEEP(133) = INFOG(11)
      CALL MPI_REDUCE( INFO(12), INFOG(12), 3, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      IF ( PROKG ) THEN
          WRITE(MPG,99984) RINFOG(2),RINFOG(3),INFOG(9),INFOG(10),
     *                    INFOG(11)
          IF (id%SYM == 0) THEN
            WRITE(MPG, 99982) 'OFF DIAGONAL PIVOTS',INFOG(12)
          END IF
          IF (id%SYM .NE. 1) THEN
            WRITE(MPG, 99982) 'DELAYED PIVOTS     ',INFOG(13)
          END IF
          WRITE(MPG, 99981) INFOG(14)
          IF ( id%Deficiency .ne. 0 )
     *    WRITE(MPG, 99983) id%Deficiency
      END IF
99983 FORMAT( '            DEFICIENCY  DETECTED                =',
     *          I12)
99984 FORMAT(/' GLOBAL STATISTICS '/
     1  ' RINFOG(2)  OPERATIONS DURING NODE ASSEMBLY     =',1PD10.3/
     2  ' ------(3)  OPERATIONS DURING NODE ELIMINATION  =',1PD10.3/
     3  ' INFOG (9)  REAL SPACE FOR FACTORS              =',I12/
     4  ' INFOG(10)  INTEGER SPACE FOR FACTORS           =',I12/
     5  ' INFOG(11)  MAXIMUM FRONT SIZE                  =',I12)
99982 FORMAT( ' INFOG(12)  NUMBER OF ',A19,'       =',I12)
99981 FORMAT(
     1  ' INFOG(14)  NUMBER OF MEMORY COMPRESS           =',I12)
 500  IF ( I_AM_SLAVE ) THEN
        CALL CMUMPS_183( INFO(1), IERR )
      END IF
 550  id%KEEP(13) = KEEP13_SAVE
      RETURN
 120  FORMAT(/' LOCAL REDISTRIB: DATA LOCAL/SENT     =',I12,I12)
 125  FORMAT(/' REDISTRIB: TOTAL DATA LOCAL/SENT     =',I12,I12)
 130  FORMAT(/' ****** FACTORIZATION STEP ********'/)
 160  FORMAT(' GLOBAL TIME FOR MATRIX DISTRIBUTION  =',F12.4)
 170  FORMAT(/' STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     *        ' MAXS                                 =',I12/
     *        ' MAXIS                                =',I12/
     *        ' MINIMUM (ICNTL(14)=0) VALUE OF MAXS  =',I12/
     *        ' MINIMUM (ICNTL(14)=0) VALUE OF MAXIS =',I12/
     *        ' REAL SPACE FOR ORIGINAL MATRIX       =',I12/
     *        ' INTEGER SPACE FOR ORIGINAL MATRIX    =',I12/
     *        ' REAL SPACE FOR FACTORS               =',I12/
     *        ' INTEGER SPACE FOR FACTORS            =',I12/
     *        ' MAXIMUM FRONTAL SIZE (ESTIMATED)     =',I12)
 172  FORMAT(/' GLOBAL STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     *        ' NUMBER OF WORKING PROCESSES          =',I12/
     *        ' REAL SPACE FOR FACTORS               =',I12/
     *        ' INTEGER SPACE FOR FACTORS            =',I12/
     *        ' MAXIMUM FRONTAL SIZE (ESTIMATED)     =',I12/
     *        ' NUMBER OF NODES IN THE TREE          =',I12)
 175  FORMAT(/' NUMBER OF ENTRIES FOR // ROOT                 =',I12)
 180  FORMAT(/' ELAPSED TIME FOR FACTORIZATION       =',F12.4)
      END SUBROUTINE CMUMPS_142
      SUBROUTINE CMUMPS_83
     * ( N, MAPPING, NZ, IRN, JCN, PROCNODE, STEP,
     *   SLAVEF, PERM, FILS,
     *   RG2L, KEEP, MBLOCK, NBLOCK, NPROW, NPCOL )
      IMPLICIT NONE
      INTEGER N, NZ, SLAVEF, MBLOCK, NBLOCK, NPROW, NPCOL
      INTEGER KEEP(150)
      INTEGER IRN( NZ ), JCN( NZ ) 
      INTEGER MAPPING( NZ ), STEP( N )
      INTEGER PROCNODE( KEEP(28) ), PERM( N ), FILS( N ), RG2L( N )
      INTEGER CMUMPS_275, CMUMPS_330
      EXTERNAL CMUMPS_275, CMUMPS_330
      INTEGER K, IOLD, JOLD, INEW, JNEW, ISEND, JSEND, IARR, INODE
      INTEGER TYPE_NODE, DEST
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INODE = KEEP(38)
      K = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = K
        INODE = FILS( INODE )
        K = K + 1
      END DO
      DO K = 1, NZ
        IOLD = IRN( K )
        JOLD = JCN( K )
        IF ( IOLD .GT. N .OR. IOLD .LT. 1 .OR.
     *       JOLD .GT. N .OR. JOLD .LT. 1 ) THEN
           MAPPING( K ) = -1
           CYCLE
        END IF
        IF ( IOLD .eq. JOLD ) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM( IOLD )
          JNEW = PERM( JOLD )
          IF ( INEW .LT. JNEW ) THEN
            ISEND = IOLD
            IF ( KEEP(50) .ne. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          END IF
        END IF
        IARR = ABS( ISEND )
        TYPE_NODE = CMUMPS_330( ABS(STEP(IARR)),
     *                              PROCNODE, SLAVEF )
        IF ( TYPE_NODE .eq. 1 .or. TYPE_NODE .eq. 2 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = CMUMPS_275( ABS(STEP(IARR)),
     *                             PROCNODE, SLAVEF ) + 1
          ELSE
            DEST = CMUMPS_275( ABS(STEP(IARR)),
     *                             PROCNODE, SLAVEF )
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L( JSEND )
            JPOSROOT = RG2L( IARR  )
          ELSE
            IPOSROOT = RG2L( IARR  )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = MOD( ( IPOSROOT - 1 )/MBLOCK, NPROW )
          JCOL_GRID = MOD( ( JPOSROOT - 1 )/NBLOCK, NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * NPCOL + JCOL_GRID
          END IF
        END IF
        MAPPING( K ) = DEST
      END DO
      RETURN
      END SUBROUTINE CMUMPS_83
      SUBROUTINE CMUMPS_282(
     * N, NZ_loc, id,
     * DBLARR, LDBLARR, INTARR, LINTARR,
     * PTRAIW, PTRARW, KEEP, MYID, COMM, NBRECORDS,
     *
     * A, LA, root, PROCNODE_STEPS, SLAVEF, PERM, STEP,
     * ICNTL, INFO, NSEND, NLOCAL )
      IMPLICIT NONE
      INCLUDE 'cmumps_struc.h'
      INTEGER N, NZ_loc
      TYPE (CMUMPS_STRUC) :: id
      INTEGER LDBLARR, LINTARR
      COMPLEX DBLARR( LDBLARR )
      INTEGER INTARR( LINTARR )
      INTEGER PTRAIW( N ), PTRARW( N )
      INTEGER KEEP(150)
      INTEGER MYID, COMM, NBRECORDS
      INTEGER LA, SLAVEF
      COMPLEX A( LA )
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER PROCNODE_STEPS(KEEP(28)), PERM( N ), STEP( N )
      INTEGER INFO( 40 ), ICNTL(40)
      INTEGER CMUMPS_275, CMUMPS_330, NUMROC
      EXTERNAL CMUMPS_275, CMUMPS_330, NUMROC
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE ), MSGSOU
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      INTEGER END_MSG_2_RECV
      INTEGER I, K, I1, IA
      INTEGER TYPE_NODE, DEST
      INTEGER IOLD, JOLD, IARR, ISEND, JSEND, INEW, JNEW
      INTEGER allocok 
      COMPLEX VAL
      INTEGER PTR_ROOT, LOCAL_M, LOCAL_N, ARROW_ROOT
      INTEGER IROW_GRID, JCOL_GRID, IPOSROOT, JPOSROOT
      INTEGER MP,LP
      INTEGER KPROBE, FREQPROBE
      PARAMETER( FREQPROBE = 1000 )
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFI
      COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: BUFR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI
      COMPLEX, ALLOCATABLE, DIMENSION(:) :: BUFRECR
      INTEGER IACT( SLAVEF ), IREQI( SLAVEF ), IREQR( SLAVEF )
      LOGICAL SEND_ACTIVE( SLAVEF )
      LOGICAL FLAG
      INTEGER NSEND, NLOCAL
      NSEND = 0
      NLOCAL = 0
      LP = ICNTL(1)
      MP = ICNTL(2)
      END_MSG_2_RECV = SLAVEF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, 2, SLAVEF ), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *     '** Error allocating int buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = ( NBRECORDS * 2 + 1 ) * SLAVEF * 2
      END IF
      ALLOCATE( BUFR( NBRECORDS, 2, SLAVEF), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *     '** Error allocating real buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFRECI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * 2 + 1
        GOTO 20
      END IF
      ALLOCATE( BUFRECR( NBRECORDS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS
        GOTO 20
      END IF
      ALLOCATE( IW4( N, 2 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(LP,*) '** Error allocating IW4 for matrix distribution'
        INFO(1) = -13
        INFO(2) = N * 2
      END IF
 20   CONTINUE
      CALL CMUMPS_276( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      ARROW_ROOT = 0
      DO I = 1, N
          I1 = PTRAIW( I )
          IA = PTRARW( I )
          IF ( IA .GT. 0 ) THEN
            DBLARR( IA ) = CMPLX(ZERO)
            IW4( I, 1 ) = INTARR( I1 )
            IW4( I, 2 ) = -INTARR( I1 + 1 )
            INTARR( I1 + 2 ) = I
          END IF
      END DO
      IF ( KEEP(38) .NE. 0 ) THEN
          LOCAL_M = NUMROC( root%ROOT_SIZE, root%MBLOCK,
     *               root%MYROW, 0, root%NPROW )
          LOCAL_M = MAX( 1, LOCAL_M )
          LOCAL_N = NUMROC( root%ROOT_SIZE, root%NBLOCK,
     *               root%MYCOL, 0, root%NPCOL )
          PTR_ROOT = LA - LOCAL_M * LOCAL_N + 1
          IF ( PTR_ROOT .LE. LA ) THEN
            A( PTR_ROOT:LA ) = CMPLX(ZERO)
          END IF
      END IF
      DO I = 1, SLAVEF
        BUFI( 1, 1, I ) = 0
      END DO
      DO I = 1, SLAVEF
        BUFI( 1, 2, I ) = 0
      END DO
      DO I = 1, SLAVEF
        SEND_ACTIVE( I ) = .FALSE.
        IACT( I ) = 1
      END DO
      KPROBE = 0
      DO K = 1, NZ_loc
        KPROBE = KPROBE + 1
        IF ( KPROBE .eq. FREQPROBE ) THEN
          KPROBE = 0
          CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     *                     FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            MSGSOU = STATUS( MPI_SOURCE )
            CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, 
     *                 MPI_INTEGER,
     *                 MSGSOU, ARR_INT, COMM, STATUS, IERR )
            CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_COMPLEX,
     *                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
            CALL CMUMPS_102(
     *             BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     *             KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *             A, LA,
     *             END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *             ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *             INTARR, LINTARR, DBLARR, LDBLARR
     *             )
          END IF
        END IF
        IOLD = id%IRN_loc(K)
        JOLD = id%JCN_loc(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     *                 .OR.(JOLD.LT.1) ) CYCLE
        VAL = id%A_loc(K)
        IF (IOLD.EQ.JOLD) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM(IOLD)
          JNEW = PERM(JOLD)
          IF (INEW.LT.JNEW) THEN
            ISEND = IOLD
            IF ( KEEP(50) .NE. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          ENDIF
        ENDIF
        IARR = ABS( ISEND )
        TYPE_NODE = CMUMPS_330( ABS(STEP(IARR)),
     *              PROCNODE_STEPS, SLAVEF )
        IF ( TYPE_NODE .eq. 1 ) THEN
          DEST = CMUMPS_275( ABS(STEP(IARR)),
     *           PROCNODE_STEPS, SLAVEF )
        ELSE IF ( TYPE_NODE .eq. 2 ) THEN
          IF ( ISEND .LT. 0 ) THEN
            DEST = -1
          ELSE
            DEST = CMUMPS_275( ABS(STEP(IARR)),
     *             PROCNODE_STEPS, SLAVEF )
          END IF
        ELSE
          IF ( ISEND < 0 ) THEN
            IPOSROOT = root%RG2L_ROW(JSEND)
            JPOSROOT = root%RG2L_ROW(IARR )
          ELSE
            IPOSROOT = root%RG2L_ROW(IARR )
            JPOSROOT = root%RG2L_ROW(JSEND)
          END IF
          IROW_GRID = MOD( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = MOD( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          DEST = IROW_GRID * root%NPCOL + JCOL_GRID
        END IF
        if (dest .eq. -1) then
          nlocal = nlocal + 1
          nsend = nsend + slavef -1
        else
          if (dest .eq.myid ) then
            nlocal = nlocal + 1
          else
            nsend = nsend + 1
          endif
        end if
        CALL CMUMPS_101( DEST, ISEND, JSEND, VAL,
     *  BUFI, BUFR, BUFRECI, BUFRECR,
     *  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     *  N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1), 
     *  root, KEEP )
      END DO
      DEST = -2
        CALL CMUMPS_101( DEST, ISEND, JSEND, VAL,
     *  BUFI, BUFR, BUFRECI, BUFRECR,
     *  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     *  N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, 
     *  IW4(1,1), root, KEEP )
      DO WHILE ( END_MSG_2_RECV .NE. 0 )
        CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, MPI_INTEGER,
     *                 MPI_ANY_SOURCE, ARR_INT, COMM, STATUS, IERR )
        MSGSOU = STATUS( MPI_SOURCE )
        CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_COMPLEX,
     *                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
        CALL CMUMPS_102(
     *           BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     *           KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *           A, LA,
     *           END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *           ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *           INTARR, LINTARR, DBLARR, LDBLARR
     *           )
      END DO
      DO I = 1, SLAVEF
        IF ( SEND_ACTIVE( I ) ) THEN
          CALL MPI_WAIT( IREQI( I ), STATUS, IERR )
          CALL MPI_WAIT( IREQR( I ), STATUS, IERR )
        END IF
      END DO
      KEEP(49) = ARROW_ROOT
      DEALLOCATE( IW4 )
      DEALLOCATE( BUFI )
      DEALLOCATE( BUFR )
      DEALLOCATE( BUFRECI )
      DEALLOCATE( BUFRECR )
      RETURN
      END SUBROUTINE CMUMPS_282
      SUBROUTINE CMUMPS_192
     *( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     *  LDLT, MTYPE)
      IMPLICIT NONE
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      COMPLEX A_loc( NZ_loc ), X( N ), Y_loc( N )
      INTEGER LDLT, MTYPE
      INTEGER I, J, K
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      Y_loc = CMPLX(ZERO)
      IF ( LDLT .eq. 0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *          (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
        ENDDO
        ELSE
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     *        .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
        ENDDO
        END IF
      ELSE
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *        (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE CMUMPS_192
      SUBROUTINE CMUMPS_101( DEST, ISEND, JSEND, VAL,
     *  BUFI, BUFR, BUFRECI, BUFRECR,
     *  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR, N,
     *  PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4, root,
     *  KEEP )
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER ISEND, JSEND, DEST, NBRECORDS, SLAVEF, COMM, MYID, N
      INTEGER KEEP(150)
      INTEGER ARROW_ROOT, END_MSG_2_RECV, LOCAL_M, LOCAL_N
      INTEGER LINTARR, LDBLARR, LA, PTR_ROOT
      INTEGER BUFI( NBRECORDS * 2 + 1, 2, SLAVEF )
      INTEGER BUFRECI( NBRECORDS * 2 + 1 )
      INTEGER IREQI(SLAVEF), IREQR(SLAVEF), IACT(SLAVEF)
      INTEGER IW4( N, 2 )
      INTEGER PTRAIW( N ), PTRARW( N ), PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR( LINTARR )
      COMPLEX DBLARR( LDBLARR ), A( LA )
      LOGICAL SEND_ACTIVE(SLAVEF)
      COMPLEX BUFR( NBRECORDS, 2, SLAVEF )
      COMPLEX BUFRECR( NBRECORDS )
      COMPLEX VAL
      INTEGER ISLAVE, IBEG, IEND, NBREC, IREQ
      INTEGER TAILLE_SEND_I, TAILLE_SEND_R, MSGSOU
      LOGICAL FLAG, SEND_LOCAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, STATUS(MPI_STATUS_SIZE)
      IF ( DEST .eq. -1 .or. DEST .eq. -2 ) THEN
        IBEG = 1
        IEND = SLAVEF
      ELSE
        IBEG = DEST + 1
        IEND = DEST + 1
      END IF
      SEND_LOCAL = .FALSE.
      DO ISLAVE = IBEG, IEND
        NBREC = BUFI(1,IACT(ISLAVE),ISLAVE)
        IF ( DEST .eq. -2 ) THEN
          BUFI(1,IACT(ISLAVE),ISLAVE) = - NBREC
        END IF
        IF ( DEST .eq. -2 .or. NBREC + 1 > NBRECORDS ) THEN
          DO WHILE ( SEND_ACTIVE( ISLAVE ) )
            CALL MPI_TEST( IREQI( ISLAVE ), FLAG, STATUS, IERR )
            IF ( .NOT. FLAG ) THEN
                CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     *                           FLAG, STATUS, IERR )
                IF ( FLAG ) THEN
                  MSGSOU = STATUS(MPI_SOURCE)
                  CALL MPI_RECV( BUFRECI(1), 2*NBRECORDS+1,
     *                  MPI_INTEGER, MSGSOU, ARR_INT, COMM,
     *                  STATUS, IERR )
                  CALL MPI_RECV( BUFRECR(1), NBRECORDS,
     *                  MPI_COMPLEX, MSGSOU,
     *                  ARR_REAL, COMM, STATUS, IERR )
                  CALL CMUMPS_102(
     *              BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     *              KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *              A, LA,
     *              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *              INTARR, LINTARR, DBLARR, LDBLARR
     *              )
                END IF
            ELSE
                CALL MPI_WAIT( IREQR( ISLAVE ), STATUS, IERR )
                SEND_ACTIVE( ISLAVE ) = .FALSE.
            END IF
          END DO
          IF ( ISLAVE - 1 .ne. MYID ) THEN
            TAILLE_SEND_I = NBREC * 2 + 1
            TAILLE_SEND_R = NBREC
            CALL MPI_ISEND( BUFI(1, IACT(ISLAVE), ISLAVE ),
     *          TAILLE_SEND_I,
     *          MPI_INTEGER, ISLAVE - 1, ARR_INT, COMM,
     *          IREQI( ISLAVE ), IERR )
            CALL MPI_ISEND( BUFR(1, IACT(ISLAVE), ISLAVE ),
     *          TAILLE_SEND_R,
     *          MPI_COMPLEX, ISLAVE - 1, ARR_REAL, COMM,
     *          IREQR( ISLAVE ), IERR )
            SEND_ACTIVE( ISLAVE ) = .TRUE.
          ELSE
            SEND_LOCAL = .TRUE.
          END IF
          IACT( ISLAVE ) = 3 - IACT( ISLAVE )
          BUFI( 1, IACT( ISLAVE ), ISLAVE ) = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          IREQ = BUFI(1,IACT(ISLAVE),ISLAVE) + 1
          BUFI(1,IACT(ISLAVE),ISLAVE) = IREQ
          BUFI(IREQ*2,IACT(ISLAVE),ISLAVE)  = ISEND
          BUFI(IREQ*2+1,IACT(ISLAVE),ISLAVE) = JSEND
          BUFR(IREQ,IACT(ISLAVE),ISLAVE )    = VAL
        END IF
      END DO
      IF ( SEND_LOCAL ) THEN
            ISLAVE = MYID + 1
            CALL CMUMPS_102(
     *              BUFI(1,3-IACT(ISLAVE),ISLAVE),
     *              BUFR(1,3-IACT(ISLAVE),ISLAVE),
     *              NBRECORDS, N, IW4(1,1),
     *              KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *              A, LA,
     *              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *              INTARR, LINTARR, DBLARR, LDBLARR
     *              )
      END IF
      RETURN
      END SUBROUTINE CMUMPS_101
      SUBROUTINE CMUMPS_102
     *           ( BUFI, BUFR, NBRECORDS, N, IW4,
     *             KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT, A, LA,
     *             END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     *             SLAVEF, ARROW_ROOT,
     *             PTRAIW, PTRARW, PERM, STEP,
     *             INTARR, LINTARR, DBLARR, LDBLARR )
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER NBRECORDS, N, ARROW_ROOT, MYID, SLAVEF
      INTEGER BUFI( NBRECORDS * 2 + 1 )
      COMPLEX BUFR( NBRECORDS )
      INTEGER IW4( N, 2 )
      INTEGER KEEP(150)
      INTEGER END_MSG_2_RECV
      INTEGER PTRAIW( N ), PTRARW( N ), PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER LINTARR, LDBLARR
      INTEGER INTARR( LINTARR )
      INTEGER LOCAL_M, LOCAL_N, PTR_ROOT, LA
      COMPLEX A( LA ), DBLARR( LDBLARR )
      INTEGER CMUMPS_330, CMUMPS_275
      EXTERNAL CMUMPS_330, CMUMPS_275
      INTEGER IREC, NB_REC, NODE_TYPE, IPROC
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID,
     *        ILOCROOT, JLOCROOT
      INTEGER IA, IS1, ISHIFT, IIW, IS, IAS, IARR, JARR
      INTEGER TAILLE
      COMPLEX VAL
      NB_REC = BUFI( 1 )
      IF ( NB_REC .LE. 0 ) THEN
        END_MSG_2_RECV = END_MSG_2_RECV - 1
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) GOTO 100
      DO IREC = 1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        NODE_TYPE = CMUMPS_330( ABS(STEP(ABS( IARR ))),
     *              PROCNODE_STEPS, SLAVEF )
        IF ( NODE_TYPE .eq. 3 ) THEN
          ARROW_ROOT = ARROW_ROOT + 1
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
          END IF
          IROW_GRID = MOD( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = MOD( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          IF ( IROW_GRID .NE. root%MYROW .OR.
     *       JCOL_GRID .NE. root%MYCOL ) THEN
            WRITE(*,*) MYID,':INTERNAL Error: recvd root arrowhead '
            WRITE(*,*) MYID,':not belonging to me. IARR,JARR=',IARR,JARR
            WRITE(*,*) MYID,':IROW_GRID,JCOL_GRID=',IROW_GRID,JCOL_GRID
            WRITE(*,*) MYID,':MYROW, MYCOL=', root%MYROW, root%MYCOL
            WRITE(*,*) MYID,':IPOSROOT,JPOSROOT=', IPOSROOT, JPOSROOT
            CALL CMUMPS_ABORT()
          END IF
          ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
          JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
          A( PTR_ROOT + ( JLOCROOT - 1 ) * LOCAL_M
     *      + ILOCROOT - 1 ) =  A( PTR_ROOT + (JLOCROOT - 1)
     *      * LOCAL_M + ILOCROOT - 1 ) + VAL
        ELSE IF (IARR.GE.0) THEN
         IF (IARR.EQ.JARR) THEN
          IA = PTRARW(IARR)
          DBLARR(IA) = DBLARR(IA) + VAL
         ELSE
          IS1 =  PTRAIW(IARR)
          ISHIFT      = INTARR(IS1) + IW4(IARR,2)
          IW4(IARR,2) = IW4(IARR,2) - 1
          IIW         = IS1 + ISHIFT + 2
          INTARR(IIW)     = JARR
          IS          = PTRARW(IARR)
          IAS         = IS + ISHIFT
          DBLARR(IAS) = VAL
         ENDIF
        ELSE
           IARR = -IARR
           ISHIFT      = PTRAIW(IARR)+IW4(IARR,1)+2
           INTARR(ISHIFT)  = JARR
           IAS         = PTRARW(IARR)+IW4(IARR,1)
           IW4(IARR,1) = IW4(IARR,1) - 1
           DBLARR(IAS)      = VAL
           IPROC = CMUMPS_275( ABS(STEP(IARR)),
     *             PROCNODE_STEPS,SLAVEF )
           IF ( KEEP(50) .NE. 0 .AND.
     *          IW4(IARR,1) .EQ. 0 .AND.
     *          IPROC .EQ. MYID
     *          .AND. STEP(IARR) > 0 ) THEN
             TAILLE = INTARR( PTRAIW(IARR) )
             CALL CMUMPS_310( N, PERM,
     *            INTARR( PTRAIW(IARR) + 3 ),
     *            DBLARR( PTRARW(IARR) + 1 ),
     *            TAILLE, MYID )
           END IF
        ENDIF
      ENDDO
 100  CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_102
      SUBROUTINE CMUMPS_193
     *( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     *  LDLT, MTYPE)
      IMPLICIT NONE
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      COMPLEX A_loc( NZ_loc ), X( N ), Y_loc( N )
      INTEGER LDLT, MTYPE
      INTEGER I, J, K
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      Y_loc = CMPLX(ZERO)
      IF ( LDLT .eq. 0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *          (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + ABS( A_loc(K) * X(J) )
        ENDDO
        ELSE
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     *        .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(J) = Y_loc(J) + ABS( A_loc(K) * X(I) )
        ENDDO
        END IF
      ELSE
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *        (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + ABS( A_loc(K) * X(J) )
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + ABS( A_loc(K) * X(I) )
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE CMUMPS_193
