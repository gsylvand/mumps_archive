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
      INTEGER FUNCTION SMUMPS_50 
     *         ( SLAVEF, K48, K10, K50,
     *         NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K10, K50, NFRONT, NCB
      INTEGER NSLAVESMIN, NASS, KMAX
      REAL Wmaster, Wtotal, Wmax
      INTEGER I,ACC,X
      REAL SMUMPS_45
      INTEGER SMUMPS_497    
      EXTERNAL  SMUMPS_45, SMUMPS_497
      KMAX = SMUMPS_497( K10, NCB )
      NASS = NFRONT - NCB
      NSLAVESMIN  = 1
      IF ( K48 .EQ.0 ) THEN 
         NSLAVESMIN = MAX(NCB/MAX(1,KMAX),1)
      ELSE IF (K48 .EQ. 3 ) THEN
        Wmax    = SMUMPS_45(KMAX,NFRONT,NASS)
        Wtotal  = SMUMPS_45(NCB,NFRONT,NASS)
        Wmaster = REAL(NASS*NASS)*REAL(NASS)/(3.0)
        IF ( Wmaster .GT.  Wmax ) THEN 
         NSLAVESMIN = MAX ( NINT ( Wtotal / Wmaster ), 1 )
        ELSE
         NSLAVESMIN = MAX ( NINT ( Wtotal / Wmax ), 1 )
        ENDIF
      ELSE IF (K48 .EQ. 4 ) THEN
         IF ( K10 > 0 ) THEN
           WRITE(*,*) 'Internal Error 1 in SMUMPS_50'
           STOP
         ENDIF
         KMAX=ABS(K10)
         IF(K50.EQ.0)THEN
            NSLAVESMIN = MAX(INT((NFRONT*NCB)/KMAX),1)
         ELSE
            ACC=0
            NSLAVESMIN=0
            DO WHILE (ACC.NE.NCB)
               x=INT((-(NFRONT-NCB+ACC)
     $              +SQRT(REAL(((NFRONT-NCB+ACC)*
     $              (NFRONT-NCB+ACC))+REAL(4)*
     *              (KMAX))))/
     $              REAL(2))
               ACC=ACC+x
               NSLAVESMIN=NSLAVESMIN+1
               IF (((NCB-(ACC+NFRONT-NCB))*NFRONT).LT.KMAX)THEN
                  ACC=NCB 
                  NSLAVESMIN=NSLAVESMIN+1
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      NSLAVESMIN = MIN ( NSLAVESMIN,(SLAVEF-1) )
      SMUMPS_50 = 
     *               MIN ( NSLAVESMIN, NCB )
      RETURN
      END FUNCTION SMUMPS_50 
      INTEGER FUNCTION SMUMPS_52 
     *        ( SLAVEF, K48, K10, K50,
     *          NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K10, K50,NFRONT, NCB
      INTEGER NSLAVESMAX, KMAX, KMIN
      INTEGER NSLAVESMIN
      INTEGER SMUMPS_497,SMUMPS_442,
     *        SMUMPS_50,
     *        SMUMPS_46
      EXTERNAL SMUMPS_497,SMUMPS_442,
     *        SMUMPS_50,
     *        SMUMPS_46
      IF (K48 .eq. 0 .OR. K48.eq.3) THEN
         KMAX = SMUMPS_497( K10, NCB )
         KMIN = SMUMPS_442( K10, K50, KMAX, NCB)
         NSLAVESMAX = SMUMPS_46(
     &                SLAVEF, K48, KMIN, NFRONT, NCB )
      ELSE
         NSLAVESMAX = SLAVEF-1
      ENDIF
      NSLAVESMIN = SMUMPS_50(
     &     SLAVEF, K48, K10, K50, NFRONT, NCB )
      NSLAVESMAX = MAX ( NSLAVESMAX, NSLAVESMIN )
      SMUMPS_52 = 
     &               MIN ( NSLAVESMAX, NCB )
      RETURN
      END FUNCTION SMUMPS_52
      SUBROUTINE SMUMPS_503( WHAT, KEEP,
     *           NCB, NFR, SLAVEF, NBROWMAX, MAXSURFCB )
      IMPLICIT NONE
      INTEGER, intent(in) :: WHAT, NCB, NFR, SLAVEF
      INTEGER, intent(in) :: KEEP(150)
      INTEGER, intent(out) :: NBROWMAX, MAXSURFCB
      INTEGER KMAX, KMIN, NSLAVES, SIZEDUMMY, TABDUMMY
      EXTERNAL SMUMPS_497, SMUMPS_442,
     *         SMUMPS_50
      INTEGER SMUMPS_497, SMUMPS_442,
     *        SMUMPS_50
      IF ( WHAT .NE. 1 .and. WHAT .NE. 2 ) THEN
        WRITE(*,*) "Internal error 1 in SMUMPS_503"
        CALL SMUMPS_ABORT()
      ENDIF
      KMAX    = SMUMPS_497( KEEP(10), NCB )
      NSLAVES = SMUMPS_50( SLAVEF, KEEP(48),
     *          KEEP(10), KEEP(50),
     *          NFR, NCB )
      IF ( KEEP(48) == 0 ) THEN
        NBROWMAX = NCB / NSLAVES + MOD( NCB, NSLAVES )
        IF ( WHAT == 2 ) MAXSURFCB = NBROWMAX * NCB
      ELSE IF ( KEEP(48) == 3 ) THEN
        KMIN = SMUMPS_442( KEEP(10), KEEP(50), KMAX, NCB )  
        SIZEDUMMY        = 1
        CALL  SMUMPS_440(
     *    WHAT, NSLAVES, NFR, NCB, 
     *    KMIN, KMAX, SLAVEF,
     *    NBROWMAX, MAXSURFCB, TABDUMMY, SIZEDUMMY) 
      ELSE IF ( KEEP(48) == 4 ) THEN
         IF (KEEP(10) > 0) THEN
            WRITE(*,*) "Internal error 2 in SMUMPS_503"
            CALL SMUMPS_ABORT()
         END IF
         IF(KEEP(50).EQ.0)THEN
            IF ( ABS(KEEP(10)) * ( SLAVEF - 1 ) > NCB *NFR ) THEN
              NBROWMAX = (NCB + SLAVEF -2 ) / ( SLAVEF - 1 )
              IF ( WHAT == 2 ) MAXSURFCB = NBROWMAX *NFR
            ELSE
              NBROWMAX=INT (( ABS(KEEP(10)) + NFR - 1 )/ NFR)
              IF ( WHAT == 2 ) MAXSURFCB = ABS(KEEP(10))
            ENDIF
         ELSE
            NBROWMAX=INT((-(NFR-NCB)
     $              +SQRT(REAL(((NFR-NCB)*
     $              (NFR-NCB))+REAL(4)*
     *              (ABS(KEEP(10))))))/
     $              REAL(2))
            IF ( WHAT == 2 ) MAXSURFCB = ABS(KEEP(10))
         ENDIF
      ELSE
        NBROWMAX = NCB
        IF (WHAT == 2) MAXSURFCB = NCB * NCB
      ENDIF 
      NBROWMAX = MIN ( MAX(NBROWMAX, 1), NCB)
      RETURN
      END SUBROUTINE SMUMPS_503
      INTEGER FUNCTION SMUMPS_46( SLAVEF, K48, 
     *         BLSIZE, NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, BLSIZE, NFRONT, NCB
      INTEGER NSLAVES, NASS
      REAL Wtotal, Wblsize
      REAL SMUMPS_45
      EXTERNAL          SMUMPS_45
      NASS = NFRONT - NCB
      NSLAVES  = SLAVEF-1
      IF ( K48 .EQ.0 ) THEN 
         NSLAVES = MAX(NCB/MAX(1,BLSIZE),1)
      ELSE IF (K48 .EQ. 3 ) THEN
        Wblsize = SMUMPS_45(BLSIZE,NFRONT,NASS)
        Wtotal  = SMUMPS_45(NCB,NFRONT,NASS)
        NSLAVES = MAX(NINT ( Wtotal / Wblsize ), 1)
      ENDIF
      SMUMPS_46 = 
     *               MIN ( NSLAVES,(SLAVEF-1) )
      RETURN
      END FUNCTION SMUMPS_46
      INTEGER FUNCTION SMUMPS_439( SLAVEF, K48, K10, 
     *         K50, NFRONT, NCB, NSLAVES)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K10, K50, NFRONT, 
     *                        NCB, NSLAVES
      INTEGER NBROWMAX, KMIN, KMAX, 
     *        GETPOSITIONS, TABDUMMY(1), SIZEDUMMY, 
     *        MAXSURFDUMMY
      INTEGER SMUMPS_442, SMUMPS_497
      EXTERNAL SMUMPS_442, SMUMPS_497, 
     *        SMUMPS_440
      IF ( K48 .EQ.0 ) THEN 
       NBROWMAX =  NCB/NSLAVES + MOD(NCB,NSLAVES)
      ELSE IF (K48 .EQ. 3 ) THEN
        KMAX = SMUMPS_497(K10, NCB)
        KMIN = SMUMPS_442(K10, K50, KMAX, NCB)
        GETPOSITIONS = 1
        SIZEDUMMY        = 1
        CALL  SMUMPS_440(
     *    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     *    KMIN, KMAX, SLAVEF,
     *    NBROWMAX, MAXSURFDUMMY, TABDUMMY, SIZEDUMMY) 
      ELSE IF(K48.EQ.4)THEN
         NBROWMAX =NCB
      ELSE
        NBROWMAX = NCB
      ENDIF
      SMUMPS_439 = MIN(MAX(NBROWMAX,1), NCB)
      RETURN
      END FUNCTION SMUMPS_439
      SUBROUTINE  SMUMPS_440( 
     *    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     *    KMIN, KMAX, SLAVEF,
     *    NBROWMAX, MAXSURFCB, TABPOS, SIZETABPOS) 
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: GETPOSITIONS, 
     *    NSLAVES, NFRONT, NCB, 
     *    KMIN, KMAX, SLAVEF, SIZETABPOS 
      INTEGER, INTENT (OUT) :: NBROWMAX, MAXSURFCB
      INTEGER, INTENT (OUT) :: TABPOS(SIZETABPOS) 
      REAL Wtotal, W, COSTni
      REAL delta
      INTEGER  SumNi, NCOLim1, I, BLSIZE, NASS
      LOGICAL GETROW, GETSURF, GETPOS
      REAL SMUMPS_45
      EXTERNAL          SMUMPS_45
      GETROW =((GETPOSITIONS.EQ.1).OR.(GETPOSITIONS.EQ.2))
      GETSURF= (GETPOSITIONS.EQ.2)
      GETPOS = (GETPOSITIONS.GE.3)
      NBROWMAX  = 0
      MAXSURFCB = 0
      IF (GETPOS) THEN 
        TABPOS (1) = 1
        TABPOS (NSLAVES+1)= NCB+1
        TABPOS (SLAVEF+2) = NSLAVES
      ENDIF
      IF (NSLAVES.EQ.1) THEN
       IF ( GETSURF ) THEN
         NBROWMAX  = NCB
         MAXSURFCB = NCB * NCB
       ELSEIF ( GETROW ) THEN
         NBROWMAX  = NCB
       ENDIF
      ELSE
        NASS    = NFRONT - NCB
        W       = SMUMPS_45(NCB,NFRONT,NASS)
        SumNi   = 0
        NCOLim1 = NASS
        DO I = 1, NSLAVES-1
          delta   = REAL(2*NCOLim1-NASS)**2 + 
     &                  (REAL(4)*W)/REAL(NASS*(NSLAVES-I+1))
          delta   = SQRT(delta)
          delta   = (REAL(-2*NCOLim1+NASS) + delta )/REAL(2)
          BLSIZE  = MAX(INT(delta), 1)
          IF ( (NFRONT-NCOLim1-BLSIZE) .LE. NSLAVES-I ) THEN
            BLSIZE = 1
          ENDIF
          NCOLim1 = NCOLim1+BLSIZE
          COSTni  = SMUMPS_45(BLSIZE,NCOLim1,NASS)
          W       = W - COSTni
          IF (GETPOS) TABPOS(I) = SumNi + 1
          IF (GETSURF) THEN
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
            MAXSURFCB = MAX ( MAXSURFCB, 
     *       BLSIZE* (SumNi+BLSIZE) )
          ELSEIF ( GETROW ) THEN         
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
             RETURN
          ENDIF
          SumNi   = SumNi + BLSIZE
        ENDDO
        BLSIZE = NCB - SumNi
        IF (BLSIZE.LE.0) THEN
          write(*,*) ' Error in SMUMPS_BLOC2_init: ', 
     &     ' size lastbloc ', BLSIZE
          CALL SMUMPS_ABORT()
        ENDIF
        if (NCOLim1+BLSIZE.NE.NFRONT) then
          write(*,*) ' Error in SMUMPS_BLOC2_init: ', 
     &     ' NCOLim1, BLSIZE, NFRONT=', 
     *       NCOLim1, BLSIZE, NFRONT
          CALL SMUMPS_ABORT()
        endif
        IF (GETPOS) TABPOS(NSLAVES) = SumNi + 1
        IF (GETSURF) THEN
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
            MAXSURFCB = MAX ( MAXSURFCB, 
     *       BLSIZE* (SumNi+BLSIZE) )
        ELSEIF ( GETROW ) THEN         
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_440
      SUBROUTINE SMUMPS_441(      
     *            KEEP, SLAVEF,
     *            TAB_POS_IN_PERE,
     *            NSLAVES, NFRONT, NCB
     *             )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NSLAVES, SLAVEF, NFRONT,
     *                         KEEP(150)
      INTEGER TAB_POS_IN_PERE(SLAVEF+2)
      INTEGER :: I, BLSIZE
      INTEGER KMIN, KMAX, NBROWDUMMY, MAXSURFDUMMY,
     *        GETPOSITIONS, SIZECOLTAB
      INTEGER SMUMPS_442, SMUMPS_497 
      EXTERNAL SMUMPS_442, SMUMPS_497, 
     *        SMUMPS_440
       IF (KEEP(48).EQ.0) THEN
        BLSIZE = NCB / NSLAVES
        TAB_POS_IN_PERE( 1 ) = 1
        DO I = 1, NSLAVES-1
          TAB_POS_IN_PERE( I+1 ) = TAB_POS_IN_PERE(I) +
     *    BLSIZE
        ENDDO
        TAB_POS_IN_PERE(NSLAVES+1) = NCB+1
        TAB_POS_IN_PERE(SLAVEF+2)  = NSLAVES
        RETURN
      ELSE IF (KEEP(48).EQ.3 ) THEN
        KMAX = SMUMPS_497(KEEP(10), NCB)
        KMIN = SMUMPS_442(KEEP(10), KEEP(50), KMAX, NCB)
        GETPOSITIONS = 3
        SIZECOLTAB       = SLAVEF+2
        CALL  SMUMPS_440(
     *    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     *    KMIN, KMAX, SLAVEF,
     *    NBROWDUMMY, MAXSURFDUMMY, 
     *    TAB_POS_IN_PERE(1), SIZECOLTAB)
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_441
      SUBROUTINE SMUMPS_49(
     *            KEEP, INODE, STEP, N, SLAVEF,
     *            ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *            ISLAVE, NCB, NSLAVES, SIZE, FIRST_INDEX )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: ISLAVE, NCB, NSLAVES, SLAVEF,
     *                         KEEP(150), INODE, N 
      INTEGER, INTENT( IN ) :: STEP(N), 
     *          ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER, INTENT( OUT ):: SIZE, FIRST_INDEX
      INTEGER BLSIZE, J
      IF (KEEP(48).EQ.0) THEN
       BLSIZE = NCB / NSLAVES
       IF ( ISLAVE .NE. NSLAVES ) THEN
        SIZE = BLSIZE
       ELSE
        SIZE = BLSIZE + MOD( NCB, NSLAVES )
       END IF
       FIRST_INDEX = ( ISLAVE - 1 ) * BLSIZE + 1
      ELSEIF (KEEP(48).EQ.3) THEN
       J = ISTEP_TO_INIV2 ( STEP(INODE) )
       FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
       SIZE        = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSEIF (KEEP(48).EQ.4) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
         SIZE = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSE
       WRITE(*,*) 'Error in SMUMPS_BLOC2 undef strat'  
       CALL SMUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_49
      REAL FUNCTION SMUMPS_45(NROW,NCOL,NASS)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: NROW,NCOL,NASS
      SMUMPS_45 = REAL(NASS*NROW)*
     &                 REAL(2*NCOL - NASS - NROW + 1)
      RETURN
      END FUNCTION SMUMPS_45
      INTEGER FUNCTION SMUMPS_12 
     &      (K10, K48, K50, SLAVEF, 
     &      NCB, NFRONT, NSLAVES_less, NMB_OF_CAND )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NFRONT, NSLAVES_less, 
     *                  K10, K48, K50, SLAVEF, NMB_OF_CAND
      INTEGER J, NSLAVES, NSLAVES_REAL, ISLAVE
      INTEGER KMIN, KMAX, BLSIZE, NPIV,
     *        NSLAVES_ref, KMINSURF, NSLAVES_max
      REAL WK_MASTER, WK_SLAVE
      INTEGER  SMUMPS_497, SMUMPS_50, 
     *         SMUMPS_52
      REAL  SMUMPS_45
      EXTERNAL SMUMPS_497, SMUMPS_50, 
     *         SMUMPS_52
      EXTERNAL SMUMPS_45
      IF (NMB_OF_CAND.LE.0) THEN
      ENDIF
      IF ( (K48.EQ.0).OR. (K48.EQ.3) ) THEN
         KMAX = SMUMPS_497( K10, NCB )
         NSLAVES_ref = SMUMPS_50(
     &     SLAVEF, K48, K10, K50, NFRONT, NCB )
         NSLAVES = NSLAVES_ref
         IF ( NSLAVES_ref.LT.SLAVEF ) THEN
           NSLAVES_max = SMUMPS_52(
     *       SLAVEF, K48, K10, K50, NFRONT, NCB )
           IF ( NSLAVES_max .LT. NSLAVES_less ) THEN
            NSLAVES =  NSLAVES_max
           ELSE 
            NSLAVES =  NSLAVES_less
           ENDIF
           NSLAVES = MAX(NSLAVES_ref,NSLAVES)
         ENDIF
         IF ( NSLAVES.GT.NSLAVES_ref) THEN
          NPIV = NFRONT - NCB
          IF ( K50.EQ.0 ) THEN
           WK_SLAVE = REAL( NPIV ) * REAL( NCB ) *
     *         ( 2.0E0 * NFRONT - NPIV ) / NSLAVES
           WK_MASTER = 0.6667E0 * NPIV*NPIV*NPIV+
     *                NPIV*NPIV*NCB
          ELSE
           WK_SLAVE = SMUMPS_45(NCB,NFRONT,NPIV)
     *         / NSLAVES
           WK_MASTER =  REAL(NPIV*NPIV)*REAL(NPIV)/(3.0)
          ENDIF
          IF ( (WK_MASTER.GT.WK_SLAVE).AND.
     *     (WK_SLAVE.GT.1.0) ) THEN
           NSLAVES = 
     *         INT( REAL(NSLAVES) * (WK_SLAVE/WK_MASTER))
           NSLAVES = MAX(NSLAVES_ref, NSLAVES)
          ENDIF
         ENDIF
      ELSE
       NSLAVES = NSLAVES_LESS
      ENDIF
      NSLAVES = MIN (NSLAVES, NCB)
      NSLAVES = MIN (NSLAVES, NMB_OF_CAND)
      SMUMPS_12 = NSLAVES
      RETURN
      END FUNCTION SMUMPS_12
      SUBROUTINE SMUMPS_47( 
     *   KEEP, INODE, STEP, N, SLAVEF, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *   NASS, NCB, 
     *   NSLAVES, POSITION, ISLAVE, IPOSSLAVE )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: KEEP(150),INODE,N,SLAVEF 
      INTEGER, INTENT( IN ) :: STEP(N), 
     *          ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER, INTENT( IN  ) :: NASS, NCB, 
     &                          NSLAVES, POSITION
      INTEGER, INTENT( OUT ) :: ISLAVE, IPOSSLAVE
      INTEGER BLSIZE, J, ISHIFT
      IF ((NSLAVES.LE.0).OR.(POSITION.LE.NASS)) THEN
       ISLAVE = 0
       IPOSSLAVE = POSITION
       RETURN
      ENDIF
      IF ( KEEP(48).EQ.0) THEN
       BLSIZE = NCB / NSLAVES
       ISLAVE    = MIN( NSLAVES,
     *               ( POSITION - NASS - 1 ) / BLSIZE + 1 )
       IPOSSLAVE = POSITION - NASS - ( ISLAVE - 1 ) * BLSIZE
      ELSEIF (KEEP(48).EQ.3) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSEIF (KEEP(48).EQ.4) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSE
       WRITE(*,*) 'Error in SMUMPS_47: undef strat'  
       CALL SMUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_47
      INTEGER FUNCTION SMUMPS_442( K10, K50, KMAX, NCB )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: K10, KMAX, NCB, K50
      INTEGER KMIN, MINGRAN, KMINSURF, MINSURF
      IF ( ( NCB .LE.0 ).OR. (KMAX.LE.0) ) THEN 
        SMUMPS_442 = 1
        RETURN
      ENDIF
      IF (K50.EQ.0) THEN
       KMINSURF = 60000
#if defined(t3e) || defined(sgi)
       MINGRAN = 40 
#else
       MINGRAN = 50 
#endif
      ELSE
       KMINSURF = 30000
#if defined(t3e) || defined(sgi)
       MINGRAN = 10 
#else
       MINGRAN = 20 
#endif
      ENDIF
      IF (K10.GT.0) THEN
#if defined(t3e) || defined(sgi)
           KMIN = MAX(MINGRAN,KMAX/10)
#else
           KMIN = MAX(MINGRAN,KMAX/20)
#endif
      ELSE
           KMINSURF = MAX ( 
     &                INT(IABS(K10)/500), KMINSURF)
           KMIN     = MAX(INT(KMINSURF/MAX(NCB,1)), 1)
      ENDIF
      KMIN = min(KMIN,KMAX)
      KMIN = max(KMIN,1)  
      SMUMPS_442 = KMIN
      RETURN
      END FUNCTION SMUMPS_442
      INTEGER FUNCTION SMUMPS_497( KEEP10, NCB )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: KEEP10, NCB
      INTEGER KMAX 
      IF ( NCB .LE.0 ) THEN 
        SMUMPS_497 = 1
        RETURN
      ENDIF
      IF ( KEEP10.GT.0 ) THEN
       KMAX = KEEP10
      ELSE
       KMAX =  -KEEP10/NCB 
      ENDIF
      KMAX = MIN (NCB, KMAX)
      SMUMPS_497 = MAX ( KMAX, 1 )
      RETURN
      END FUNCTION SMUMPS_497
