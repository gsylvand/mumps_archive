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
      SUBROUTINE SMUMPS_128(N, NELT, ELTPTR, ELTVAR, LIW,
     *            IKEEP, PTRAR,
     *            IORD, NFSIZ, FILS, FRERE, 
     *            LISTVAR_SCHUR, SIZE_SCHUR,
     *            ICNTL, INFO, KEEP,
     *            ELTNOD, NSLAVES, 
     *            XNODEL, NODEL)
      IMPLICIT NONE
      INTEGER N,NELT,LIW,IORD, SIZE_SCHUR, NSLAVES
      INTEGER PTRAR(N,3), NFSIZ(N), FILS(N), FRERE(N)
      INTEGER ELTPTR(NELT+1) 
      INTEGER XNODEL(N+1), NODEL(ELTPTR(NELT+1)-1)
      INTEGER ELTVAR(ELTPTR(NELT+1)-1)
      INTEGER IKEEP(N,3)
      INTEGER LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER INFO(40), ICNTL(40), KEEP(150)
      INTEGER ELTNOD(NELT)
      INTEGER K,I,L1,L2,IWFR,NCMPA,LLIW,IFSON,IN
      INTEGER NEMIN, MPRINT, LP, MP, LDIAG
      INTEGER NZ, allocok, ITEMP
      LOGICAL PROK, NOSUPERVAR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW2
      INTEGER IDUM 
      INTEGER IDUMARR(1)
      INTEGER IDUMARR4(1)
      EXTERNAL SMUMPS_197, SMUMPS_130, SMUMPS_131,
     *         SMUMPS_129, SMUMPS_132, 
     *         SMUMPS_133, SMUMPS_134,
     *         SMUMPS_199,
     *         SMUMPS_200, SMUMPS_201
        ALLOCATE( IW ( LIW ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO( 1 ) = -13
          INFO( 2 ) = LIW
          RETURN
        ENDIF
      MPRINT= ICNTL(3)
      PROK  = (MPRINT.GT.0)
      LP    = ICNTL(1)
      MP    = ICNTL(3)
      LDIAG = ICNTL(4)
      IF (KEEP(60).NE.0) THEN
       NOSUPERVAR=.TRUE.
      ELSE
       NOSUPERVAR=.FALSE.
      ENDIF
      IF (KEEP(1).LT.1) KEEP(1) = 1
      NEMIN = KEEP(1)
      IF (LDIAG.LE.2 .OR. MP.LE.0) GO TO 10
      WRITE (MP,99999) N, NELT, LIW, INFO(1)
      K = MIN0(10,NELT+1)
      IF (LDIAG.EQ.4) K = NELT+1
      IF (K.GT.0) WRITE (MP,99998) (ELTPTR(I),I=1,K)
      K = MIN0(10,ELTPTR(NELT+1)-1)
      IF (LDIAG.EQ.4) K = ELTPTR(NELT+1)-1
      IF (K.GT.0) WRITE (MP,99995) (ELTVAR(I),I=1,K)
      K = MIN0(10,N)
      IF (LDIAG.EQ.4) K = N
      IF (IORD.EQ.1 .AND. K.GT.0) THEN
        WRITE (MP,99997) (IKEEP(I,1),I=1,K)
      ENDIF
   10 L1 = 1
      L2 = L1 + N
      IF (NOSUPERVAR) THEN
        IF ( LIW .LT. 2*N ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        END IF
      ELSE
        IF ( LIW .LT.  4*N-4 ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        END IF
      ENDIF
      CALL SMUMPS_258(NELT, N, ELTPTR(NELT+1)-1, ELTPTR, ELTVAR,
     *           XNODEL, NODEL, IW(L1), IDUM, ICNTL)
      IF (IORD.NE.1) THEN
        IF (NOSUPERVAR) THEN
          CALL SMUMPS_129(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              PTRAR(1,2), IW(L1))
        ELSE
         CALL SMUMPS_130(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              PTRAR(1,2), 4*N+4, IW(L1))
        ENDIF
        LLIW = NZ
        ALLOCATE( IW2(LLIW), stat = allocok )
        IF (allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = LLIW
          RETURN
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL SMUMPS_132(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              IW2, LLIW, PTRAR, PTRAR(1,2),
     *              IW(L1), IWFR)
        ELSE
         CALL SMUMPS_131(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              IW2, LLIW, PTRAR, PTRAR(1,2),
     *              IW(L1), IWFR)
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL SMUMPS_162(N, LLIW, PTRAR, IWFR, PTRAR(1,2), IW2,
     *   IW(L1), IKEEP,
     *   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), PTRAR(1,3),
     *   LISTVAR_SCHUR, SIZE_SCHUR)
         KEEP(20) = LISTVAR_SCHUR(1)
        ELSE
         CALL SMUMPS_23(N, LLIW, PTRAR, IWFR, PTRAR(1,2), IW2, 
     *   IW(L1), IKEEP, 
     *   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), PTRAR(1,3))
        ENDIF
      ELSE
       DO K=1,N
         IW(L1+K) = 0
       ENDDO
       DO K=1,N
         IF ((IKEEP(K,1).LE.0).OR.(IKEEP(K,1).GT.N)) 
     *    GO TO 40
         IF (IW(L1+IKEEP(K,1)).EQ.1) THEN
          GOTO 40
         ELSE
          IW(L1+IKEEP(K,1)) = 1
         ENDIF
       ENDDO
       CALL SMUMPS_133(N, NZ, NELT, ELTPTR(NELT+1)-1,
     *             ELTPTR, ELTVAR, XNODEL, NODEL, 
     *             IKEEP, PTRAR(1,2), IW(L1))
       LLIW = NZ+N
       ALLOCATE( IW2(LLIW), stat = allocok )
       IF (allocok.GT.0) THEN
         INFO(1) = -13
         INFO(2) = LLIW
         RETURN
       ENDIF
       CALL SMUMPS_134(N, NZ, NELT, ELTPTR(NELT+1)-1,
     *             ELTPTR, ELTVAR, XNODEL, NODEL, 
     *             IKEEP, IW2, LLIW, PTRAR, PTRAR(1,2),
     *             IW(L1), IWFR)
       ITEMP = 0
       CALL SMUMPS_199(N, PTRAR, IW2, LLIW, IWFR, IKEEP, 
     *    IKEEP(1,2), IW(L1),
     *    IW(L2), NCMPA, ITEMP)
      ENDIF
      DEALLOCATE(IW2)
      CALL SMUMPS_200(N, PTRAR, IW(L1), IKEEP, IKEEP(1,2),
     *     IKEEP(1,3),
     *     NFSIZ, INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN, KEEP(60))
      IF (KEEP(60).NE.0) THEN
         IN = KEEP(20)
         DO WHILE (IN.GT.0)
          IN = FILS (IN)
         END DO
         IFSON = -IN
         IN = KEEP(20)
         DO I=2,SIZE_SCHUR
          FILS(IN) = LISTVAR_SCHUR (I)
          IN       = FILS(IN)
          FRERE (IN) = N+1
         ENDDO
         FILS(IN) = -IFSON
      ENDIF
      CALL SMUMPS_201(IKEEP(1,2),
     *  PTRAR(1,3), INFO(6),
     *  INFO(5), KEEP(2),KEEP(50),INFO(20))
      IF ( KEEP(53) .NE. 0 ) THEN
        CALL SMUMPS_209( N, FRERE, FILS, NFSIZ )
      END IF
      IF ( KEEP(48) == 4 .OR.
     *   ( (KEEP(24).NE.0).AND.(KEEP(10).GT.0) ) ) THEN
          CALL SMUMPS_510(KEEP(10), KEEP(2),
     *    KEEP(50), NSLAVES)
      END IF
      IF (.TRUE.) THEN
       IF ( KEEP(62) .GE. 1 )
     *   CALL SMUMPS_97(N, FRERE, FILS, NFSIZ , INFO(6),
     *                        NSLAVES, KEEP )
       CALL SMUMPS_202(N, IDUM, IKEEP(1,1), FILS, FRERE,
     *      NFSIZ, IKEEP(1,3), IKEEP(1,2),
     *      IDUMARR4, IDUMARR4, IDUMARR, IDUMARR, PTRAR,
     *      ICNTL, INFO, KEEP)
      ELSE
      CALL SMUMPS_348(N, IDUM, IKEEP(1,1), FILS, FRERE,
     *      NFSIZ, IKEEP(1,3), IKEEP(1,2),
     *      IDUMARR4, IDUMARR4, PTRAR, ICNTL, INFO, KEEP)
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
       K = MIN0(10,N)
       IF (LDIAG.EQ.4) K = N
       IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
       IF (K.GT.0) WRITE (MP,99991) (IKEEP(I,2),I=1,K)
       IF (K.GT.0) WRITE (MP,99990) (IKEEP(I,3),I=1,K)
       IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
   40 INFO(1) = -4
      INFO(2) = K
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996) INFO(1)
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
   90 CONTINUE
      DEALLOCATE(IW)
      RETURN
99999 FORMAT (/'Entering analysis phase with ...'/
     * '                N         NELT       LIW       INFO(1)'/,
     * 9X, I8, I11, I12, I14)
99998 FORMAT ('Element pointers:  ELTPTR()   '/(9X, 7I10))
99995 FORMAT ('Element variables: ELTVAR()   '/(9X, 7I10))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT (/'** Error return ** from Analysis   *  INFO(1)=', I3)
99991 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99990 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99989 FORMAT ('FILS (.)  =', 10I6/(12X, 10I6))
99988 FORMAT ('FRERE(.)  =', 10I6/(12X, 10I6))
99987 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99982 FORMAT ('Error in permutation array KEEP   INFO(2)=', I3)
      END SUBROUTINE SMUMPS_128
      SUBROUTINE SMUMPS_258( NELT, N, NELNOD, XELNOD, ELNOD,
     $                        XNODEL, NODEL, FLAG, IERROR, ICNTL ) 
      IMPLICIT NONE
      INTEGER NELT, N, NELNOD, IERROR, ICNTL(40)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     $        FLAG(N)
      INTEGER I, J, K, MP, NBERR
      MP = ICNTL(2)
      FLAG(1:N) = 0
      XNODEL(1:N) = 0
      IERROR = 0
      DO I = 1, NELT
        DO K = XELNOD(I), XELNOD(I+1)-1
          J = ELNOD(K)
          IF ( J.LT.1 .OR. J.GT.N ) THEN
            IERROR = IERROR + 1
          ELSE
            IF ( FLAG(J).NE.I ) THEN
              XNODEL(J) = XNODEL(J) + 1
              FLAG(J) = I
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF ( IERROR.GT.0 .AND. MP.GT.0 .AND. ICNTL(4).GE.2 ) THEN
        NBERR = 0
        WRITE(MP,99999)
        DO I = 1, NELT
          DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF ( J.LT.1 .OR. J.GT.N ) THEN
              NBERR = NBERR + 1
              IF (NBERR.LE.10) THEN
                WRITE(MP,'(A,I8,A,I8,A)')
     *          'Element ',I,' variable ',J,' ignored.' 
              ELSE
                GO TO 100
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
  100 CONTINUE
      K = 1
      DO I = 1, N
         K = K + XNODEL(I)
         XNODEL(I) = K
      ENDDO
      XNODEL(N+1) = XNODEL(N)
      FLAG(1:N) = 0
      DO I = 1, NELT
         DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF (FLAG(J).NE.I) THEN   
              XNODEL(J) = XNODEL(J) - 1
              NODEL(XNODEL(J)) = I
              FLAG(J) = I
            ENDIF
         ENDDO
      ENDDO
      RETURN
99999 FORMAT (/'*** Warning message from subroutine SMUMPS_258 ***')
      END SUBROUTINE SMUMPS_258
      SUBROUTINE SMUMPS_129(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  LEN, FLAG)
      IMPLICIT NONE
      INTEGER N, NELT, NELNOD, NZ
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER  XNODEL(N+1), NODEL(NELNOD),
     $        FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN 
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                LEN(I) = LEN(I) + 1
                LEN(J) = LEN(J) + 1
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0
      DO I = 1,N
        NZ = NZ + LEN(I)
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_129
      SUBROUTINE SMUMPS_132(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER LEN(N)
      INTEGER IPE(N)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER  XNODEL(N+1), NODEL(NELNOD), 
     *          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IWFR = IWFR + LEN(I)
        IF (LEN(I).GT.0) THEN
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                IPE(I) = IPE(I) - 1
                IW(IPE(I)) = J
                IPE(J) = IPE(J) - 1
                IW(IPE(J)) = I
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_132
      SUBROUTINE SMUMPS_133(N, NZ, NELT, NELNOD,
     * XELNOD, ELNOD, XNODEL, NODEL, 
     * PERM, LEN, FLAG)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  LEN(I) = LEN(I) + 1
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0  
      DO I = 1,N 
        NZ = NZ + LEN(I)
      ENDDO   
      RETURN  
      END SUBROUTINE SMUMPS_133
      SUBROUTINE SMUMPS_134(N, NZ, NELT, NELNOD,
     * XELNOD, ELNOD, XNODEL, NODEL, 
     * PERM, IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER  XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER IPE(N), LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), IW(LW), 
     *          FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 0
      DO I = 1,N
        IWFR = IWFR + LEN(I) + 1
        IPE(I) = IWFR 
      ENDDO
      IWFR = IWFR + 1
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  IW(IPE(I)) = J
                  IPE(I) = IPE(I) - 1
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DO I = 1,N
        J = IPE(I)
        IW(J) = LEN(I)
        IF (LEN(I).EQ.0) IPE(I) = 0
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_134
      SUBROUTINE SMUMPS_25( MYID, SLAVEF, N, 
     *           PROCNODE, STEP, PTRAIW, PTRARW, 
     *           NELT, FRTPTR, FRTELT, 
     *           KEEP, ICNTL, SYM )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, N, NELT, SYM
      INTEGER KEEP( 150 ), ICNTL( 40 )
      INTEGER PTRAIW( NELT+1 ), PTRARW( NELT+1 )
      INTEGER STEP( N )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PROCNODE( KEEP(28) )
      INTEGER MASTER
      PARAMETER (MASTER=0)
      INTEGER SMUMPS_330, SMUMPS_275
      EXTERNAL SMUMPS_330, SMUMPS_275
      INTEGER ELT, I, K, IPTRI, IPTRR, NVAR
      INTEGER TYPE_PARALL, ITYPE, IRANK
      TYPE_PARALL = KEEP(46)
      PTRAIW( 1:NELT ) = 0
      DO I = 1, N
        IF (STEP(I).LT.0) CYCLE
        ITYPE = SMUMPS_330( ABS(STEP(I)), PROCNODE, SLAVEF )
        IRANK = SMUMPS_275( ABS(STEP(I)), PROCNODE, SLAVEF )
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK = IRANK + 1
        END IF
        IF ( (ITYPE .EQ. 2) .OR.
     *       (ITYPE .EQ. 1 .AND. IRANK .EQ. MYID) ) THEN
          DO K = FRTPTR(I),FRTPTR(I+1)-1
            ELT = FRTELT(K)
            PTRAIW( ELT ) = PTRARW(ELT+1) - PTRARW(ELT)
          ENDDO
        ELSE 
        END IF
      END DO
      IPTRI = 1
      DO ELT = 1,NELT
        NVAR = PTRAIW( ELT )
        PTRAIW( ELT ) = IPTRI
        IPTRI = IPTRI + NVAR
      ENDDO
      PTRAIW( NELT+1 ) = IPTRI
      KEEP( 14 ) = IPTRI - 1
      IF ( .TRUE. ) THEN  
        IF (SYM .EQ. 0) THEN
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + NVAR*NVAR
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ELSE
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + (NVAR*(NVAR+1))/2
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ENDIF
      ELSE
        IF (SYM .EQ. 0) THEN
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + NVAR*NVAR 
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ELSE
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + (NVAR*(NVAR+1))/2  
          ENDDO 
          PTRARW( NELT+1 ) = IPTRR
        ENDIF
      ENDIF 
      KEEP( 13 ) = IPTRR - 1
      RETURN
      END SUBROUTINE SMUMPS_25
      SUBROUTINE SMUMPS_120( N, NELT, ELTPROC, SLAVEF, PROCNODE )
      IMPLICIT NONE
      INTEGER N, NELT, SLAVEF
      INTEGER PROCNODE( N ), ELTPROC( NELT )
      INTEGER ELT, I, ITYPE, SMUMPS_330, SMUMPS_275
      EXTERNAL SMUMPS_330, SMUMPS_275
      DO ELT = 1, NELT
          I = ELTPROC(ELT)
          ITYPE = SMUMPS_330(I,PROCNODE,SLAVEF)
          IF (ITYPE.EQ.1) THEN
            ELTPROC(ELT) = SMUMPS_275(I,PROCNODE,SLAVEF)
          ELSE IF (ITYPE.EQ.2) THEN
            ELTPROC(ELT) = -1
          ELSE
            ELTPROC(ELT) = -2
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_120
      SUBROUTINE SMUMPS_153(N, NELT, NELNOD, FRERE, FILS, NA, NE,
     *           XNODEL, NODEL, FRTPTR, FRTELT, ELTNOD) 
      IMPLICIT NONE
      INTEGER  N, NELT, NELNOD
      INTEGER  FRERE(N), FILS(N), NA(N), NE(N)
      INTEGER  FRTPTR(N+1), FRTELT(NELT), ELTNOD(NELT)
      INTEGER  XNODEL(N+1), NODEL(NELNOD) 
      INTEGER TNSTK( N ), IPOOL( N )
      INTEGER I, K, IFATH
      INTEGER INODE, LEAF, NBLEAF, NBROOT, III, IN
      TNSTK = NE
      LEAF = 1
      IF (N.EQ.1) THEN
        NBROOT = 1
        NBLEAF = 1
        IPOOL(1) = 1
        LEAF = LEAF + 1
      ELSEIF (NA(N).LT.0) THEN
        NBLEAF = N
        NBROOT = N
        DO 20 I=1,NBLEAF-1
           INODE = NA(I)
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
 20     CONTINUE
        INODE = -NA(N)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSEIF (NA(N-1).LT.0) THEN
        NBLEAF = N-1
        NBROOT = NA(N)
        IF (NBLEAF-1.GT.0) THEN
         DO 30 I=1,NBLEAF-1
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 30      CONTINUE
        ENDIF
        INODE = -NA(N-1)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSE
        NBLEAF = NA(N-1)
        NBROOT = NA(N)
        DO 40 I = 1,NBLEAF
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 40     CONTINUE
      ENDIF
      ELTNOD(1:NELT) = 0
      III = 1
 90   CONTINUE
        IF (III.NE.LEAF) THEN
           INODE=IPOOL(III)
           III = III + 1
        ELSE 
           WRITE(6,*) ' ERROR 1 in file SMUMPS_153 '
           CALL SMUMPS_ABORT()
        ENDIF
 95     CONTINUE 
        IN = INODE
 100    CONTINUE
        DO K = XNODEL(IN),XNODEL(IN+1)-1
          I = NODEL(K)
          IF (ELTNOD(I).EQ.0) ELTNOD(I) = INODE
        ENDDO
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
        IN = INODE
 110    IN = FRERE(IN)
        IF (IN.GT.0) GO TO 110
        IF (IN.EQ.0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
        ELSE
         IFATH = -IN
        ENDIF
        TNSTK(IFATH) = TNSTK(IFATH) - 1
        IF ( TNSTK(IFATH) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
        ELSE
            GOTO 90
        ENDIF
  115 CONTINUE
      FRTPTR(1:N) = 0
      DO I = 1,NELT
        FRTPTR(ELTNOD(I)) = FRTPTR(ELTNOD(I)) + 1
      ENDDO
      K = 1
      DO I = 1,N
        K = K + FRTPTR(I)
        FRTPTR(I) = K
      ENDDO
      FRTPTR(N+1) = FRTPTR(N)
      DO K = 1,NELT
        INODE = ELTNOD(K)
        FRTPTR(INODE) = FRTPTR(INODE) - 1
        FRTELT(FRTPTR(INODE)) = K
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_153
      SUBROUTINE SMUMPS_130(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  LEN, LW, IW)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     $        IW(LW)
      INTEGER I,J,K1,K2,K3,LP,NSUP,SUPVAR
      INTEGER INFO44(6)
      EXTERNAL SMUMPS_315
      LP = 6
      CALL SMUMPS_315(N,NELT,XELNOD(NELT+1)-1,ELNOD,XELNOD,
     *           NSUP,IW(3*N+3+1),3*N+3,IW,LP,INFO44)
      IF (INFO44(1) .LT. 0) THEN
        IF (LP.GE.0) WRITE(LP,*) 
     *     'Error return from SMUMPS_315. INFO(1) = ',INFO44(1) 
      ENDIF
      IW(1:NSUP) = 0
      LEN(1:N) = 0
      DO I = 1,N
        SUPVAR = IW(3*N+3+1+I)
        IF (SUPVAR .EQ. 0) CYCLE
        IF (IW(SUPVAR).NE.0) THEN
          LEN(I) = -IW(SUPVAR)
        ELSE
          IW(SUPVAR) = I
        ENDIF
      ENDDO
      IW(N+1:2*N) = 0
      NZ = 0
      DO SUPVAR = 1,NSUP
        I = IW(SUPVAR)
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J).GE.0) THEN
                IF ((I.NE.J) .AND. (IW(N+J).NE.I)) THEN
                  IW(N+J) = I
                  LEN(I) = LEN(I) + 1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        NZ = NZ + LEN(I)
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_130
      SUBROUTINE SMUMPS_131(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER IPE(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     *          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IF (LEN(I).GT.0) THEN
          IWFR = IWFR + LEN(I)
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        IF (LEN(I).LE.0) CYCLE
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J) .GT. 0) THEN
                IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                  IPE(I) = IPE(I) - 1
                  IW(IPE(I)) = J
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_131
      SUBROUTINE SMUMPS_315(N,NELT,NZ,ELTVAR,ELTPTR,NSUP,SVAR,
     +                 LIW,IW,LP,INFO)
      INTEGER LIW,LP,N,NELT,NSUP,NZ
      INTEGER INFO(6)
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER IW(LIW),SVAR(0:N)
      INTEGER FLAG,NEW,VARS
      EXTERNAL SMUMPS_316
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      INFO(4) = 0
      IF (N.LT.1) GO TO 10
      IF (NELT.LT.1) GO TO 20
      IF (NZ.LT.ELTPTR(NELT+1)-1) GO TO 30
      IF (LIW.LT.6) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      END IF
      NEW = 1
      VARS = NEW + LIW/3
      FLAG = VARS + LIW/3
      CALL SMUMPS_316(N,NELT,ELTPTR,NZ,ELTVAR,SVAR,NSUP,LIW/3-1,
     +           IW(NEW),IW(VARS),IW(FLAG),INFO)
      IF (INFO(1).EQ.-4) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      ELSE
         INFO(4) = 3*NSUP + 3
      END IF
      GO TO 50
   10 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   20 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   30 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   40 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) INFO(4)
      END IF
   50 RETURN
 9000 FORMAT (/3X,'Error message from SMUMPS_315: INFO(1) = ',I2)
 9010 FORMAT (3X,'LIW is insufficient. Upper bound on required work',
     +       'space is ',I8)
      END
      SUBROUTINE SMUMPS_316( N, NELT, ELTPTR, NZ, ELTVAR,
     *           SVAR, NSUP, MAXSUP, NEW, VARS, FLAG, INFO )
      INTEGER MAXSUP,N,NELT,NSUP,NZ
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER INFO(6)
      INTEGER FLAG(0:MAXSUP), NEW(0:MAXSUP),SVAR(0:N),
     *          VARS(0:MAXSUP)
      INTEGER I,IS,J,JS,K,K1,K2
      DO 10 I = 0,N
         SVAR(I) = 0
   10 CONTINUE
      VARS(0) = N + 1
      NEW(0) = -1
      FLAG(0) = 0
      NSUP = 0
      DO 40 J = 1,NELT
         K1 = ELTPTR(J)
         K2 = ELTPTR(J+1) - 1
         DO 20 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) THEN
               INFO(2) = INFO(2) + 1
               GO TO 20
            END IF
            IS = SVAR(I)
            IF (IS.LT.0) THEN
               ELTVAR(K) = 0
               INFO(3) = INFO(3) + 1
               GO TO 20
            END IF
            SVAR(I) = SVAR(I) - N - 2
            VARS(IS) = VARS(IS) - 1
   20    CONTINUE
         DO 30 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) GO TO 30
            IS = SVAR(I) + N + 2
            IF (FLAG(IS).LT.J) THEN
               FLAG(IS) = J
               IF (VARS(IS).GT.0) THEN
                  NSUP = NSUP + 1
                  IF (NSUP.GT.MAXSUP) THEN
                     INFO(1) = -4
                     RETURN
                  END IF
                  VARS(NSUP) = 1
                  FLAG(NSUP) = J
                  NEW(IS) = NSUP
                  SVAR(I) = NSUP
               ELSE
                  VARS(IS) = 1
                  NEW(IS) = IS
                  SVAR(I) = IS
               END IF
            ELSE
               JS = NEW(IS)
               VARS(JS) = VARS(JS) + 1
               SVAR(I) = JS
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SMUMPS_125(
     *           NUMELT, LIST_ELT, 
     *           MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF, PERM, 
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, NELT, 
     *           IW, LIW, 
     *           INTARR, LINTARR, ITLOC, FILS, FRERE_STEPS, 
     *           KEEP28, SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      IMPLICIT NONE
      INTEGER NELT, INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG,
     *        LINTARR, NUMELT, NFRONT_EFF, KEEP28
      INTEGER LIST_ELT(*)
      INTEGER STEP(N), PIMASTER(KEEP28),
     *  PTRAIW(NELT+1), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE_STEPS(KEEP28),
     *        NBPROCFILS(KEEP28), PERM(N)
      INTEGER INTARR(LINTARR)
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K, JDEBROW, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     *        LAST_J_ASS, JMIN, MIN_PERM
      INTEGER ELTI
      LOGICAL LEVEL1
      INTEGER allocok
      INTEGER , ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      IF (NUMSTK.GT.0) THEN
        ALLOCATE(PTTRI(NUMSTK), stat=allocok)
        IF (allocok .GT. 0) THEN
         IFLAG = -13
         GOTO 800
        ENDIF
        ALLOCATE(PTLAST(NUMSTK), stat=allocok)
        IF (allocok .GT. 0) THEN
         IFLAG = -13
         GOTO 800
        ENDIF
      ENDIF
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1
      ICT11  = IOLDP2 + NFRONT
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      DO WHILE (IN.GT.0)
       ITLOC(IN)        = INEW1
       IW(INEW)         = IN
       IW(INEW+NFRONT)  = IN
       INEW1     = INEW1 + 1
       INEW      = INEW + 1
       IN = FILS(IN)
      END DO
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
      MIN_PERM = N + 1
      JMIN     = -1
      DO IELL = 1, NUMSTK 
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
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
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT) THEN 
       DO IELL = 1,NUMELT
        ELTI = LIST_ELT(IELL)
         J1= PTRAIW(ELTI)
         J2= PTRAIW(ELTI+1)-1
         DO JJ=J1,J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
       ENDDO
      END IF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
       IW( IOLDP2+NASS1+1 : IOLDP2+NFRONT_EFF ) = 
     *  IW( ICT11+NASS1+1 : ICT11+NFRONT_EFF ) 
      ELSE
          CALL SMUMPS_308( N, PERM, 
     *            IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
          CALL SMUMPS_309( N, NASS1, PERM, ITLOC,
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
      DO IELL=1,NUMELT
       ELTI = LIST_ELT(IELL)
       J1= PTRAIW(ELTI)
       J2= PTRAIW(ELTI+1)-1
       DO JJ=J1,J2
          J = INTARR(JJ)
          INTARR(JJ) = ITLOC(J)
       END DO
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
      END SUBROUTINE SMUMPS_125
      SUBROUTINE SMUMPS_124(
     *           NUMELT, LIST_ELT,
     *           MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF,
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, NELT,
     *           IW, LIW, 
     *           INTARR, LINTARR, ITLOC, FILS, FRERE_STEPS, 
     *           KEEP28,
     *           SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      IMPLICIT NONE
      INTEGER NELT, INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG,
     *        LINTARR, NUMELT, NFRONT_EFF, KEEP28
      INTEGER LIST_ELT(*)
      INTEGER STEP(N), 
     * PIMASTER(KEEP28),
     *  PTRAIW(NELT+1), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE_STEPS(KEEP28),
     *        NBPROCFILS(KEEP28)
      INTEGER INTARR(LINTARR)
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, INEW, IOLDP2, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, 
     *        I, K, ELTI
      LOGICAL LEVEL1
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1
      NEWEL = IOLDP2 + NASS1
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      DO WHILE (IN.GT.0)
       ITLOC(IN) = INEW1
       IW(INEW)  = IN
       INEW1     = INEW1 + 1
       INEW      = INEW + 1
       IN = FILS(IN)
      END DO
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
          ISON = FRERE_STEPS(STEP(ISON))
  100   CONTINUE
      ENDIF
      DO IELL=1,NUMELT 
       ELTI = LIST_ELT(IELL)
       J1= PTRAIW(ELTI)
       J2= PTRAIW(ELTI+1)-1
       DO JJ=J1,J2
          J = INTARR(JJ)
          IF (ITLOC(J) .EQ. 0) THEN
              NEWEL = NEWEL + 1
              NFRONT_EFF = NFRONT_EFF + 1
              IW(NEWEL)  = J
              INTARR(JJ) = NFRONT_EFF
              ITLOC(J)   = NFRONT_EFF
          ELSE
              INTARR(JJ) = ITLOC(J)
          ENDIF
       END DO
      ENDDO
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
      END SUBROUTINE SMUMPS_124
      SUBROUTINE SMUMPS_36( COMM_LOAD, ASS_IRECV,
     *    NELT, FRT_PTR, FRT_ELT,
     *    N, INODE, IW, LIW, A, LA, IFLAG,
     *    IERROR, ND, 
     *    NIRBDU, FILS, FRERE_STEPS, MAXFRW, root,
     *    OPASSW, OPELIW, PTRIST, PTLUST_S, PTRAST, 
     *    STEP, PIMASTER, PAMASTER,PTRARW, 
     *    PTRAIW, ITLOC, NSTEPS, SON_LEVEL2,
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, 
     *    ICNTL, KEEP,INTARR,DBLARR, 
     *
     *    NSTK_S,NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE SMUMPS_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER IZERO 
      PARAMETER (IZERO=0)
      INTEGER NELT,N,LIW,LA,NIRBDU,NSTEPS
      INTEGER KEEP(150), ICNTL(40)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     *        LRLU, IPTRLU,IWPOS, LRLUS,
     *        POSFAC, IWPOSCB, COMP, IERR_MPI
      INTEGER IDUMMY(1)
      INTEGER IW(LIW), ITLOC(N),
     *        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)), PERM(N), 
     *        FILS(N), FRERE_STEPS(KEEP(28)),
     *        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28)),
     *        STEP(N), 
     * PIMASTER(KEEP(28)),
     * PAMASTER(KEEP(28))
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      LOGICAL SON_LEVEL2
      REAL A(LA)
      REAL  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER IPOOL( LPOOL )
      INTEGER NBPROCFILS(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INTEGER SMUMPS_275, SMUMPS_330
      EXTERNAL SMUMPS_275, SMUMPS_330
      INTEGER LP, HS, HF
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ,LAELL
      INTEGER SIZFI, SIZFR, NCB
      INTEGER LAPOS2,J1,J2
      INTEGER NCOL, NROW, NCOLS, NROWS, LDA_SON
      INTEGER NELIM,JJ,JJ1,JJ2,J3,
     *        IORG
      INTEGER IACHK,JPOS,ICT11
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4,
     *        NUMELT, ELBEG
      INTEGER APOS, APOS2, AINPUT, POSELT, POSEL1, ICT12,
     *        AII, ICT21, J
      INTEGER NSLAVES, NSLSON, NPIVS
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, K, I
      LOGICAL FLAG, BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER NCBSON
      LOGICAL FREE, SAME_PROC
      INTRINSIC REAL
      REAL ZERO
      DATA ZERO /0.0E0/
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      LEVEL = SMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF (LEVEL.NE.1) THEN 
       write(6,*) 'Error1 in mpi51f_niv1 '
       CALL SMUMPS_ABORT()
      END IF
      NSLAVES = 0
      HF = 6 + NSLAVES
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .ne. 0 ) THEN
        ELBEG  = FRT_PTR(INODE)
      ELSE
        ELBEG  = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NSTEPS = NSTEPS + 1
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      DO WHILE (ISON .GT. 0)
           NUMSTK = NUMSTK + 1
           NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1)
           ISON = FRERE_STEPS(STEP(ISON))
      END DO
      NFRONT = ND(STEP(INODE)) + NASS
      NASS1 = NASS + NUMORG
      LREQ = 2 * NFRONT + HF
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL SMUMPS_94(N, KEEP(28),
     *        IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
          COMP = COMP+1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          END IF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      END IF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .TRUE.
      IF (KEEP(50).EQ.0) THEN
        CALL  SMUMPS_124(
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF, NFRONT, NFRONT_EFF,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS,
     *        KEEP(28),
     *        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      ELSE
        CALL SMUMPS_125( 
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF,
     *        NFRONT, NFRONT_EFF, PERM,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS,
     *        KEEP(28),
     *        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
        IF (IFLAG.LT.0) GOTO 300
      END IF
      IF (NFRONT_EFF.NE.NFRONT) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
        IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF))
        NFRONT = NFRONT_EFF
        ELSE
         Write(*,*) ' ERROR 1 during ass_niv1_ELT'
         GOTO 270
        ENDIF
      ENDIF
      NCB   = NFRONT - NASS1
      MAXFRW = MAX0(MAXFRW, NFRONT)
      ICT11 = IOLDPS + HF - 1 + NFRONT
      LAELL = NFRONT * NFRONT
      IF (LRLU .LT. LAELL) THEN
        IF (LRLUS .LT. LAELL) THEN
          GOTO 280
        ELSE
          CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *         NIRBDU, LRLU, IPTRLU,
     *         IWPOS, IWPOSCB, PTRIST, PTRAST,
     *         STEP, PIMASTER, PAMASTER, ITLOC)
          COMP = COMP + 1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          END IF
        END IF
      END IF
      LRLU = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,LAELL,KEEP)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL - 1
      A(POSELT:LAPOS2) = REAL(ZERO)
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL - 1
        A(POSELT:LAPOS2) = REAL(ZERO)
      ELSE
        APOS = POSELT
        DO JJ = 0, NFRONT -1
          A(APOS:APOS+JJ) = REAL(ZERO)
          APOS = APOS + NFRONT
        END DO
      END IF
#endif
      NASS = NASS1
      PTRAST(STEP(INODE)) = POSELT
      PTLUST_S(STEP(INODE)) = IOLDPS
      IW(IOLDPS)   = NFRONT
      IW(IOLDPS + 1) = 0
      IW(IOLDPS + 2) = -NASS1
      IW(IOLDPS + 3) = -NASS1
      IW(IOLDPS + 4) = 1
      IW(IOLDPS+5)   = NSLAVES
      IF (NUMSTK.NE.0) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK    = PIMASTER(STEP(ISON))
          LSTK      = IW(ISTCHK)
          NELIM     = IW(ISTCHK + 1)
          NPIVS     = IW(ISTCHK + 3)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON    = IW(ISTCHK + 5)
          HS        = 6 + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          LEVEL1    = NSLSON.EQ.0
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW( ISTCHK + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          SIZFI   = HS + NROWS + NCOLS 
          J1 = ISTCHK + HS + NROWS + NPIVS
          IF ( .NOT. LEVEL1 .AND. NELIM.EQ.0 ) GOTO 205
          IF (LEVEL1) THEN
           J2 = J1 + LSTK - 1
           SIZFR  = LSTK*LSTK
          ELSE
           IF ( KEEP(50).eq.0 ) THEN
             SIZFR = NELIM * LSTK
           ELSE
             SIZFR = NELIM * NELIM
           END IF
           J2 = J1 + NELIM - 1
          ENDIF
          OPASSW = OPASSW + REAL(SIZFR)
          IACHK = PAMASTER(STEP(ISON))
          IF ( KEEP(50) .eq. 0 ) THEN
            POSEL1 = PTRAST(STEP(INODE)) - NFRONT
            IF (J2.GE.J1) THEN
              DO 170 JJ = J1, J2
                APOS = POSEL1 + IW(JJ) * NFRONT
                DO 160 JJ1 = 1, LSTK
                  JJ2 = APOS + IW(J1 + JJ1 - 1) - 1
                  A(JJ2) = A(JJ2) + A(IACHK + JJ1 - 1)
  160           CONTINUE
                IACHK = IACHK + LSTK
  170         CONTINUE
            END IF
          ELSE
            IF (LEVEL1) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            CALL SMUMPS_178(
     *           A( PTRAST(STEP(INODE)) ), NFRONT, NASS1,
     *           A( IACHK ), LDA_SON,
     *           IW( J1 ), J2 - J1 + 1, NELIM
     *          )
          ENDIF
  205     IF (LEVEL1) THEN 
           IF (SAME_PROC) ISTCHK = PTRIST(STEP(ISON))
           IF (SAME_PROC) THEN
             IF (KEEP(50).NE.0) THEN
              J2 = J1 + LSTK - 1
              DO JJ = J1, J2
               IW(JJ) = IW(JJ - NROWS)
              END DO
             ELSE
              J2 = J1 + LSTK - 1
              J3 = J1 + NELIM
              DO JJ = J3, J2
               IW(JJ) = IW(JJ - NROWS)
              END DO
              IF (NELIM .NE. 0) THEN
                J3 = J3 - 1
                DO JJ = J1, J3
                 JPOS = IW(JJ) + ICT11
                 IW(JJ) = IW(JPOS)
                END DO
             ENDIF
           ENDIF
          ENDIF
           IF ( SAME_PROC ) THEN
               PTRIST(STEP( ISON )) = -99999999
           ELSE
               PIMASTER(STEP( ISON )) = -99999999
           ENDIF
           IF (ISTCHK .EQ. (IWPOSCB + 1)) THEN
             IF (SAME_PROC) SIZFI = 2
             IPTRLU = IPTRLU + SIZFR
             IWPOSCB = IWPOSCB + SIZFI
             LRLU = LRLU + SIZFR
             LRLUS = LRLUS + SIZFR
       CALL SMUMPS_471(.FALSE.,LA-LRLUS,0, -SIZFR,KEEP)
  200        IF (IWPOSCB .EQ. NIRBDU) GOTO 210
             SAME_PROC = IW(IWPOSCB+1).LT.0
             FREE = .FALSE.
             IF (SAME_PROC) THEN
               FREE  = IW(IWPOSCB+1).LT.-N
               SIZFI = 2
               IF (FREE) THEN
                 NROW = -(IW(IWPOSCB+1) + N)
                 NCOL = IW(IWPOSCB+2)
               ENDIF
             ELSEIF (IW(IWPOSCB + 5) .EQ. 0) THEN
               FREE = .TRUE.
               NCOL = IW(IWPOSCB + 1) + IW(IWPOSCB + 4)
               NROW = IW(IWPOSCB + 3)
               NSLSON  = IW(IWPOSCB + 6)
               HS      = NSLSON + 6
               IF ( IW( IWPOSCB + 4 ) .LT. 0 ) THEN
                 SIZFI = NROW + IW(IWPOSCB + 1) + HS
               ELSE
                 SIZFI = NROW + NCOL + HS
               END IF
             ENDIF
             IF (FREE) THEN
               SIZFR = NROW * NCOL
               IPTRLU = IPTRLU + SIZFR
               IWPOSCB = IWPOSCB + SIZFI
               LRLU = LRLU + SIZFR
               GOTO 200
             ENDIF
           ELSE
            IF (SAME_PROC) THEN
             IW(ISTCHK) = IW(ISTCHK) - N
            ELSE
              IW(ISTCHK + 4) = 0
            ENDIF
            LRLUS = LRLUS + SIZFR
             CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,
     *                                  -SIZFR,KEEP)
           ENDIF
          ELSE
           PDEST = ISTCHK + 6
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL SMUMPS_49( 
     &                KEEP, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
               CALL SMUMPS_210( COMM_LOAD, ASS_IRECV, 
     *                           BUFR, LBUFR, LBUFR_BYTES,
     *                           INODE, ISON, NSLAVES, IDUMMY,
     *                           NFRONT, NASS1,
     *         TROW_SIZE, IW( INDX ),
     *         PROCNODE_STEPS,
     *         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *         LRLUS, N, IW,
     *         LIW, A, LA, NIRBDU,
     *         PTRIST, PTLUST_S, PTRAST, STEP,
     *         PIMASTER, PAMASTER, NSTK_S, COMP,
     *         IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL,
     *         LEAF, ICNTL, KEEP, root,
     *         OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *         INTARR, DBLARR, ND, FRERE_STEPS,
     *         NELT+1, NELT, FRT_PTR, FRT_ELT,
     *   
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *         )
               IF ( IFLAG .LT. 0 ) GOTO 500
             ENDIF
           END DO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6
            CALL  SMUMPS_71( INODE, NFRONT, 
     *       NASS1, ISON, MYID,
     *       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     *       COMM, IERR, IW(PDEST), NSLSON, 
     *       SLAVEF, 
     *       KEEP, STEP, N, 
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *         MPI_ANY_SOURCE, MPI_ANY_TAG,
     *         STATUS, 
     *         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *         IWPOS, IWPOSCB, IPTRLU,
     *         LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *         PTLUST_S,
     *         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *         IFLAG, IERROR, COMM,
     *         NBPROCFILS,
     *         IPOOL, LPOOL, LEAF,
     *         NBFIN, MYID, SLAVEF,
     *         root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *         INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *         NELT+1, NELT, FRT_PTR, FRT_ELT, 
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
               IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
           ENDIF
          ENDIF
  210   ISON = FRERE_STEPS(STEP(ISON))
  220 CONTINUE
      END IF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          AINPUT    = AII + II - J1
          ICT12 = POSELT + (I-1) * NFRONT
          DO JJ=J1,J2
           APOS2 = ICT12 + INTARR(JJ) - 1
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
           AINPUT = AINPUT + SIZE_ELTI
          END DO
         ELSE
          ICT12 = POSELT - NFRONT + I - 1 
          ICT21 = POSELT + (I-1)*NFRONT - 1 
          DO JJ=II,J2
           J =  INTARR(JJ)
           IF (I.LT.J) THEN
              APOS2 = ICT12 + J*NFRONT
           ELSE
              APOS2 = ICT21 + J 
           ENDIF
           A(APOS2) = A(APOS2) + DBLARR(AII)
           AII = AII + 1
          END DO
         END IF
        END DO
      END DO
      GOTO 500
  270 CONTINUE
      IFLAG = -8
      IERROR = LREQ
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER ALLOCATION during assembly'
      ENDIF
      GOTO 490
  280 CONTINUE
      IFLAG = -9
      IERROR = LAELL - LRLUS
      GOTO 500
  290 CONTINUE
      IFLAG = -17
      LREQ = NCBSON + 6*NSLSON
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION during assembly'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK 
  490 CALL  SMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_36
      SUBROUTINE SMUMPS_37( COMM_LOAD, ASS_IRECV,
     *    NELT, FRT_PTR, FRT_ELT,
     *    N, INODE, IW, LIW, A, LA, IFLAG,
     *    IERROR, ND, NIRBDU, FILS, FRERE_STEPS,
     *    CAND, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    MAXFRW, root,
     *    OPASSW, OPELIW, PTRIST, PTLUST_S, 
     *    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_S,
     *    PTRAIW, ITLOC, NSTEPS, 
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, 
     *    ICNTL, KEEP,INTARR,DBLARR, 
     *    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, LEAF, IPOOL, LPOOL,
     *    PERM,
     *    MEM_DISTRIB)
      USE SMUMPS_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER NELT, N,LIW,LA,NIRBDU,NSTEPS, NBFIN
      INTEGER KEEP(150), ICNTL(40)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     *        LPOOL, LEAF, 
     *        LRLU, IPTRLU,IWPOS, LRLUS,
     *        POSFAC, IWPOSCB, COMP, SLAVEF
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER IW(LIW), ITLOC(N),
     *        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)),
     *        FILS(N), FRERE_STEPS(KEEP(28)),
     *        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28)),
     *        STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), PERM(N)
      INTEGER   CAND(SLAVEF+1,MAX(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      REAL A(LA)
      REAL  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
      INTEGER MYID, COMM
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER NBPROCFILS(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INTEGER LP, HS, HF, HF_OLD, NSLAVES_OLD,NCBSON
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ,LAELL
      INTEGER NCB, IERR_MPI
      INTEGER LAPOS2,J1,J2,J3
      INTEGER NELIM,JJ,JJ1,JJ2,NPIVS,NCOLS,NROWS,
     *        IORG
      INTEGER IACHK, LDAFS, LDA_SON
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4
      INTEGER APOS, APOS2, AINPUT, POSELT, POSEL1, ICT12, ICT21
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, ISLAVE, PDEST
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, K, ELBEG, NUMELT, I, J, AII
      LOGICAL FLAG, SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER NSLAVES_less, ITEMP, NMB_OF_CAND
      logical :: force_cand
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST
      INTEGER IZERO 
      INTEGER IDUMMY(1)
      INTEGER PDEST1(1)
      PARAMETER( IZERO = 0 )
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      INTRINSIC REAL
      REAL DATA_SIZE
      REAL ZERO
      DATA ZERO /0.0E0/
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      NSTEPS = NSTEPS + 1
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .NE. 0 ) THEN
        ELBEG = FRT_PTR(INODE)
      ELSE
        ELBEG = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1)
         ISON = FRERE_STEPS(STEP(ISON))
      END DO
      NFRONT = ND(STEP(INODE)) + NASS
      MAXFRW = MAX0(MAXFRW, NFRONT)
      NASS1 = NASS + NUMORG
      NCB   = NFRONT - NASS1
      DATA_SIZE = NASS1 * NCB
      IF((KEEP(24).eq.0).or.(KEEP(24).eq.1)) then
         force_cand=.FALSE.
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      end if
      IF (force_cand) THEN
         INIV2 = ISTEP_TO_INIV2( STEP( INODE ))
         SIZE_TMP_SLAVES_LIST = CAND( SLAVEF+1, INIV2 )
      ELSE
         INIV2 = 1
         SIZE_TMP_SLAVES_LIST = SLAVEF - 1
      ENDIF
      ALLOCATE(TMP_SLAVES_LIST(SIZE_TMP_SLAVES_LIST),stat=allocok)
      IF (allocok > 0 ) THEN
        GOTO 265
      ENDIF
      CALL SMUMPS_472( SLAVEF, KEEP,
     *     CAND(1,INIV2),
     *     MEM_DISTRIB(0), NCB, NFRONT, NSLAVES,
     *     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *     TMP_SLAVES_LIST,
     *     SIZE_TMP_SLAVES_LIST )
      HF   = NSLAVES + 6
      LREQ = 2 * NFRONT + HF
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL SMUMPS_94(N, KEEP(28),
     *        IW, LIW, A, LA, NIRBDU,
     *        LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
          COMP = COMP+1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv2'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          ENDIF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      ENDIF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .FALSE.
      IF (KEEP(50).EQ.0) THEN
        CALL  SMUMPS_124(
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF, NFRONT,NFRONT_EFF,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS, KEEP(28),
     *        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      ELSE
        CALL SMUMPS_125(
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF,
     *        NFRONT, NFRONT_EFF, PERM,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS,
     *        KEEP(28), SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
        IF (IFLAG.LT.0) GOTO 250
      ENDIF
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
            NCB    = NFRONT_EFF - NASS1
            NSLAVES_OLD = NSLAVES
            HF_OLD      = HF
            CALL SMUMPS_472( SLAVEF, KEEP,
     *      CAND(1,INIV2),
     *      MEM_DISTRIB(0), NCB, NFRONT_EFF, NSLAVES,
     *      TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *      TMP_SLAVES_LIST, SIZE_TMP_SLAVES_LIST )
            HF = NSLAVES + 6
            IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF)) -
     &                   (NSLAVES_OLD - NSLAVES)
            IF (NSLAVES_OLD .NE. NSLAVES) THEN
              IF (NSLAVES_OLD > NSLAVES) THEN
               IW(IOLDPS+HF: IOLDPS+HF+2*NFRONT_EFF-1) =
     &         IW(IOLDPS+HF_OLD: IOLDPS+HF_OLD+2*NFRONT_EFF-1)
              ELSE
               IF (IWPOS - 1 > IWPOSCB ) GOTO 270
               DO JJ=2*NFRONT_EFF-1, 0, -1
                 IW(IOLDPS+HF+JJ) = IW(IOLDPS+HF_OLD+JJ)
               ENDDO
              END IF
            END IF
            NFRONT = NFRONT_EFF
        ELSE
          Write(*,*) ' ERROR 2 during ass_niv2'
          GOTO 270
        ENDIF
      ENDIF
      MAXFRW = MAX0(MAXFRW, NFRONT)
      PTLUST_S(STEP(INODE)) = IOLDPS
      IW(IOLDPS + 1) = 0
      IW(IOLDPS + 2) = -NASS1
      IW(IOLDPS + 3) = -NASS1
      IW(IOLDPS + 4) = 1
      IW(IOLDPS)   = NFRONT
      IW(IOLDPS+5) = NSLAVES
      IW(IOLDPS+6:IOLDPS+5+NSLAVES)=TMP_SLAVES_LIST(1:NSLAVES)
      DEALLOCATE(TMP_SLAVES_LIST)
#if defined(OLD_LOAD_MECHANISM)
      IF (KEEP(73) .EQ. 0) THEN
#endif
        CALL SMUMPS_461(MYID, SLAVEF, COMM_LOAD,
     *     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *     NASS1, KEEP, IW(IOLDPS+6), NSLAVES)
#if defined(OLD_LOAD_MECHANISM)
      ENDIF
#endif
      IF (KEEP(50).EQ.0) THEN
        LAELL = NASS1 * NFRONT
        LDAFS = NFRONT
      ELSE
        LAELL = NASS1 * NASS1
        LDAFS = NASS1
      ENDIF
      IF (LRLU .LT. LAELL) THEN
        IF (LRLUS .LT. LAELL) THEN
          GOTO 280
        ELSE
          CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv2'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          ENDIF
        ENDIF
      ENDIF
      LRLU = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      POSELT = POSFAC
      PTRAST(STEP(INODE)) = POSELT
      POSFAC = POSFAC + LAELL
      CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,LAELL,KEEP)
      POSEL1 = POSELT - LDAFS
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL - 1
      A(POSELT:LAPOS2) = REAL(ZERO)
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL - 1
        A(POSELT:LAPOS2) = REAL(ZERO)
      ELSE
        APOS = POSELT
        DO JJ = 0, LDAFS - 1
          A(APOS:APOS+JJ) = REAL(ZERO)
          APOS = APOS + LDAFS
        END DO
      END IF
#endif
      IF ((NUMSTK.NE.0).AND.(NASS.NE.0)) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK = PIMASTER(STEP(ISON))
          NELIM = IW(ISTCHK + 1)
          IF (NELIM.EQ.0) GOTO 210
          LSTK    = IW(ISTCHK)
          NPIVS   = IW(ISTCHK + 3)
          IF (NPIVS.LT.0) NPIVS=0
          NSLSON  = IW(ISTCHK + 5)
          HS      = 6 + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW(ISTCHK + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          OPASSW = OPASSW + REAL(NELIM*LSTK)
          J1 = ISTCHK + HS + NROWS + NPIVS
          J2 = J1 + NELIM - 1
          IACHK = PAMASTER(STEP(ISON))
          IF (KEEP(50).eq.0) THEN
           DO 170 JJ = J1, J2
            APOS = POSEL1 + IW(JJ) * LDAFS
            DO 160 JJ1 = 1, LSTK
              JJ2 = APOS + IW(J1 + JJ1 - 1) - 1
              A(JJ2) = A(JJ2) + A(IACHK + JJ1 - 1)
  160       CONTINUE
            IACHK = IACHK + LSTK
  170      CONTINUE
          ELSE
            IF (NSLSON.EQ.0) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            CALL SMUMPS_178(
     *           A( POSELT ), LDAFS, NASS1,
     *           A( IACHK ), LDA_SON,
     *           IW( J1 ), NELIM, NELIM
     *          )
          ENDIF
  210     ISON = FRERE_STEPS(STEP(ISON))
  220   CONTINUE
      ENDIF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          IF (I.LE.NASS1) THEN
           AINPUT    = AII + II - J1
           ICT12 = POSELT + (I-1) * LDAFS
           DO JJ=J1,J2
            APOS2 = ICT12 + INTARR(JJ) - 1
            A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            AINPUT = AINPUT + SIZE_ELTI
           END DO
          ENDIF
         ELSE
          ICT12 = POSELT - LDAFS + I - 1
          ICT21 = POSELT + (I-1)*LDAFS - 1
          IF ( I .GT. NASS1 ) THEN
           AII = AII + J2 - II + 1
           CYCLE
          ENDIF
          DO JJ=II,J2
           J =  INTARR(JJ)
           IF ( J .LE. NASS1) THEN
              IF (I.LT.J) THEN
                  APOS2 = ICT12 + J*LDAFS
              ELSE
                  APOS2 = ICT21 + J
              ENDIF
              A(APOS2) = A(APOS2) + DBLARR(AII)
           ENDIF
           AII = AII + 1
          END DO
         END IF
        END DO
      END DO
      PTRCOL = IOLDPS + HF + NFRONT 
      PTRROW = IOLDPS + HF + NASS1 
      PDEST  = IOLDPS + 6
      DO ISLAVE = 1, NSLAVES
              CALL SMUMPS_49( 
     &                KEEP, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB,
     &                NSLAVES, 
     &                NBLIG, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
        IERR = -1
        DO WHILE (IERR .EQ.-1)
         IF ( KEEP(50) .eq. 0 ) THEN
           NBCOL =  NFRONT
           CALL SMUMPS_68( INODE,
     *      NBPROCFILS(STEP(INODE)),
     *      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     *      IZERO, IDUMMY,
     *      IW(PDEST), NFRONT, COMM, IERR)
         ELSE
           NBCOL = NASS1+SHIFT_INDEX+NBLIG
           CALL SMUMPS_68( INODE,
     *      NBPROCFILS(STEP(INODE)),
     *      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     *      NSLAVES-ISLAVE, IW( PTLUST_S(STEP(INODE))+6+ISLAVE),
     *      IW(PDEST), NFRONT, COMM, IERR)
         ENDIF
         IF (IERR.EQ.-1) THEN
          BLOCKING  = .FALSE.
          SET_IRECV = .TRUE.
          MESSAGE_RECEIVED = .FALSE.
          CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *     MPI_ANY_SOURCE, MPI_ANY_TAG,
     *     STATUS, BUFR, LBUFR,
     *     LBUFR_BYTES,
     *     PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *     LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *     IERROR, COMM,
     *     NBPROCFILS,
     *     IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *     root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *     INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *     NELT+1, NELT, FRT_PTR, FRT_ELT,
     *     ISTEP_TO_INIV2, TAB_POS_IN_PERE)
           IF ( IFLAG .LT. 0 ) GOTO 500
          IF (MESSAGE_RECEIVED) THEN
           IOLDPS = PTLUST_S(STEP(INODE))
           PTRCOL = IOLDPS + HF + NFRONT
           PTRROW = IOLDPS + HF + NASS1 + SHIFT_INDEX
          ENDIF
         ENDIF
        END DO
        IF (IERR .EQ. -2) GOTO 300
        PTRROW = PTRROW + NBLIG
        PDEST  = PDEST + 1
      END DO
      IF (NUMSTK.EQ.0) GOTO 500
      ISON = IFSON
      DO IELL = 1, NUMSTK
        ISTCHK = PIMASTER(STEP(ISON))
        NELIM = IW(ISTCHK + 1)
        LSTK    = IW(ISTCHK)
        NPIVS   = IW(ISTCHK + 3)
        IF ( NPIVS .LT. 0 ) NPIVS = 0
        NSLSON  = IW(ISTCHK + 5)
        HS      = 6 + NSLSON
        NCOLS     = NPIVS + LSTK
        SAME_PROC     = (ISTCHK.LE.IWPOS)
        IF (.NOT.SAME_PROC) THEN
         NROWS = IW(ISTCHK + 2)
        ELSE
         NROWS = NCOLS
        ENDIF
        PDEST   = ISTCHK + 6
        NCBSON  = LSTK - NELIM
        PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
        IF (NSLSON.EQ.0) THEN
          NSLSON = 1
          PDEST1(1)  = SMUMPS_275(STEP(ISON),
     *                 PROCNODE_STEPS, SLAVEF)
          IF (PDEST1(1).EQ.MYID) THEN
            CALL SMUMPS_211( COMM_LOAD, ASS_IRECV, 
     *      BUFR, LBUFR, LBUFR_BYTES,
     *      INODE, ISON, NSLAVES, IW( PTLUST_S(STEP(INODE)) + 6 ),
     *      NFRONT, NASS1, NCBSON, IW( PTRCOL ),
     *      PROCNODE_STEPS,
     *      SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *      LRLUS, N, IW,
     *      LIW, A, LA, NIRBDU,
     *      PTRIST, PTLUST_S, PTRAST, STEP,
     *      PIMASTER, PAMASTER, NSTK_S, COMP,
     *      IFLAG, IERROR, MYID, COMM, NBPROCFILS, 
     *      IPOOL, LPOOL, LEAF,
     *      ICNTL, KEEP, root,
     *      OPASSW, OPELIW,
     *      ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *      ND, FRERE_STEPS, NELT+1, NELT, 
     *      FRT_PTR, FRT_ELT, 
     *       
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *      )
           IF ( IFLAG .LT. 0 ) GOTO 500
          ELSE
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            CALL  SMUMPS_71( INODE, NFRONT,
     *      NASS1, ISON, MYID,
     *      NSLAVES, IW( PTLUST_S(STEP(INODE)) + 6 ),
     *      IW(PTRCOL), NCBSON,
     *      COMM, IERR, PDEST1, NSLSON, SLAVEF, 
     *      KEEP, STEP, N, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *        MPI_ANY_SOURCE, MPI_ANY_TAG,
     *        STATUS, BUFR, LBUFR, LBUFR_BYTES,
     *        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *        LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *        IERROR, COMM,
     *        NBPROCFILS,
     *        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *        root,OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE)
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
          ENDIF
        ELSE
          DO ISLAVE = 0, NSLSON-1
            IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
               CALL SMUMPS_49( 
     &                KEEP, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON,
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX        = PTRCOL + SHIFT_INDEX
              CALL SMUMPS_210( COMM_LOAD, ASS_IRECV, 
     *        BUFR, LBUFR, LBUFR_BYTES,
     *        INODE, ISON, NSLAVES, IW( PTLUST_S(STEP(INODE))+6),
     *        NFRONT, NASS1,
     *        TROW_SIZE, IW( INDX ),
     *        PROCNODE_STEPS,
     *        SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *        LRLUS, N, IW,
     *        LIW, A, LA, NIRBDU,
     *        PTRIST, PTLUST_S, PTRAST, STEP,
     *        PIMASTER, PAMASTER, NSTK_S, COMP,
     *        IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     *        ICNTL, KEEP, root,
     *        OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ND, FRERE_STEPS,
     *        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     * 
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *        )
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
          END DO
          IF (PIMASTER(STEP(ISON)).GT.0) THEN
          IERR = -1
          DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6
            CALL  SMUMPS_71( INODE, NFRONT,
     *      NASS1, ISON, MYID,
     *      NSLAVES, IW(PTLUST_S(STEP(INODE))+6),
     *      IW(PTRCOL), NCBSON,
     *      COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     *      KEEP, STEP, N, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *        MPI_ANY_SOURCE, MPI_ANY_TAG,
     *        STATUS, BUFR, LBUFR,
     *        LBUFR_BYTES,
     *        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *        LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *        IERROR, COMM,
     *        NBPROCFILS,
     *        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *        root,OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *        NELT+1, NELT, FRT_PTR, FRT_ELT,
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE )
             IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
          END DO
          IF (IERR .EQ. -2) GOTO 290
          ENDIF
        ENDIF
       ISON = FRERE_STEPS(STEP(ISON))
      END DO
      GOTO 500
  250 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION during assembly'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK + 1
      GOTO 490
  265 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE ALLOCATING TMP_SLAVES_LIST',
     &                 ' DYNAMIC ALLOCATION during assembly'
      ENDIF
      IFLAG  = -13
      IERROR = SIZE_TMP_SLAVES_LIST
      GOTO 490
  270 CONTINUE
      IFLAG = -8
      IERROR = LREQ
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER ALLOCATION during assembly'
      ENDIF
      GOTO 490
  280 CONTINUE
      IFLAG = -9
      IERROR = LAELL - LRLUS
      GOTO 490
  290 CONTINUE
      IFLAG = -17
      LREQ = NCBSON + 6*NSLSON
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IFLAG = -17
      LREQ = NBLIG + NBCOL + 4
      IERROR =  LREQ  * KEEP( 34 ) 
  490 CALL SMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_37
      SUBROUTINE SMUMPS_123(
     *    NELT, FRT_PTR, FRT_ELT,
     *    N, INODE, IW, LIW, A, LA, 
     *    NBROWS, NBCOLS, ROWLIST, COLLIST, VALSON, 
     *    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *    FILS, PTRARW, PTRAIW, INTARR, DBLARR, 
     *    ICNTL, KEEP, MYID)
      IMPLICIT NONE
      INTEGER NELT, N,LIW,LA
      INTEGER KEEP(150), ICNTL(40)
      INTEGER INODE, MYID
      INTEGER NBROWS, NBCOLS 
      INTEGER ROWLIST(NBROWS), COLLIST(NBCOLS)
      INTEGER IW(LIW), ITLOC(N), STEP(N),
     *        PTRIST(KEEP(28)), PTRAST(KEEP(28)),
     *        FILS(N), PTRARW(NELT+1), 
     *        PTRAIW(NELT+1)
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      REAL A(LA), VALSON(NBCOLS,NBROWS),
     *        DBLARR(MAX(1,KEEP(13)))
      REAL OPASSW, OPELIW
      INTEGER IOLDPS, POSELT, NBCOLF, NBROWF, NSLAVES, HF,
     *        K1,K2,K,I,J,POSEL1,APOS,JPOS,NASS,JJ,
     *        IN,AINPUT,JK,J1,J2,IJROW,ICT12, ILOC, 
     *        ELBEG, NUMELT, ELTI, SIZE_ELTI, IPOS, 
     *        IPOS1, IPOS2, AII, APOS2, II, IELL
      REAL  ZERO
      PARAMETER (ZERO=0.0E0)
      INTRINSIC REAL
      IOLDPS  = PTRIST(STEP(INODE))
      POSELT  = PTRAST(STEP(INODE))
      NBCOLF  = IW(IOLDPS)
      NBROWF  = IW(IOLDPS+2)
      NASS    = IW(IOLDPS+1)
       IF ( NBROWS .GT. NBROWF ) THEN
          WRITE(*,*) ' ERR: ERROR : NBROWS > NBROWF'
          WRITE(*,*) ' ERR: INODE =', INODE
          WRITE(*,*) ' ERR: NBROW=',NBROWS,'NBROWF=',NBROWF
          WRITE(*,*) ' ERR: ROW_LIST=', ROWLIST
          CALL SMUMPS_ABORT()
       END IF
      NSLAVES = IW(IOLDPS+5)
      HF      = 6 + NSLAVES
      IF ((NASS.GE.0).AND.(NBROWS.LE.0)) GOTO 500
      IF (NASS.LT.0) THEN
          NASS         = -NASS
          IW(IOLDPS+1) = NASS
          A(POSELT:POSELT+NBROWF*NBCOLF-1) = REAL(ZERO)
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -JPOS
           JPOS     = JPOS + 1
          END DO
          K1 = IOLDPS + HF 
          K2 = K1 + NBROWF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
           JPOS     = JPOS + 1
          END DO
          ELBEG  = FRT_PTR(INODE)
          NUMELT = FRT_PTR(INODE+1) - ELBEG
          DO IELL=ELBEG,ELBEG+NUMELT-1
           ELTI = FRT_ELT(IELL)
           J1= PTRAIW(ELTI)
           J2= PTRAIW(ELTI+1)-1
           AII = PTRARW(ELTI)
           SIZE_ELTI = J2 - J1 + 1
           DO II=J1,J2
            I = ITLOC(INTARR(II))
            IF (KEEP(50).EQ.0) THEN
             IF (I.LE.0) CYCLE
             AINPUT    = AII + II - J1
             IPOS = MOD(I,NBCOLF)
             ICT12 = POSELT + (IPOS-1) * NBCOLF
             DO JJ = J1, J2
              JPOS = ITLOC(INTARR(JJ))
              IF (JPOS.LE.0) THEN 
                   JPOS = -JPOS
              ELSE
                   JPOS = JPOS/NBCOLF
              END IF
              APOS2    = ICT12 + JPOS - 1
              A(APOS2) = A(APOS2) +  DBLARR(AINPUT)
              AINPUT   = AINPUT + SIZE_ELTI
             END DO
            ELSE
              IF ( I .EQ. 0 ) THEN 
               AII = AII + J2 - II + 1
               CYCLE
              ENDIF
              IF ( I .LE. 0 ) THEN 
               IPOS1 = -I
               IPOS2 = 0
              ELSE 
               IPOS1 = I/NBCOLF
               IPOS2 = MOD(I,NBCOLF)
              END IF
              ICT12 =  POSELT + (IPOS2-1)*NBCOLF
              DO JJ=II,J2
               AII = AII + 1
               J = ITLOC(INTARR(JJ))
               IF ( J .EQ. 0 ) CYCLE
               IF ( IPOS2.EQ.0 .AND. J.LE.0) CYCLE
               IF ( J .LE. 0 ) THEN
                JPOS = -J
               ELSE
                JPOS = J/NBCOLF
               END IF
               IF ( (IPOS1.GE.JPOS) .AND. (IPOS2.GT.0) ) THEN
                 APOS2 = ICT12  + JPOS - 1
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
               IF ( (IPOS1.LT.JPOS) .AND. (J.GT.0) ) THEN
                 IPOS = MOD(J,NBCOLF)
                 JPOS = IPOS1
                 APOS2 = POSELT + (IPOS-1)*NBCOLF + JPOS - 1
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
              END DO
            END IF
           END DO
          END DO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          END DO
      END IF
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          END DO
          POSEL1 = POSELT - NBCOLF
          IF (KEEP(50).EQ.0) THEN
           DO I=1,NBROWS
            APOS = POSEL1 + ROWLIST(I) * NBCOLF
            DO J=1,NBCOLS
             K = APOS + ITLOC(COLLIST(J)) - 1
             A(K) = A(K) + VALSON(J,I)
            END DO
           END DO
          ELSE
           DO I=1,NBROWS
            APOS = POSEL1 + ROWLIST(I) * NBCOLF
            DO J=1,NBCOLS
             IF (ITLOC(COLLIST(J)) .EQ. 0) EXIT
             K = APOS + ITLOC(COLLIST(J)) - 1
             A(K) = A(K) + VALSON(J,I)
            END DO
           END DO
          END IF
          OPASSW = OPASSW + NBROWS*NBCOLS
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          END DO
      END IF
 500  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_123
      SUBROUTINE SMUMPS_126(
     *            N, NELT,
     *            COMM, MYID, SLAVEF,
     *            IELPTR_LOC, RELPTR_LOC,
     *            ELTVAR_LOC, ELTVAL_LOC,
     *            KEEP, MAXELT_SIZE,
     *            FRTPTR, FRTELT, A, LA, FILS,
     *            id, root )
      IMPLICIT NONE
      INCLUDE 'smumps_struc.h'
      INTEGER N, NELT, COMM, MYID, SLAVEF, MAXELT_SIZE, MSGLEN, LA
      INTEGER FRTPTR( N+1 )
      INTEGER FRTELT( NELT ), FILS ( N )
      INTEGER KEEP(150)
      INTEGER IELPTR_LOC( NELT + 1 ), RELPTR_LOC( NELT + 1 )
      INTEGER ELTVAR_LOC( MAX(1,KEEP(14)) )
      REAL ELTVAL_LOC( MAX(1,KEEP(13)) )
      REAL A( LA )
      TYPE(SMUMPS_STRUC)     :: id
      TYPE(SMUMPS_ROOT_STRUC) :: root
      INTEGER NUMROC
      EXTERNAL NUMROC
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR_MPI
      INTEGER MSGTAG
      INTEGER allocok
      INTEGER I, DEST, MAXELT_REAL_SIZE, MPG, IEL, SIZEI, SIZER
      INTEGER NBRECORDS, NBUF
      INTEGER RECV_IELTPTR, RECV_RELTPTR
      INTEGER IELTPTR, RELTPTR, INODE
      LOGICAL FINI, PROKG, I_AM_SLAVE
      INTEGER PTR_ROOT, LOCAL_M, LOCAL_N, LP, IBEG, IGLOB, JGLOB
      INTEGER ARROW_ROOT
      INTEGER IELT, J, K, NB_REC, IREC
      INTEGER ILOCROOT, JLOCROOT, IPOSROOT, JPOSROOT, IPTR
      INTEGER JCOL_GRID, IROW_GRID
      INTEGER IVALPTR
      INTEGER NBELROOT
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      REAL  VAL
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER, DIMENSION( :, : ), ALLOCATABLE :: BUFI
      REAL, DIMENSION( :, : ), ALLOCATABLE :: BUFR
      REAL, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_R
      INTEGER, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_I
      INTEGER, DIMENSION( : ), ALLOCATABLE :: ELROOTPOS
      INTEGER, DIMENSION( : ), ALLOCATABLE, TARGET :: RG2LALLOC
      INTEGER, DIMENSION( : ), POINTER     :: RG2L
      MPG = id%ICNTL(3)
      LP  = id%ICNTL(1)
      I_AM_SLAVE = ( KEEP(46) .eq. 1 .or. MYID .ne.MASTER )
      PROKG = ( MPG > 0 .and. MYID .eq. MASTER )
      KEEP(49) = 0
      ARROW_ROOT = 0
      IF ( MYID .eq. MASTER ) THEN
        IF ( KEEP(46) .eq. 0 ) THEN
          NBUF = SLAVEF
        ELSE
          NBUF = SLAVEF - 1
        END IF
        NBRECORDS = KEEP(39)
        IF ( KEEP(50) .eq. 0 ) THEN
          MAXELT_REAL_SIZE = MAXELT_SIZE * MAXELT_SIZE
        ELSE
          MAXELT_REAL_SIZE = MAXELT_SIZE * (MAXELT_SIZE+1)/2
        END IF
        IF ( MAXELT_REAL_SIZE .GT. KEEP(39) ) THEN
          NBRECORDS = MAXELT_REAL_SIZE
          IF ( MPG .GT. 0 ) THEN
            WRITE(MPG,*)
     * ' ** Warning : For element distrib NBRECORDS set to ',
     * MAXELT_REAL_SIZE,' because one element is large'
          END IF
        END IF
        ALLOCATE( BUFI( 2*NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = 2*NBRECORDS + 1
          GOTO 100
        END IF
        ALLOCATE( BUFR( NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = NBRECORDS + 1
          GOTO 100
        END IF
        IF ( KEEP(52) .ne. 0 ) THEN
          ALLOCATE( TEMP_ELT_R( MAXELT_REAL_SIZE ), stat =allocok )
          IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_REAL_SIZE
            GOTO 100
          END IF
        END IF
        ALLOCATE( TEMP_ELT_I( MAXELT_SIZE ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_SIZE
            GOTO 100
        END IF
        IF ( KEEP(38) .ne. 0 ) THEN
          NBELROOT = FRTPTR(KEEP(38)+1)-FRTPTR(KEEP(38))
          WRITE(*,*) 'NBELROOT = ', NBELROOT
          ALLOCATE( ELROOTPOS( NBELROOT ), stat = allocok )
          IF ( allocok .gt. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = NBELROOT
              GOTO 100
          END IF
          IF (KEEP(46) .eq. 0 ) THEN
           ALLOCATE( RG2LALLOC( N ), stat = allocok )
           IF ( allocok .gt. 0 ) THEN
               id%INFO(1) = -13
               id%INFO(2) = N
               GOTO 100
           END IF
           INODE = KEEP(38)
           I     = 1
           DO WHILE ( INODE .GT. 0 )
             RG2LALLOC( INODE ) = I
             INODE = FILS( INODE )
             I = I + 1
           END DO
           RG2L => RG2LALLOC
          ELSE 
           RG2L => root%RG2L_ROW
          END IF
        END IF
        DO I = 1, NBUF
          BUFI( 1, I ) = 0
          BUFR( 1, I ) = REAL(0)
        END DO
      END IF
 100  CONTINUE
      CALL SMUMPS_276( id%ICNTL, id%INFO, COMM, MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( NBRECORDS, 1, MPI_INTEGER, MASTER,
     *                COMM, IERR_MPI )
      RECV_IELTPTR = 1
      RECV_RELTPTR = 1
      IF ( MYID .eq. MASTER ) THEN
        NBELROOT = 0
        RELTPTR = 1
        RELPTR_LOC(1) = 1
        DO IEL = 1, NELT
          IELTPTR = id%ELTPTR( IEL )
          SIZEI   = id%ELTPTR( IEL + 1 ) - IELTPTR
          IF ( KEEP( 50 ) .eq. 0 ) THEN
            SIZER = SIZEI * SIZEI
          ELSE
            SIZER = SIZEI * ( SIZEI + 1 ) / 2
          END IF
          DEST = id%ELTPROC( IEL )
          IF ( DEST .eq. -2 ) THEN
            NBELROOT = NBELROOT + 1
            FRTELT( FRTPTR(KEEP(38)) + NBELROOT - 1 ) = IEL
            ELROOTPOS( NBELROOT ) = RELTPTR
            GOTO 200
          END IF
          IF ( DEST .ge. 0 .and. KEEP(46) .eq. 0 ) DEST = DEST + 1
          IF ( KEEP(52) .ne. 0 ) THEN
            CALL SMUMPS_288( N, SIZEI, SIZER,
     *               id%ELTVAR( IELTPTR ), id%A_ELT( RELTPTR ),
     *               TEMP_ELT_R(1), MAXELT_REAL_SIZE,
     *               id%ROWSCA, id%COLSCA, KEEP(50) )
          END IF
          IF ( DEST .eq. 0 .or. ( DEST .eq. -1 .and. KEEP(46) .ne. 0 ) )
     *      THEN
            ELTVAR_LOC( RECV_IELTPTR: RECV_IELTPTR + SIZEI - 1 )
     *      = id%ELTVAR( IELTPTR: IELTPTR + SIZEI - 1 )
            RECV_IELTPTR = RECV_IELTPTR + SIZEI
            IF ( KEEP(52) .ne. 0 ) THEN
              ELTVAL_LOC( RECV_RELTPTR: RECV_RELTPTR + SIZER - 1)
     *        = TEMP_ELT_R( 1: SIZER )
              RECV_RELTPTR = RECV_RELTPTR + SIZER
            END IF
          END IF
          IF ( DEST .NE. 0 ) THEN
            IF ( KEEP(52) .eq. 0 ) THEN
              CALL SMUMPS_127(
     *           id%ELTVAR(IELTPTR),
     *           id%A_ELT (RELTPTR),
     *           SIZEI, SIZER,
     *
     *           DEST, NBUF, NBRECORDS,
     *           BUFI, BUFR, COMM )
            ELSE
              CALL SMUMPS_127(
     *           id%ELTVAR(IELTPTR),
     *           TEMP_ELT_R( 1 ),
     *           SIZEI, SIZER,
     *
     *           DEST, NBUF, NBRECORDS,
     *           BUFI, BUFR, COMM )
            END IF
          END IF
 200      CONTINUE
          RELTPTR = RELTPTR + SIZER
          IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
            RELPTR_LOC( IEL + 1 ) = RELTPTR
          ELSE
            RELPTR_LOC( IEL + 1 ) = RECV_RELTPTR
          ENDIF
        END DO
        IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
          KEEP(13) = RELTPTR - 1
        ELSE
          KEEP(13) = RECV_RELTPTR - 1
        ENDIF
        IF ( RELTPTR - 1 .ne. id%NA_ELT ) THEN
          WRITE(*,*) ' ** ERROR ELT DIST: RELPTR - 1 / id%NA_ELT=',
     *               RELTPTR - 1,id%NA_ELT
          CALL SMUMPS_ABORT()
        END IF
        DEST = -2
        IELTPTR = 1
        RELTPTR = 1
        SIZEI   = 1
        SIZER   = 1
        CALL SMUMPS_127(
     *           id%ELTVAR(IELTPTR),
     *           id%A_ELT (RELTPTR),
     *           SIZEI, SIZER,
     *
     *           DEST, NBUF, NBRECORDS,
     *           BUFI, BUFR, COMM )
        IF ( KEEP(52) .NE. 0 ) DEALLOCATE( TEMP_ELT_R )
      ELSE
        FINI = ( RECV_IELTPTR .eq. IELPTR_LOC( NELT+1 )
     *     .and. RECV_RELTPTR .eq. RELPTR_LOC( NELT+1 ) )
        DO WHILE ( .not. FINI )
          CALL MPI_PROBE( MASTER, MPI_ANY_TAG,
     *                    COMM, STATUS, IERR_MPI )
          MSGTAG = STATUS( MPI_TAG    )
          SELECT CASE ( MSGTAG )
             CASE( ELT_INT )
               CALL MPI_GET_COUNT( STATUS, MPI_INTEGER,
     *                             MSGLEN, IERR_MPI )
               CALL MPI_RECV( ELTVAR_LOC( RECV_IELTPTR ), MSGLEN,
     *            MPI_INTEGER, MASTER, ELT_INT,
     *            COMM, STATUS, IERR_MPI )
               RECV_IELTPTR = RECV_IELTPTR + MSGLEN
             CASE( ELT_REAL )
                CALL MPI_GET_COUNT( STATUS, MPI_REAL,
     *                              MSGLEN, IERR_MPI )
                CALL MPI_RECV( ELTVAL_LOC( RECV_RELTPTR ), MSGLEN,
     *            MPI_REAL, MASTER, ELT_REAL,
     *            COMM, STATUS, IERR_MPI )
                RECV_RELTPTR = RECV_RELTPTR + MSGLEN
          END SELECT
          FINI = ( RECV_IELTPTR .eq. IELPTR_LOC( NELT+1 )
     *       .and. RECV_RELTPTR .eq. RELPTR_LOC( NELT+1 ) )
        END DO
      END IF
      IF ( KEEP(38) .NE. 0 ) THEN
         IF ( I_AM_SLAVE .and. root%yes ) THEN
           LOCAL_M = NUMROC( root%ROOT_SIZE, root%MBLOCK,
     *             root%MYROW, 0, root%NPROW )
           LOCAL_M = MAX( 1, LOCAL_M )
           LOCAL_N = NUMROC( root%ROOT_SIZE, root%NBLOCK,
     *               root%MYCOL, 0, root%NPCOL )
           PTR_ROOT = LA - LOCAL_M * LOCAL_N + 1
           IF ( PTR_ROOT .LE. LA ) THEN
             A( PTR_ROOT:LA ) = REAL(ZERO)
           END IF
         END IF
        IF ( MYID .NE. MASTER ) THEN
          ALLOCATE( BUFI( NBRECORDS * 2 + 1, 1 ), stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS * 2 + 1
            GOTO 250
          END IF
          ALLOCATE( BUFR( NBRECORDS, 1 )        , stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS
          END IF
        END IF
 250    CONTINUE
        CALL SMUMPS_276( id%ICNTL, id%INFO, COMM, MYID )
        IF ( id%INFO(1) .LT. 0 ) RETURN
        IF ( MYID .eq. MASTER ) THEN
        DO IPTR = FRTPTR(KEEP(38)), FRTPTR(KEEP(38)+1) - 1
          IELT = FRTELT( IPTR )
          SIZEI = id%ELTPTR( IELT + 1 ) - id%ELTPTR( IELT )
          DO I = 1, SIZEI
            TEMP_ELT_I( I ) = RG2L
     *              ( id%ELTVAR( id%ELTPTR(IELT) + I - 1 ) )
          END DO
          IVALPTR = ELROOTPOS( IPTR - FRTPTR(KEEP(38)) + 1 ) - 1
          K = 1
          DO J = 1, SIZEI
            JGLOB = id%ELTVAR( id%ELTPTR( IELT ) + J - 1 )
            IF ( KEEP(50).eq. 0 ) THEN
              IBEG = 1
            ELSE
              IBEG = J
            END IF
            DO I = IBEG, SIZEI
              IGLOB = id%ELTVAR( id%ELTPTR( IELT ) + I - 1 )
              IF ( KEEP(52) .eq. 0 ) THEN
                VAL = id%A_ELT( IVALPTR + K )
              ELSE
                VAL = id%A_ELT( IVALPTR + K ) *
     *                id%ROWSCA( IGLOB ) * id%COLSCA( JGLOB )
              END IF
              IF ( KEEP(50).eq.0 ) THEN
                IPOSROOT = TEMP_ELT_I( I )
                JPOSROOT = TEMP_ELT_I( J )
              ELSE
                IF ( TEMP_ELT_I(I) .GT. TEMP_ELT_I(J) ) THEN
                  IPOSROOT = TEMP_ELT_I(I)
                  JPOSROOT = TEMP_ELT_I(J)
                ELSE
                  IPOSROOT = TEMP_ELT_I(J)
                  JPOSROOT = TEMP_ELT_I(I)
                END IF
              END IF
              IROW_GRID = MOD( ( IPOSROOT - 1 )/root%MBLOCK,
     *                           root%NPROW )
              JCOL_GRID = MOD( ( JPOSROOT - 1 )/root%NBLOCK,
     *                           root%NPCOL )
              IF ( KEEP(46) .eq. 0 ) THEN
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID + 1
              ELSE
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID
              END IF
              IF ( DEST .eq. MASTER ) THEN
                ARROW_ROOT = ARROW_ROOT + 1
                ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
                JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
                A( PTR_ROOT + ( JLOCROOT - 1 ) * LOCAL_M
     *        + ILOCROOT - 1 ) =  A( PTR_ROOT + (JLOCROOT - 1)
     *        * LOCAL_M + ILOCROOT - 1 ) + VAL
              ELSE
                CALL SMUMPS_34(
     *          IPOSROOT, JPOSROOT, VAL, DEST, BUFI, BUFR, NBRECORDS,
     *          NBUF, LP, COMM, KEEP(46) )
              END IF
              K = K + 1
            END DO
          END DO
        END DO
        DEST = -2
        CALL SMUMPS_34(
     *          IPOSROOT, JPOSROOT, VAL,
     *          DEST, BUFI, BUFR, NBRECORDS,
     *          NBUF, LP, COMM, KEEP(46) )
        ELSE
          FINI = .FALSE.
          DO WHILE ( .not. FINI )
            CALL MPI_RECV( BUFI(1,1), 2*NBRECORDS+1,
     *                MPI_INTEGER, MASTER,
     *                ARROWHEAD,
     *                COMM, STATUS, IERR_MPI )
            NB_REC = BUFI(1,1)
            IF (NB_REC.LE.0) THEN
              FINI = .TRUE.
              NB_REC = -NB_REC
            ENDIF
            IF (NB_REC.EQ.0) EXIT
            CALL MPI_RECV( BUFR(1,1), NBRECORDS, MPI_REAL,
     *                     MASTER, ARROWHEAD,
     *                     COMM, STATUS, IERR_MPI )
            ARROW_ROOT = ARROW_ROOT + NB_REC
            DO IREC = 1, NB_REC
              IPOSROOT = BUFI( IREC * 2, 1 )
              JPOSROOT = BUFI( IREC * 2 + 1, 1 )
              VAL      = BUFR( IREC, 1 )
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
              A( PTR_ROOT + ( JLOCROOT - 1 ) * LOCAL_M
     *      + ILOCROOT - 1 ) =  A( PTR_ROOT + (JLOCROOT - 1)
     *      * LOCAL_M + ILOCROOT - 1 ) + VAL
            END DO
          END DO
          DEALLOCATE( BUFI )
          DEALLOCATE( BUFR )
        END IF
      END IF
      IF ( MYID .eq. MASTER ) THEN
        DEALLOCATE( BUFI )
        DEALLOCATE( BUFR )
        IF (KEEP(38).ne.0) THEN 
          DEALLOCATE(ELROOTPOS)
          IF (KEEP(46) .eq. 0 ) THEN
             DEALLOCATE(RG2LALLOC)
          ENDIF
        ENDIF
        DEALLOCATE( TEMP_ELT_I )
      END IF
      KEEP(49) = ARROW_ROOT
      RETURN
      END SUBROUTINE SMUMPS_126
      SUBROUTINE SMUMPS_127(
     *       ELNODES, ELVAL, SIZEI, SIZER,
     *       DEST, NBUF, NBRECORDS, BUFI, BUFR, COMM )
      IMPLICIT NONE
      INTEGER SIZEI, SIZER, DEST, NBUF, NBRECORDS, COMM
      INTEGER ELNODES( SIZEI ), BUFI( 2*NBRECORDS + 1, NBUF )
      REAL ELVAL( SIZER ), BUFR( NBRECORDS + 1, NBUF )
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER I, IBEG, IEND, IERR_MPI, NBRECR
      INTEGER NBRECI
      IF ( DEST .lt. 0 ) THEN
        IBEG = 1
        IEND = NBUF
      ELSE
        IBEG = DEST
        IEND = DEST
      END IF
      DO I = IBEG, IEND
        NBRECI = BUFI(1,I)
        IF ( NBRECI .ne.0  .and.
     *       ( DEST.eq.-2 .or.
     *         NBRECI + SIZEI .GT. 2*NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFI(2, I), NBRECI, MPI_INTEGER,
     *                    I, ELT_INT, COMM, IERR_MPI )
           BUFI(1,I) = 0
           NBRECI    = 0
        END IF
        NBRECR = INT(BUFR(1,I)+0.5E0)
        IF ( NBRECR .ne.0  .and.
     *       ( DEST.eq.-2 .or.
     *         NBRECR + SIZER .GT. NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFR(2, I), NBRECR, MPI_REAL,
     *                    I, ELT_REAL, COMM, IERR_MPI )
           BUFR(1,I) = 0.0E0
           NBRECR    = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          BUFI( 2 + NBRECI : 2 + NBRECI + SIZEI - 1, I ) =
     *    ELNODES( 1: SIZEI )
          BUFR( 2 + NBRECR : 2 + NBRECR + SIZER - 1, I ) =
     *    ELVAL( 1: SIZER )
          BUFI(1,I) = NBRECI + SIZEI
          BUFR(1,I) = REAL( NBRECR + SIZER )
        END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_127
      SUBROUTINE SMUMPS_213( ELTPTR, NELT, MAXELT_SIZE )
      INTEGER NELT, MAXELT_SIZE
      INTEGER ELTPTR( NELT + 1 )
      INTEGER I, S
      MAXELT_SIZE = 0
      DO I = 1, NELT
        S = ELTPTR( I + 1 ) - ELTPTR( I )
        MAXELT_SIZE = MAX( S, MAXELT_SIZE )
      END DO
      RETURN
      END SUBROUTINE SMUMPS_213
      SUBROUTINE SMUMPS_288( N, SIZEI, SIZER,
     *               ELTVAR, ELTVAL,
     *               SELTVAL, LSELTVAL,
     *               ROWSCA, COLSCA, K50 )
      INTEGER N, SIZEI, SIZER, LSELTVAL, K50
      INTEGER ELTVAR( SIZEI )
      REAL ELTVAL( SIZER )
      REAL SELTVAL( LSELTVAL )
      REAL ROWSCA( N ), COLSCA( N )
      INTEGER I, J, K
      K = 1
      IF ( K50 .eq. 0 ) THEN
        DO J = 1, SIZEI
          DO I = 1, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     *                   ROWSCA(ELTVAR(I)) *
     *                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      ELSE
        DO J = 1, SIZEI
          DO I = J, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     *                   ROWSCA(ELTVAR(I)) *
     *                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE SMUMPS_288
      SUBROUTINE SMUMPS_257( N, NELT, ELTPTR, ELTVAR, A_ELT,
     *                         X, Y, K50, MTYPE )
      IMPLICIT NONE
      INTEGER N, NELT, K50, MTYPE
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( * )
      REAL A_ELT( * ), X( N ), Y( N )
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      REAL TEMP
      REAL ZERO 
      PARAMETER( ZERO = 0.0E0 )
      Y = REAL(ZERO)
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
     *          Y( ELTVAR( IELPTR + I ) ) +
     *             A_ELT( K ) * TEMP
                K = K + 1
              END DO
            END DO
          ELSE
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP + 
     *          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
            END DO
          END IF
        ELSE
          DO J = 1, SIZEI
            Y( ELTVAR( IELPTR + J ) ) =
     *      Y( ELTVAR( IELPTR + J ) ) +
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            K = K + 1
            DO I = J+1, SIZEI
              Y( ELTVAR( IELPTR + I ) ) =
     *        Y( ELTVAR( IELPTR + I ) ) +
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     *        Y( ELTVAR( IELPTR + J ) ) +
     *           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_257
      SUBROUTINE SMUMPS_F77( JOB, SYM, PAR, COMM_F77, N, ICNTL, CNTL,
     *                      NZ, IRN, IRNhere, JCN, JCNhere, A, Ahere,
     *                      NZ_loc, IRN_loc, IRN_lochere,
     *                      JCN_loc, JCN_lochere,
     *                      A_loc, A_lochere,
     *                      NELT, ELTPTR, ELTPTRhere,  ELTVAR,
     *                      ELTVARhere, A_ELT, A_ELThere,
     *                      PERM_IN, PERM_INhere,
     *                      RHS, RHShere, INFO, RINFO, INFOG, RINFOG,
     *                      DEFICIENCY, SIZE_SCHUR, LISTVAR_SCHUR,
     *                      LISTVAR_SCHURhere, SCHUR, SCHURhere,
     *                      COLSCA, COLSCAhere, ROWSCA, ROWSCAhere,
     *                      INSTANCE_NUMBER )
      IMPLICIT NONE
      INTEGER JOB, SYM, PAR, COMM_F77, N, NZ, NZ_LOC, NELT,
     *        DEFICIENCY, SIZE_SCHUR, INSTANCE_NUMBER
      INTEGER ICNTL(40), INFO(40), INFOG(40)
      REAL CNTL(5), RINFO(20), RINFOG(20)
      INTEGER, TARGET :: IRN(*), JCN(*), ELTPTR(*), ELTVAR(*)
      INTEGER, TARGET :: PERM_IN(*), IRN_loc(*), JCN_loc(*)
      INTEGER, TARGET :: LISTVAR_SCHUR(*)
      REAL, TARGET :: A(*), A_ELT(*), A_LOC(*), RHS(*)
      REAL, TARGET :: ROWSCA(*), COLSCA(*)
      REAL, TARGET :: SCHUR(*)
      INTEGER IRNhere, JCNhere, Ahere, ELTPTRhere, ELTVARhere,
     *        A_ELThere, PERM_INhere, RHShere, IRN_lochere,
     *        JCN_lochere, A_lochere, LISTVAR_SCHURhere,
     *        SCHURhere, COLSCAhere, ROWSCAhere
      INCLUDE 'smumps_struc.h'
      INCLUDE 'mpif.h'
      TYPE SMUMPS_STRUC_PTR
          TYPE (SMUMPS_STRUC), POINTER :: ptr
      END TYPE SMUMPS_STRUC_PTR
      TYPE (SMUMPS_STRUC), POINTER :: mumps_par
      TYPE (SMUMPS_STRUC_PTR), DIMENSION (:), POINTER, SAVE ::
     *  mumps_par_array
      TYPE (SMUMPS_STRUC_PTR), DIMENSION (:), POINTER ::
     *  mumps_par_array_bis
      INTEGER, SAVE :: SMUMPS_STRUC_ARRAY_SIZE = 0
      INTEGER, SAVE :: N_INSTANCES = 0
      INTEGER A_ELT_SIZE, I, Np, IERR
      INTEGER SMUMPS_STRUC_ARRAY_SIZE_INIT
      PARAMETER (SMUMPS_STRUC_ARRAY_SIZE_INIT=10)
      EXTERNAL SMUMPS_AFFECT_MAPPING,
     *         SMUMPS_AFFECT_NULLSPACE,
     *         SMUMPS_AFFECT_SYM_PERM,
     *         SMUMPS_AFFECT_UNS_PERM
      IF (JOB == -1) THEN
        DO I = 1, SMUMPS_STRUC_ARRAY_SIZE
          IF ( .NOT. ASSOCIATED(mumps_par_array(I)%PTR) ) GOTO 10
        END DO
        ALLOCATE( mumps_par_array_bis(SMUMPS_STRUC_ARRAY_SIZE +
     *  SMUMPS_STRUC_ARRAY_SIZE_INIT), stat=ierr)
        IF (ierr /= 0) THEN
          WRITE(*,*) ' ** Allocation Error 1 in SMUMPS_F77.'
          CALL SMUMPS_ABORT()
        END IF
        DO I = 1, SMUMPS_STRUC_ARRAY_SIZE
          mumps_par_array_bis(I)%PTR=>mumps_par_array(I)%PTR
        ENDDO
        IF (ASSOCIATED(mumps_par_array)) DEALLOCATE(mumps_par_array)
        mumps_par_array=>mumps_par_array_bis
        NULLIFY(mumps_par_array_bis)
        DO I = SMUMPS_STRUC_ARRAY_SIZE+1, SMUMPS_STRUC_ARRAY_SIZE +
     *  SMUMPS_STRUC_ARRAY_SIZE_INIT
          NULLIFY(mumps_par_array(I)%PTR)
        ENDDO
        I = SMUMPS_STRUC_ARRAY_SIZE+1
        SMUMPS_STRUC_ARRAY_SIZE = SMUMPS_STRUC_ARRAY_SIZE +
     *  SMUMPS_STRUC_ARRAY_SIZE_INIT
 10     CONTINUE
        INSTANCE_NUMBER = I
        N_INSTANCES = N_INSTANCES+1
        ALLOCATE( mumps_par_array(INSTANCE_NUMBER)%PTR,stat=IERR )
        IF (IERR /= 0) THEN
          WRITE(*,*) '** Allocation Error 2 in SMUMPS_F77.'
          CALL SMUMPS_ABORT()
        ENDIF
        mumps_par_array(INSTANCE_NUMBER)%PTR%INSTANCE_NUMBER =
     *  INSTANCE_NUMBER
      END IF
      IF ( INSTANCE_NUMBER .LE. 0 .OR. INSTANCE_NUMBER .GT.
     *     SMUMPS_STRUC_ARRAY_SIZE ) THEN
        WRITE(*,*) ' ** Instance Error 1 in SMUMPS_F77',
     *             INSTANCE_NUMBER
        CALL SMUMPS_ABORT()
      END IF
      IF ( .NOT. ASSOCIATED ( mumps_par_array(INSTANCE_NUMBER)%PTR ) )
     *  THEN
        WRITE(*,*) ' Instance Error 2 in SMUMPS_F77',
     *             INSTANCE_NUMBER
        CALL SMUMPS_ABORT()
      END IF
      mumps_par => mumps_par_array(INSTANCE_NUMBER)%PTR
      mumps_par%SYM = SYM
      mumps_par%PAR = PAR
      mumps_par%JOB = JOB
      mumps_par%N   = N
      mumps_par%NZ  = NZ
      mumps_par%NZ_loc  = NZ_loc
      mumps_par%SIZE_SCHUR  = SIZE_SCHUR
      mumps_par%NELT= NELT
      mumps_par%ICNTL(1:40)=ICNTL(1:40)
      mumps_par%CNTL(1:5)=CNTL(1:5)
      IF ( COMM_F77 .NE. -987654 ) THEN
        mumps_par%COMM = COMM_F77
      ELSE
        mumps_par%COMM = MPI_COMM_WORLD
      ENDIF
      IF ( IRNhere /= 0 ) mumps_par%IRN => IRN(1:NZ)
      IF ( JCNhere /= 0 ) mumps_par%JCN => JCN(1:NZ)
      IF ( Ahere /= 0 )   mumps_par%A   => A(1:NZ)
      IF ( IRN_lochere /= 0 ) mumps_par%IRN_loc => IRN_loc(1:NZ_loc)
      IF ( JCN_lochere /= 0 ) mumps_par%JCN_loc => JCN_loc(1:NZ_loc)
      IF ( A_lochere /= 0 )   mumps_par%A_loc   => A_loc(1:NZ_loc)
      IF ( ELTPTRhere /= 0 ) mumps_par%ELTPTR => ELTPTR(1:NELT+1)
      IF ( ELTVARhere /= 0 ) mumps_par%ELTVAR =>
     *   ELTVAR(1:ELTPTR(NELT+1)-1)
      IF ( A_ELThere /= 0 ) THEN
        A_ELT_SIZE = 0
        DO I = 1, NELT
          Np = ELTPTR(I+1) -ELTPTR(I)
          IF (SYM == 0) THEN
            A_ELT_SIZE = A_ELT_SIZE + Np * Np
          ELSE
            A_ELT_SIZE = A_ELT_SIZE + Np * ( Np + 1 ) / 2
          END IF
        END DO
        mumps_par%A_ELT => A_ELT(1:A_ELT_SIZE)
      END IF
      IF ( PERM_INhere /= 0) mumps_par%PERM_IN => PERM_IN(1:N)
      IF ( LISTVAR_SCHURhere /= 0)
     *   mumps_par%LISTVAR_SCHUR =>LISTVAR_SCHUR(1:SIZE_SCHUR)
      IF ( SCHURhere /= 0 )
     *   mumps_par%SCHUR => SCHUR(1:SIZE_SCHUR*SIZE_SCHUR)
      IF ( RHShere /= 0 ) mumps_par%RHS => RHS(1:N)
      IF ( COLSCAhere /= 0) mumps_par%COLSCA => COLSCA(1:N)
      IF ( ROWSCAhere /= 0) mumps_par%ROWSCA => ROWSCA(1:N)
      CALL SMUMPS( mumps_par )
      INFO(1:40)=mumps_par%INFO(1:40)
      INFOG(1:40)=mumps_par%INFOG(1:40)
      RINFO(1:20)=mumps_par%RINFO(1:20)
      RINFOG(1:20)=mumps_par%RINFOG(1:20)
      ICNTL(1:40) = mumps_par%ICNTL(1:40)
      CNTL(1:5) = mumps_par%CNTL(1:5)
      SYM = mumps_par%SYM
      PAR = mumps_par%PAR
      JOB = mumps_par%JOB
      N   = mumps_par%N
      NZ  = mumps_par%NZ
      NZ_loc  = mumps_par%NZ_loc
      SIZE_SCHUR  = mumps_par%SIZE_SCHUR
      NELT= mumps_par%NELT
      DEFICIENCY = mumps_par%Deficiency
      IF ( ASSOCIATED (mumps_par%MAPPING) ) THEN
         CALL SMUMPS_AFFECT_MAPPING(mumps_par%MAPPING(1))
      ELSE
         CALL SMUMPS_NULLIFY_C_MAPPING()
      ENDIF
      IF ( ASSOCIATED (mumps_par%NULL_SPACE) ) THEN
         CALL SMUMPS_AFFECT_NULLSPACE(mumps_par%NULL_SPACE(1,1))
      ELSE
         CALL SMUMPS_NULLIFY_C_NULLSPACE()
      ENDIF
      IF ( ASSOCIATED (mumps_par%SYM_PERM) ) THEN
         CALL SMUMPS_AFFECT_SYM_PERM(mumps_par%SYM_PERM(1))
      ELSE
         CALL SMUMPS_NULLIFY_C_SYM_PERM()
      ENDIF
      IF ( ASSOCIATED (mumps_par%UNS_PERM) ) THEN
         CALL SMUMPS_AFFECT_UNS_PERM(mumps_par%UNS_PERM(1))
      ELSE
         CALL SMUMPS_NULLIFY_C_UNS_PERM()
      ENDIF
      IF ( JOB == -2 ) THEN
         IF (ASSOCIATED(mumps_par_array(INSTANCE_NUMBER)%PTR))THEN
           DEALLOCATE(mumps_par_array(INSTANCE_NUMBER)%PTR)
           NULLIFY   (mumps_par_array(INSTANCE_NUMBER)%PTR)
           N_INSTANCES = N_INSTANCES - 1
           IF ( N_INSTANCES == 0 ) THEN
             DEALLOCATE(mumps_par_array)
             SMUMPS_STRUC_ARRAY_SIZE = 0
           END IF
         ELSE
           WRITE(*,*) "** Warning: instance already freed"
           WRITE(*,*) "            this should normally not happen."
         ENDIF
      END IF
      RETURN
      END SUBROUTINE SMUMPS_F77
