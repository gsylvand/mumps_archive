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
      SUBROUTINE DMUMPS( id )
      IMPLICIT NONE
      INTERFACE
      SUBROUTINE DMUMPS_26( id )
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC), TARGET :: id
      END SUBROUTINE DMUMPS_26
      SUBROUTINE DMUMPS_142( id )
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC), TARGET :: id
      END SUBROUTINE DMUMPS_142
      SUBROUTINE DMUMPS_301( id )
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC), TARGET :: id
      END SUBROUTINE DMUMPS_301
      SUBROUTINE DMUMPS_349(id, LP)
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC), TARGET :: id
      INTEGER  :: LP
      END SUBROUTINE DMUMPS_349
      END INTERFACE
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) :: id
      INTEGER JOBMIN, JOBMAX, OLDJOB
      INTEGER I, J, MP, LP, MPG
      LOGICAL LANAL, LFACTO, LSOLVE, PROK, FLAG, PROKG
      INTEGER COMM_SAVE
      INTEGER JOB, N, NZ, NELT
      JOB  = id%JOB
      N    = id%N
      NZ   = id%NZ
      NELT = id%NELT
      id%INFO(1) = 0
      id%INFO(2) = 0
      IF ( JOB .NE. -1 ) THEN
        LP      = id%ICNTL(1)
        MP      = id%ICNTL(2)
        MPG     = id%ICNTL(3)
        PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.3))
        PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
        IF ((id%MYID .eq. MASTER) .AND. PROK .AND. (id%ICNTL(5).EQ.0 ) ) 
     *       WRITE(MP,'(A,I4,I12,I15)') 
     *       'Entering driver (DMUMPS) WITH JOB, N, NZ =', JOB,N,NZ
        IF ((id%MYID .eq. MASTER).AND. PROK .AND. (id%ICNTL(5).EQ.1 ) ) 
     *       WRITE(MP,'(A,I4,I12,I15)') 
     *      'Entering driver (DMUMPS) WITH JOB, N, NELT =', JOB,N,NELT
      ELSE
        MPG = 0
        PROK = .FALSE.
        PROKG = .FALSE.
        LP = 6
        MP = 6
      END IF
      CALL MPI_INITIALIZED( FLAG, IERR )
      IF ( .NOT. FLAG ) THEN
        WRITE(LP,990)
 990  FORMAT(' Error in DMUMPS initialization: MPI is not running.')
        id%INFO(1) = -23
        id%INFO(2) =   0
        GOTO 500
      END IF
       COMM_SAVE = id%COMM
       CALL MPI_COMM_DUP( COMM_SAVE, id%COMM, IERR )
      CALL MPI_ALLREDUCE(JOB,JOBMIN,1,MPI_INTEGER,MPI_MAX,
     *                   id%COMM,IERR)
      CALL MPI_ALLREDUCE(JOB,JOBMAX,1,MPI_INTEGER,MPI_MIN,
     *                   id%COMM,IERR)
      IF ( JOBMIN .NE. JOBMAX ) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
      IF ( JOB .EQ. -1 ) THEN
        CALL DMUMPS_163( id )
        GOTO 500
      END IF
      IF ( JOB .EQ. -2 ) THEN
        id%KEEP(40)=-2
        CALL DMUMPS_136( id )
        GOTO 500
      END IF
      IF ((JOB.LT.1).OR.(JOB.GT.6)) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        IF ( id%ICNTL(18) .eq. 0 ) THEN
        IF ((N.LE.0).OR.((N+N+N)/3.NE.N)) THEN
          id%INFO(1) = -16
          id%INFO(2) = N
        END IF
        IF (id%ICNTL(5).EQ.0) THEN
          IF (NZ.LE.0) THEN
            id%INFO(1) = -2
            id%INFO(2) = NZ
          END IF
        ELSE
          IF (NELT.LE.0) THEN
            id%INFO(1) = -24
            id%INFO(2) = NELT
          END IF
        ENDIF
        END IF
        IF ( (id%KEEP(46).EQ.0).AND.(id%NPROCS.LE.1) ) 
     &     THEN
          id%INFO(1) = -21
          id%INFO(2) = id%NPROCS
        ENDIF
      END IF
      CALL DMUMPS_276( id%ICNTL,
     &                    id%INFO,
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 499
      LANAL  = .FALSE.
      LFACTO = .FALSE.
      LSOLVE = .FALSE.
      IF ((JOB.EQ.1).OR.(JOB.EQ.4).OR.
     *    (JOB.EQ.6))               LANAL  = .TRUE.
      IF ((JOB.EQ.2).OR.(JOB.EQ.4).OR.
     *    (JOB.EQ.5).OR.(JOB.EQ.6)) LFACTO = .TRUE.
      IF ((JOB.EQ.3).OR.(JOB.EQ.5).OR.
     *    (JOB.EQ.6))               LSOLVE = .TRUE.
      IF (MP.GT.0) CALL DMUMPS_349(id, MP)
      OLDJOB = id%KEEP( 40 )
      IF ( LANAL ) THEN
        IF ( OLDJOB .EQ. 0 .OR. OLDJOB .GT. 3 .OR. OLDJOB .LT. -1 ) THEN
          id%INFO(1) = -3
          id%INFO(2) = JOB
          GOTO 499
        END IF
        IF ( OLDJOB .GE. 2 ) THEN
          IF (ASSOCIATED(id%IS)) THEN
            DEALLOCATE  (id%IS)
            NULLIFY     (id%IS)
          END IF
          IF (ASSOCIATED(id%S)) THEN
            DEALLOCATE  (id%S)
            NULLIFY     (id%S)
          END IF
        END IF   
      END IF
      IF ( LFACTO ) THEN
         IF ( OLDJOB .LT. 1 .and. .NOT. LANAL ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
      IF ( LSOLVE ) THEN
         IF ( OLDJOB .LT. 2 .AND. .NOT. LFACTO ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
      IF (LANAL) THEN
        id%KEEP(40)=-1
        IF (id%MYID.EQ.MASTER) THEN
          IF (associated(id%IS1)) DEALLOCATE(id%IS1)
          IF ( id%ICNTL(5) .EQ. 0 ) THEN 
            IF ( id%ICNTL(6) .NE. 0 .and.
     *           id%ICNTL(7) .ne.1 ) THEN
              id%MAXIS1 = 11 * N
            ELSE
              id%MAXIS1 = 10 * N
            END IF
          ELSE
            id%MAXIS1 = 6 * N + 2 * NELT + 2
          ENDIF
          ALLOCATE( id%IS1(id%MAXIS1), stat=ierr )
          IF (ierr.gt.0) THEN
            id%INFO(1) = -13
            id%INFO(2) = id%MAXIS1
            IF ( LP .GT.0 ) 
     *      WRITE(LP,*) 'Problem in allocating work array for analysis.'
            GO TO 100
          END IF
          IF ( associated( id%PROCNODE ) )
     *          DEALLOCATE( id%PROCNODE )
          ALLOCATE( id%PROCNODE(id%N), stat=ierr )
          IF (ierr.gt.0) THEN
            id%INFO(1) = -13
            id%INFO(2) = id%N
            IF ( LP .GT. 0 ) THEN
              WRITE(LP,*) 'Problem in allocating work array PROCNODE'
            END IF
            GOTO 100
          END IF
          id%PROCNODE(1:id%N) = 0
          IF ( id%ICNTL(5) .NE. 0 ) THEN
            IF ( associated( id%ELTPROC ) )
     *            DEALLOCATE( id%ELTPROC )
            ALLOCATE( id%ELTPROC(id%NELT), stat=ierr )
            IF (ierr.gt.0) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%NELT
              IF ( LP .GT. 0 ) THEN
                WRITE(LP,*) 'Problem in allocating work array ELTPROC'
              END IF
              GOTO 100
            END IF
          END IF
          IF ( id%ICNTL(5) .EQ. 0 ) THEN
            IF ( id%ICNTL(18).EQ.0 ) THEN
             IF ( .not. ASSOCIATED( id%IRN ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
             ELSE IF ( SIZE( id%IRN ) < id%NZ ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
             ELSE IF ( .not. ASSOCIATED( id%JCN ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
             ELSE IF ( SIZE( id%JCN ) < id%NZ ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
             END IF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF (LP.GT.0) WRITE(LP,*)
     *           'Error in analysis: IRN/JCN badly allocated.'
            END IF
          ELSE
            IF ( .not. ASSOCIATED( id%ELTPTR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( SIZE( id%ELTPTR ) < id%NELT+1 ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( .not. ASSOCIATED( id%ELTVAR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
            ELSE 
              id%LELTVAR = id%ELTPTR( id%NELT+1 ) - 1
              IF ( SIZE( id%ELTVAR ) < id%LELTVAR ) THEN 
                id%INFO(1) = -22
                id%INFO(2) = 2
              ENDIF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF (LP.GT.0) WRITE(LP,*) 
     *           'Error in analysis: ELTPTR/ELTVAR badly allocated.'
            END IF
          ENDIF
 100    CONTINUE
        END IF
        CALL DMUMPS_276( id%ICNTL,
     *                    id%INFO,
     *                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
        CALL DMUMPS_26( id )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
        id%KEEP(40) = 1
      END IF
      IF (LFACTO) THEN
         id%KEEP(40) = 1
        IF ( id%MYID .EQ. MASTER ) THEN
          IF (id%KEEP(60).ne. 0) THEN
            IF ( .NOT. ASSOCIATED (id%SCHUR)) THEN
              id%INFO(1)=-22
              id%INFO(2)=9
            ELSE IF ( size(id%SCHUR) .NE.
     *                id%SIZE_SCHUR * id%SIZE_SCHUR ) THEN
              id%INFO(1)=-22
              id%INFO(2)=9
            END IF
          END IF
          IF ( id%KEEP(55) .EQ. 0 ) THEN
           IF ( id%KEEP(54).eq.0 ) THEN
            IF ( .not. ASSOCIATED( id%A ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            ELSE IF ( size( id%A ) < id%NZ ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            END IF
           END IF
          ELSE
            IF ( .not. ASSOCIATED( id%A_ELT ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 4
            ELSE 
              id%NA_ELT = 0
              IF ( id%KEEP(50) .EQ. 0 ) THEN
                DO I = 1,NELT
                  J = id%ELTPTR(I+1) - id%ELTPTR(I)
                  J = (J * J)
                  id%NA_ELT = id%NA_ELT + J
                ENDDO
              ELSE
                DO I = 1,NELT
                  J = id%ELTPTR(I+1) - id%ELTPTR(I)
                  J = (J * (J+1))/2
                  id%NA_ELT = id%NA_ELT + J
                ENDDO
              ENDIF
              IF ( size( id%A_ELT ) < id%NA_ELT ) THEN
                id%INFO( 1 ) = -22
                id%INFO( 2 ) = 4
              ENDIF
            END IF
          ENDIF
          id%KEEP(12) = id%ICNTL(14) 
          IF ( (id%KEEP(54).NE.0).AND. (id%KEEP(50).NE.1)
     *    .AND. (id%KEEP(12) .GT. 0) ) id%KEEP(12)= id%KEEP(12)*2
#ifndef try_null_space
          id%KEEP(19) = id%ICNTL(16)
          IF ( id%KEEP(19) .LT. 0
     *    .OR. id%KEEP(19).GT.10) THEN
#endif
              id%KEEP(19) = 0
#ifndef try_null_space
          END IF
#endif
#if defined(try_null_space_not_tested)
#else
          IF ( id%KEEP(19) > 2 .and.
     *         id%KEEP(19) < 9 ) THEN
            id%KEEP(19) = 0
          END IF
#endif
          IF ( id%KEEP(53) .LE. 0 .and.
     *         id%KEEP(19) .NE. 0 ) THEN
             id%KEEP(19) = 0
             IF ( MPG .GT. 0 ) THEN
             WRITE( MPG,'(A)') '** Warning:  null space disabled.'
             WRITE( MPG,'(A)') '** (incompatibility with analysis)'
             END IF
          END IF
          IF ( id%KEEP(19) .ne. 0 .and. id%KEEP(60) .ne. 0 ) THEN
             IF ( MPG .GT. 0 ) THEN
             WRITE( MPG,'(A)') '** Warning:  null space disabled.'
             WRITE( MPG,'(A)') '** (incompatibility with Schur)'
             END IF
            id%KEEP(19) = 0
          END IF
          id%KEEP(21) = MIN(id%ICNTL(17),N)
          id%KEEP( 52 ) = id%ICNTL( 8 )
          IF ( id%KEEP(52) .GT. 6 .or. id%KEEP(52).LT.-1)
     *    id%KEEP(52) = 0
          IF ( id%KEEP( 19 ) .ne. 0 .and. id%KEEP( 52 ).ne. 0 ) THEN
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with scaling)'
            END IF
            id%KEEP(52) = 0
          END IF
          IF ( id%KEEP(60) .ne. 0 .and. id%KEEP(52) .ne. 0 ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 .and. id%ICNTL(8) .ne. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with Schur)'
            END IF
          END IF
          IF (id%KEEP(54) .NE. 0 ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 .and. id%ICNTL(8) .ne. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with dist entry)'
            END IF
          END IF
          IF ( id%KEEP(50) .NE. 0 ) THEN
            IF ( id%KEEP(52).gt.1 .or.
     *           id%KEEP(52).lt.-1 ) THEN
              IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)')
     *  ' ** Warning: Scaling option n.a. for symmetric matrix'
              END IF
              id%KEEP(52) = 0
            END IF
          END IF
          IF (id%KEEP(55) .NE. 0 .AND. 
     *        ( id%KEEP(52) .gt. 0 ) ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)')
     *        ' ** (only user scaling av. for elt. entry)'
            END IF
          END IF
          IF ( id%KEEP(52) .eq. -1 ) THEN
            IF ( .not. ASSOCIATED( id%ROWSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( SIZE( id%ROWSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( .not. ASSOCIATED( id%COLSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            ELSE IF ( SIZE( id%COLSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            END IF
          END IF
          IF (id%KEEP(52).GT.0 .AND.
     *        id%KEEP(52) .LE.6) THEN
            IF ( associated(id%COLSCA))
     *             DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA))
     *             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=ierr)
            IF (ierr .GT.0) id%INFO(1)=-13
            ALLOCATE( id%ROWSCA(N), stat=ierr)
            IF (ierr .GT.0) id%INFO(1)=-13
          END IF
          IF (.NOT. associated(id%COLSCA)) THEN
            ALLOCATE( id%COLSCA(1), stat=ierr)
          END IF
          IF (ierr .GT.0) id%INFO(1)=-13
          IF (.NOT. associated(id%ROWSCA))
     *    ALLOCATE( id%ROWSCA(1), stat=ierr)
          IF (ierr .GT.0) id%INFO(1)=-13
          IF ( id%INFO(1) .eq. -13 ) THEN
            IF ( LP .GT. 0 )
     *      WRITE(LP,*) 'Problems in allocations before facto'
            GOTO 200
          END IF
 200      CONTINUE
        END IF
        CALL DMUMPS_276( id%ICNTL,
     *                      id%INFO,
     *                      id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
        CALL MPI_BCAST( id%KEEP(12), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( id%KEEP(19), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( id%KEEP(21), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        CALL DMUMPS_142(id)
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
        id%KEEP(40) = 2
      END IF
      IF (LSOLVE) THEN
        id%KEEP(40) = 2
        IF ( id%MYID .EQ. MASTER ) THEN
           IF ( .not. ASSOCIATED( id%RHS ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 7
           ELSE IF ( size( id%RHS ) < id%N ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 7
           END IF
        ELSE
           IF ( ASSOCIATED( id%RHS ) ) THEN
             id%INFO( 1 ) = -22
             id%INFO( 2 ) = 7
           END IF
        END IF
        CALL DMUMPS_276( id%ICNTL,
     *                      id%INFO,
     *                      id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
        CALL DMUMPS_301(id)
        IF (id%INFO(1).LT.0) GOTO 499
        id%KEEP(40) = 3
      ENDIF
      IF (MP.GT.0) CALL DMUMPS_349(id, MP)
      GOTO 500
  499 PROK  = ((id%ICNTL(1).GT.0).AND.
     *         (id%ICNTL(4).GE.1))
      IF (PROK) WRITE (id%ICNTL(1),99995) id%INFO(1)
      IF (PROK) WRITE (id%ICNTL(1),99994) id%INFO(2)
500   CONTINUE
      CALL DMUMPS_300( id%INFO, id%INFOG, id%COMM, id%MYID )
      IF (id%MYID.EQ.MASTER.and.MPG.GT.0.and.
     * id%INFOG(1).lt.0) THEN
        WRITE(MPG,'(A,I12)') ' On return from DMUMPS, INFOG(1)=',
     *      id%INFOG(1)
        WRITE(MPG,'(A,I12)') ' On return from DMUMPS, INFOG(2)=',
     *      id%INFOG(2)
      END IF
       CALL MPI_COMM_FREE( id%COMM, IERR )
       id%COMM = COMM_SAVE
      RETURN
99995 FORMAT (' ** ERROR RETURN ** FROM DMUMPS INFO(1)=', I3)
99994 FORMAT (' ** INFO(2)=', I10)
      END SUBROUTINE DMUMPS
      SUBROUTINE DMUMPS_300( INFO, INFOG, COMM, MYID )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INFO(2), INFOG(2), COMM, MYID
      INTEGER TMP(2), ROOT, IERR
      IF ( INFO(1) .ge. 0  .and. INFO(2) .ge. 0 ) THEN
        INFOG(1) = INFO(1)
        INFOG(2) = INFO(2)
      ELSE
        INFOG(1) = INFO(1)
        INFOG(2) = MYID
        CALL MPI_ALLREDUCE(INFOG,TMP,1,MPI_2INTEGER,
     *                     MPI_MINLOC,COMM,ierr )
        INFOG(2) = INFO(2)
        ROOT = TMP(2)
        CALL MPI_BCAST( INFOG(1), 1, MPI_INTEGER, ROOT, COMM, IERR )
        CALL MPI_BCAST( INFOG(2), 1, MPI_INTEGER, ROOT, COMM, IERR )
      END IF
      RETURN
      END SUBROUTINE DMUMPS_300
      SUBROUTINE DMUMPS_349(id, LP)
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER  :: LP
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LT.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      IF (ID%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19)
           IF ((ICNTL(6).EQ.5).OR.(ICNTL(6).EQ.6))THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           WRITE (LP,993) ICNTL(14)
#if defined(try_null_space)
           WRITE (LP,996) ICNTL(16)
#endif
         CASE(2);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (ICNTL(6).EQ.0)THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           WRITE (LP,993) ICNTL(14)
#if defined(try_null_space)
           WRITE (LP,996) ICNTL(16)
           WRITE (LP,994) ICNTL(17)
#endif
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,995) ICNTL(9),ICNTL(10),ICNTL(11)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (ICNTL(6).NE.0)THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF  
           WRITE (LP,995) ICNTL(9),ICNTL(10),ICNTL(11)
           WRITE (LP,993) ICNTL(14)
#if defined(try_null_space)
           WRITE (LP,996) ICNTL(16)
           WRITE (LP,994) ICNTL(17)
#endif
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19)
           IF ((ICNTL(6).EQ.5).OR.(ICNTL(6).EQ.6))THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF              
           IF (ICNTL(6).EQ.0)THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           WRITE (LP,993) ICNTL(14)
#if defined(try_null_space)
           WRITE (LP,996) ICNTL(16)
           WRITE (LP,994) ICNTL(17)
#endif
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),ICNTL(18),ICNTL(19)
           IF ((ICNTL(6).EQ.5).OR.(ICNTL(6).EQ.6))THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           IF (ICNTL(6).EQ.0)THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           WRITE (LP,995) ICNTL(9),ICNTL(10),ICNTL(11)
           WRITE (LP,993) ICNTL(14)
#if defined(try_null_space)
           WRITE (LP,996) ICNTL(16)
           WRITE (LP,994) ICNTL(17)
#endif
        END SELECT
      ENDIF
 980  FORMAT (/'***********CONTROL PARAMETERS (ICNTL)**************'/)
 990  FORMAT (
     1     'ICNTL(1)   Output stream for error messages        =',I10/
     2     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     3     'ICNTL(3)   Output stream for global information    =',I10/
     4     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     1     'ICNTL(5)   Matrix format                           =',I10/
     2     'ICNTL(6)   Maximum transversal                     =',I10/
     3     'ICNTL(7)   Ordering                                =',I10/
     4     'ICNTL(12)  Node level parallelism (0=on, 1=off)    =',I10/
     5     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     7     'ICNTL(18)  Distributed matrix                      =',I10/
     8     'ICNTL(19)  Schur option (  0=off,else=on )         =',I10)
 992  FORMAT (
     1     'ICNTL(8)   Scaling strategy                        =',I10)
 993  FORMAT (
     1     'ICNTL(14)  Percent of memory increase              =',I10)
#if defined(try_null_space)
 996  FORMAT (
     1     'ICNTL(16)  Null space functionality                =',I10) 
 994  FORMAT (
     1     'ICNTL(17)  Estimate of null space size             =',I10)
#endif
 995  FORMAT (
     1     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)       =',I10/    
     2     'ICNTL(10)  Max steps iterative refinement          =',I10/
     3     'ICNTL(11)  Error analysis ( 0= off, else=on)       =',I10)
      END SUBROUTINE DMUMPS_349
      SUBROUTINE DMUMPS_350(id, LP)
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER ::LP
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL, KEEP
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LT.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      KEEP=>id%KEEP
      IF (id%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),ICNTL(12),
     &          ICNTL(13),KEEP(54),KEEP(60)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6))THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
#if defined(try_null_space)
           WRITE (LP,996) KEEP(19)
#endif
         CASE(2);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
#if defined(try_null_space)
           WRITE (LP,996) KEEP(19)
           WRITE (LP,994) KEEP(21)
#endif
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4) 
           WRITE (LP,995) ICNTL(9),ICNTL(10),ICNTL(11)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).NE.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF  
           WRITE (LP,995) ICNTL(9),ICNTL(10),ICNTL(11)
           WRITE (LP,993) KEEP(12)
#if defined(try_null_space)
           WRITE (LP,996) KEEP(19)
           WRITE (LP,994) KEEP(21)
#endif
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),ICNTL(12),
     &          ICNTL(13),KEEP(54),KEEP(60)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF              
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
#if defined(try_null_space)
           WRITE (LP,996) KEEP(19)
           WRITE (LP,994) KEEP(21)
#endif
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),ICNTL(12),
     &          ICNTL(13),KEEP(54),KEEP(60)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,995) ICNTL(9),ICNTL(10),ICNTL(11)
           WRITE (LP,993) KEEP(12)
#if defined(try_null_space)
           WRITE (LP,996) KEEP(19)
           WRITE (LP,994) KEEP(21)
#endif
        END SELECT
      ENDIF
 980  FORMAT (/'******INTERNAL VALUE OF PARAMETERS (ICNTL/KEEP)****'/)
 990  FORMAT (
     1     'ICNTL(1)   Output stream for error messages        =',I10/
     2     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     3     'ICNTL(3)   Output stream for global information    =',I10/
     4     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     1     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     2     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     3     'ICNTL(7)   Ordering                                =',I10/
     4     'ICNTL(12)  Node level parallelism (0=on, 1=off)    =',I10/
     5     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     7     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     8     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10)
 992  FORMAT (
     1     'ICNTL(8)   Scaling strategy ( keep(52) )           =',I10)
 993  FORMAT (
     1     'ICNTL(14)  Percent of memory increase ( keep(12) ) =',I10)
#if defined(try_null_space)
 996  FORMAT (
     1     'ICNTL(16)  Null space functionality ( keep(19) )   =',I10) 
 994  FORMAT (
     1     'ICNTL(17)  Estimate of null space size ( keep(21) )=',I10)
#endif
 995  FORMAT (
     1     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)       =',I10/
     2     'ICNTL(10)  Max steps iterative refinement          =',I10/
     3     'ICNTL(11)  Error analysis ( 0= off, else=on)       =',I10)
      END SUBROUTINE DMUMPS_350
      SUBROUTINE DMUMPS_24( MYID, SLAVEF, N,
     *           PROCNODE, STEP, PTRAIW, PTRARW,
     *           KEEP, ICNTL, ma41_par )
      IMPLICIT NONE
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) :: ma41_par
      INTEGER MYID, N, SLAVEF
      INTEGER KEEP( 150 ), ICNTL( 40 )
      INTEGER PROCNODE( KEEP(28) ), STEP( N ),
     *        PTRAIW( N ), PTRARW( N )
      INTEGER DMUMPS_330, DMUMPS_275
      EXTERNAL DMUMPS_330, DMUMPS_275
      INTEGER I, IPTRI, IPTRR, NCOL, NROW, allocok
      INTEGER TYPE_PARALL, ITYPE, IRANK
      TYPE_PARALL = KEEP(46)
      KEEP(14) = 0
      KEEP(13) = 0
      DO I = 1, N
        ITYPE = DMUMPS_330( ABS(STEP(I)), PROCNODE, SLAVEF )
        IRANK = DMUMPS_275( ABS(STEP(I)), PROCNODE, SLAVEF )
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK = IRANK + 1
        END IF
        IF ( ITYPE .EQ. 2 ) THEN
          IF ( IRANK .NE. MYID ) THEN
            PTRARW( I ) = 0
          END IF
          KEEP(14) = KEEP(14) + 3 + PTRAIW( I ) + PTRARW( I )
          KEEP(13) = KEEP(13) + 1 + PTRAIW( I ) + PTRARW( I )
        ELSE IF ( ITYPE .EQ. 1 .AND.
     *            IRANK .EQ. MYID ) THEN
          KEEP( 14 ) = KEEP( 14 ) + 3 + PTRAIW( I ) + PTRARW( I )
          KEEP( 13 ) = KEEP( 13 ) + 1 + PTRAIW( I ) + PTRARW( I )
        ELSE IF ( ITYPE .EQ. 3 ) THEN
        END IF
      END DO
      IF ( ASSOCIATED( ma41_par%INTARR ) ) THEN
        DEALLOCATE( ma41_par%INTARR )
        NULLIFY( ma41_par%INTARR )
      END IF
      IF ( KEEP(14) > 0 ) THEN
      ALLOCATE( ma41_par%INTARR( KEEP(14) ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        ma41_par%INFO(1) = -13
        ma41_par%INFO(2) = KEEP(14)
        RETURN
      END IF
      ELSE
      ALLOCATE( ma41_par%INTARR( 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        ma41_par%INFO(1) = -13
        ma41_par%INFO(2) = 1
        RETURN
      END IF
      END IF
      IPTRI = 1
      IPTRR = 1
      DO I = 1, N
        ITYPE = DMUMPS_330( ABS(STEP(I)), PROCNODE, SLAVEF )
        IRANK = DMUMPS_275( ABS(STEP(I)), PROCNODE, SLAVEF )
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK =IRANK + 1
        END IF
        IF (
     *      ( ITYPE .eq. 2 .and.
     *        IRANK .eq. MYID )
     * .or.
     *      ( ITYPE .eq. 1 .and.
     *        IRANK .eq. MYID )
     *     )  THEN
          NCOL = PTRAIW( I )
          NROW = PTRARW( I )
          ma41_par%INTARR( IPTRI     ) = NCOL
          ma41_par%INTARR( IPTRI + 1 ) = -NROW
          ma41_par%INTARR( IPTRI + 2 ) = I
          PTRAIW( I ) = IPTRI
          PTRARW( I ) = IPTRR
          IPTRI = IPTRI + NCOL + NROW + 3
          IPTRR = IPTRR + NCOL + NROW + 1
        ELSE IF ( ITYPE .eq. 2 .and.
     *            IRANK .ne. MYID )
     *    THEN
          NCOL = PTRAIW( I )
          NROW = 0
          ma41_par%INTARR( IPTRI     ) = NCOL
          ma41_par%INTARR( IPTRI + 1 ) = -NROW
          ma41_par%INTARR( IPTRI + 2 ) = I
          PTRAIW( I ) = IPTRI
          PTRARW( I ) = IPTRR
          IPTRI = IPTRI + NCOL + NROW + 3
          IPTRR = IPTRR + NCOL + NROW + 1
        ELSE
          PTRAIW(I) = 0
          PTRARW(I) = 0
        END IF
      END DO
      IF ( IPTRI - 1 .NE. KEEP(14) ) THEN
        WRITE(*,*) 'Error 1 in anal_arrowheads'
        CALL DMUMPS_ABORT()
      END IF
      IF ( IPTRR - 1 .NE. KEEP(13) ) THEN
        WRITE(*,*) 'Error 2 in anal_arrowheads'
        CALL DMUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE DMUMPS_24
      SUBROUTINE DMUMPS_148(N, NZ, ASPK, 
     *   IRN, ICN, PERM,
     *   LSCAL,COLSCA,ROWSCA,
     *   MYID, SLAVEF, PROCNODE_STEPS, NBRECORDS,
     *   LP, COMM, root, KEEP, FILS, RG2L,
     *   INTARR, DBLARR, PTRAIW, PTRARW, FRERE_STEPS,
     *   STEP, A, LA )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER N,NZ, COMM, NBRECORDS
      INTEGER KEEP( 150 )
      DOUBLE PRECISION ASPK(NZ), COLSCA(*), ROWSCA(*)
      INTEGER IRN(NZ), ICN(NZ) 
      INTEGER PERM(N), PROCNODE_STEPS(KEEP(28))
      INTEGER RG2L( N ), FILS( N )
      INTEGER LP, SLAVEF, MYID
      LOGICAL LSCAL
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER LA
      INTEGER PTRAIW( N ), PTRARW( N ), FRERE_STEPS( KEEP(28) )
      INTEGER STEP(N)
      INTEGER INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION A( LA ), DBLARR(MAX(1,KEEP(13)))
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: BUFI
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BUFR
      INTEGER DMUMPS_275, DMUMPS_330, NUMROC
      EXTERNAL DMUMPS_275, DMUMPS_330, NUMROC
      DOUBLE PRECISION VAL
      INTEGER IOLD,JOLD,INEW,JNEW,ISEND,JSEND,DEST,I,K,IARR
      INTEGER IPOSROOT, JPOSROOT
      INTEGER IROW_GRID, JCOL_GRID
      INTEGER INODE
      INTEGER NBUFS
      INTEGER ARROW_ROOT, TAILLE
      INTEGER LOCAL_M, LOCAL_N, PTR_ROOT
      INTEGER I1, IA, JARR, ILOCROOT, JLOCROOT
      INTEGER IS1, ISHIFT, IIW, IS, IAS, IPROC
      INTEGER allocok
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      ARROW_ROOT = 0
      IF ( KEEP(46) .eq. 0 ) THEN
        NBUFS = SLAVEF
      ELSE
        NBUFS = SLAVEF - 1
        ALLOCATE( IW4( N, 2 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          WRITE(*,*) 'Error allocating IW4'
          CALL DMUMPS_ABORT()
        END IF
        DO I = 1, N
          I1 = PTRAIW( I )
          IA = PTRARW( I )
          IF ( IA .GT. 0 ) THEN
            DBLARR( IA ) = DBLE(ZERO)
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
            A( PTR_ROOT:LA ) = DBLE(ZERO)
          END IF
        END IF
      END IF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, NBUFS), stat=allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) 'Error allocating BUFI'
        CALL DMUMPS_ABORT()
      END IF
      ALLOCATE( BUFR( NBRECORDS, NBUFS ), stat=allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) 'Error allocating BUFR'
        CALL DMUMPS_ABORT()
      END IF
      DO I = 1, NBUFS
        BUFI( 1, I ) = 0
      ENDDO
      INODE = KEEP(38)
      I     = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = I
        INODE = FILS( INODE )
        I = I + 1
      END DO
      DO 120 K=1,NZ
        IOLD = IRN(K)
        JOLD = ICN(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     *                 .OR.(JOLD.LT.1) ) THEN
           GOTO 120
        END IF
        IF (LSCAL) THEN
          VAL = ASPK(K)*ROWSCA(IOLD)*COLSCA(JOLD)
        ELSE
          VAL = ASPK(K)
        ENDIF
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
        IF ( DMUMPS_330( ABS(STEP(IARR)),
     *       PROCNODE_STEPS, SLAVEF ) .EQ. 1 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = DMUMPS_275( ABS(STEP(IARR)),
     *             PROCNODE_STEPS, SLAVEF ) + 1
          ELSE
            DEST = DMUMPS_275( ABS(STEP(IARR)),
     *             PROCNODE_STEPS, SLAVEF )
          END IF
        ELSE IF ( DMUMPS_330( ABS(STEP(IARR)),
     *            PROCNODE_STEPS, SLAVEF )
     *            .EQ. 2 ) THEN
          IF ( ISEND .LT. 0  ) THEN
            DEST = -1
          ELSE
            IF ( KEEP( 46 ) .eq. 0 ) THEN
              DEST = DMUMPS_275( ABS(STEP(IARR)),
     *               PROCNODE_STEPS, SLAVEF ) + 1
            ELSE 
              DEST = DMUMPS_275( ABS(STEP(IARR)),
     *               PROCNODE_STEPS, SLAVEF )
            END IF
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L(JSEND)
            JPOSROOT = RG2L(IARR)
          ELSE
            IPOSROOT = RG2L( IARR )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = MOD( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = MOD( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * root%NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * root%NPCOL + JCOL_GRID
          END IF
        END IF
        IF ( DEST .eq. 0 .or.
     *     ( DEST .eq. -1 .and. KEEP( 46 ) .eq. 1 ) ) THEN
          IARR = ISEND
          JARR = JSEND
          IF ( DMUMPS_330( ABS(STEP(ABS(IARR))),
     *         PROCNODE_STEPS, SLAVEF ) .eq. 3 ) THEN
            ARROW_ROOT = ARROW_ROOT + 1
            IF ( IROW_GRID .EQ. root%MYROW .AND.
     *         JCOL_GRID .EQ. root%MYCOL ) THEN
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
              WRITE(*,*) MYID,':INTERNAL Error: root arrowhead '
              WRITE(*,*) MYID,':is not belonging to me. IARR,JARR='
     *        ,IARR,JARR
              CALL DMUMPS_ABORT()
            END IF
          ELSE IF ( IARR .GE. 0 ) THEN
            IF ( IARR .eq. JARR ) THEN
              IA = PTRARW( IARR )
              DBLARR( IA ) = DBLARR( IA ) + VAL
            ELSE
              IS1 =  PTRAIW(IARR)
              ISHIFT      = INTARR(IS1) + IW4(IARR,2)
              IW4(IARR,2) = IW4(IARR,2) - 1
              IIW         = IS1 + ISHIFT + 2
              INTARR(IIW)     = JARR
              IS          = PTRARW(IARR)
              IAS         = IS + ISHIFT
              DBLARR(IAS) = VAL
            END IF
          ELSE
            IARR = -IARR
            ISHIFT      = PTRAIW(IARR)+IW4(IARR,1)+2
            INTARR(ISHIFT)  = JARR
            IAS         = PTRARW(IARR)+IW4(IARR,1)
            IW4(IARR,1) = IW4(IARR,1) - 1
            DBLARR(IAS)      = VAL
            IPROC = DMUMPS_275( ABS(STEP(IARR)),
     *              PROCNODE_STEPS,SLAVEF )
            IF ( KEEP(50) .NE. 0 .AND.
     *           IW4(IARR,1) .EQ. 0 .AND.
     *           IPROC .EQ. MYID
     *           .AND. STEP( IARR) > 0 ) THEN
              TAILLE = INTARR( PTRAIW(IARR) )
              CALL DMUMPS_310( N, PERM,
     *             INTARR( PTRAIW(IARR) + 3 ),
     *             DBLARR( PTRARW(IARR) + 1 ),
     *             TAILLE, MYID )
            END IF
          ENDIF
        END IF
        IF ( DEST .NE. 0 ) THEN
         CALL DMUMPS_34( ISEND, JSEND, VAL,
     *    DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     *    LP, COMM, KEEP(46) )
        END IF
  120 CONTINUE
      DEST = -2
      KEEP(49) = ARROW_ROOT
      CALL DMUMPS_34( ISEND, JSEND, VAL,
     *   DEST, BUFI, BUFR, NBRECORDS, NBUFS,
     *   LP, COMM, KEEP( 46 ) )
      IF ( KEEP( 46 ) .NE. 0 ) DEALLOCATE( IW4 )
      DEALLOCATE( BUFI )
      DEALLOCATE( BUFR )
      RETURN
      END SUBROUTINE DMUMPS_148
      SUBROUTINE DMUMPS_34(ISEND, JSEND, VAL,
     *   DEST, BUFI, BUFR, NBRECORDS, NBUFS, LP, COMM,
     *   TYPE_PARALL )
      IMPLICIT NONE
      INTEGER ISEND, JSEND, DEST, NBUFS, NBRECORDS, TYPE_PARALL
      INTEGER BUFI( NBRECORDS * 2 + 1, NBUFS )
      DOUBLE PRECISION BUFR( NBRECORDS, NBUFS )
      INTEGER COMM
      INTEGER LP
      DOUBLE PRECISION VAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER ISLAVE, IBEG,IEND, IERR
      INTEGER TAILLE_SENDI, TAILLE_SENDR, IREQ
      IF (DEST.EQ.-2) THEN
        DO ISLAVE = 1,NBUFS 
          TAILLE_SENDI = BUFI(1,ISLAVE) * 2 + 1
          TAILLE_SENDR = BUFI(1,ISLAVE)
          BUFI(1,ISLAVE) = - BUFI(1,ISLAVE)
          CALL MPI_SEND(BUFI(1,ISLAVE),TAILLE_SENDI,
     *                   MPI_INTEGER,
     *                   ISLAVE, ARROWHEAD, COMM, IERR )
          IF ( TAILLE_SENDR .NE. 0 ) THEN
            CALL MPI_SEND( BUFR(1,ISLAVE), TAILLE_SENDR,
     *                     MPI_DOUBLE_PRECISION, ISLAVE,
     *                     ARROWHEAD, COMM, IERR )
          END IF
        ENDDO
        GOTO 500
      ENDIF
      IF (DEST.EQ.-1) THEN
       IBEG = 1
       IEND = NBUFS
      ELSE
       IBEG = DEST
       IEND = DEST
      ENDIF
      DO ISLAVE = IBEG, IEND
         IF (BUFI(1,ISLAVE)+1.GT.NBRECORDS) THEN
          TAILLE_SENDI = BUFI(1,ISLAVE) * 2 + 1
          TAILLE_SENDR = BUFI(1,ISLAVE)
          CALL MPI_SEND(BUFI(1,ISLAVE),TAILLE_SENDI,
     *                   MPI_INTEGER,
     *                   ISLAVE, ARROWHEAD, COMM, IERR )
          CALL MPI_SEND( BUFR(1,ISLAVE), TAILLE_SENDR,
     *                   MPI_DOUBLE_PRECISION, ISLAVE,
     *                   ARROWHEAD, COMM, IERR )
          BUFI(1,ISLAVE) = 0
         ENDIF
         IREQ = BUFI(1,ISLAVE) + 1
         BUFI(1,ISLAVE) = IREQ
         BUFI( IREQ * 2, ISLAVE )     = ISEND
         BUFI( IREQ * 2 + 1, ISLAVE ) = JSEND
         BUFR( IREQ, ISLAVE )         = VAL
      ENDDO
500   CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_34
      SUBROUTINE DMUMPS_310( N, PERM, 
     *            INTLIST, DBLLIST, TAILLE, MYID )
      IMPLICIT NONE
      INTEGER N, TAILLE, MYID
      INTEGER PERM( N ) 
      INTEGER INTLIST( TAILLE )
      DOUBLE PRECISION DBLLIST( TAILLE )
      INTEGER I, SWAPI
      DOUBLE PRECISION SWAPR
      LOGICAL SORTED
      SORTED = .FALSE.
      DO WHILE ( .NOT. SORTED )
        SORTED = .TRUE.
        DO I = 1, TAILLE - 1
          IF ( PERM( INTLIST( I ) ) .GT. PERM( INTLIST( I + 1 ) ) ) THEN
            SORTED = .FALSE.
            SWAPI            = INTLIST( I )
            INTLIST( I )     = INTLIST( I + 1 )
            INTLIST( I + 1 ) = SWAPI
            SWAPR            = DBLLIST( I )
            DBLLIST( I )     = DBLLIST( I + 1 )
            DBLLIST( I + 1 ) = SWAPR
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_310
      SUBROUTINE DMUMPS_145(  N,
     *    DBLARR, LDBLARR, INTARR, LINTARR, PTRAIW, PTRARW, 
     *    KEEP, MYID,  COMM, NBRECORDS,
     *    A, LA, root,
     *    PROCNODE_STEPS,
     *    SLAVEF, PERM, FRERE_STEPS, STEP, INFO1, INFO2
     *   )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INTEGER N, MYID, LDBLARR, LINTARR,
     *        COMM
      INTEGER INTARR(LINTARR) 
      INTEGER PTRAIW(N), PTRARW(N) 
      INTEGER   KEEP(150)
      INTEGER LA
      INTEGER PROCNODE_STEPS( KEEP(28) ), PERM( N )
      INTEGER SLAVEF, NBRECORDS
      DOUBLE PRECISION A( LA )
      INTEGER INFO1, INFO2
      DOUBLE PRECISION DBLARR(LDBLARR)
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER, POINTER, DIMENSION(:) :: BUFI
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: BUFR
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      LOGICAL FINI 
      INTEGER IREC, NB_REC, IARR, JARR, IA, I1, I, ALLOCOK
      INTEGER IS, IS1, ISHIFT, IIW, IAS
      INTEGER LOCAL_M, LOCAL_N, ILOCROOT, JLOCROOT, 
     *        IPOSROOT, JPOSROOT, TAILLE,
     *        IPROC
      INTEGER FRERE_STEPS( KEEP(28) ), STEP(N)
      INTEGER ARROW_ROOT, TYPE_PARALL, PTR_ROOT
      INTEGER DMUMPS_330, DMUMPS_275
      EXTERNAL DMUMPS_330, DMUMPS_275
      DOUBLE PRECISION VAL
      DOUBLE PRECISION  ZERO
      PARAMETER( ZERO = 0.0D0 )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MASTER
      PARAMETER(MASTER=0)
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR
      INTEGER NUMROC
      EXTERNAL NUMROC
      TYPE_PARALL = KEEP(46)
      ARROW_ROOT=0
      ALLOCATE( BUFI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = NBRECORDS * 2 + 1
        WRITE(*,*) MYID,': Could not allocate BUFI: goto 500'
        GOTO 500
      END IF
      ALLOCATE( BUFR( NBRECORDS )        , stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = NBRECORDS
        WRITE(*,*) MYID,': Could not allocate BUFR: goto 500'
        GOTO 500
      END IF
      ALLOCATE( IW4(N,2), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = 2 * N
        WRITE(*,*) MYID,': Could not allocate IW4: goto 500'
        GOTO 500
      END IF
      IF ( KEEP(38).NE.0) THEN
         LOCAL_M = NUMROC( root%ROOT_SIZE, root%MBLOCK,
     *             root%MYROW, 0, root%NPROW )
         LOCAL_M = MAX( 1, LOCAL_M )
         LOCAL_N = NUMROC( root%ROOT_SIZE, root%NBLOCK,
     *             root%MYCOL, 0, root%NPCOL )
         PTR_ROOT = LA - LOCAL_M * LOCAL_N + 1
         IF ( PTR_ROOT .LE. LA ) THEN
           A( PTR_ROOT:LA ) = DBLE(ZERO)
         END IF
      END IF
      FINI = .FALSE.
      DO I=1,N
       I1 = PTRAIW(I)
       IA = PTRARW(I)
       IF (IA.GT.0) THEN
        DBLARR(IA) = DBLE(ZERO)
        IW4(I,1) = INTARR(I1)
        IW4(I,2) = -INTARR(I1+1)
        INTARR(I1+2)=I
       ENDIF
      ENDDO
      DO WHILE (.NOT.FINI) 
       CALL MPI_RECV( BUFI(1), 2*NBRECORDS+1, 
     *                MPI_INTEGER, MASTER, 
     *                ARROWHEAD,
     *                COMM, STATUS, IERR )
       NB_REC = BUFI(1)
       IF (NB_REC.LE.0) THEN
         FINI = .TRUE.
         NB_REC = -NB_REC 
       ENDIF
       IF (NB_REC.EQ.0) EXIT
       CALL MPI_RECV( BUFR(1), NBRECORDS, MPI_DOUBLE_PRECISION,
     *                  MASTER, ARROWHEAD,
     *                COMM, STATUS, IERR )
       DO IREC=1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        IF ( DMUMPS_330( ABS(STEP(ABS(IARR))),
     *       PROCNODE_STEPS, SLAVEF ) .eq. 3 ) THEN
          ARROW_ROOT = ARROW_ROOT + 1
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
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
           IPROC = DMUMPS_275( ABS(STEP(IARR)),
     *             PROCNODE_STEPS,SLAVEF )
           IF ( TYPE_PARALL .eq. 0 ) THEN
             IPROC = IPROC + 1
           END IF 
           IF ( KEEP(50) .NE. 0 .AND.
     *          IW4(IARR,1) .EQ. 0 .AND.
     *          IPROC .EQ. MYID
     *          .AND. STEP(IARR) > 0 ) THEN
             TAILLE = INTARR( PTRAIW(IARR) )
             CALL DMUMPS_310( N, PERM,
     *            INTARR( PTRAIW(IARR) + 3 ),
     *            DBLARR( PTRARW(IARR) + 1 ),
     *            TAILLE, MYID )
           END IF
        ENDIF
       ENDDO
      END DO
      DEALLOCATE( BUFI )
      DEALLOCATE( BUFR )
      DEALLOCATE( IW4 )
 500  CONTINUE
      KEEP(49) = ARROW_ROOT
      RETURN 
      END SUBROUTINE DMUMPS_145
      SUBROUTINE DMUMPS_266( MYID, BUFR, LBUFR, 
     *     LBUFR_BYTES,
     *     IWPOS, IWPOSCB,
     *     IPTRLU, LRLU, LRLUS,
     *     TNBPROCFILS, N, IW, LIW, A, LA, NIRBDU,
     *     PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, COMP,
     *     KEEP, ITLOC,
     *     IFLAG, IERROR )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER MYID
      INTEGER KEEP(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, LIW, LA
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER NIRBDU
      INTEGER PTRIST(KEEP(28)), STEP(N), 
     * PIMASTER(KEEP(28)), PAMASTER(KEEP(28)),
     * PTRAST(KEEP(28)), TNBPROCFILS( KEEP(28) ), ITLOC( N )
      INTEGER COMP, IFLAG, IERROR
      INTEGER INODE, NBPROCFILS, NCOL, NROW, NASS, NSLAVES
      INTEGER NSLAVES_RECU, NFRONT
      INTEGER LREQ, LREQCB
      DOUBLE PRECISION FLOP1
      INODE = BUFR( 1 )
      NBPROCFILS = BUFR( 2 )
      NROW = BUFR( 3 )
      NCOL = BUFR( 4 )
      NASS = BUFR( 5 )
      NFRONT = BUFR( 6 )
      NSLAVES_RECU = BUFR( 7 )
      IF ( KEEP(50) .eq. 0 ) THEN
         FLOP1 = DBLE( NASS * NROW ) +
     *     DBLE(NROW*NASS)*DBLE(2*NCOL-NASS-1)
      ELSE
         FLOP1 = DBLE( NASS ) * DBLE( NROW )
     *            * ( 2 * NCOL - NROW - NASS + 1)
      END IF
      CALL DMUMPS_190(.TRUE.,FLOP1, KEEP)
      IF ( KEEP(50) .eq. 0 ) THEN
        NSLAVES = NSLAVES_RECU
      ELSE
        NSLAVES = NSLAVES_RECU + 2
      END IF
      LREQ   = NROW + NCOL + 6 + NSLAVES
      LREQCB = NCOL * NROW
      CALL DMUMPS_22(.TRUE.,
     *   MYID,N, KEEP, IW, LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST, STEP, PIMASTER,PAMASTER, ITLOC,
     *   LREQ, LREQCB,
     *   COMP, LRLUS, IFLAG, IERROR )
      IF ( IFLAG .LT. 0 ) RETURN
      PTRIST(STEP(INODE)) = IWPOSCB + 1
      PTRAST(STEP(INODE)) = IPTRLU  + 1
      IW( IWPOSCB + 1 ) = NCOL
      IW( IWPOSCB + 2 ) = - NASS
      IW( IWPOSCB + 3 ) = NROW
      IW( IWPOSCB + 4 ) = 0
      IW( IWPOSCB + 5 ) = NASS
      IW( IWPOSCB + 6 ) = NSLAVES
      IW( IWPOSCB + 7+NSLAVES : IWPOSCB + 6+NSLAVES + NROW + NCOL )
     *= BUFR( 8 + NSLAVES_RECU : 7 + NSLAVES_RECU + NROW + NCOL )
      IF ( KEEP(50) .eq. 0 ) THEN
        IF( NSLAVES.NE. NSLAVES_RECU )THEN
          WRITE(*,*) 'error in process_bande'
          stop
        END IF
        IW( IWPOSCB + 7: IWPOSCB+6+NSLAVES ) =
     *       BUFR( 8: 7 + NSLAVES_RECU )
      ELSE
        IW( IWPOSCB + 7 ) = 0
        IW( IWPOSCB + 8 ) = NFRONT
        IW( IWPOSCB + 9: IWPOSCB+8+NSLAVES_RECU ) =
     *       BUFR( 8: 7 + NSLAVES_RECU )
      END IF
      TNBPROCFILS(STEP( INODE )) = NBPROCFILS
      RETURN
      END SUBROUTINE DMUMPS_266
      SUBROUTINE DMUMPS_163( id )
      USE DMUMPS_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) id
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      INTEGER color
      CALL MPI_COMM_SIZE(id%COMM, id%NPROCS, IERR )
      CALL MPI_COMM_RANK(id%COMM, id%MYID, IERR )
      IF ( id%PAR .eq. 0 ) THEN
        IF ( id%MYID .eq. MASTER ) THEN
          color = MPI_UNDEFINED
        ELSE
          color = 0
        END IF
        CALL MPI_COMM_SPLIT( id%COMM, color, 0,
     *                       id%COMM_NODES, ierr )
        id%NSLAVES = id%NPROCS - 1
      ELSE
        CALL MPI_COMM_DUP( id%COMM, id%COMM_NODES, IERR )
        id%NSLAVES = id%NPROCS
      END IF
      IF (id%PAR .ne. 0 .or. id%MYID .NE. MASTER) THEN
        CALL MPI_COMM_DUP( id%COMM_NODES, id%COMM_LOAD, IERR )
      ENDIF
      CALL DMUMPS_20( id%NSLAVES,
     *    id%CNTL, id%ICNTL,
     *    id%KEEP, id%RINFO, id%RINFOG,
     *    id%SYM, id%PAR )
      CALL DMUMPS_61( id%KEEP( 34 ) )
      NULLIFY(id%BUFR)
      id%MAXS  = 0
      id%MAXIS = 0
      id%MAXIS1 = 0
      id%INST_Number = -1
      NULLIFY(id%IRN)
      NULLIFY(id%JCN)
      NULLIFY(id%A)
      NULLIFY(id%IRN_LOC)
      NULLIFY(id%JCN_LOC)
      NULLIFY(id%A_loc)
      NULLIFY(id%MAPPING)
      NULLIFY(id%RHS)
      NULLIFY(id%COLSCA)
      NULLIFY(id%ROWSCA)
      NULLIFY(id%PERM_IN)
      NULLIFY(id%IS)
      NULLIFY(id%IS1)
      NULLIFY(id%STEP)
      NULLIFY(id%NE_STEPS)
      NULLIFY(id%ND_STEPS)
      NULLIFY(id%FRERE_STEPS)
      NULLIFY(id%SYM_PERM)
      NULLIFY(id%UNS_PERM)
      NULLIFY(id%FILS)
      NULLIFY(id%PTRAR)
      NULLIFY(id%FRTPTR)
      NULLIFY(id%FRTELT)
      NULLIFY(id%NA)
      NULLIFY(id%PROCNODE_STEPS)
      NULLIFY(id%MEM_SUBTREE)
      id%LNA=0
      NULLIFY(id%S)
      NULLIFY(id%PROCNODE)
      NULLIFY(id%POIDS)
      NULLIFY(id%PTLUST_S)
      NULLIFY(id%INTARR) 
      NULLIFY(id%DBLARR)
      NULLIFY(id%NULL_SPACE)
      NULLIFY(id%root%RG2L_ROW)
      NULLIFY(id%root%RG2L_COL)
      NULLIFY(id%root%IPIV)
      NULLIFY(id%ELTPTR)
      NULLIFY(id%ELTVAR)
      NULLIFY(id%A_ELT)
      NULLIFY(id%ELTPROC)
      id%SIZE_SCHUR = 0
      NULLIFY( id%LISTVAR_SCHUR )
      NULLIFY( id%SCHUR )
      NULLIFY(id%ISTEP_TO_INIV2)
      NULLIFY(id%FUTURE_NIV2)
      NULLIFY(id%TAB_POS_IN_PERE)
      NULLIFY(id%CANDIDATES)
      NULLIFY(id%root%QR_TAU)
      NULLIFY(id%root%GROW)
      NULLIFY(id%root%GSIN)
      NULLIFY(id%root%GCOS)
      NULLIFY(id%root%PERM_ROW)
      NULLIFY(id%root%PERM_COL)
      NULLIFY(id%root%ELROW)
      NULLIFY(id%root%EUROW)
      NULLIFY(id%root%PTREL)
      NULLIFY(id%root%PTREU)
      NULLIFY(id%root%ELELG)
      NULLIFY(id%root%EUELG)
      NULLIFY(id%root%DL)
      NULLIFY(id%MEM_DIST)
      id%Deficiency = 0
      id%root%LPIV = -1
      id%root%yes  = .FALSE.
      id%root%gridinit_done  = .FALSE.
        IF ( id%KEEP( 46 ) .ne. 0  .OR.
     *     id%MYID .ne. MASTER ) THEN
          CALL MPI_COMM_RANK
     *         (id%COMM_NODES, id%MYID_NODES, IERR )
        ENDIF
      RETURN
      END SUBROUTINE DMUMPS_163
      SUBROUTINE DMUMPS_252( COMM_LOAD, ASS_IRECV,
     *    N, INODE, IW, LIW, A, LA, IFLAG,
     *    IERROR, ND, 
     *    NIRBDU, FILS, FRERE, MAXFRW, root,
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
      USE DMUMPS_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER IZERO 
      PARAMETER (IZERO=0)
      INTEGER ASS_IRECV, COMM_LOAD
      INTEGER N,LIW,LA,NIRBDU,NSTEPS
      INTEGER KEEP(150), ICNTL(40)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     *        LRLU, IPTRLU,IWPOS, LRLUS,
     *        POSFAC, IWPOSCB, COMP, IERR_MPI
      LOGICAL SON_LEVEL2
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER NBPROCFILS(KEEP(28)) 
      INTEGER NSTK_S(KEEP(28)),PROCNODE_STEPS(KEEP(28))
      INTEGER IPOOL( LPOOL )
      INTEGER BUFR( LBUFR )
      INTEGER IDUMMY(1)
      INTEGER IW(LIW), ITLOC(N),
     *        PTRARW(N), PTRAIW(N), ND(KEEP(28)), PERM(N), 
     *        FILS(N), FRERE(KEEP(28)),
     *        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28)),
     *        STEP(N), 
     * PIMASTER(KEEP(28)),
     * PAMASTER(KEEP(28))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER      INTARR(MAX(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(MAX(1,KEEP(13)))
      INTEGER DMUMPS_275, DMUMPS_330
      EXTERNAL DMUMPS_275, DMUMPS_330
      INTEGER LP, HS, HF
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ,LAELL
      INTEGER SIZFI, SIZFR, NCB
      INTEGER LAPOS2,J1,J2
      INTEGER NCOL, NROW, NCOLS, NROWS, LDA_SON
      INTEGER NELIM,JJ,JJ1,JJ2,J3,
     *        IBROT,IORG
      INTEGER IACHK,JPOS,ICT11
      INTEGER JK,IJROW,NBCOL,ICT13,NUMORG,IOLDPS,J4
      INTEGER APOS, APOS2, AINPUT, POSELT, POSEL1, ICT12
      INTEGER NSLAVES, NSLSON, NPIVS
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE
      INTEGER INDX, FIRST_INDEX, SHIFT_INDEX
      LOGICAL FLAG, BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER NCBSON
      LOGICAL FREE, SAME_PROC
      INTRINSIC REAL
      DOUBLE PRECISION ZERO
      DATA ZERO /0.0D0/
      INTEGER NELT, LPTRAR
      NELT = 1
      LPTRAR = N
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      LEVEL = DMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF (LEVEL.NE.1) THEN 
       write(6,*) 'Error1 in mpi51f_niv1 '
       CALL DMUMPS_ABORT()
      ENDIF
      NSLAVES = 0
      HF = 6 + NSLAVES
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      ENDDO
      NSTEPS = NSTEPS + 1
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .NE. 0) THEN
         DO WHILE (ISON .GT. 0)
           NUMSTK = NUMSTK + 1
           NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1)
           ISON = FRERE(STEP(ISON))
         ENDDO
      ENDIF
      NFRONT = ND(STEP(INODE)) + NASS
      NASS1 = NASS + NUMORG
      LREQ = 2 * NFRONT + HF
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL DMUMPS_94(N, KEEP(28),
     *        IW, LIW, A, LA, NIRBDU,
     *        LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
          COMP = COMP+1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          ENDIF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      ENDIF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .TRUE.
      IF (KEEP(50).EQ.0) THEN
        CALL  DMUMPS_81(MYID, INODE, N, IOLDPS, HF, NFRONT, 
     *        NFRONT_EFF,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     *        INTARR, ITLOC, FILS, FRERE,
     *        SON_LEVEL2, NIV1, NBPROCFILS, KEEP, IFLAG)
      ELSE
        CALL DMUMPS_86( MYID, INODE, N, IOLDPS, HF,
     *        NFRONT, NFRONT_EFF, PERM,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     *        INTARR, ITLOC, FILS, FRERE,
     *        SON_LEVEL2, NIV1, NBPROCFILS, KEEP, IFLAG)
        IF (IFLAG.LT.0) GOTO 300
      ENDIF
      IF (NFRONT_EFF.NE.NFRONT) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
        IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF))
        NFRONT = NFRONT_EFF
        ELSE
         Write(*,*) ' ERROR 1 during ass_niv1'
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
          CALL DMUMPS_94
     *        (N, KEEP(28), IW, LIW, A, LA, NIRBDU, LRLU, IPTRLU,
     *         IWPOS, IWPOSCB, PTRIST, PTRAST, STEP, PIMASTER,
     *         PAMASTER,ITLOC)
          COMP = COMP + 1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          ENDIF
        ENDIF
      ENDIF
      LRLU = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL
      CALL DMUMPS_471(.FALSE.,LA-LRLUS,0,LAELL,KEEP)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT +LAELL - 1
      A(POSELT:LAPOS2) = DBLE(ZERO)
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL - 1
        A(POSELT:LAPOS2) = DBLE(ZERO)
      ELSE
        APOS = POSELT
        DO JJ = 0, NFRONT - 1
          A(APOS:APOS+JJ) = DBLE(ZERO)
          APOS = APOS + NFRONT
        END DO
      END IF
#endif
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
          OPASSW = OPASSW + DBLE(SIZFR)
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
            CALL DMUMPS_178(
     *           A( PTRAST(STEP( INODE )) ), NFRONT, NASS1,
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
              ENDDO
             ELSE
              J2 = J1 + LSTK - 1
              J3 = J1 + NELIM
              DO JJ = J3, J2
               IW(JJ) = IW(JJ - NROWS)
              ENDDO
              IF (NELIM .NE. 0) THEN
                J3 = J3 - 1
                DO JJ = J1, J3
                 JPOS = IW(JJ) + ICT11
                 IW(JJ) = IW(JPOS)
                ENDDO
              ENDIF
             ENDIF
           ENDIF
           IF ( SAME_PROC ) THEN 
               PTRIST(STEP(ISON))   = -99999999
           ELSE
               PIMASTER(STEP( ISON )) = -99999999
           ENDIF
      CALL DMUMPS_152( MYID, N, ISTCHK,
     *     PAMASTER(STEP(ISON)),
     *     IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *     IWPOSCB, LA, KEEP
     *     )
          ELSE
           PDEST = ISTCHK + 6
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL DMUMPS_49( 
     &                KEEP, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
              CALL DMUMPS_210( COMM_LOAD, ASS_IRECV, 
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
     *         INTARR, DBLARR, ND, FRERE,
     *         LPTRAR, NELT, IW, IW, 
     *
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 500
             ENDIF
           ENDDO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6
            CALL  DMUMPS_71( INODE, NFRONT, 
     *       NASS1, ISON, MYID,
     *       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     *       COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     *       KEEP, STEP, N, 
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *        )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL DMUMPS_329( 
     *         COMM_LOAD, ASS_IRECV,
     *         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *         MPI_ANY_SOURCE, MPI_ANY_TAG,
     *         STATUS,
     *         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *         IWPOS, IWPOSCB, IPTRLU,
     *         LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *         PTRIST, PTLUST_S,
     *         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *         IFLAG, IERROR, COMM,
     *         NBPROCFILS,
     *         IPOOL, LPOOL, LEAF,
     *         NBFIN, MYID, SLAVEF,
     *         root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *         INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *         LPTRAR, NELT, IW, IW,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
               IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           ENDDO
           IF (IERR .EQ. -2) GOTO 290
           ENDIF
          ENDIF
  210   ISON = FRERE(STEP(ISON))
  220 CONTINUE
      END IF
      POSELT = PTRAST(STEP(INODE))
      IBROT = INODE
      DO 260 IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        IBROT = FILS(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        IJROW = INTARR(J1)
        ICT12 = POSELT - NFRONT + IJROW - 1
        DO 240 JJ = J1, J2
           APOS2 = ICT12 + INTARR(JJ) * NFRONT
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
          AINPUT = AINPUT + 1
  240   CONTINUE
        IF (J3 .LE. J4) THEN
          ICT13 = POSELT + (IJROW - 1) * NFRONT
          NBCOL = J4 - J3 + 1
          DO 250 JJ = 1, NBCOL
            JJ1 = ICT13 + INTARR(J3 + JJ - 1) - 1
            A(JJ1) = A(JJ1) + DBLARR(AINPUT + JJ - 1)
  250     CONTINUE
        ENDIF
  260 CONTINUE
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
      GOTO 490
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
      IERROR  = NUMSTK + 1
  490 CALL  DMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_252
      SUBROUTINE DMUMPS_253(COMM_LOAD, ASS_IRECV,
     *    N, INODE, IW, LIW, A, LA, IFLAG,
     *    IERROR, ND, NIRBDU, FILS, FRERE,
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
     *    PERM , MEM_DISTRIB)
      USE DMUMPS_BUFFER
      USE DMUMPS_LOAD
#if defined (PAGING_)
      USE DMUMPS_PAGING
#endif
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N,LIW,LA,NIRBDU,NSTEPS, NBFIN
      INTEGER KEEP(150), ICNTL(40)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     *        LPOOL, LEAF, 
     *        LRLU, IPTRLU,IWPOS, LRLUS,
     *        POSFAC, IWPOSCB, COMP
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER COMM, SLAVEF, MYID,  LBUFR, LBUFR_BYTES
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER IW(LIW), ITLOC(N),
     *        PTRARW(N), PTRAIW(N), ND(KEEP(28)),
     *        FILS(N), FRERE(KEEP(28)),
     *        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28)),
     *        STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), PERM(N)
      INTEGER   CAND(SLAVEF+1, MAX(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER NBPROCFILS(KEEP(28)),
     *        PROCNODE_STEPS(KEEP(28)), BUFR(LBUFR)
      INTEGER      INTARR(MAX(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(MAX(1,KEEP(13)))
      INTEGER LP, HS, HF, HF_OLD,NCBSON, NSLAVES_OLD
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ,LAELL
      INTEGER NCB, IERR_MPI
      INTEGER LAPOS2,J1,J2,J3
      INTEGER NELIM,JJ,JJ1,JJ2,NPIVS,NCOLS,NROWS,
     *        IBROT,IORG
      INTEGER IACHK, LDAFS, LDA_SON
      INTEGER JK,IJROW,NBCOL,ICT13,NUMORG,IOLDPS,J4
      INTEGER APOS, APOS2, AINPUT, POSELT, POSEL1, ICT12
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, ISLAVE, PDEST
      INTEGER PDEST1(1)
      INTEGER NSLAVES_less, ITEMP, NMB_OF_CAND
      LOGICAL FLAG, SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER IZERO
      INTEGER IDUMMY(1)
      PARAMETER( IZERO = 0 )
      INTEGER DMUMPS_275 
      EXTERNAL DMUMPS_275
      INTRINSIC REAL
      DOUBLE PRECISION DATA_SIZE
      DOUBLE PRECISION ZERO
      DATA ZERO /0.0D0/
      INTEGER NELT, LPTRAR
      logical :: force_cand
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      NSTEPS = NSTEPS + 1
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      ENDDO
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .NE. 0) THEN
         DO WHILE (ISON .GT. 0)
           NUMSTK = NUMSTK + 1
           NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1)
           ISON = FRERE(STEP(ISON))
         ENDDO
      ENDIF
      NFRONT = ND(STEP(INODE)) + NASS
      NASS1 = NASS + NUMORG
      NCB   = NFRONT - NASS1
      DATA_SIZE = NASS1 * NCB
      if((KEEP(24).eq.0).or.(KEEP(24).eq.1)) then
         force_cand=.FALSE.
      else
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
      CALL DMUMPS_472( SLAVEF, KEEP,
     *     CAND(1,INIV2),
     *     MEM_DISTRIB(0), NCB, NFRONT, NSLAVES,
     *     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *     TMP_SLAVES_LIST,
     *     SIZE_TMP_SLAVES_LIST )
      HF   = NSLAVES + 6
      LREQ = 2 * NFRONT + HF
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL DMUMPS_94(N, KEEP(28),
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
        CALL  DMUMPS_81(MYID, INODE, N, IOLDPS, HF, NFRONT,
     *        NFRONT_EFF,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     *        INTARR, ITLOC, FILS, FRERE,
     *        SON_LEVEL2, NIV1, NBPROCFILS, KEEP, IFLAG)
      ELSE
        CALL DMUMPS_86( MYID, INODE, N, IOLDPS, HF,
     *        NFRONT, NFRONT_EFF, PERM,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, IW, LIW,
     *        INTARR, ITLOC, FILS, FRERE,
     *        SON_LEVEL2, NIV1, NBPROCFILS, KEEP, IFLAG)
        IF (IFLAG.LT.0) GOTO 250
      ENDIF
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
            NCB    = NFRONT_EFF - NASS1
            NSLAVES_OLD = NSLAVES
            HF_OLD      = HF
            CALL DMUMPS_472( SLAVEF, KEEP,
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
      IF ( KEEP(73) .EQ. 0 ) THEN
#endif
      CALL DMUMPS_461(MYID, SLAVEF, COMM_LOAD,
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
         CALL DMUMPS_94(N, KEEP(28),
     *      IW, LIW, A, LA, NIRBDU,
     *      LRLU, IPTRLU,
     *      IWPOS, IWPOSCB, PTRIST, PTRAST,
     *      STEP, PIMASTER, PAMASTER, ITLOC)
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
      CALL DMUMPS_471(.FALSE.,LA-LRLUS,0,LAELL,KEEP)
      POSEL1 = POSELT - LDAFS
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL - 1
      A(POSELT:LAPOS2) = DBLE(ZERO)
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL - 1
        A(POSELT:LAPOS2) = DBLE(ZERO)
      ELSE
        APOS = POSELT
        DO JJ = 0, LDAFS - 1
          A(APOS:APOS+JJ) = DBLE(ZERO)
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
          OPASSW = OPASSW + DBLE(NELIM*LSTK)
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
            CALL DMUMPS_178(
     *           A( POSELT ), LDAFS, NASS1,
     *           A( IACHK ), LDA_SON,
     *           IW( J1 ), NELIM, NELIM
     *          )
          ENDIF
  210     ISON = FRERE(STEP(ISON))
  220   CONTINUE
      ENDIF
      IBROT = INODE
      DO 260 IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        IBROT = FILS(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        IJROW = INTARR(J1)
        ICT12 = POSELT - LDAFS + IJROW - 1
        DO 240 JJ = J1, J2
          IF (INTARR(JJ).LE.NASS1) THEN
           APOS2 = ICT12 + INTARR(JJ) * LDAFS
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
          ENDIF
          AINPUT = AINPUT + 1
  240   CONTINUE
        IF (J3 .GT. J4) GOTO 260
        ICT13 = POSELT + (IJROW - 1) * LDAFS
        NBCOL = J4 - J3 + 1
        DO JJ = 1, NBCOL
          JJ1 = ICT13 + INTARR(J3 + JJ - 1) - 1
          A(JJ1) = A(JJ1) + DBLARR(AINPUT + JJ - 1)
        ENDDO
  260 CONTINUE
      PTRCOL = IOLDPS + HF + NFRONT 
      PTRROW = IOLDPS + HF + NASS1 
      PDEST  = IOLDPS + 6
      DO ISLAVE = 1, NSLAVES
              CALL DMUMPS_49( 
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
           CALL DMUMPS_68( INODE,
     *      NBPROCFILS(STEP(INODE)),
     *      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     *      IZERO, IDUMMY,
     *      IW(PDEST), NFRONT, COMM, IERR)
         ELSE
           NBCOL = NASS1+SHIFT_INDEX+NBLIG
           CALL DMUMPS_68( INODE,
     *      NBPROCFILS(STEP(INODE)),
     *      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     *      NSLAVES-ISLAVE, IW( PTLUST_S(STEP(INODE))+6+ISLAVE),
     *      IW(PDEST), NFRONT, COMM, IERR)
         ENDIF
         IF (IERR.EQ.-1) THEN
          BLOCKING  = .FALSE.
          SET_IRECV = .TRUE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     *     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *     MPI_ANY_SOURCE, MPI_ANY_TAG,
     *     STATUS, BUFR, LBUFR,
     *     LBUFR_BYTES,
     *     PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *     LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *     PTRIST, PTLUST_S,
     *     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *     IERROR, COMM,
     *     NBPROCFILS,
     *     IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *     root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *     INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *     LPTRAR, NELT, IW, IW,
     *     ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
          IF (MESSAGE_RECEIVED) THEN
           IOLDPS = PTLUST_S(STEP(INODE))
           PTRCOL = IOLDPS + HF + NFRONT
           PTRROW = IOLDPS + HF + NASS1 + SHIFT_INDEX
          ENDIF
         ENDIF
        ENDDO
        IF (IERR .EQ. -2) GOTO 300
        PTRROW = PTRROW + NBLIG
        PDEST  = PDEST + 1
      ENDDO
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
          PDEST1(1)  = DMUMPS_275(STEP(ISON),
     *                 PROCNODE_STEPS, SLAVEF)
          IF (PDEST1(1).EQ.MYID) THEN
            CALL DMUMPS_211( COMM_LOAD, ASS_IRECV, 
     *      BUFR, LBUFR, LBUFR_BYTES,
     *      INODE, ISON, NSLAVES, IW( PTLUST_S(STEP(INODE)) + 6 ),
     *      NFRONT, NASS1, NCBSON, IW( PTRCOL ),
     *      PROCNODE_STEPS,
     *      SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *      LRLUS, N, IW,
     *      LIW, A, LA, NIRBDU,
     *      PTRIST, PTLUST_S, PTRAST, STEP,
     *      PIMASTER, PAMASTER, NSTK_S, COMP,
     *      IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     *      ICNTL, KEEP, root,
     *      OPASSW, OPELIW,
     *      ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *      ND, FRERE, LPTRAR, NELT, IW, IW,
     *
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *      )
           IF ( IFLAG .LT. 0 ) GOTO 500
          ELSE
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON))+HS+NROWS+NPIVS+NELIM
            CALL  DMUMPS_71( INODE, NFRONT,
     *      NASS1, ISON, MYID,
     *      NSLAVES, IW( PTLUST_S(STEP(INODE)) + 6 ),
     *      IW(PTRCOL), NCBSON,
     *      COMM, IERR, PDEST1, NSLSON, SLAVEF, 
     *      KEEP, STEP, N, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *      )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     *        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *        MPI_ANY_SOURCE, MPI_ANY_TAG,
     *        STATUS, BUFR, LBUFR, LBUFR_BYTES,
     *        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *        LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *        PTLUST_S,
     *        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *        IERROR, COMM,
     *        NBPROCFILS,
     *        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *        root,OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ICNTL, KEEP, ND, FRERE, LPTRAR,
     *        NELT, IW, IW, 
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           ENDDO
           IF (IERR .EQ. -2) GOTO 290
          ENDIF
        ELSE
          DO ISLAVE = 0, NSLSON-1
            IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
            CALL DMUMPS_49( 
     &                KEEP, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON,
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX        = PTRCOL + SHIFT_INDEX
              CALL DMUMPS_210( COMM_LOAD, ASS_IRECV, 
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
     *        INTARR, DBLARR, ND, FRERE, LPTRAR, NELT, IW,
     *        IW, 
     *        
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *        )
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
          ENDDO
          IF (PIMASTER(STEP(ISON)).GT.0) THEN
          IERR = -1
          DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  =  PIMASTER(STEP(ISON)) + 6
            CALL  DMUMPS_71( INODE, NFRONT,
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
             CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     *        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *        MPI_ANY_SOURCE, MPI_ANY_TAG,
     *        STATUS, BUFR, LBUFR,
     *        LBUFR_BYTES,
     *        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *        LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *        PTLUST_S,
     *        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *        IERROR, COMM,
     *        NBPROCFILS,
     *        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *        root,OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *        LPTRAR, NELT, IW, IW, 
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
          ENDDO
          IF (IERR .EQ. -2) GOTO 290
          ENDIF
        ENDIF
       ISON = FRERE(STEP(ISON))
      ENDDO
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
  490 CALL DMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_253
      SUBROUTINE DMUMPS_39(N, INODE, IW, LIW, A, LA, 
     *    ISON, NBROWS, NBCOLS, ROWLIST, 
     *    VALSON, PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *    NIRBDU, OPASSW,  OPELIW,
     *    POSFAC, LRLU, IPTRLU, IWPOSCB, LRLUS, 
     *    NBPROCFILS, IPOOL, LPOOL, LEAF, FREEMASTER, PROCNODE_STEPS,
     *    SLAVEF, MYID, KEEP, FILS, ND, COMM_LOAD )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER SLAVEF, KEEP(150), COMM_LOAD
      INTEGER N,LIW,LA,NIRBDU, MYID
      INTEGER INODE,ISON,
     *        LRLU, IPTRLU,LRLUS,
     *        IWPOSCB, POSFAC
      INTEGER NBROWS, NBCOLS
      INTEGER IW(LIW), STEP(N), 
     * PIMASTER(KEEP(28)), PAMASTER(KEEP(28)),
     * NBPROCFILS(KEEP(28)),
     *        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRAST(KEEP(28)),
     *        ROWLIST(NBROWS), FILS(N), ND(KEEP(28))
      DOUBLE PRECISION A(LA), VALSON(NBCOLS,NBROWS)
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER LPOOL, LEAF 
      INTEGER IPOOL( LPOOL ), PROCNODE_STEPS( KEEP(28) )
      LOGICAL FREEMASTER
      INTEGER HF,HS, NSLAVES, NFRONT, NASS1,
     *        IOLDPS, POSELT, POSEL1, ISTCHK,
     *        LSTK, NSLSON,NELIM,NPIVS,NCOLS, J1,J2,J3, JJ,APOS,
     *        JJ1,JJ2, ICT11, JPOS, SIZFI, SIZFR, NCOL, NROW,
     *        NROWS, LDAFS_PERE, IBEG
      LOGICAL SAME_PROC, FREE
      INTRINSIC REAL
      EXTERNAL DMUMPS_275
      INTEGER DMUMPS_275
      IOLDPS = PTLUST_S(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      NFRONT = IW(IOLDPS)
      NASS1  = IABS(IW(IOLDPS + 2))
      NSLAVES= IW(IOLDPS+5)
      IF (KEEP(50).EQ.0) THEN
        LDAFS_PERE = NFRONT
      ELSE
        IF ( NSLAVES .eq. 0 ) THEN
          LDAFS_PERE = NFRONT
        ELSE
          LDAFS_PERE = NASS1
        ENDIF
      ENDIF
      HF      = 6 + NSLAVES
      POSEL1 = POSELT - LDAFS_PERE
      ISTCHK = PIMASTER(STEP(ISON))
      LSTK = IW(ISTCHK)
      NSLSON  = IW(ISTCHK + 5)
      HS      = 6 + NSLSON
      OPASSW = OPASSW + DBLE(NBROWS*NBCOLS)
      NELIM = IW(ISTCHK + 1)
      NPIVS = IW(ISTCHK + 3)
      IF (NPIVS.LT.0) NPIVS = 0
      NCOLS = NPIVS + LSTK
      SAME_PROC = (ISTCHK.LT.IWPOSCB)
      IF (SAME_PROC) THEN
       NROWS = NCOLS
      ELSE
       NROWS = IW(ISTCHK+2)
      ENDIF
      J1 = ISTCHK + NROWS + HS + NPIVS
      IF (KEEP(50).EQ.0) THEN
       DO 170 JJ = 1, NBROWS
        APOS = POSEL1 + ROWLIST(JJ) * LDAFS_PERE
        DO 160 JJ1 = 1, NBCOLS
          JJ2 = APOS + IW(J1 + JJ1 - 1) - 1
          A(JJ2) = A(JJ2) + VALSON(JJ1,JJ) 
  160   CONTINUE
  170  CONTINUE
      ELSE
       DO JJ = 1, NBROWS
        IF (ROWLIST(JJ).LE.NASS1) THEN
         APOS = POSEL1 + ROWLIST(JJ) - 1
         DO JJ1 = 1, NELIM
          JJ2 = APOS + IW(J1 + JJ1 - 1)*LDAFS_PERE
           A(JJ2) = A(JJ2) + VALSON(JJ1,JJ)
         ENDDO
         IBEG = NELIM+1
        ELSE
         IBEG = 1
        ENDIF
        APOS = POSEL1 + ROWLIST(JJ) * LDAFS_PERE
        DO JJ1 = IBEG, NBCOLS
          IF (ROWLIST(JJ).LT.IW(J1 + JJ1 - 1)) EXIT
          JJ2 = APOS + IW(J1 + JJ1 - 1) - 1
          A(JJ2) = A(JJ2) + VALSON(JJ1,JJ)
        ENDDO
       ENDDO
      ENDIF
      IF (NBPROCFILS(STEP(ISON)).EQ.0) THEN
          IF (SAME_PROC) THEN
           IF (KEEP(50).NE.0) THEN
            J2 = J1 +  LSTK - 1
             DO JJ = J1, J2
              IW(JJ) = IW(JJ - NROWS)
             ENDDO
           ELSE
            J2 = J1 + LSTK - 1
            J3 = J1 + NELIM
            DO JJ = J3, J2
             IW(JJ) = IW(JJ - NROWS)
            ENDDO
            IF (NELIM .NE. 0) THEN
              ICT11 = IOLDPS + HF - 1 + NFRONT
              J3 = J3 - 1
              DO 190 JJ = J1, J3
               JPOS = IW(JJ) + ICT11
               IW(JJ) = IW(JPOS)
  190         CONTINUE
            ENDIF
           ENDIF
          ENDIF
          IF (.NOT.FREEMASTER) GOTO 210
          IF (SAME_PROC) THEN
            ISTCHK = PTRIST(STEP(ISON))
            PTRIST(STEP( ISON) ) = -99999999
          ELSE 
            PIMASTER(STEP( ISON )) = -99999999
          END IF
      CALL DMUMPS_152( MYID, N, ISTCHK,
     *     PAMASTER(STEP(ISON)),
     *     IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *     IWPOSCB, LA, KEEP
     *     )
      ENDIF
  210 CONTINUE
      IF ( NBPROCFILS(STEP(INODE)) .EQ. 0 ) THEN
        CALL DMUMPS_507( N, IPOOL, LPOOL, PROCNODE_STEPS,
     *SLAVEF, KEEP(28), KEEP(76), STEP, INODE + N )
        IF (KEEP(47) .GE. 3) THEN
           CALL DMUMPS_500(
     $          IPOOL, LPOOL, 
     *          PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *          MYID, STEP, N, ND, FILS )
        ENDIF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_39
      SUBROUTINE DMUMPS_40(N, INODE, IW, LIW, A, LA, 
     *    NBROWS, NBCOLS, ROWLIST, COLLIST, VALSON, 
     *    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *    FILS, PTRARW, PTRAIW, INTARR, DBLARR, 
     *    ICNTL, KEEP, MYID)
      IMPLICIT NONE
      INTEGER N,LIW,LA
      INTEGER KEEP(150), ICNTL(40)
      INTEGER INODE, MYID
      INTEGER NBROWS, NBCOLS 
      INTEGER ROWLIST(NBROWS), COLLIST(NBCOLS)
      INTEGER IW(LIW), ITLOC(N), STEP(N),
     *        PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28)) , FILS(N), PTRARW(N), PTRAIW(N)
      INTEGER INTARR(MAX(1,KEEP(14)))
      DOUBLE PRECISION A(LA), VALSON(NBCOLS,NBROWS),
     *        DBLARR(MAX(1,KEEP(13)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER IOLDPS, POSELT, NBCOLF, NBROWF, NSLAVES, HF,
     *        K1,K2,K,I,J,POSEL1,APOS,JPOS,NASS,JJ,
     *        IN,AINPUT,JK,J1,J2,IJROW,ICT12, ILOC
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
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
          CALL DMUMPS_ABORT()
       END IF
      NSLAVES = IW(IOLDPS+5)
      HF      = 6 + NSLAVES
      IF ((NASS.GE.0).AND.(NBROWS.LE.0)) GOTO 500
      IF (NASS.LT.0) THEN
          NASS         = -NASS
          IW(IOLDPS+1) = NASS
          A(POSELT:POSELT+NBROWF*NBCOLF-1) = DBLE(ZERO)
          K1 = IOLDPS + HF 
          K2 = K1 + NBROWF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          ENDDO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NASS - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -JPOS
           JPOS     = JPOS + 1
          ENDDO
          IN = INODE
          DO WHILE (IN.GT.0) 
           AINPUT = PTRARW(IN)
           JK     = PTRAIW(IN)
           JJ     = JK + 1
           J1     = JJ + 1
           J2 = J1 + INTARR(JK)
           IJROW = -ITLOC(INTARR(J1))
           ICT12 = POSELT - NBCOLF + IJROW - 1
           DO JJ= J1,J2
            ILOC = ITLOC(INTARR(JJ))
            IF (ILOC.GT.0) THEN
              APOS = ICT12 + ILOC*NBCOLF
              A(APOS) = A(APOS) + DBLARR(AINPUT)
            ENDIF
            AINPUT  = AINPUT + 1
           ENDDO
           IN = FILS(IN)
          ENDDO
          K1 = IOLDPS + HF
          K2 = K1 + NBROWF + NASS - 1
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          ENDDO
      ENDIF
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          ENDDO
          POSEL1 = POSELT - NBCOLF
          IF (KEEP(50).EQ.0) THEN
           DO I=1,NBROWS
            APOS = POSEL1 + ROWLIST(I) * NBCOLF
            DO J=1,NBCOLS
             K = APOS + ITLOC(COLLIST(J)) - 1
             A(K) = A(K) + VALSON(J,I)
            ENDDO
           ENDDO
          ELSE
           DO I=1,NBROWS
            APOS = POSEL1 + ROWLIST(I) * NBCOLF
            DO J=1,NBCOLS
             IF (ITLOC(COLLIST(J)) .EQ. 0) EXIT
             K = APOS + ITLOC(COLLIST(J)) - 1
             A(K) = A(K) + VALSON(J,I)
            ENDDO
           ENDDO
          ENDIF
          OPASSW = OPASSW + NBROWS*NBCOLS
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          ENDDO
      ENDIF
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_40
      SUBROUTINE DMUMPS_178(
     *             AFATH, NFRONT, NASS1,
     *             ACB, NCOLS,
     *             IW, NROWS, NELIM )
      IMPLICIT NONE
      INTEGER NFRONT, NASS1
      DOUBLE PRECISION AFATH( NFRONT * NFRONT )
      INTEGER NCOLS, NROWS, NELIM
      INTEGER IW( NCOLS )
      DOUBLE PRECISION ACB( NCOLS, NROWS )
      INTEGER APOS, POSELT, I, J
      DO I = 1, NELIM
        POSELT = ( IW( I ) - 1 ) * NFRONT
        DO J = 1, I
          APOS = POSELT + IW( J )
          AFATH( APOS ) = AFATH( APOS ) + ACB(J, I)
        END DO
      END DO
      DO I = NELIM + 1, NROWS
        POSELT = IW( I )
        IF (POSELT.LE.NASS1) THEN 
         DO J = 1, NELIM
           APOS = POSELT + ( IW( J ) - 1 ) * NFRONT
           AFATH( APOS ) = AFATH( APOS ) + ACB( J, I )
         END DO
        ELSE
         POSELT = ( IW( I ) - 1 ) * NFRONT
         DO J = 1, NELIM
           APOS = POSELT + IW( J )
           AFATH( APOS ) = AFATH( APOS ) + ACB(J, I)
         END DO
        ENDIF
        POSELT = ( IW( I ) - 1 ) * NFRONT
        DO J = NELIM + 1, I
          APOS = POSELT + IW( J )
          AFATH( APOS ) = AFATH( APOS ) + ACB( J, I )
        END DO
      END DO
      RETURN
      END SUBROUTINE DMUMPS_178
      RECURSIVE SUBROUTINE DMUMPS_264(
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
     *    )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
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
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER COMM, MYID
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTLUST_S(KEEP(28)), ITLOC(N), FILS(N), ND(KEEP(28))
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION  DBLARR(MAX(1,KEEP(13)))
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER INODE, POSITION, NPIV, IERR
      INTEGER LAELL, NCOL, POSBLOCFACTO, NROW
      INTEGER IOLDPS, POSELT, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW
      INTEGER ICT11, LPOS, LPOS1, LPOS2
      INTEGER I, IPOS, KPOS, IPIV, FPERE
      INTEGER LCONT,NELIM,NASS, LDA, NCOL_TO_SEND,
     *        SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON, SHIFT_VAL_SON
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION ONE,ALPHA
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      FPERE    = -1
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LAELL = NPIV * NCOL
      IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          IERROR = LAELL - LRLU
          GOTO 700
        END IF
        CALL DMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *      NIRBDU, LRLU, IPTRLU,
     *      IWPOS, IWPOSCB, PTRIST, PTRAST,
     *      STEP, PIMASTER, PAMASTER, ITLOC )
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
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL DMUMPS_471(.FALSE.,LA-LRLUS,0,LAELL,KEEP)
      IPIV = IWPOS
      IWPOS = IWPOS + NPIV
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 IW( IPIV ), NPIV,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 A(POSBLOCFACTO), NPIV*NCOL, 
     *                 MPI_DOUBLE_PRECISION,
     *                 COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
        MSGSOU = DMUMPS_275( STEP(INODE), PROCNODE_STEPS,
     *           SLAVEF )
        WRITE(*,*) MYID,
     *   ': Internal ERROR 1  in DMUMPS_264',
     *   ' INODE =', INODE,
     *   ' MAITRE_DESC_BANDE not yet received from ', MSGSOU
        CALL DMUMPS_ABORT()
      ENDIF
      DO WHILE ( NBPROCFILS( STEP(INODE)) .NE. 0 ) 
        BLOCKING = .TRUE.
        SET_IRECV = .FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL DMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
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
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END  DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
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
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IOLDPS = PTRIST(STEP(INODE))
      IF (NPIV.EQ.0) GOTO 200
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS )
      NASS1  = IW( IOLDPS + 1 )
      NROW1  = IW( IOLDPS + 2 )
      NPIV1  = IW( IOLDPS + 3 )
      NSLAV1 = IW( IOLDPS + 5 )
      HS     = 6 + NSLAV1
      NCOL1  = LCONT1 + NPIV1
      ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
      DO 100 I = 1, NPIV
        IF (IW(IPIV+I-1).EQ.I) GOTO 100
        ISW = IW(ICT11+I)
        IW(ICT11+I) = IW(ICT11+IW(IPIV+I-1))
        IW(ICT11+IW(IPIV+I-1)) = ISW
        IPOS = POSELT + NPIV1 + I - 1
        KPOS = POSELT + NPIV1 + IW(IPIV+I-1) - 1
        CALL DSWAP(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
 100  CONTINUE
      LPOS2 = POSELT + NPIV1
      CALL DTRSM('L','L','N','N',NPIV, NROW1, ONE, 
     *           A(POSBLOCFACTO), NCOL, A(LPOS2), NCOL1)
      LPOS1 = POSBLOCFACTO+NPIV
      LPOS  = LPOS2 + NPIV
      CALL DGEMM('N','N', NCOL-NPIV,NROW1,NPIV,ALPHA,A(LPOS1),NCOL,
     *           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
 200  CONTINUE
      IW(IOLDPS) = IW(IOLDPS) - NPIV
      IW(IOLDPS + 3) = IW(IOLDPS+3) + NPIV
      IF (LASTBL) IW(IOLDPS+1) = IW(IOLDPS + 3)
      IF ( .not. LASTBL .AND. (IW(IOLDPS+1) .EQ. IW(IOLDPS + 3)) ) THEN
        write(*,*) ' ERROR 1 **** IN BLACFACTO '
        CALL DMUMPS_ABORT()
      ENDIF
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      POSFAC = POSFAC - LAELL
      CALL DMUMPS_471(.FALSE.,LA-LRLUS,0,-LAELL,KEEP)
      IWPOS = IWPOS - NPIV
      IF ( LASTBL  .AND. (KEEP(38).EQ.FPERE) ) THEN
       LCONT  =IW(IOLDPS)
       NROW   = IW(IOLDPS+2)
       NPIV   = IW(IOLDPS+3)
       NASS   = IW(IOLDPS+4)
       NELIM  = NASS-NPIV
       NCOL_TO_SEND =  LCONT-NELIM
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5) 
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NASS
       SHIFT_VAL_SON      = NASS
       LDA                = LCONT + NPIV
       CALL DMUMPS_80(
     *    COMM_LOAD, ASS_IRECV, 
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
        CALL DMUMPS_314( N, INODE,  
     *    PTRIST, PTRAST, PTLUST_S, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP)
       ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_264
      SUBROUTINE DMUMPS_301( id)
      USE DMUMPS_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC), TARGET :: id
      INTEGER MP,LP, MPG
      LOGICAL PROK, PROKG
      INTEGER MTYPE
      LOGICAL LSCAL, ERANAL, GIVSOL
      INTEGER ICNTL10, ICNTL11
      INTEGER I,K,JPERM
      INTEGER LIW,LIWW
      INTEGER LA, LWC, ISRW2, LIW_PASSED, LA_PASSED
      INTEGER DMUMPS_LBUF, DMUMPS_LBUF_INT
      INTEGER NITREF, NOITER, SOLVET, KASE, JOBIREF
      DOUBLE PRECISION RSOL(1), ONE
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0 )
      DOUBLE PRECISION, ALLOCATABLE :: SAVERHS(:), RW1(:),
     *                                 RW2(:), SRW1(:),
     *                                 Y(:), W(:), D(:)
      INTEGER LIWK_SOLVE
      INTEGER, ALLOCATABLE :: IW1(:), IWK_SOLVE(:)
      DOUBLE PRECISION RHSSOL( id%N )
      INTEGER, POINTER :: JOB, N, NZ, MAXIS, MAXS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CNTL
      INTEGER, DIMENSION (:), POINTER :: KEEP,ICNTL,INFO
      INTEGER, DIMENSION (:), POINTER :: IS
      DOUBLE PRECISION, DIMENSION(:),POINTER::ASPK,COLSCA,ROWSCA,S
      DOUBLE PRECISION, DIMENSION(:),POINTER::   RINFOG
      DOUBLE PRECISION ARRET
      DOUBLE PRECISION DUMMY(1)
      INTEGER allocok, PERLU
      LOGICAL I_AM_SLAVE
      JOB  =>id%JOB
      N    =>id%N
      NZ   =>id%NZ
      MAXIS=>id%MAXIS
      MAXS =>id%MAXS
      CNTL =>id%CNTL
      KEEP =>id%KEEP
      IS   =>id%IS
      S    =>id%S
      ICNTL=>id%ICNTL
      INFO =>id%INFO
      ASPK =>id%A
      COLSCA =>id%COLSCA
      ROWSCA =>id%ROWSCA
      RINFOG =>id%RINFOG
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               KEEP(46) .eq. 1 ) )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LP  = ICNTL( 1 )
      PROK  = (MP.GT.0)
      PROKG = (MPG.GT.0 .and. id%MYID.eq.MASTER)
      IF ( PROK  ) WRITE(MP,100)
      IF ( PROKG ) WRITE(MPG,100)
      IF (id%MYID .eq. MASTER) THEN
        MTYPE = ICNTL(  9 )
        IF ( PROKG )  WRITE( MPG, 150 )
     *                ICNTL(9), ICNTL(10), ICNTL(11)
        LSCAL = (((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 6)) .OR. (
     *    KEEP(52) .EQ. -1))
        ERANAL = ((ICNTL(11) .GT. 0) .OR. (ICNTL(10) .GT. 0))
        IF (KEEP(23) .NE. 0) THEN
          IF (MTYPE .NE. 1) THEN
            DO 10 I = 1, N
              RHSSOL(I)=id%RHS(I)
   10       CONTINUE
            DO 20 I = 1, N
              JPERM = id%UNS_PERM(I)
              id%RHS(I) = RHSSOL(JPERM)
   20       CONTINUE
          ENDIF
        ENDIF
        IF (ERANAL) THEN
          ALLOCATE(SAVERHS(N),stat = allocok)
            IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem in solve: error allocating SAVERHS'
              INFO(1) = -13
              INFO(2) = N
              GOTO 444
            END IF
          DO 30 K = 1, N
            SAVERHS(K) = id%RHS(K)
   30     CONTINUE
        ENDIF
        IF (LSCAL) THEN
          IF (MTYPE .EQ. 1) THEN
            DO 40 K = 1, N
              id%RHS(K) = id%RHS(K) * ROWSCA(K)
   40       CONTINUE
          ELSE
            DO 50 K = 1, N
              id%RHS(K) = id%RHS(K) * COLSCA(K)
   50       CONTINUE
          ENDIF
        END IF
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        PERLU      = KEEP(12)
        IF ( KEEP( 50 ) .EQ. 1 ) PERLU = 0
        DMUMPS_LBUF = INT( 2.01D0 * KEEP( 43 ) *KEEP( 35 )  )
        DMUMPS_LBUF = MAX( DMUMPS_LBUF, 50000 )
        DMUMPS_LBUF = DMUMPS_LBUF
     *                 + INT( 2.0D0 * DBLE(MAX(PERLU,0))*
     *                   DBLE(DMUMPS_LBUF)/100D0)
        DMUMPS_LBUF_INT = id%NSLAVES * id%NSLAVES * 4
     *                 * KEEP(34)
        CALL DMUMPS_55( DMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = DMUMPS_LBUF_INT
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     *      ':Error allocating small Send buffer:IERR=',IERR
          END IF
          GOTO 444
        END IF
        CALL DMUMPS_53( DMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = DMUMPS_LBUF
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     *      ':Error allocating Send buffer:IERR=', IERR
          END IF
          GOTO 444
        END IF
      END IF
 444  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF ( INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST(MTYPE,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      IF ( id%MYID .NE. MASTER ) THEN
        IF ( associated(id%RHS) ) THEN
          WRITE(*,*) id%MYID,
     *    ": Warning: RHS associated in DMUMPS_301"
          DEALLOCATE(id%RHS)
        END IF
        ALLOCATE(id%RHS(N),stat=ierr)
        IF ( ierr .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N
          IF (LP > 0)
     *      WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
        END IF
      ELSE IF (.NOT. associated(id%RHS) ) THEN
        WRITE(LP,*) 'Error in solve: RHS not initialized'
        INFO(1) = -3
        INFO(2) = id%JOB
      END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF ( INFO( 1 ) .LT. 0 ) GOTO 90
      CALL MPI_BCAST(id%RHS, N, MPI_DOUBLE_PRECISION,
     *               MASTER,id%COMM,IERR)
      IF ( I_AM_SLAVE ) THEN
        LIW = KEEP(32)
        LA  = KEEP(31)
        LIWK_SOLVE = 5 * KEEP(28) + 2
        ALLOCATE ( IWK_SOLVE( LIWK_SOLVE), stat = allocok )
        IF (allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=LIWK_SOLVE
          GOTO 555
        END IF
        LIWW  = MAXIS - LIW
        IF ( LIWW .LT. 0 ) THEN
          WRITE(LP,*) 'ERROR in DMUMPS_301: LIWW<0'
          CALL DMUMPS_ABORT()
        END IF
        LWC     = MAXS - LA
        ISRW2   = LA + 1
        ALLOCATE(SRW1(KEEP(133)), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=KEEP(133)
          GOTO 555
        END IF
      END IF
 555  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF ( INFO(1) .LT. 0 ) GOTO 90
      IF ( I_AM_SLAVE ) THEN
        LIW_PASSED = MAX( LIW, 1 )
        LA_PASSED  = MAX( LA, 1 )
        CALL DMUMPS_245(id%root, N, S(1), LA_PASSED,
     *    IS(1), LIW_PASSED,
     *    S(ISRW2), LWC, IS( LIW + 1 ),
     *    LIWW, id%RHS,RHSSOL,
     *    id%NA,id%LNA,id%NE_STEPS(1), SRW1, MTYPE,
     *    ICNTL, id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     *    id%PTLUST_S,
     *    IWK_SOLVE, LIWK_SOLVE,
     *    id%PROCNODE_STEPS,
     *    id%NSLAVES, INFO, KEEP,
     *    id%COMM, id%COMM_NODES, id%MYID,
     *    id%MYID_NODES,
     *    id%BUFR,
     *    id%LBUFR, id%LBUFR_BYTES, 
     *
     *    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1))
        DEALLOCATE(SRW1)
      END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1).eq.-2) then
        INFO(1)=-11
        write(*,*) ' WARNING : -11 error code obtained in solve'
      END IF
      IF (INFO(1).eq.-3) then
        INFO(1)=-14
        write(*,*) ' WARNING : -14 error code obtained in solve'
      END IF
      IF (INFO(1).LT.0) GO TO 90
      IF ( .NOT. I_AM_SLAVE ) THEN
        DO 95 I = 1, N
          RHSSOL(I)= DBLE(ZERO)
 95     CONTINUE
      END IF
      CALL MPI_REDUCE(RHSSOL,id%RHS,N,
     *                MPI_DOUBLE_PRECISION,
     *                MPI_SUM,MASTER,id%COMM, IERR)
      IF ( id%MYID.eq.MASTER) THEN
        IF (LSCAL) THEN
          IF (MTYPE .EQ. 1) THEN
              DO K = 1, N
                id%RHS(K) = id%RHS(K) * COLSCA(K)
              END DO
          ELSE
              DO K = 1, N
                id%RHS(K) = id%RHS(K) * ROWSCA(K)
              END DO
          ENDIF
        END IF
      END IF
      IF ( id%MYID.eq.MASTER) THEN
        IF ( (KEEP(55).eq.0) .AND. KEEP(54).eq.0 .AND. 
     &      .NOT.ASSOCIATED(id%A) ) THEN
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        ELSE
          ICNTL10 = ICNTL(10)
          ICNTL11 = ICNTL(11)
        ENDIF
      END IF
      CALL MPI_BCAST(ICNTL10,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL11,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      IF ( id%MYID.ne.MASTER ) THEN
        ERANAL = ((ICNTL11 .GT. 0) .OR. (ICNTL10 .GT. 0))
      END IF
      IF (ERANAL) THEN
        IF ((ICNTL10 .GT. 0) .AND. (ICNTL11 .GT. 0)) THEN
          IF (id%MYID .EQ. MASTER) THEN
            GIVSOL = .FALSE.
            IF (MP .GT. 0) WRITE( MP, 170 )
            ALLOCATE(RW1(N),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
              GOTO 776
            endif
            ALLOCATE(RW2(N),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
              GOTO 776
            endif
          END IF
          IF ( KEEP(54) .ne. 0 ) THEN
            ALLOCATE( SRW1( N ), stat =allocok )
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
            endif
          END IF
 776      CONTINUE
          CALL DMUMPS_276( ICNTL, INFO,
     *                  id%COMM,id%MYID)
          IF ( INFO(1) .LT. 0 ) GOTO 90
          IF ( KEEP(54) .eq. 0 ) THEN
            IF (id%MYID .EQ. MASTER) THEN
              IF (KEEP(55).EQ.0) THEN
                CALL DMUMPS_278( ICNTL(9), N, NZ, id%A(1),
     *             id%IRN(1), id%JCN(1),
     *             id%RHS, SAVERHS, RW1, RW2, KEEP )
              ELSE
                CALL DMUMPS_121( ICNTL(9), N, 
     *          id%NELT, id%ELTPTR, 
     *          id%LELTVAR, id%ELTVAR,
     *          id%NA_ELT, id%A_ELT,
     *          id%RHS, SAVERHS, RW1, RW2, KEEP )
              ENDIF
            END IF
          ELSE
            CALL MPI_BCAST( id%RHS, id%N,
     *              MPI_DOUBLE_PRECISION, MASTER,
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        id%RHS, SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, RW2,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
              RW2 = SAVERHS - RW2
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, RW1,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            DEALLOCATE( SRW1 )
          END IF
          IF ( id%MYID .EQ. MASTER )  THEN
            CALL DMUMPS_205(ICNTL(9),INFO(1),N,NZ,
     *        id%RHS, SAVERHS,RW1,RW2,GIVSOL,
     *        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MP,ICNTL,KEEP)
            DEALLOCATE(RW1)
            DEALLOCATE(RW2)
          END IF
        END IF
      IF ( PROK  .AND. ICNTL10 .GT. 0 ) WRITE( MP, 270 )
      IF ( PROKG .AND. ICNTL10 .GT. 0 ) WRITE( MPG, 270 )
      ALLOCATE(Y(N), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        GOTO 777
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
        ALLOCATE( IW1( 2 * N ),stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=2 * N
          GOTO 777
        ENDIF
        ALLOCATE( D(N),stat =allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N
          GOTO 777
        ENDIF
        ALLOCATE( W(3*N), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N
          GOTO 777
        ENDIF
        NITREF = ICNTL10
        JOBIREF= ICNTL11
        IF ( PROKG .AND. ICNTL10 .GT. 0 )
     *    WRITE( MPG, 240) 'MAXIMUM NUMBER OF STEPS =', NITREF
        DO I = 1, N
          D( I ) = ONE
        END DO
      END IF
      ALLOCATE(SRW1(N),stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        GOTO 777
      ENDIF
      KASE = 0
 777  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF ( INFO(1) .LT. 0 ) GOTO 90
  22    CONTINUE
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 0 ) THEN
              IF (KEEP(55).NE.0) THEN 
               CALL DMUMPS_119(MTYPE, N, 
     *           id%NELT, id%ELTPTR(1), 
     *           id%LELTVAR, id%ELTVAR(1),
     *           id%NA_ELT, id%A_ELT,
     *           W(N+1), KEEP )
              ELSE
               IF ( MTYPE .eq. 1 ) THEN
                 CALL DMUMPS_207
     *       ( id%A(1), NZ, N, id%IRN(1), id%JCN(1), W(N+1), KEEP)
               ELSE
                 CALL DMUMPS_207
     *       ( id%A(1), NZ, N, id%JCN(1), id%IRN(1), W(N+1), KEEP)
               END IF
              ENDIF
            ENDIF
          END IF
        ELSE
          IF ( KASE .eq. 0 ) THEN
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              IF ( MTYPE .eq. 1 ) THEN
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP )
              ELSE
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%JCN_loc, id%IRN_loc,
     *          SRW1, id%KEEP )
              END IF
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, W( N + 1 ),
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
          END IF
        END IF
        IF ( id%MYID .eq. MASTER ) THEN
            ARRET = CNTL(2)
            IF (ARRET .LT. 0.0D0) THEN
              ARRET = SQRT(EPSILON(0.0D0))
            END IF
            CALL DMUMPS_206(NZ,N,SAVERHS,id%RHS,
     *      Y, D, W,
     *      IW1, KASE,RINFOG(7),
     *      RINFOG(9), JOBIREF, RINFOG(10), NITREF, NOITER, MP,
     *      KEEP, ARRET )
        END IF
        IF ( KEEP(54) .ne. 0 ) THEN
          CALL MPI_BCAST( KASE, 1, MPI_INTEGER, MASTER,
     *    id%COMM, IERR )
        END IF
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 14 ) THEN
              IF (KEEP(55).NE.0) THEN
               CALL DMUMPS_122( MTYPE, N, 
     *            id%NELT, id%ELTPTR, id%LELTVAR,
     *            id%ELTVAR, id%NA_ELT, id%A_ELT,
     *            SAVERHS, id%RHS,
     *            Y, W, KEEP(50))
              ELSE
                 IF ( MTYPE .eq. 1 ) THEN
                   CALL DMUMPS_208
     *    (id%A(1), NZ, N, id%IRN(1), id%JCN(1), SAVERHS, id%RHS,
     *              Y, W, KEEP)
                 ELSE
                   CALL DMUMPS_208
     *    (id%A(1), NZ, N, id%JCN(1), id%IRN(1), SAVERHS, id%RHS,
     *              Y, W, KEEP)
                 END IF
              ENDIF
              GOTO 22
            END IF
          END IF
        ELSE
          IF ( KASE.eq.14 ) THEN
            CALL MPI_BCAST( id%RHS, id%N,
     *              MPI_DOUBLE_PRECISION, MASTER,
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        id%RHS, SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, Y,
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
              Y = SAVERHS - Y
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and. id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_193( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        id%RHS, SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, W,
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            GOTO 22
          END IF
        END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF ( KASE .GT. 0 ) THEN
          IF ( MTYPE .EQ. 1 ) THEN
            SOLVET = KASE - 1
          ELSE
            SOLVET = KASE
          END IF
          IF ( LSCAL ) THEN
            IF ( SOLVET .EQ. 1 ) THEN
              DO K = 1, N
                Y( K ) = Y( K ) * ROWSCA( K )
              END DO
            ELSE
              DO K = 1, N
                Y( K ) = Y( K ) * COLSCA( K )
              END DO
            END IF
          END IF
        END IF
      END IF
      CALL MPI_BCAST( KASE , 1, MPI_INTEGER, MASTER,
     *                id%COMM, IERR)
      CALL MPI_BCAST( SOLVET, 1, MPI_INTEGER, MASTER,
     *                id%COMM, IERR)
      IF ( KASE .GT. 0 ) THEN
        CALL MPI_BCAST( Y, N, MPI_DOUBLE_PRECISION, MASTER,
     *                id%COMM, IERR )
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED = MAX( LIW, 1 )
          LA_PASSED = MAX( LA, 1 )
          CALL DMUMPS_245( id%root, N,
     *    S(1), LA_PASSED, id%IS( 1 ),
     *    LIW_PASSED, S( ISRW2 ), LWC, id%IS( LIW + 1 ),
     *    LIWW,
     *    Y, RHSSOL, id%NA, id%LNA, id%NE_STEPS,
     *    SRW1, SOLVET, ICNTL,
     *    id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     *    id%PTLUST_S,
     *    IWK_SOLVE, LIWK_SOLVE,
     *    id%PROCNODE_STEPS, id%NSLAVES, INFO, KEEP,
     *    id%COMM,
     *    id%COMM_NODES,
     *    id%MYID, id%MYID_NODES,
     *    id%BUFR, id%LBUFR, id%LBUFR_BYTES , 
     *
     *    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1))
        END IF
        CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
        IF (INFO(1).eq.-2) INFO(1)=-12
        IF (INFO(1).eq.-3) INFO(1)=-15
        IF (INFO(1).LT.0) GO TO 90
        IF ( .NOT. I_AM_SLAVE ) THEN
          DO I = 1, N
            RHSSOL(I)= DBLE(ZERO)
            Y     (I)= DBLE(ZERO)
          END DO
        END IF
        CALL MPI_REDUCE(RHSSOL,Y,N,
     *                  MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,MASTER,id%COMM, IERR)
        IF ( id%MYID.eq.MASTER) THEN
          IF (LSCAL) THEN
            IF (SOLVET .EQ. 1) THEN
               DO K = 1, N
                 Y(K) = Y(K) * COLSCA(K)
               END DO
            ELSE
               DO K = 1, N
                 Y(K) = Y(K) * ROWSCA(K)
               END DO
            ENDIF
          END IF
        END IF
        GO TO 22
      ELSEIF ( KASE .LT. 0 ) THEN
           INFO( 1 ) = INFO( 1 ) + 8
      END IF
      IF ( id%MYID .eq. MASTER ) DEALLOCATE(W)
      IF ( id%MYID .eq. MASTER ) DEALLOCATE(IW1)
      IF ( PROKG .AND. NITREF .GT. 0 .AND.
     *id%MYID .EQ. MASTER ) THEN
        WRITE( MPG, 81 ) 
        WRITE( MPG, 141 ) 'NUMBER OF STEPS OF ITERATIVE REFINEMENTS  
     *=', NOITER
      ENDIF
      IF ( id%MYID .EQ. MASTER .AND. NITREF .GT. 0 ) THEN
        id%INFOG(15) = NOITER
      END IF
      IF ( PROK .AND. NITREF .GT.0 ) WRITE( MP, 131 ) 
      IF (ICNTL11 .GT. 0) THEN
        IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
              CALL DMUMPS_278( MTYPE, N, NZ, id%A(1),
     *          id%IRN(1), id%JCN(1),
     *          id%RHS, SAVERHS, Y, D, KEEP )
            ELSE
              CALL DMUMPS_121( MTYPE, N, 
     *          id%NELT, id%ELTPTR, 
     *          id%LELTVAR, id%ELTVAR,
     *          id%NA_ELT, id%A_ELT,
     *          id%RHS, SAVERHS, Y, D, KEEP )
            ENDIF
          END IF
        ELSE
            CALL MPI_BCAST( id%RHS, id%N,
     *              MPI_DOUBLE_PRECISION, MASTER, 
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        id%RHS, SRW1, KEEP(50), MTYPE ) 
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, D,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
              D = SAVERHS - D
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, Y,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
        END IF
        IF (id%MYID .EQ. MASTER) THEN
         IF ((MPG .GT. 0) .AND. (NITREF .GT. 0)) WRITE( MPG, 65 ) 
         IF ((MPG .GT. 0) .AND. (NITREF .LE. 0)) WRITE( MPG, 170 ) 
         GIVSOL = .FALSE.
         CALL DMUMPS_205(MTYPE,INFO(1),N,NZ,id%RHS,
     *        SAVERHS,Y,D,GIVSOL,
     *        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MPG,ICNTL, KEEP)
         IF ( MPG .GT. 0 ) THEN
          WRITE( MPG, 115 )
     *'RINFOG(7):COMPONENTWISE SCALED RESIDUAL(W1)=', RINFOG(7)
          WRITE( MPG, 115 )
     *'------(8):---------------------------- (W2)=', RINFOG(8)
          WRITE( MPG, 115 )
     *'------(9):Upper bound ERROR ...............=', RINFOG(9)
          WRITE( MPG, 115 )
     *'-----(10):CONDITION NUMBER (1) ............=', RINFOG(10)
          WRITE( MPG, 115 )
     *'-----(11):CONDITION NUMBER (2) ............=', RINFOG(11)
         END IF
        END IF ! MASTER
      END IF ! ICNTL11 > 0
      IF (id%myid == MASTER) DEALLOCATE(D)
      DEALLOCATE(Y)
      DEALLOCATE(SRW1)
      END IF
      IF ( id%MYID .EQ. MASTER ) THEN
        IF ((ICNTL(9) .EQ. 1) .AND. (KEEP(23) .NE. 0)) THEN
          ALLOCATE( RW1( N ),stat =allocok )
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=N
            WRITE(*,*) 'could not allocate ', N, 'integers.'
            CALL DMUMPS_ABORT()
          END IF
          DO 70 I = 1, N
            RW1(I) = id%RHS(I)
 70       CONTINUE
          DO 80 I = 1, N
            JPERM = id%UNS_PERM(I)
            IF ( JPERM .eq. I ) GO TO 80
            id%RHS( JPERM ) = RW1( I )
 80       CONTINUE
          DEALLOCATE( RW1 )
        END IF
        IF (PROK) WRITE( MP, 120 ) INFO(1), INFO(2)
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        IF ( INFO(1) .GE. 0 .AND. ICNTL(4).GE.3 .AND. ICNTL(3).GT.0)
     *    THEN
          K = min0(10, N)
          IF (ICNTL(4) .eq. 4 ) K = N
          WRITE(ICNTL(3),110)
          IF ( N .GT. 0 ) WRITE(ICNTL(3),160)
     *    (id%RHS(I),I=1,K)
        END IF
      END IF
 90   CONTINUE
      IF ( I_AM_SLAVE ) THEN
        IF (allocated(IWK_SOLVE)) DEALLOCATE( IWK_SOLVE )
        CALL DMUMPS_57( IERR )
        CALL DMUMPS_59( IERR )
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF (allocated(SAVERHS)) DEALLOCATE( SAVERHS)
      ELSE
        IF (ASSOCIATED(id%RHS)) THEN
          DEALLOCATE(id%RHS)
          NULLIFY(id%RHS)
        END IF
      END IF
      RETURN
 65   FORMAT (//' ERROR ANALYSIS AFTER ITERATIVE REFINEMENT')
 100  FORMAT(//' ****** SOLVE & CHECK STEP ********'/)
 110  FORMAT (//' LEAVING SOLVER WITH SOLUTION')
 115  FORMAT(1X, A44,1P,D9.2)
 120  FORMAT(//' LEAVING SOLVER WITH:  INFO(1) ............ =',I12/
     *         '                       INFO(2) ............ =',I12)
 150  FORMAT (/' STATISTICS PRIOR SOLVE PHASE     ...........'/
     *        ' ICNTL (9)                                     =',I12/
     *        '  --- (10)                                     =',I12/
     *        '  --- (11)                                     =',I12)
 160  FORMAT (' RHS'/(1X,1P,5D14.6))
 170  FORMAT (//' ERROR ANALYSIS' )
 240  FORMAT (1X, A42,I4)
 270  FORMAT (//' BEGIN ITERATIVE REFINEMENT' )
  81  FORMAT (/' STATISTICS AFTER ITERATIVE REFINEMENT ')
 131  FORMAT (/' END   ITERATIVE REFINEMENT ')
 141  FORMAT(1X, A42,I4)
      END SUBROUTINE DMUMPS_301
      SUBROUTINE DMUMPS_265( COMM_LOAD, ASS_IRECV, 
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, POSFAC,
     *   N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S, PTRAST,
     *   STEP, PIMASTER, PAMASTER, NBPROCFILS,
     *   COMP, root, OPASSW, OPELIW, ITLOC, NSTK_S,
     *   FILS, PTRARW, PTRAIW, INTARR, DBLARR, NBFIN,
     *   MYID, COMM, ICNTL, KEEP, IFLAG, IERROR,
     *   IPOOL, LPOOL, LEAF, ND, FRERE_STEPS, LPTRAR, NELT,
     *   FRTPTR, FRTELT, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 150 )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, LIW, LA
      INTEGER POSFAC, NBFIN
      INTEGER COMP
      INTEGER NIRBDU
      INTEGER NELT, LPTRAR
      INTEGER PROCNODE_STEPS( KEEP(28) ), PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER PTLUST_S( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER ITLOC( N ), NSTK_S( KEEP(28) ), FILS( N )
      INTEGER ND(KEEP(28)), FRERE_STEPS( KEEP(28) )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( MAX( 1,KEEP(13)) )
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER COMM, MYID, IFLAG, IERROR
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER FRTPTR(N+1), FRTELT( NELT )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      INTEGER IERR
      INTEGER I, INODE, ISON, POSITION, NBROW, LROW, IROW, INDCOL
      INTEGER LREQI, LREQA
      INTEGER POSCONTRIB
      INTEGER MASTER
      LOGICAL SLAVE_NODE, FREEMASTER
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, ISON, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NBROW, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, LROW, 1,
     *                 MPI_INTEGER, COMM, IERR )
      SLAVE_NODE = DMUMPS_275( STEP(INODE),
     *             PROCNODE_STEPS, SLAVEF ) .NE. MYID
      IF ( SLAVE_NODE ) THEN
        LREQI = LROW + NBROW
      ELSE
        LREQI = NBROW
      END IF
      LREQA = LROW * NBROW
      IF ( LRLU .LT. LREQA .OR. IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LREQA ) THEN
          IFLAG = -9
          IERROR = LREQA - LRLUS
          CALL DMUMPS_44( MYID, SLAVEF, COMM )
          RETURN
        END IF
        CALL DMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *       NIRBDU, LRLU, IPTRLU,
     *       IWPOS, IWPOSCB, PTRIST, PTRAST,
     *       STEP, PIMASTER, PAMASTER, ITLOC)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..process_contrib'
             WRITE(*,*) 'LRLU,LRLUS=',LRLU,LRLUS
             IFLAG = -9
             IERROR = LREQA - LRLUS
             CALL DMUMPS_44( MYID, SLAVEF, COMM )
             RETURN
        END IF
        IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
          IFLAG  = -8
          IERROR = IWPOS + LREQI - 1 - IWPOSCB
          CALL DMUMPS_44( MYID, SLAVEF, COMM )
          RETURN
        END IF
      END IF
      LRLU  = LRLU - LREQA
      LRLUS = LRLUS - LREQA
      POSCONTRIB = POSFAC
      POSFAC = POSFAC + LREQA
      CALL DMUMPS_471(.FALSE.,LA-LRLUS,0,LREQA,KEEP)
      IF  ( SLAVE_NODE ) THEN
         IROW   = IWPOS
         INDCOL = IWPOS + NBROW
      ELSE
         IROW   = IWPOS
         INDCOL = -1
      END IF
      IWPOS = IWPOS + LREQI
      IF ( SLAVE_NODE ) THEN
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 IW( INDCOL ), LROW, MPI_INTEGER,
     *                 COMM, IERR )
      END IF
      DO I = 1, NBROW
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   IW( IROW + I - 1 ), 1, MPI_INTEGER,
     *                   COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   A( POSCONTRIB + ( I - 1 ) * LROW ), LROW,
     *                   MPI_DOUBLE_PRECISION,
     *                   COMM, IERR )
      END DO
      IF ( SLAVE_NODE ) THEN
        DO WHILE ( PTRIST( STEP(INODE) ) .EQ. 0 )
          MASTER = DMUMPS_275(STEP(INODE),PROCNODE_STEPS,SLAVEF)
          BLOCKING = .TRUE.
          SET_IRECV = .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_329( COMM_LOAD, ASS_IRECV,
     *     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *     MASTER, MAITRE_DESC_BANDE,
     *     STATUS, 
     *     BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *     IWPOS, IWPOSCB, IPTRLU,
     *     LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, 
     *     PTRIST, PTLUST_S,
     *     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *     IFLAG, IERROR, COMM,
     *     NBPROCFILS, IPOOL, LPOOL, LEAF,
     *     NBFIN, MYID, SLAVEF,
     *
     *     root, OPASSW, OPELIW, ITLOC, FILS, 
     *     PTRARW, PTRAIW,
     *     INTARR, DBLARR, ICNTL, KEEP, ND, FRERE_STEPS,
     *     LPTRAR, NELT, FRTPTR, FRTELT, 
     *     ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF (IFLAG.LT.0) RETURN
        END DO
        NBPROCFILS(STEP(INODE))=NBPROCFILS(STEP(INODE))-1
        IF ( KEEP(55) .eq. 0 ) THEN
        CALL DMUMPS_40(N, INODE, IW, LIW, A, LA,
     *  NBROW, LROW, IW( IROW ), IW( INDCOL ), A( POSCONTRIB ),
     *  OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *  FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL, KEEP, MYID )
        ELSE
        CALL DMUMPS_123(
     *  NELT, FRTPTR, FRTELT,
     *  N, INODE, IW, LIW, A, LA,
     *  NBROW, LROW, IW( IROW ), IW( INDCOL ), A( POSCONTRIB ),
     *  OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *  FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL, KEEP, MYID )
        END IF
      ELSE
        NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) - 1
        NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
        FREEMASTER = .TRUE.
        CALL DMUMPS_39(N, INODE, IW, LIW, A, LA,
     *  ISON, NBROW, LROW, IW( IROW ),
     *  A(POSCONTRIB), PTRIST, PTLUST_S, PTRAST,
     *  STEP, PIMASTER, PAMASTER, NIRBDU, OPASSW, OPELIW,
     *  POSFAC, LRLU, IPTRLU, IWPOSCB, LRLUS, NBPROCFILS,
     *  IPOOL, LPOOL, LEAF, FREEMASTER, 
     *  PROCNODE_STEPS, SLAVEF, MYID, KEEP, FILS, ND,
     *  COMM_LOAD)
      END IF
      LRLU = LRLU + LREQA
      LRLUS = LRLUS + LREQA
      POSFAC = POSFAC - LREQA
      CALL DMUMPS_471(.FALSE.,LA-LRLUS,0,-LREQA,KEEP)
      IWPOS = IWPOS - LREQI
      RETURN
      END SUBROUTINE DMUMPS_265
      SUBROUTINE DMUMPS_143( N, INODE, IW, LIW, A, LA,
     *                           IOLDPS, POSELT, IFLAG, UU, NOFFW,
     *                           NPVW,
     *                           KEEP, STEP,
     *                           PROCNODE_STEPS, MYID, SLAVEF, SEUIL )
      IMPLICIT NONE
      INTEGER N, INODE, LIW, LA, IFLAG, NOFFW, NPVW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER MYID, SLAVEF, IOLDPS, POSELT
      INTEGER KEEP( 150 )
      INTEGER PROCNODE_STEPS( KEEP(28) ), STEP(N)
      DOUBLE PRECISION UU, SEUIL
      INTEGER INOPV, IFINB, NFRONT, NPIV
      INTEGER NASS, NEL1, NPIVB, NPIVE, NBTLKJ
      EXTERNAL DMUMPS_330, DMUMPS_221, DMUMPS_233, 
     *         DMUMPS_229,
     *         DMUMPS_225, DMUMPS_232, DMUMPS_231,
     *         DMUMPS_220,
     *         DMUMPS_228, DMUMPS_236
      INTEGER  DMUMPS_330
      NFRONT = IW(IOLDPS)
      NASS   = IABS(IW(IOLDPS+2))
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5), NASS )
      ENDIF
 50   CONTINUE
      CALL DMUMPS_221(NFRONT,NASS,N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                IFLAG,IOLDPS,POSELT,UU,SEUIL)
      IF (IFLAG.LT.0) GOTO 500
      IF (INOPV.EQ.1) GO TO 80
      IF (INOPV.EQ.2) THEN
         CALL DMUMPS_233(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4))
         GOTO 50
      ENDIF
      NPVW = NPVW + 1
      IF (NASS.LE.1) THEN
       CALL DMUMPS_229(NFRONT,N,INODE,IW,LIW,A,LA,
     *                 IOLDPS,POSELT)
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       GO TO 500
      ENDIF
       CALL DMUMPS_225(NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     *             IOLDPS,POSELT,IFINB,
     *             NBTLKJ,KEEP(4))
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       IF (IFINB.EQ.0) GOTO 50
       IF (IFINB.EQ.(-1)) GOTO 80
       NPIV   = IW(IOLDPS+1)
       NEL1   = NASS - NPIV
      CALL DMUMPS_232(A,LA,
     *           NFRONT,NPIV,NASS,POSELT,NBTLKJ)
      GO TO 50
 80   CONTINUE
      NPIV   = IW(IOLDPS+1)
      IF (NPIV.LE.0) GO TO 110
      NEL1   = NFRONT - NASS
      IF (NEL1.LE.0) GO TO 110
        CALL DMUMPS_231(A,LA,NFRONT, NPIV,NASS,POSELT)
 110  CONTINUE
      IF (DMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
     *                   .EQ.1) THEN
        NPIV   = IW(IOLDPS+1)
        IW(IOLDPS+4) = NPIV
        IF (NASS.EQ.NPIV) GOTO 500
 120    CALL DMUMPS_220(NFRONT,NASS,N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                IOLDPS,POSELT,UU)
        IF (INOPV.NE.1) THEN
         NPVW = NPVW + 1
         CALL DMUMPS_228(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *                 IOLDPS,POSELT,IFINB)
         IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       IF (IFINB.EQ.0) GOTO 120
        ENDIF
        NPIV   = IW(IOLDPS+1)
        NPIVB  = IW(IOLDPS+4)
        NPIVE  = NPIV - NPIVB
        NEL1   = NFRONT - NASS
        IF ((NPIVE.LE.0).OR.(NEL1.EQ.0)) GO TO 500
        CALL DMUMPS_236(A,LA,NPIVB,
     *                NFRONT,NPIV,NASS,POSELT)
      ENDIF
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_143
