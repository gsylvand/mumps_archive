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
        MODULE DMUMPS_BUFFER
        PRIVATE
        PUBLIC :: DMUMPS_61,
     *   DMUMPS_53 ,       DMUMPS_57 ,
     *   DMUMPS_55, DMUMPS_59,
     *   DMUMPS_54,DMUMPS_58,
     *   DMUMPS_66,     DMUMPS_78,
     *   DMUMPS_62,       DMUMPS_68,
     *   DMUMPS_71, DMUMPS_70,
     *   DMUMPS_67,
     *   DMUMPS_65,
     *   DMUMPS_64,
     *   DMUMPS_72,
     *   DMUMPS_75,
     *   DMUMPS_76,
     *   DMUMPS_73, DMUMPS_74,
     *   DMUMPS_63,DMUMPS_77, 
     *   DMUMPS_60,
     *   DMUMPS_460, DMUMPS_469,
     *   DMUMPS_498, DMUMPS_502
        INTEGER NEXT, REQ, CONTENT, OVHSIZE
        PARAMETER( NEXT = 0, REQ = 1, CONTENT = 2, OVHSIZE = 2 )
        INTEGER, SAVE :: SIZEofINT
        TYPE DMUMPS_BUFFER_TYPE
          INTEGER LBUF, HEAD, TAIL,LBUF_INT, ILASTMSG
          INTEGER, DIMENSION(:),POINTER :: CONTENT
        END TYPE DMUMPS_BUFFER_TYPE
        TYPE ( DMUMPS_BUFFER_TYPE ), SAVE :: BUF_CB
        TYPE ( DMUMPS_BUFFER_TYPE ), SAVE :: BUF_SMALL
        TYPE ( DMUMPS_BUFFER_TYPE ), SAVE :: BUF_LOAD
        INTEGER, SAVE :: SIZE_RBUF_BYTES
      CONTAINS
        SUBROUTINE DMUMPS_61( IntSize )
        IMPLICIT NONE
        INTEGER IntSize
        SIZEofINT = IntSize
        NULLIFY(BUF_CB  %CONTENT)
        NULLIFY(BUF_SMALL%CONTENT)
        NULLIFY(BUF_LOAD%CONTENT)
        BUF_CB%LBUF     = 0
        BUF_CB%LBUF_INT = 0
        BUF_CB%HEAD     = 1
        BUF_CB%TAIL     = 1
        BUF_CB%ILASTMSG = 1
        BUF_SMALL%LBUF     = 0
        BUF_SMALL%LBUF_INT = 0
        BUF_SMALL%HEAD     = 1
        BUF_SMALL%TAIL     = 1
        BUF_SMALL%ILASTMSG = 1
        BUF_LOAD%LBUF     = 0
        BUF_LOAD%LBUF_INT = 0
        BUF_LOAD%HEAD     = 1
        BUF_LOAD%TAIL     = 1
        BUF_LOAD%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_61
        SUBROUTINE DMUMPS_53( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL DMUMPS_2( BUF_CB, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_53
        SUBROUTINE DMUMPS_55( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL DMUMPS_2( BUF_SMALL, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_55
        SUBROUTINE DMUMPS_54( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        INTEGER stat
        CALL DMUMPS_2( BUF_LOAD, SIZE, IERR )        
        RETURN
        END SUBROUTINE DMUMPS_54
        SUBROUTINE DMUMPS_58( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL DMUMPS_3( BUF_LOAD, IERR )
        END SUBROUTINE DMUMPS_58
        SUBROUTINE DMUMPS_57( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL DMUMPS_3( BUF_CB, IERR )
        RETURN
        END SUBROUTINE DMUMPS_57
        SUBROUTINE DMUMPS_59( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL DMUMPS_3( BUF_SMALL, IERR )
        RETURN
        END SUBROUTINE DMUMPS_59
        SUBROUTINE DMUMPS_2( BUF, SIZE, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: BUF
        INTEGER SIZE, IERR
        IERR         = 0
        BUF%LBUF     = SIZE
        BUF%LBUF_INT = ( SIZE + SIZEOFINT - 1 ) / SIZEOFINT
        IF ( ASSOCIATED ( BUF%CONTENT ) ) DEALLOCATE( BUF%CONTENT )
        ALLOCATE( BUF%CONTENT( BUF%LBUF_INT ), stat = ierr )
        IF (IERR .NE. 0) THEN
          NULLIFY( BUF%CONTENT )
          IERR         = -1
          BUF%LBUF     =  0
          BUF%LBUF_INT =  0
        END IF
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_2
        SUBROUTINE DMUMPS_3( BUF, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: BUF
        INCLUDE 'mpif.h'
        INTEGER IERR
        INTEGER STATUS( MPI_STATUS_SIZE )
        LOGICAL FLAG
        IF ( .NOT. ASSOCIATED ( BUF%CONTENT ) ) THEN
          BUF%HEAD     = 1
          BUF%LBUF     = 0
          BUF%LBUF_INT = 0
          BUF%TAIL     = 1
          BUF%ILASTMSG = 1
          RETURN
        END IF
        DO WHILE ( BUF%HEAD.NE.0 .AND. BUF%HEAD .NE. BUF%TAIL )
          CALL MPI_TEST(BUF%CONTENT( BUF%HEAD + REQ ), FLAG,
     *                  STATUS, IERR)
          IF ( .not. FLAG ) THEN
            WRITE(*,*) '** Warning: trying to cancel a request.'
            WRITE(*,*) '** This might be problematic on SGI'
            CALL MPI_CANCEL( BUF%CONTENT( BUF%HEAD + REQ ), IERR )
            CALL MPI_REQUEST_FREE( BUF%CONTENT( BUF%HEAD + REQ ), IERR )
          END IF
          BUF%HEAD = BUF%CONTENT( BUF%HEAD + NEXT )
        END DO
        DEALLOCATE( BUF%CONTENT )
        NULLIFY( BUF%CONTENT )
        BUF%LBUF     = 0
        BUF%LBUF_INT = 0
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_3
        SUBROUTINE DMUMPS_66( INODE, FPERE, NFRONT, LCONT,
     *                               NASS, NPIV,
     *                               IWROW, IWCOL, A,
     *       DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER DEST, TAG, COMM, IERR
        INTEGER INODE, FPERE, NFRONT, LCONT, NASS, NPIV
        INTEGER IWROW( LCONT ), IWCOL( LCONT )
        DOUBLE PRECISION A( * )
        INCLUDE 'mpif.h'
        INTEGER POSITION, IREQ, IPOS, I, J1
        INTEGER SIZE1, SIZE2, SIZE
        INTEGER IZERO, IONE
        INTEGER DEST2(1)
        PARAMETER( IZERO = 0, IONE = 1 )
        DEST2(1) = DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 9 + LCONT + LCONT, MPI_INTEGER,
     *                      COMM, SIZE1,  IERR)
        CALL MPI_PACK_SIZE( LCONT * LCONT, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2,  IERR)
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -2
             RETURN
          ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LCONT, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LCONT, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( NASS-NPIV, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LCONT , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IONE,  1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IWROW, LCONT, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IWCOL, LCONT, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( LCONT .NE. 0 ) THEN
          J1 = 1
          DO I = 1, LCONT
            CALL MPI_PACK( A( J1 ), LCONT, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
            J1 = J1 + NFRONT
          END DO
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_send_cb: SIZE, POSITION = ', SIZE, POSITION
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_66
        SUBROUTINE DMUMPS_72( INODE, IFATH,
     *             LCB, LSOL, CB, SOL,
     *             DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER INODE, IFATH, LCB, LSOL 
        INTEGER DEST, COMM, IERR
        DOUBLE PRECISION CB( LCB ), SOL( MAX( 1, LSOL ) )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, SIZE1, SIZE2
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 4, MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( LCB + LSOL, MPI_DOUBLE_PRECISION, COMM,
     *                      SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IFATH, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LCB  , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LSOL , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( CB,    LCB, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( LSOL .GT. 0 ) THEN
          CALL MPI_PACK( SOL  , LSOL, MPI_DOUBLE_PRECISION,
     *                          BUF_CB%CONTENT( IPOS ), SIZE,
     *                          POSITION, COMM, IERR )
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, MASTER2SLAVE, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_send_master2slave: SIZE, POSITION = ',
     *               SIZE, POSITION
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_72
        SUBROUTINE DMUMPS_78( NODE1, NODE2, NCB, LONG,
     *             IW, W,
     *             DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER DEST, TAG, COMM, IERR
        INTEGER NODE1, NODE2, NCB, LONG
        INTEGER IW( MAX( 1, LONG ) )
        DOUBLE PRECISION W( MAX( 1, LONG ) )
        INCLUDE 'mpif.h'
        INTEGER POSITION, IREQ, IPOS
        INTEGER SIZE1, SIZE2, SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        IF ( NODE2 .EQ. 0 ) THEN
         CALL MPI_PACK_SIZE( 2+LONG, MPI_INTEGER, COMM, SIZE1, IERR )
        ELSE
         CALL MPI_PACK_SIZE( 4+LONG, MPI_INTEGER, COMM, SIZE1, IERR )
        END IF
        SIZE2 = 0
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK_SIZE( LONG, MPI_DOUBLE_PRECISION,
     *                        COMM, SIZE2, IERR )
        END IF
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
        POSITION = 0
        CALL MPI_PACK( NODE1, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( NODE2 .NE. 0 ) THEN
          CALL MPI_PACK( NODE2, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
          CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        END IF
        CALL MPI_PACK( LONG,  1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK( IW, LONG, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
          CALL MPI_PACK( W,  LONG, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_78
        SUBROUTINE DMUMPS_62( I, DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER I
        INTEGER DEST, TAG, COMM, IERR
        INCLUDE 'mpif.h'
        INTEGER IPOS, IREQ, MSG_SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        MSG_SIZE = SIZEofINT
        CALL DMUMPS_4( BUF_SMALL, IPOS, IREQ, MSG_SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN 
         write(6,*) ' Internal error in DMUMPS_62',
     *       ' Buf size (bytes)= ',BUF_SMALL%LBUF
         RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS ) = I
        CALL MPI_ISEND( BUF_SMALL%CONTENT(IPOS), MSG_SIZE,
     *                  MPI_PACKED, DEST, TAG, COMM,
     *                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_62
        SUBROUTINE DMUMPS_469(FLAG)
        LOGICAL FLAG
        LOGICAL FLAG1, FLAG2, FLAG3
        CALL DMUMPS_468( BUF_SMALL, FLAG1 )
        CALL DMUMPS_468( BUF_CB, FLAG2 )
        CALL DMUMPS_468( BUF_LOAD, FLAG3 )
        FLAG = FLAG1 .AND. FLAG2 .AND. FLAG3
        RETURN
        END SUBROUTINE DMUMPS_469
        SUBROUTINE DMUMPS_468( B, FLAG )
        TYPE ( DMUMPS_BUFFER_TYPE ) :: B
        LOGICAL :: FLAG
        INTEGER SIZE_AVAIL
        CALL DMUMPS_79(B, SIZE_AVAIL)
        FLAG = ( B%HEAD == B%TAIL )
        RETURN
        END SUBROUTINE DMUMPS_468
        SUBROUTINE DMUMPS_79( B, SIZE_AV )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: B
        INTEGER SIZE_AV
        INCLUDE 'mpif.h'
        INTEGER IERR
        INTEGER STATUS( MPI_STATUS_SIZE )
        LOGICAL FLAG
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
          SIZE_AV = MAX( B%LBUF_INT - B%TAIL + 1, B%HEAD - 1 )
        ELSE
          SIZE_AV = B%HEAD - B%TAIL
        END IF
        SIZE_AV = SIZE_AV * SIZEofINT
        RETURN
        END SUBROUTINE DMUMPS_79
        SUBROUTINE DMUMPS_4( B, IPOS, IREQ, MSG_SIZE, IERR, 
     &    NDEST , PDEST
     &         )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: B
        INTEGER, INTENT(IN)        :: MSG_SIZE
        INTEGER, INTENT(OUT)       :: IPOS, IREQ, IERR
        INTEGER NDEST
        INTEGER, INTENT(IN)        :: PDEST(MAX(1,NDEST))
        INCLUDE 'mpif.h'
        INTEGER MSG_SIZE_INT
        INTEGER IBUF
        LOGICAL FLAG
        INTEGER STATUS( MPI_STATUS_SIZE )
        IERR = 0
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END iF
        MSG_SIZE_INT = ( MSG_SIZE + ( SIZEofINT - 1 ) ) / SIZEofINT
        MSG_SIZE_INT = MSG_SIZE_INT + OVHSIZE
        FLAG = (     ( B%HEAD .LE. B%TAIL )
     *               .AND. (
     *                 ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL )
     *                 .OR. ( MSG_SIZE_INT .LE. B%HEAD - 2 ) ) )
     *         .OR.
     *               ( ( B%HEAD .GT. B%TAIL )
     *               .AND. ( MSG_SIZE_INT .LE. B%HEAD - B%TAIL - 1 ) )
        IF ( .NOT. FLAG
     *    ) THEN
          IERR = -1
          IF ( MSG_SIZE_INT .GT. B%LBUF_INT ) IERR = -2
          IPOS = -1
          IREQ = -1
          RETURN
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
          IF ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL + 1 ) THEN
            IBUF = B%TAIL
          ELSE IF ( MSG_SIZE_INT .LE. B%HEAD - 1 ) THEN
            IBUF = 1
          END IF
        ELSE
          IBUF = B%TAIL
        END IF
        B%CONTENT( B%ILASTMSG + NEXT ) = IBUF
        B%ILASTMSG = IBUF
        B%TAIL = IBUF + MSG_SIZE_INT
        B%CONTENT( IBUF + NEXT ) = 0
        IPOS = IBUF + CONTENT
        IREQ = IBUF + REQ
        RETURN
        END SUBROUTINE DMUMPS_4
        SUBROUTINE DMUMPS_1( BUF, SIZE )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: BUF
        INTEGER SIZE
        INTEGER SIZE_INT
        SIZE_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        SIZE_INT = SIZE_INT + OVHSIZE
        BUF%TAIL = BUF%ILASTMSG + SIZE_INT
        RETURN
        END SUBROUTINE DMUMPS_1
      SUBROUTINE DMUMPS_68(
     *             INODE, NBPROCFILS, NLIG, ILIG, NCOL, ICOL,
     *             NASS, NSLAVES, LIST_SLAVES,
     *             DEST, NFRONT, COMM, IERR )
      IMPLICIT NONE
        INTEGER COMM, IERR, NFRONT
        INTEGER INODE
        INTEGER NLIG, NCOL, NASS, NSLAVES
        INTEGER NBPROCFILS, DEST
        INTEGER ILIG( NLIG )
        INTEGER ICOL( NCOL )
        INTEGER LIST_SLAVES( NSLAVES )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE = ( 6 + NLIG + NCOL + NSLAVES + 1 ) * SIZEofINT
        IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
         IERR = -2
         RETURN
        END IF
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = INODE
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NBPROCFILS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NLIG
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NCOL
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NASS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NFRONT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = 
     *  LIST_SLAVES( 1: NSLAVES )
        POSITION = POSITION + NSLAVES
        BUF_CB%CONTENT( POSITION:POSITION + NLIG - 1 ) = ILIG
        POSITION = POSITION + NLIG
        BUF_CB%CONTENT( POSITION:POSITION + NCOL - 1 ) = ICOL
        POSITION = POSITION + NCOL
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
          WRITE(*,*) 'Error in DMUMPS_68 :',
     *               ' wrong estimated size'
          CALL DMUMPS_ABORT()
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE, MPI_PACKED,
     *                  DEST, MAITRE_DESC_BANDE, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_68
        SUBROUTINE DMUMPS_70( IPERE, ISON, NROW,
     *  IROW, NCOL, ICOL, VAL, LDA, NELIM, TYPE_SON,
     *  NSLAVES, SLAVES, DEST, COMM, IERR, 
     * 
     *  SLAVEF, KEEP, INIV2, TAB_POS_IN_PERE )
        IMPLICIT NONE
        INTEGER LDA, NELIM, TYPE_SON
        INTEGER IPERE, ISON, NROW, NCOL, NSLAVES
        INTEGER IROW( NROW )
        INTEGER ICOL( NCOL )
        INTEGER SLAVES( NSLAVES )
        DOUBLE PRECISION VAL(LDA, *)
        INTEGER IPOS, IREQ, DEST, COMM, IERR
        INTEGER SLAVEF, KEEP(150), INIV2
        INTEGER TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE1, SIZE2, SIZE3, SIZE, POSITION, I
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        IF ( NELIM .NE. NROW ) THEN
          WRITE(*,*) 'Error in TRY_SEND_MAITRE2:',NELIM, NROW
          CALL DMUMPS_ABORT()
        END IF
        CALL MPI_PACK_SIZE( NROW+NCOL+5+NSLAVES, MPI_INTEGER,
     *                      COMM, SIZE1, IERR )
        IF ( ( KEEP(48).NE. 0 ).AND.(TYPE_SON .eq. 2 )) THEN
          CALL MPI_PACK_SIZE( NSLAVES+1, MPI_INTEGER,
     *                        COMM, SIZE3, IERR )
        ELSE
          SIZE3 = 0
        ENDIF
        IF ( KEEP(50).ne.0  .AND. TYPE_SON .eq. 2 ) THEN
          CALL MPI_PACK_SIZE( NROW * NROW, MPI_DOUBLE_PRECISION,
     *                        COMM, SIZE2, IERR )
        ELSE
          CALL MPI_PACK_SIZE( NROW * NCOL, MPI_DOUBLE_PRECISION,
     *                        COMM, SIZE2, IERR )
        END IF
        SIZE = SIZE1 + SIZE3 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
        IF (NELIM.NE.0) THEN
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -2 
             RETURN
          ENDIF
        ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( ISON,  1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        IF (NSLAVES.GT.0) THEN
         CALL MPI_PACK( SLAVES, NSLAVES, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        ENDIF
        CALL MPI_PACK( IROW, NROW, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( ICOL, NCOL, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        IF ( ( KEEP(48).NE. 0 ).AND.(TYPE_SON .eq. 2 ) ) THEN
         CALL MPI_PACK( TAB_POS_IN_PERE(1,INIV2), NSLAVES+1, 
     *                 MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        ENDIF
        IF (NROW.GE.1) THEN
          IF ( KEEP(50).ne.0  .AND. TYPE_SON .eq. 2 ) THEN
            DO I=1,NROW
              CALL MPI_PACK( VAL(1,I), NROW, 
     *                 MPI_DOUBLE_PRECISION,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
            ENDDO
          ELSE
            DO I=1,NROW
              CALL MPI_PACK( VAL(1,I), NCOL, 
     *                 MPI_DOUBLE_PRECISION,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
            ENDDO
          ENDIF
        ENDIF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, MAITRE2, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          write(*,*) 'Try_send_maitre2, SIZE,POSITION=',SIZE,POSITION
          Stop
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_70
        SUBROUTINE DMUMPS_67(DESC_IN_LU,
     *  IPERE, NFRONT_PERE, NASS_PERE, NSLAVES_PERE,
     *  ISON, NBROW, LMAP, MAPROW, IW_CBSON, A_CBSON,
     *  ISLAVE, PDEST, PDEST_MASTER, COMM, IERR, 
     *  
     * KEEP, STEP, N, SLAVEF,
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *   )
        IMPLICIT NONE
        INTEGER IPERE, ISON, NBROW 
        INTEGER PDEST, ISLAVE, COMM, IERR
        INTEGER PDEST_MASTER, NASS_PERE, NSLAVES_PERE,
     *          NFRONT_PERE, LMAP
        INTEGER MAPROW( LMAP )
        INTEGER IW_CBSON( * )
        DOUBLE PRECISION A_CBSON( * )
        LOGICAL DESC_IN_LU
       INTEGER   KEEP(150), N , SLAVEF
       INTEGER   STEP(N), 
     *          ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        INTEGER LROW
        INTEGER NPIV, NFRONT, HS
        INTEGER SIZE, SIZE1, SIZE2, POSITION,I
        INTEGER APOS, IPOS_IN_SLAVE
        INTEGER INDICE_PERE, NROW, IPOS, IREQ, NOSLA
        INTEGER IONE
        INTEGER PDEST2(1)
        PARAMETER ( IONE=1 )
        PDEST2(1) = PDEST
        IERR   = 0
        LROW   = IW_CBSON( 1 )
        NPIV   = IW_CBSON( 4 )
        IF ( NPIV .LT. 0 ) THEN
          NPIV = 0
        END IF
        NROW   = IW_CBSON( 3 )
        NFRONT = LROW + NPIV
        HS     = 6 + IW_CBSON( 6 )
        IF ( PDEST .EQ. PDEST_MASTER ) THEN
          CALL MPI_PACK_SIZE( NBROW + 4, MPI_INTEGER,
     *                        COMM, SIZE1, IERR )
        ELSE
          CALL MPI_PACK_SIZE( NBROW + 4 + LROW, MPI_INTEGER,
     *                        COMM, SIZE1, IERR )
        END IF
        CALL MPI_PACK_SIZE( LROW * NBROW, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , PDEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -2
             RETURN
          ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NBROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( LROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        IF ( PDEST .NE. PDEST_MASTER ) THEN
          CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ), LROW,
     *                 MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        END IF
        DO I = 1, LMAP
          INDICE_PERE = MAPROW( I )
          CALL DMUMPS_47(
     *         KEEP, IPERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *         NASS_PERE,
     *         NFRONT_PERE - NASS_PERE,
     *         NSLAVES_PERE,
     *         INDICE_PERE,
     *         NOSLA,
     *         IPOS_IN_SLAVE )
          INDICE_PERE = IPOS_IN_SLAVE
          IF ( NOSLA .EQ. ISLAVE .OR. NBROW .EQ. LMAP ) THEN
            CALL MPI_PACK( INDICE_PERE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
            IF (DESC_IN_LU) THEN 
               APOS = ( I - 1 ) * LROW + 1
            ELSE
               APOS = ( I - 1 ) * NFRONT + NPIV + 1
            ENDIF
            CALL MPI_PACK( A_CBSON( APOS ), LROW, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
          END IF
        END DO
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  PDEST, CONTSIPERENIV2, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' contniv2: SIZE, POSITION =',SIZE, POSITION
          WRITE(*,*) ' NBROW, LROW = ', NBROW, LROW
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_67
        SUBROUTINE DMUMPS_71( 
     *                INODE, NFRONT, NASS1,
     *                ISON, MYID, NSLAVES, SLAVES_PERE,
     *                TROW, NCBSON,
     *                COMM, IERR,
     *                DEST, NDEST, SLAVEF, 
     * 
     *                KEEP, STEP, N, 
     *                ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *
     *                                  )
        IMPLICIT NONE
      INTEGER INODE, NFRONT, NASS1, NCBSON, NSLAVES, 
     *          NDEST
      INTEGER SLAVEF, MYID, ISON
      INTEGER TROW( NCBSON )
      INTEGER DEST( NDEST )
      INTEGER SLAVES_PERE( NSLAVES )
      INTEGER COMM, IERR
      INTEGER KEEP(150), N
      INTEGER STEP(N), 
     *        ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        INTEGER SIZE_AV, IDEST, NSEND, SIZE
        INTEGER TROW_SIZE, POSITION, INDX, INIV2
        INTEGER IPOS, IREQ
        INTEGER IONE
        PARAMETER ( IONE=1 )
        INTEGER NASS_SON
        NASS_SON = -99998
        IERR = 0
        IF ( NDEST .eq. 1 ) THEN
          IF ( DEST(1).EQ.MYID )  GOTO 500
          SIZE = SIZEofINT * ( 6 + NSLAVES + NCBSON )
          IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
           SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
          ENDIF
          CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST
     &                 )
          IF (IERR .LT. 0 ) RETURN
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -2
             RETURN
          END IF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NCBSON
              POSITION = POSITION + 1
              IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     *          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     *          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+NCBSON-1 ) =
     *        TROW( 1: NCBSON )
              POSITION = POSITION + NCBSON
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
                WRITE(*,*) 'Error in DMUMPS_71 :',
     *                     ' wrong estimated size'
                CALL DMUMPS_ABORT()
              END IF
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     *                        MPI_PACKED,
     *                        DEST( NDEST ), MAPLIG, COMM,
     *                        BUF_CB%CONTENT( IREQ ),
     *                        IERR )
        ELSE
          NSEND = 0
          DO IDEST = 1, NDEST
            IF ( DEST( IDEST ) .ne. MYID ) NSEND = NSEND + 1
          END DO
          SIZE = SIZEofINT * ( ( OVHSIZE + 6 + NSLAVES )
     *                           * NSEND + NCBSON )
          IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
           SIZE = SIZE + SIZEofINT * NSEND*( NSLAVES + 1 )
          ENDIF
          CALL DMUMPS_79( BUF_CB, SIZE_AV )
          IF ( SIZE_AV .LT. SIZE ) THEN
            IERR = -1
            RETURN
          END IF
          DO IDEST= 1, NDEST
            CALL DMUMPS_49( 
     &                KEEP, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                IDEST, NCBSON, 
     &                NDEST, 
     &                TROW_SIZE, INDX  )
            SIZE = SIZEofINT * ( NSLAVES + TROW_SIZE + 6 )
            IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
            ENDIF
            IF ( MYID .NE. DEST( IDEST ) ) THEN
              CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR,
     &                       IONE, DEST(IDEST)
     &                     )
              IF ( IERR .LT. 0 )  THEN
                 WRITE(*,*) 'Problem in DMUMPS_4: IERR<0'
                 CALL DMUMPS_ABORT()
              END IF
              IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
                 IERR = -2
                 RETURN
              ENDIF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = TROW_SIZE
              POSITION = POSITION + 1
              IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     *          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     *          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+TROW_SIZE-1 ) =
     *        TROW( INDX: INDX + TROW_SIZE - 1 )
              POSITION = POSITION + TROW_SIZE
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
               WRITE(*,*) ' ERROR 1 in TRY_SEND_MAPLIG:',
     *          'Wrong estimated size'
               CALL DMUMPS_ABORT()
              END IF
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     *                        MPI_PACKED,
     *                        DEST( IDEST ), MAPLIG, COMM,
     *                        BUF_CB%CONTENT( IREQ ),
     *                        IERR )
            END IF
          END DO
        END IF
 500    RETURN
        END SUBROUTINE DMUMPS_71
        SUBROUTINE DMUMPS_65( INODE, NFRONT,
     *             NCOL, NPIV, FPERE, LASTBL, IPIV, VAL,
     *             PDEST, NDEST, KEEP50, NB_BLOC_FAC, COMM, IERR )
        IMPLICIT NONE
        INTEGER INODE, NCOL, NPIV, FPERE, NFRONT, NDEST
        INTEGER IPIV( NPIV )
        DOUBLE PRECISION VAL( NFRONT, * )
        INTEGER PDEST( NDEST ) 
        INTEGER   KEEP50, NB_BLOC_FAC, COMM, IERR
        LOGICAL LASTBL
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE,
     *          IDEST, IPOSMSG, I
        INTEGER NPIVSENT
        INTEGER SSS, SS2
        IERR = 0
        IF ( LASTBL ) THEN
          IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                          MPI_INTEGER, COMM, SIZE1, IERR )
          ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                          MPI_INTEGER, COMM, SIZE1, IERR )
          END IF
        ELSE
          IF ( KEEP50 .eq. 0 ) THEN
          CALL MPI_PACK_SIZE( 3 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                        MPI_INTEGER, COMM, SIZE1, IERR )
          ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                          MPI_INTEGER, COMM, SIZE1, IERR )
          END IF
        END IF
        SIZE2 = 0
        IF (NPIV.GT.0) 
     *    CALL MPI_PACK_SIZE( NPIV*NCOL, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 NDEST , PDEST
     &               )
        IF ( IERR .LT. 0 ) RETURN
        IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
          SSS = 0 
          IF ( LASTBL ) THEN
           IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV ,
     *                        MPI_INTEGER, COMM, SSS, IERR )
           ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV , 
     *                           MPI_INTEGER, COMM, SSS, IERR )
           END IF
          ELSE
           IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 3 + NPIV ,
     *                        MPI_INTEGER, COMM, SSS, IERR )
           ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV ,
     *                        MPI_INTEGER, COMM, SSS, IERR )
           END IF
          END IF
          IF (NPIV.GT.0) 
     *    CALL MPI_PACK_SIZE( NPIV*NCOL, MPI_DOUBLE_PRECISION,
     *                      COMM, SS2, IERR )
          SSS = SSS + SS2
          IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -2 
           RETURN
          ENDIF
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        NPIVSENT = NPIV
        IF (LASTBL) NPIVSENT = -NPIV
        CALL MPI_PACK( NPIVSENT, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( LASTBL .or. KEEP50.ne.0 ) THEN
          CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     *                   BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
        END IF
        IF ( LASTBL .AND. KEEP50 .NE. 0 ) THEN
            CALL MPI_PACK( NDEST, 1, MPI_INTEGER,
     *                   BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
            CALL MPI_PACK( NB_BLOC_FAC, 1, MPI_INTEGER,
     *                   BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
        END IF
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( NPIV.GT.0) THEN
          CALL MPI_PACK( IPIV, NPIV, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
          DO I = 1, NPIV
            CALL MPI_PACK( VAL(1,I), NCOL,
     *                        MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
          END DO
        ENDIF
        DO IDEST = 1, NDEST
        IF ( KEEP50.eq.0) THEN
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     *                  PDEST(IDEST), BLOC_FACTO, COMM,
     *                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     *                  IERR )
        ELSE
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     *                  PDEST(IDEST), BLOC_FACTO_SYM, COMM,
     *                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     *                  IERR )
        END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blocfacto : size < position'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_65
        SUBROUTINE DMUMPS_64( INODE,
     *             NPIV, FPERE, IPOSK, JPOSK, UIP21K, NCOLU,
     *             NDEST, PDEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER INODE, NCOLU, IPOSK, JPOSK, NPIV, NDEST, FPERE
        DOUBLE PRECISION UIP21K( NPIV, NCOLU )
        INTEGER PDEST( NDEST ) 
        INTEGER   COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE,
     *          IDEST, IPOSMSG, SSS, SS2
        IERR = 0
        CALL MPI_PACK_SIZE( 6 + ( NDEST - 1 ) * OVHSIZE,
     *                      MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( ABS(NPIV)*NCOLU, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
         CALL MPI_PACK_SIZE( 6 ,
     *                      MPI_INTEGER, COMM, SSS, IERR )
         CALL MPI_PACK_SIZE( ABS(NPIV)*NCOLU, MPI_DOUBLE_PRECISION,
     *                      COMM, SS2, IERR )
         SSS = SSS+SS2
         IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -2
           RETURN
         ENDIF
        END IF
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 NDEST, PDEST
     &               )
        IF ( IERR .LT. 0 ) RETURN
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IPOSK, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( JPOSK, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( NPIV, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( NCOLU, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( UIP21K, ABS(NPIV) * NCOLU,
     *                        MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        DO IDEST = 1, NDEST
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     *                  PDEST(IDEST), BLOC_FACTO_SYM_SLAVE, COMM,
     *                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     *                  IERR )
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blfac slave : size < position'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_64
        SUBROUTINE DMUMPS_75( N, ISON,
     *             NBCOL_SON, NBROW_SON, INDCOL_SON, INDROW_SON,
     *             LD_SON, VAL_SON, TAG, SUBSET_ROW, SUBSET_COL,
     *             NSUBSET_ROW, NSUBSET_COL,
     *             NPROW, NPCOL, MBLOCK, RG2L_ROW, RG2L_COL,
     *             NBLOCK, PDEST, COMM, IERR , 
     *             TAB, TABSIZE, INVERT, SIZE_MSG )
        IMPLICIT NONE
        INTEGER N, ISON, NBCOL_SON, NBROW_SON, NSUBSET_ROW, NSUBSET_COL
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK, LD_SON
        INTEGER PDEST, TAG, COMM, IERR
        INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
        INTEGER SUBSET_ROW( * ), SUBSET_COL( * )
        INTEGER RG2L_ROW( N ), RG2L_COL( N ) 
        INTEGER   TABSIZE, SIZE_MSG
        DOUBLE PRECISION VAL_SON( LD_SON, * ), TAB(*)
        LOGICAL INVERT
        INCLUDE 'mpif.h'
        INTEGER SIZE1, SIZE2, SIZE, POSITION
        INTEGER IREQ, IPOS, ITAB
        INTEGER ISUB, JSUB, I, J 
        INTEGER ILOC_ROOT, JLOC_ROOT
        INTEGER IPOS_ROOT, JPOS_ROOT
        INTEGER IONE
        INTEGER PDEST2(1)
        PARAMETER ( IONE=1 )
        PDEST2(1) = PDEST
        IERR = 0
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          CALL MPI_PACK_SIZE( 3 + NSUBSET_ROW + NSUBSET_COL,
     *                      MPI_INTEGER, COMM, SIZE1, IERR )
          CALL MPI_PACK_SIZE( NSUBSET_ROW * NSUBSET_COL,
     *                      MPI_DOUBLE_PRECISION, COMM, SIZE2, IERR )
          SIZE = SIZE1 + SIZE2
        ELSE
          CALL MPI_PACK_SIZE(3,MPI_INTEGER, COMM, SIZE, IERR )
        END IF
        SIZE_MSG = SIZE
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, PDEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
        IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -2
             RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        IF ( .NOT. INVERT ) THEN
        CALL MPI_PACK( NSUBSET_ROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUBSET_COL, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        ELSE
        CALL MPI_PACK( NSUBSET_COL, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUBSET_ROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        END IF
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          IF ( .NOT. INVERT ) THEN
            DO ISUB = 1, NSUBSET_ROW
              I         = SUBSET_ROW( ISUB )
              IPOS_ROOT = RG2L_ROW( INDROW_SON( I ) )
              ILOC_ROOT = MBLOCK
     *                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     *                 + MOD( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     *                      BUF_CB%CONTENT( IPOS ),
     *                      SIZE, POSITION, COMM, IERR )
            END DO
            DO JSUB = 1, NSUBSET_COL
              J         = SUBSET_COL( JSUB )
              JPOS_ROOT = RG2L_COL( INDCOL_SON( J ) )
              JLOC_ROOT = NBLOCK
     *                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     *                  + MOD( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     *                       BUF_CB%CONTENT( IPOS ),
     *                       SIZE, POSITION, COMM, IERR )
            END DO
          ELSE
            DO JSUB = 1, NSUBSET_COL
              J         = SUBSET_COL( JSUB )
              IPOS_ROOT = RG2L_COL( INDCOL_SON( J ) )
              ILOC_ROOT = MBLOCK
     *                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     *                 + MOD( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     *                       BUF_CB%CONTENT( IPOS ),
     *                       SIZE, POSITION, COMM, IERR )
            END DO
            DO ISUB = 1, NSUBSET_ROW
              I         = SUBSET_ROW( ISUB )
              JPOS_ROOT = RG2L_ROW( INDROW_SON( I ) )
              JLOC_ROOT = NBLOCK
     *                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     *                  + MOD( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     *                      BUF_CB%CONTENT( IPOS ),
     *                      SIZE, POSITION, COMM, IERR )
            END DO
          END IF
          IF ( TABSIZE.GE.NSUBSET_ROW*NSUBSET_COL ) THEN
            IF ( .NOT. INVERT ) THEN
              ITAB = 1
              DO ISUB = 1, NSUBSET_ROW
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL
                  J              = SUBSET_COL( JSUB )
                  TAB( ITAB )    = VAL_SON(J,I)
                  ITAB           = ITAB + 1
                END DO
               END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL*NSUBSET_ROW,
     *         MPI_DOUBLE_PRECISION, 
     *         BUF_CB%CONTENT( IPOS ),
     *         SIZE, POSITION, COMM, IERR )
            ELSE
              ITAB = 1
              DO JSUB = 1, NSUBSET_COL
                J = SUBSET_COL( JSUB )
                DO ISUB = 1, NSUBSET_ROW
                  I         = SUBSET_ROW( ISUB )
                  TAB( ITAB ) = VAL_SON( J, I )
                  ITAB = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL*NSUBSET_ROW,
     *         MPI_DOUBLE_PRECISION, 
     *         BUF_CB%CONTENT( IPOS ),
     *         SIZE, POSITION, COMM, IERR )
            END IF
          ELSE
            IF ( .NOT. INVERT ) THEN
              DO ISUB = 1, NSUBSET_ROW
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL
                  J         = SUBSET_COL( JSUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     *            MPI_DOUBLE_PRECISION,
     *            BUF_CB%CONTENT( IPOS ),
     *            SIZE, POSITION, COMM, IERR )
                END DO
              END DO
            ELSE
              DO JSUB = 1, NSUBSET_COL
                J = SUBSET_COL( JSUB )
                DO ISUB = 1, NSUBSET_ROW
                  I         = SUBSET_ROW( ISUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     *            MPI_DOUBLE_PRECISION,
     *            BUF_CB%CONTENT( IPOS ),
     *            SIZE, POSITION, COMM, IERR )
                END DO
              END DO
            END IF
          ENDIF
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                PDEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending contribution to root:Size<positn'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .LT. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_75
        SUBROUTINE DMUMPS_76( ISON, NELIM,
     *             NELIM_ROW, NELIM_COL, NSLAVES, SLAVES,
     *             DEST, COMM, IERR )
        INTEGER ISON, NELIM
        INTEGER NSLAVES, DEST, COMM, IERR
        INTEGER NELIM_ROW( NELIM ), NELIM_COL( NELIM )
        INTEGER SLAVES( NSLAVES )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE = ( 3 + NSLAVES + 2 * NELIM ) * SIZEofINT
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
          IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
             IERR = -2
             RETURN
          ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = ISON
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NELIM
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_ROW
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_COL
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = SLAVES
        POSITION = POSITION + NSLAVES
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
          WRITE(*,*) 'Error in DMUMPS_BUF_SEND_ROOT_NELIM_INDICES:',
     *               'wrong estimated size'
           CALL DMUMPS_ABORT()
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE, 
     *                  MPI_PACKED,
     *                  DEST, ROOT_NELIM_INDICES, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_76
        SUBROUTINE DMUMPS_74( ISON, NELIM_ROOT,
     *             DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER ISON, NELIM_ROOT, DEST, COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER IPOS, IREQ, SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        SIZE = 2 * SIZEofINT
  10    CONTINUE
        CALL DMUMPS_4( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Increase small buffer size. Trying again.'
          GOTO 10
        END IF
        IF ( IERR .LT. 0 ) RETURN
        BUF_SMALL%CONTENT( IPOS )     = ISON
        BUF_SMALL%CONTENT( IPOS + 1 ) = NELIM_ROOT
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     *                  MPI_PACKED,
     *                  DEST, ROOT_2SON, COMM,
     *                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_74
        SUBROUTINE DMUMPS_73
     *  ( TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        SIZE = 2 * SIZEofINT
        CALL DMUMPS_4( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           WRITE(*,*)  'Increase small buf. size'
           CALL DMUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) RETURN
        BUF_SMALL%CONTENT( IPOS     ) = TOT_ROOT_SIZE
        BUF_SMALL%CONTENT( IPOS + 1 ) = TOT_CONT2RECV
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     *                  MPI_PACKED,
     *                  DEST, ROOT_2SLAVE, COMM,
     *                  BUF_SMALL%CONTENT( IREQ ), IERR )
        END SUBROUTINE DMUMPS_73
        SUBROUTINE DMUMPS_63
     *             ( INODE, W, LW, DEST,MSGTAG,COMM,IERR )
        IMPLICIT NONE
        INTEGER INODE,LW,COMM,IERR,DEST,MSGTAG
        DOUBLE PRECISION W(LW)
        INCLUDE 'mpif.h'
        INTEGER SIZE, SIZE1, SIZE2
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        CALL MPI_PACK_SIZE( 2 , MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( LW, MPI_DOUBLE_PRECISION, COMM,
     *                      SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) RETURN
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LW   , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( W    , LW, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, MSGTAG, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_update: SIZE, POSITION = ',
     *               SIZE, POSITION
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_63
        SUBROUTINE DMUMPS_77
     *             ( BDC_MEM, COMM, NPROCS, LOAD,
     *               MEM,
#if ! defined(OLD_LOAD_MECHANISM)
     *               FUTURE_NIV2,
#endif
     *               MYID, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        DOUBLE PRECISION LOAD
        DOUBLE PRECISION MEM
        LOGICAL BDC_MEM
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, WHAT, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#endif
        IF ( NDEST .eq. 0 ) RETURN
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE1, IERR )
        NREALS = 1
        IF (BDC_MEM) THEN
          NREALS = 2
        ENDIf
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) RETURN
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        WHAT = 0  ! 0: absolute load, only sender makes the update
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM, 1, MPI_DOUBLE_PRECISION,
     *                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
        END IF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
        IF ( I .ne. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                      POSITION, MPI_PACKED, I,
     *                      UPDATE_LOAD, COMM,
     *                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_77'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        END SUBROUTINE DMUMPS_77
        SUBROUTINE DMUMPS_498
     *             ( WHAT, COMM, NPROCS, LOAD,
     *               MYID, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT
        DOUBLE PRECISION LOAD
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        IF (WHAT .NE. 2 .AND. WHAT .NE. 3) THEN
          WRITE(*,*)
     *  "Internal error 1 in DMUMPS_498"
        END IF
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        IF ( NDEST .eq. 0 ) RETURN
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE1, IERR )
        NREALS = 1
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) RETURN
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        IDEST = 0
        DO I = 0, NPROCS - 1
          IF ( I .ne. MYID ) THEN
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                      POSITION, MPI_PACKED, I,
     *                      UPDATE_LOAD, COMM,
     *                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_BUF_SEND_UPDATE_POOL'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_498
        SUBROUTINE DMUMPS_502( COMM, MYID, NPROCS, IERR)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER IPOS, IREQ, IDEST, IPOSMSG, POSITION, I
        INTEGER COMM, MYID, IERR, NPROCS
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        INTEGER NDEST, NINTS, SIZE
        INTEGER WHAT
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NINTS = 1 + ( NDEST-1 ) * OVHSIZE
        CALL MPI_PACK_SIZE( NINTS, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE, IERR )
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     *       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) RETURN
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        WHAT = 4
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        IDEST = 0
        DO I = 0, NPROCS - 1
          IF ( I .NE. MYID ) THEN
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                      POSITION, MPI_PACKED, I,
     *                      UPDATE_LOAD, COMM,
     *                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_460'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_502
        SUBROUTINE DMUMPS_460( BDC_MEM,
     *      COMM, MYID, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     *      FUTURE_NIV2,
#endif
     *      NSLAVES,
     *      LIST_SLAVES,
     *      MEM_INCREMENT, FLOPS_INCREMENT, IERR )
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        LOGICAL BDC_MEM
        INTEGER COMM, MYID, NPROCS, NSLAVES, IERR
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        INTEGER LIST_SLAVES(NSLAVES)
        DOUBLE PRECISION MEM_INCREMENT(NSLAVES)
        DOUBLE PRECISION FLOPS_INCREMENT(NSLAVES)
        INTEGER NDEST, NINTS, NREALS, SIZE1, SIZE2, SIZE
        INTEGER IPOS, IPOSMSG, IREQ, POSITION
        INTEGER I, IDEST, WHAT
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1)=MYID
        IERR = 0
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#else
        NDEST = NPROCS - 1
#endif
        IF ( NDEST == 0 ) RETURN
        NINTS = 2 +  NSLAVES + ( NDEST - 1 ) * OVHSIZE
        NREALS = NSLAVES
        IF (BDC_MEM) NREALS = NREALS + NSLAVES
        CALL MPI_PACK_SIZE( NINTS, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE1, IERR )
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     *       COMM, SIZE2, IERR )
        SIZE = SIZE1+SIZE2
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     *       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) RETURN
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        WHAT = 1
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( LIST_SLAVES, NSLAVES, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( FLOPS_INCREMENT, NSLAVES,
     *      MPI_DOUBLE_PRECISION,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM_INCREMENT, NSLAVES,
     *      MPI_DOUBLE_PRECISION,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        END IF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
        IF ( I .NE. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                      POSITION, MPI_PACKED, I,
     *                      UPDATE_LOAD, COMM,
     *                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_460'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_460
        SUBROUTINE DMUMPS_60
     *             ( DMUMPS_LBUFR_BYTES)
        IMPLICIT NONE
        INTEGER DMUMPS_LBUFR_BYTES 
        SIZE_RBUF_BYTES = DMUMPS_LBUFR_BYTES
        RETURN
      END SUBROUTINE DMUMPS_60
      END MODULE DMUMPS_BUFFER
