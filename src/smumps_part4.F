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
      SUBROUTINE SMUMPS_246(MYID, N, STEP, FRERE, FILS,
     *     NA, LNA, NE, ND, PROCNODE, SLAVEF,
     *     NRLADU, NIRADU, NIRNEC, NRLNEC, NRLBDU, NIRBDU,
     *     MAXFR, OPSA, PERLU, 
     *     KEEP, LOCAL_M, LOCAL_N,
     *     SBUF_SEND, SBUF_REC, OPS_SUBTREE, NSTEPS,
     *     NMB_PAR2, ISTEP_TO_INIV2, CANDIDATES, 
     *     IFLAG, IERROR )
      IMPLICIT NONE
      INTEGER  MYID, N, LNA, IFLAG, IERROR
      INTEGER  NRLADU, NIRADU, NRLNEC, NIRNEC
      INTEGER  MAXFR, PERLU, NSTEPS
      INTEGER  STEP(N)
      INTEGER  FRERE(NSTEPS),FILS(N),NA(LNA),NE(NSTEPS),
     *         ND(NSTEPS),PROCNODE(NSTEPS)
      INTEGER  SLAVEF, KEEP(150), LOCAL_M, LOCAL_N
      INTEGER  SBUF_SEND, SBUF_REC
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      INTEGER  CANDIDATES( SLAVEF+1, MAX(NMB_PAR2,1) ), INFO
      INTEGER  NRLBDU, NIRBDU
      REAL OPSA
      REAL OPS_SUBTREE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSTKR, TNSTK, IPOOL, 
     *                                      LSTKI 
      INTEGER I,ISTKR,ITOP,NELIM,NFR
      INTEGER ISTKI, STKI, ISTKIM
      INTEGER K,LSTK,NSTK, IFATH
      INTEGER INODE, LEAF, NBLEAF, NBROOT, IN
      INTEGER LEVEL, MAXTEMPCB, MAXITEMPCB
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF 
      INTEGER LEVELF, NCB, SIZECB, SIZECBI
      INTEGER NBROWMAX, NSLAVES, TEMP,
     *         NELIMF, NFRF, NCBF, CBMAXR, CBMAXS,
     *         NSLAVESF, NBROWMAXF, NCOLMAX, LKJIB, 
     *         KMAX
      INTEGER LWK_QR, LIWK_QR
      INTEGER IROOT, SIZE_ROOT
      INTEGER allocok
      LOGICAL ROOT_OWNER
      REAL OPS_NODE
      INTRINSIC MIN,INT, REAL
      INTEGER SMUMPS_275, SMUMPS_330
      LOGICAL SMUMPS_170, SMUMPS_358
      EXTERNAL SMUMPS_275, SMUMPS_330, 
     *         SMUMPS_170, SMUMPS_358
      INTEGER WHAT
      INTEGER SMUMPS_52, SMUMPS_50
      EXTERNAL SMUMPS_52,SMUMPS_503,
     *         SMUMPS_50
      INTEGER IDUMMY
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
      LKJIB = MAX(KEEP(5),KEEP(6))
      TNSTK = NE
      LEAF = NA(1)+1
      IPOOL(1:LEAF-1) = NA(3:3+LEAF-2)
      NBROOT = NA(2)
      ISTKR = 0
      ISTKI = 0
      ISTKIM= 0
      OPSA = 0.0E0
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
         IF (SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *                                        .EQ. MYID) THEN
            NIRADU = 6+2*ND(STEP(INODE))
         ELSE
            NIRADU = 6
         ENDIF
         NIRNEC = NIRADU
      ENDIF
      IF((KEEP(24).eq.0).or.(KEEP(24).eq.1)) THEN
         force_cand=.FALSE.           
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      END IF
 90   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF - 1
         INODE = IPOOL(LEAF)
      ELSE 
         WRITE(MYID+6,*) ' ERROR 1 in file mpima41md '
         CALL SMUMPS_ABORT()
      ENDIF
 95   CONTINUE 
      MASTER = SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF) .EQ. MYID
      LEVEL  = SMUMPS_330(STEP(INODE),PROCNODE,SLAVEF)
      UPDATE=.FALSE.
      if(.NOT.force_cand) then
       UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
      else
         if(master.and.(level.ne.3)) then
            UPDATE = .TRUE.
         else if(level.eq.2) then
            if ( SMUMPS_358( MYID, SLAVEF, inode,
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
 100  NELIM = NELIM + 1 
      IN = FILS(IN)
      IF (IN .GT. 0 ) GOTO 100
      IFSON = -IN
      IF (LEVEL.EQ.3) THEN
         NRLNEC = MAX(NRLNEC,NRLADU+ISTKR+
     *        LOCAL_M*LOCAL_N)
         IF (MASTER) THEN 
            IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
      ENDIF
      NCB      = NFR-NELIM
      SIZECB   =  NCB*NCB
      SIZECBI  = 2* NCB  + 6 + 3
      IF (LEVEL.EQ.2) THEN
         IF (MASTER) THEN
            SBUF_SEND  = MAX(SBUF_SEND, NFR*LKJIB+LKJIB+4)
         elseif (UPDATE) then
            SBUF_REC   = MAX(SBUF_REC, NFR*LKJIB+LKJIB+4)
         ENDIF
      ENDIF
      IN = INODE
 660  IN = FRERE(STEP(IN))
      IF (IN.GT.0) GO TO 660
      IFATH = -IN
      IF ( UPDATE ) THEN
         IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
            SIZECB = NCB*NCB
            NRLADU = NRLADU + NELIM*(2*NFR-NELIM)
            NIRADU = NIRADU + 6 + 2*NFR
            SIZECBI    = 2* NCB  + 6 + 3
         ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
               NIRADU = NIRADU + 6 + SLAVEF - 1 + 2*NFR
               NRLADU = NRLADU + NFR*NELIM
               SIZECB = 0
               SIZECBI    = NCB + 5 +  SLAVEF - 1
            ELSE
               NSLAVES = SMUMPS_50(
     &                  SLAVEF, KEEP(48), KEEP(10), KEEP(50),
     &                  NFR, NCB )
               WHAT = 2
               CALL SMUMPS_503( WHAT,
     &         KEEP, NCB, NFR, SLAVEF, NBROWMAX, SIZECB )
               NIRADU   = NIRADU+4+NELIM+NBROWMAX
               NRLADU   = NRLADU + NELIM*NBROWMAX
               SIZECBI      = 4 + NBROWMAX + NCB
            ENDIF
         ENDIF
         NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
         NRLNEC = MAX(NRLNEC,NRLADU+ISTKR+SIZECB+
     *                  MAX(SIZECB,MAXTEMPCB) )
         IF (NFR.GT.MAXFR) MAXFR = NFR
         IF (NSTK.GT.0) THEN
            DO 70 K=1,NSTK
               LSTK = LSTKR(ITOP)
               ISTKR = ISTKR - LSTK
               STKI = LSTKI( ITOP )
               ISTKI = ISTKI - STKI
               ITOP = ITOP - 1
               IF (ITOP.LT.0) THEN
                  write(*,*) MYID,': ERROR 2 in mpima41md. ITOP = ',ITOP
                  CALL SMUMPS_ABORT()
               ENDIF
 70         CONTINUE
         ENDIF
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = SMUMPS_275(STEP(IFSON),
     *                  PROCNODE,SLAVEF).EQ.MYID
            LEVELSON  = SMUMPS_330(STEP(IFSON),
     *                  PROCNODE,SLAVEF)
            if(.NOT.force_cand) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     *                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( SMUMPS_358( MYID, SLAVEF, IFSON,
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
                  CALL SMUMPS_ABORT()
               ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) ! process next son
         END DO
      ENDIF
      IF (MASTER.OR.LEVEL.EQ.3) THEN
         IF ((INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0)) THEN
            OPS_NODE = REAL(2*NFR*NELIM)*REAL(NFR-NELIM-1)+
     *           REAL(NELIM*(NELIM+1))*REAL(2*NELIM+1)/REAL(3)
            OPS_NODE = OPS_NODE + REAL(((2*NFR-NELIM-1)*NELIM)/2)
         ELSE
            OPS_NODE = 0.0E0
         ENDIF
         IF (SMUMPS_170(STEP(INODE),
     *                               PROCNODE, SLAVEF) ) THEN
            OPS_SUBTREE = OPS_SUBTREE + OPS_NODE
         ENDIF
         IF ( LEVEL .EQ. 3 ) THEN
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
         LEVELF =SMUMPS_330(STEP(IFATH),PROCNODE,SLAVEF)
         MASTERF=SMUMPS_275(STEP(IFATH),PROCNODE,SLAVEF).EQ.MYID
         UPDATEF=.FALSE.
         if(.NOT.force_cand) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(masterf.and.(levelf.ne.3)) then
               UPDATEF = .TRUE.
            else if(levelf.eq.2) then
               if ( SMUMPS_358( MYID, SLAVEF, ifath,
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
               WHAT = 1
               CALL SMUMPS_503( WHAT, KEEP,
     &         NCBF, NFRF, SLAVEF, NBROWMAXF, IDUMMY )
               NCOLMAX  = NCB
               if (update.or.updatef) then 
                  CBMAXS = MAX(MIN(NELIMF,NCB)*NCOLMAX, 
     *                 MIN(NBROWMAXF,NCB)*NCOLMAX)
               else 
                  cbmaxs = 0
               end if
               IF (MASTERF) THEN
                  CBMAXR = MIN(NELIMF,NCB)*NCOLMAX
               elseif (UPDATEF) then
                  CBMAXR = MIN(NBROWMAXF,NCB)*NCOLMAX
               else
                  cbmaxr = 0
               endif
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
     *              SBUF_SEND = MAX(SBUF_SEND, SIZECB+SIZECBI)
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
                  WHAT = 2
                  CALL SMUMPS_503(
     &                WHAT, KEEP, NCB, NFR, SLAVEF,
     &                NBROWMAX, SIZECB )
                  SIZECB     = MIN(CBMAXR,SIZECB)
                  MAXTEMPCB  = MAX(MAXTEMPCB, MIN(SIZECB,CBMAXR))
                  MAXITEMPCB = MAX(MAXITEMPCB, NBROWMAX + NCB + 6 )
                  SBUF_REC = MAX(SBUF_REC, SIZECB+NBROWMAX + NCB + 6)
                  IF (MASTER) THEN
                     SIZECBI =  NCB + 5 +  SLAVEF - 1
                     SIZECB  = 0
                  ELSE IF (UPDATE) THEN
                     SIZECBI =  NFR + 6 + SLAVEF - 1 
                     SIZECB  = NBROWMAX * NCB
                     SIZECBI = SIZECBI + NBROWMAX + NFR + 6
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
                        WHAT = 1
                        CALL SMUMPS_503( WHAT, KEEP,
     *                  NCB, NFR, SLAVEF, NBROWMAX, IDUMMY )
                        SIZECBI =  NFR + 6 + SLAVEF - 1
                        SIZECB  = NBROWMAX * NCB
                        SIZECBI = SIZECBI + NBROWMAX + NFR + 6
                     ENDIF 
                  ENDIF
               ENDIF
            ENDIF
            IF (STACKCB) THEN
               IF (FRERE(STEP(INODE)).EQ.0) THEN 
                  write(*,*) ' ERROR 3 in mpima41md'
                  CALL SMUMPS_ABORT()
               ENDIF
               ITOP = ITOP + 1
               IF ( ITOP .GT. NSTEPS ) THEN
                 WRITE(*,*) 'ERROR 4 in mpima41md '
               ENDIF
               LSTKI(ITOP) = SIZECBI
               ISTKI=ISTKI+LSTKI(ITOP)
               ISTKIM = MAX(ISTKIM,ISTKI)
               LSTKR(ITOP) = SIZECB
               ISTKR = ISTKR + LSTKR(ITOP)
               NRLNEC = MAX(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
               NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
            ENDIF
         ENDIF ! end of (concerned) 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
      ENDIF
 115  CONTINUE
      ISTKIM = ISTKIM / 2 + MAXFR + N
      IF (KEEP(38) .ne. 0 ) THEN
         IF (SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *                                        .EQ. MYID) THEN
            ISTKIM = ISTKIM + ND( STEP( KEEP(38) ))
         END IF
      END IF
      IF ((NRLNEC - NRLADU) .LT. ISTKIM) THEN
         NRLNEC = NRLNEC + (ISTKIM-(NRLNEC - NRLADU))
      ENDIF
      IF ((NIRNEC - NIRADU) .LT. ISTKIM) THEN
         NIRNEC = NIRNEC + (ISTKIM-(NIRNEC - NIRADU))
      ENDIF
      IF ( KEEP(53) .NE. 0 ) THEN
         IF ( KEEP(38) .ne. 0 ) THEN
            IROOT = KEEP( 38 ) 
         ELSE
            IROOT = KEEP( 20 )
         END IF
         ROOT_OWNER =
     *   ( MYID .eq. SMUMPS_275(STEP(IROOT), PROCNODE, SLAVEF) )
         SIZE_ROOT = ND(STEP(IROOT))
         CALL SMUMPS_279( KEEP(51), KEEP(51), SIZE_ROOT,
     *        LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP, LIWK_QR, LWK_QR )
         IF ( NRLNEC-NRLADU .LT. LWK_QR ) THEN
            NRLNEC = NRLADU + LWK_QR
         END IF
         IF ( NIRNEC-NIRADU .LT. LIWK_QR ) THEN
            NIRNEC = NIRADU + LIWK_QR
         END IF
      END IF
      SBUF_REC = MAX(SBUF_REC, MAXTEMPCB+MAXITEMPCB )
      SBUF_REC = SBUF_REC + 10
      SBUF_SEND = SBUF_SEND + 10
      NRLBDU = NRLNEC + PERLU*(NRLNEC/100+1) 
      IF (SLAVEF.EQ.1) THEN 
         SBUF_REC = 1
         SBUF_SEND= 1
      ELSE
         NRLBDU = NRLBDU + SBUF_REC +  
     &        4*MAX(PERLU,0)*(SBUF_REC/100)
      ENDIF
      NIRBDU = NIRNEC + 2 * MAX(PERLU,10)*(NIRNEC/100+1)
      DEALLOCATE( LSTKR, TNSTK, IPOOL,
     *          LSTKI )
      RETURN
      END SUBROUTINE SMUMPS_246
      SUBROUTINE SMUMPS_271( COMM_LOAD, ASS_IRECV, 
     *    INODE, NELIM_ROOT, ROOT, 
     *
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
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      TYPE (SMUMPS_ROOT_STRUC) :: ROOT
      INTEGER KEEP(150), ICNTL( 40 )
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER INODE, NELIM_ROOT
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
      INTEGER NBPROCFILS(KEEP(28))
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER MYID, SLAVEF, NBFIN
      REAL OPASSW, OPELIW
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) ) 
      INTEGER INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mumps_tags.h'
      INTEGER I, J, OPSFAC, APOS, LCONT, NCOL_TO_SEND, LDA
      INTEGER FPERE, IOLDPS, NFRONT, NPIV, NASS, NSLAVES,
     *        H_INODE, NELIM, NBCOL, LIST_NELIM_ROW, 
     *        LIST_NELIM_COL, NELIM_LOCAL, TYPE_SON, 
     *        POSELT, NROW, NCOL, NBROW, SHIFT_LIST_ROW_SON,
     *        SHIFT_LIST_COL_SON, SHIFT_VAL_SON,LDAFS, IERR,
     *        STATUS( MPI_STATUS_SIZE ), ISON, PDEST_MASTER_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER MSGSOU, MSGTAG
      LOGICAL INVERT, FLAG
      INTEGER  SMUMPS_275, SMUMPS_330
      EXTERNAL SMUMPS_275, SMUMPS_330
      FPERE = KEEP(38)
      IF ( SMUMPS_275( STEP(INODE), PROCNODE_STEPS,
     *     SLAVEF ).EQ.MYID) THEN
       IOLDPS   = PTLUST_S(STEP(INODE))
       NFRONT   = IW(IOLDPS)
       NPIV     = IW(IOLDPS+1)
       NASS     = IABS(IW(IOLDPS + 2))
       NSLAVES  =  IW(IOLDPS+5)
       H_INODE  = 6 + NSLAVES
       NELIM    = NASS - NPIV
       NBCOL = NFRONT - NPIV
       LIST_NELIM_ROW = IOLDPS + H_INODE + NPIV
       LIST_NELIM_COL = LIST_NELIM_ROW + NFRONT
           IF (NELIM.LE.0) THEN
            write(6,*) ' ERROR 1 in SMUMPS_271 ', NELIM
            write(6,*) MYID,':Process root2son: INODE=',INODE,
     * 'Header=',IW(PTLUST_S(STEP(INODE)):PTLUST_S(STEP(INODE))+5)
            CALL SMUMPS_ABORT()
           ENDIF
       NELIM_LOCAL = NELIM_ROOT
       DO I=1, NELIM
        root%RG2L_ROW(IW(LIST_NELIM_ROW)) = NELIM_LOCAL
        root%RG2L_COL(IW(LIST_NELIM_COL)) = NELIM_LOCAL
        NELIM_LOCAL = NELIM_LOCAL + 1
        LIST_NELIM_ROW = LIST_NELIM_ROW + 1
        LIST_NELIM_COL = LIST_NELIM_COL + 1
       ENDDO
       NBROW = NFRONT - NPIV
       NROW = NELIM
       IF ( KEEP( 50 ) .eq. 0 ) THEN
         NCOL = NFRONT - NPIV
       ELSE
         NCOL = NELIM
       END IF
       SHIFT_LIST_ROW_SON = H_INODE + NPIV
       SHIFT_LIST_COL_SON = H_INODE + NFRONT + NPIV
       TYPE_SON = SMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
       IF ( KEEP(50).eq.0 .OR. TYPE_SON .eq. 1 ) THEN
         LDAFS = NFRONT
       ELSE
         LDAFS = NASS
       END IF
       SHIFT_VAL_SON      = NPIV * LDAFS + NPIV
       CALL SMUMPS_80( COMM_LOAD,
     *   ASS_IRECV, 
     *   N, INODE, FPERE,
     *   PTLUST_S(1), PTRAST(1),
     *   ROOT, NROW, NCOL, SHIFT_LIST_ROW_SON,
     *   SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDAFS,
     *   ROOT_NON_ELIM_CB, MYID, COMM,
     *   BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *   IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *   NIRBDU, PTRIST, PTLUST_S(1), PTRAST(1),
     *   STEP, PIMASTER, PAMASTER,
     *   NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *   IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *   OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *   INTARR, DBLARR, ICNTL, KEEP, .FALSE., ND, FRERE,
     *   LPTRAR, NELT, FRTPTR, FRTELT, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
       IF (IFLAG.LT.0 ) RETURN
       IF (TYPE_SON.EQ.1) THEN
        NROW = NFRONT - NASS
        NCOL = NELIM
        SHIFT_LIST_ROW_SON = H_INODE + NASS
        SHIFT_LIST_COL_SON = H_INODE + NFRONT + NPIV
        SHIFT_VAL_SON      = NASS * NFRONT + NPIV
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          INVERT = .FALSE.
        ELSE
          INVERT = .TRUE.
        END IF
        CALL SMUMPS_80( COMM_LOAD, ASS_IRECV,
     *    N, INODE, FPERE,
     *    PTLUST_S, PTRAST,
     *    ROOT, NROW, NCOL, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, NFRONT,
     *    ROOT_NON_ELIM_CB, MYID, COMM,
     *
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, INVERT, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF (IFLAG.LT.0 ) RETURN
       ENDIF
       IOLDPS = PTLUST_S(STEP(INODE))
       POSELT = PTRAST(STEP(INODE))
       IW(IOLDPS + 4) = POSELT
       IF ( TYPE_SON .eq. 1 ) THEN
         NBROW = NFRONT - NPIV
       ELSE
         NBROW = NELIM
       END IF
       IF (KEEP(50).EQ.0) THEN
        IF ((NPIV .NE. 0) .AND. (NBROW .GT. 1)) THEN
           IF ( IOLDPS .NE. PTLUST_S(STEP( INODE ) ) ) THEN
             WRITE(*,*) ' ERROR 5 in process_root2son : ',
     *       ' IOLDPS,PTLUST_S(STEP(INODE))=',
     *       IOLDPS,PTLUST_S(STEP(INODE))
             CALL SMUMPS_ABORT()
           END IF
         OPSFAC = POSELT + NPIV * NFRONT + NPIV
         APOS   = POSELT + (NPIV + 1) * NFRONT
         DO I = 1, NBROW - 1
           DO J = 1, NPIV
             A(OPSFAC) = A(APOS)
             OPSFAC = OPSFAC + 1
             APOS = APOS + 1
           END DO
           APOS = APOS + NBCOL
         END DO
        ENDIF
       ELSE
         IF ( TYPE_SON .eq. 1 ) THEN
           CALL SMUMPS_324(A(POSELT), NFRONT,
     *          NPIV, NBROW,KEEP(59))
         ELSE
           CALL SMUMPS_325(A(POSELT), NPIV, NBROW)
         END IF
       ENDIF
       IW(IOLDPS)     = NBCOL
       IW(IOLDPS + 1) = NASS - NPIV
       IF (TYPE_SON.EQ.2) THEN
        IW(IOLDPS + 2) = NASS
       ELSE
        IW(IOLDPS + 2) = NFRONT
       ENDIF
       IW(IOLDPS + 3) = NPIV
      CALL SMUMPS_93(MYID,N,IOLDPS,TYPE_SON,IW,LIW,
     *    A, LA, POSFAC, LRLU, LRLUS,
     *    IWPOS, PTRAST,STEP,  KEEP)
      ELSE 
        ISON = INODE
        PDEST_MASTER_ISON = SMUMPS_275(STEP(ISON),
     *      PROCNODE_STEPS,SLAVEF)
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
          BLOCKING  = .TRUE.
          SET_IRECV = .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, MSGTAG,
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
     *    INTARR, DBLARR, ICNTL, KEEP, ND, FRERE, LPTRAR,
     *    NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) RETURN
        END DO
       IOLDPS   = PTRIST(STEP(INODE))
       LCONT  = IW(IOLDPS)
       NROW   = IW(IOLDPS+2)
       NPIV   = IW(IOLDPS+3)
       NASS   = IW(IOLDPS+4)
       NELIM  = NASS-NPIV
           IF (NELIM.LE.0) THEN
            write(6,*) MYID,': ERROR 2 in SMUMPS_271 '
            write(6,*) MYID,': INODE,LCONT, NROW, NPIV, NASS, NELIM=',
     *      INODE,LCONT, NROW, NPIV, NASS, NELIM
            write(6,*) MYID,': IOLDPS=',IOLDPS
            CALL SMUMPS_ABORT()
           ENDIF
       NSLAVES= IW(IOLDPS+5)
       H_INODE = 6 + NSLAVES
       LIST_NELIM_COL = IOLDPS + H_INODE + NROW + NPIV
       NELIM_LOCAL = NELIM_ROOT
       DO I = 1, NELIM
        root%RG2L_COL(IW(LIST_NELIM_COL)) = NELIM_LOCAL
        NELIM_LOCAL = NELIM_LOCAL + 1
        LIST_NELIM_COL = LIST_NELIM_COL + 1
       ENDDO
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5)
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NPIV
       SHIFT_VAL_SON      = NPIV
       NCOL_TO_SEND       = NELIM
       LDA                = LCONT + NPIV
       IF ( KEEP( 50 ) .eq. 0 ) THEN
         INVERT = .FALSE.
       ELSE
         INVERT = .TRUE.
       END IF
       CALL SMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE,
     *    PTRIST, PTRAST,
     *    ROOT, NROW, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA,
     *    ROOT_NON_ELIM_CB, MYID, COMM,
     *
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP, INVERT, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF (IFLAG.LT.0 ) RETURN
        CALL SMUMPS_314( N, INODE,
     *    PTRIST, PTRAST, PTLUST_S, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP)
        IF (IFLAG.LT.0) THEN
           CALL SMUMPS_44( MYID, SLAVEF, COMM )
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_271
       SUBROUTINE SMUMPS_217(N, NZ, NSCA, 
     *      ASPK, IRN, ICN, COLSCA, ROWSCA, S, MAXS, 
     *      ICNTL, INFO)
      INTEGER N, NZ, NSCA, MAXS
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER ICNTL(40), INFO(40)
      REAL    ASPK(NZ), COLSCA(*), ROWSCA(*)
      REAL    S(MAXS)
      INTEGER MPG,LP
      INTEGER ISPW1, IWNOR
      INTEGER I, K, ITOT
      LOGICAL PROK
      REAL ONE
      PARAMETER( ONE = 1.0E0 )
      LP      = ICNTL(1)
      MPG     = ICNTL(2)
      MPG    = ICNTL(3)
      PROK   = (MPG.GT.0)
      IF (PROK) WRITE(MPG,101)
 101    FORMAT(/' ****** SCALING OF ORIGINAL MATRIX '/)
        IF (NSCA.EQ.1) THEN
         IF (PROK)
     *    WRITE (MPG,*) ' DIAGONAL SCALING '
        ELSEIF (NSCA.EQ.2) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' SCALING BASED ON (MC29)'
        ELSEIF (NSCA.EQ.3) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' COLUMN SCALING'
        ELSEIF (NSCA.EQ.4) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' ROW AND COLUMN SCALING (1 Pass)'
        ELSEIF (NSCA.EQ.5) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' MC29 FOLLOWED BY ROW &COL SCALING'
        ELSEIF (NSCA.EQ.6) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' MC29 FOLLOWED BY COLUMN SCALING'
        ENDIF
        DO 10 I=1,N
            COLSCA(I) = ONE
            ROWSCA(I) = ONE
 10     CONTINUE
        IF ((NSCA.EQ.5).OR.
     &      (NSCA.EQ.6))                   THEN
          ITOT = 5*N + NZ 
          IF (ITOT.GT.MAXS) GOTO 400
          ISPW1 = MAXS - NZ + 1
          DO 15 K=1,NZ
           S(ISPW1+K-1) = ASPK(K)
  15      CONTINUE
        ELSE
          ISPW1 = MAXS + 1
          ITOT  = 5*N
          IF (ITOT.GT.MAXS) GOTO 400
        ENDIF
        IWNOR = ISPW1 - 5*N
          IF (NSCA.EQ.1) THEN
            CALL SMUMPS_238(N,NZ,ASPK,IRN,ICN,S(IWNOR),
     *        COLSCA,ROWSCA,MPG)
          ELSEIF (NSCA.EQ.2) THEN
            CALL SMUMPS_239(N,NZ,ASPK,IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPG,MPG,NSCA)
          ELSEIF (NSCA.EQ.3) THEN
            CALL SMUMPS_241(N,NZ,ASPK,IRN,ICN,S(IWNOR),COLSCA,
     *      MPG)
          ELSEIF (NSCA.EQ.4) THEN
            CALL SMUMPS_287(N,NZ,IRN,ICN,ASPK,
     *      S(IWNOR),S(IWNOR+N),COLSCA,ROWSCA,MPG)
          ELSEIF (NSCA.EQ.5) THEN
            CALL SMUMPS_239(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPG,MPG,NSCA)
            CALL SMUMPS_241(N,NZ,S(ISPW1),IRN,ICN,S(IWNOR),
     *          COLSCA, MPG)
          ELSEIF (NSCA.EQ.6) THEN
            CALL SMUMPS_239(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPG,MPG,NSCA)
            CALL SMUMPS_240(NSCA,N,NZ,IRN,ICN,S(ISPW1),
     *            S(IWNOR+N),ROWSCA,MPG)
            CALL SMUMPS_241(N,NZ,S(ISPW1),IRN,ICN,S(IWNOR),
     *          COLSCA, MPG)
          ENDIF
      GOTO 500
 400  INFO(1) = -5
      INFO(2) = ITOT-MAXS
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1))
     * WRITE(LP,*) '*** ERROR: Not enough space to scale matrix'
 500  RETURN
      END SUBROUTINE SMUMPS_217
      SUBROUTINE SMUMPS_287(N,NZ,IRN,ICN,VAL,
     *    RNOR,CNOR,COLSCA,ROWSCA,MPRINT)
      INTEGER N, NZ
      REAL    VAL(NZ),RNOR(N),CNOR(N)
      REAL    COLSCA(N),ROWSCA(N)
      REAL    CMIN,CMAX,RMIN,ARNOR,ACNOR
      INTEGER IRN(NZ), ICN(NZ)
      REAL    VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO, ONE
      PARAMETER(ZERO=0.0E0, ONE=1.0E0)
      DO 50 J=1,N
       CNOR(J)   = REAL(ZERO)
       RNOR(J)   = REAL(ZERO)
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     *        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.ABS(CNOR(J))) THEN
              CNOR(J) =     VDIAG
            ENDIF
            IF (VDIAG.GT.ABS(RNOR(I))) THEN
              RNOR(I) =     VDIAG
            ENDIF
 100   CONTINUE
      IF (MPRINT.GT.0) THEN
       CMIN = ABS( CNOR(1) )
       CMAX = ABS( CNOR(1) )
       RMIN = ABS( RNOR(1) )
       DO 111 I=1,N
        ARNOR = ABS(RNOR(I))
        ACNOR = ABS(CNOR(I))
        IF (ACNOR.GT.CMAX) CMAX=ACNOR
        IF (ACNOR.LT.CMIN) CMIN=ACNOR
        IF (ARNOR.LT.RMIN) RMIN=ARNOR
 111   CONTINUE
       WRITE(MPRINT,*) '**** STAT. OF MATRIX PRIOR ROW&COL SCALING'
       WRITE(MPRINT,*) ' MAXIMUM NORM-MAX OF COLUMNS:',CMAX
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF COLUMNS:',CMIN
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF ROWS   :',RMIN
      ENDIF
      DO 120 J=1,N
       IF (ABS(CNOR(J)).LE.ZERO) THEN
         CNOR(J)   = REAL(ONE)
       ELSE
         CNOR(J)   = REAL(ONE)/CNOR(J)
       ENDIF
 120  CONTINUE
      DO 130 J=1,N
       IF (ABS(RNOR(J)).LE.ZERO) THEN
         RNOR(J)   = REAL(ONE)
       ELSE
         RNOR(J)   = REAL(ONE)/RNOR(J)
       ENDIF
 130  CONTINUE
       DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I)* RNOR(I)
        COLSCA(I) = COLSCA(I) * CNOR(I)
 110   CONTINUE
      IF (MPRINT.GT.0)
     *  WRITE(MPRINT,*) ' END OF SCALING BY MAX IN ROW AND COL'
      RETURN
      END SUBROUTINE SMUMPS_287
      SUBROUTINE SMUMPS_239(N,NZ,VAL,ROWIND,COLIND,
     *       RNOR,CNOR,WNOR,MPRINT,MP,
     *       NSCA)
      INTEGER N, NZ
      REAL    VAL(NZ),RNOR(N),CNOR(N),
     *        WNOR(5*N)
      INTEGER COLIND(NZ),ROWIND(NZ)
      INTEGER J,I,K
      INTEGER MPRINT,MP,NSCA
      INTEGER IFAIL9
      REAL ZERO, ONE
      PARAMETER( ZERO = 0.0E0, ONE = 1.0E0 )
      DO 15 I=1,N
       RNOR(I)   = REAL(ZERO)
       CNOR(I)   = REAL(ZERO)
  15  CONTINUE
      CALL SMUMPS_216(N,N,NZ,VAL,ROWIND,COLIND,
     *   RNOR,CNOR,WNOR, MP,IFAIL9)
      DO 30 I=1,N
       CNOR(I) = EXP(CNOR(I))
       RNOR(I) = EXP(RNOR(I))
  30  CONTINUE
      IF ((NSCA.EQ.5).OR.(NSCA.EQ.6)) THEN
        DO 100 K=1,NZ
          I   = ROWIND(K)
          J   = COLIND(K)
          IF (MIN(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 100
          VAL(K) = VAL(K) * CNOR(J) * RNOR(I)
 100    CONTINUE
      ENDIF
      IF (MPRINT.GT.0) 
     *   WRITE(MPRINT,*) ' END OF SCALING USING MC29'
      RETURN
      END SUBROUTINE SMUMPS_239
      SUBROUTINE SMUMPS_241(N,NZ,VAL,IRN,ICN,
     *       CNOR,COLSCA,MPRINT)
      INTEGER N,NZ
      REAL VAL(NZ),CNOR(N),COLSCA(N)
      INTEGER IRN(NZ), ICN(NZ)
      REAL VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO, ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      DO 10 J=1,N
       CNOR(J)   = REAL(ZERO)
  10  CONTINUE
      DO 100 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.LE.0).OR.(I.GT.N).OR.
     *      (J.LE.0).OR.(J.GT.N)) GOTO 100
        VDIAG = ABS(VAL(K))
        IF (VDIAG.GT.ABS(CNOR(J))) THEN
           CNOR(J) =     VDIAG
        ENDIF
 100  CONTINUE
      DO 110 J=1,N
       IF (ABS(CNOR(J)).LE.ZERO) THEN
         CNOR(J)   = REAL(ONE)
       ELSE
         CNOR(J)   = REAL(ONE)/CNOR(J)
       ENDIF
 110  CONTINUE
       DO 215 I=1,N
        COLSCA(I) = COLSCA(I) * CNOR(I)
 215   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF COLUMN SCALING'
      RETURN
      END SUBROUTINE SMUMPS_241
      SUBROUTINE SMUMPS_238(N,NZ,VAL,IRN,ICN,RNOR,
     *      COLSCA,ROWSCA,MPRINT)
      INTEGER   N, NZ
      REAL RNOR(N)
      REAL  VAL(NZ), ROWSCA(N),COLSCA(N)
      INTEGER   IRN(NZ),ICN(NZ)
      REAL      VDIAG
      INTEGER   MPRINT,I,J,K
      INTRINSIC SQRT
      REAL ZERO, ONE
      PARAMETER(ZERO=0.0E0, ONE=1.0E0)
      DO 10 I=1,N
       RNOR(I)   = REAL(ONE)
  10  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          IF ((I.GT.N).OR.(I.LE.0)) GOTO 100
          J = ICN(K)
          IF (I.EQ.J) THEN
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.ZERO) THEN
              RNOR(J) = REAL(ONE/(SQRT(VDIAG)))
            ENDIF
          ENDIF
 100   CONTINUE
       DO 110 I=1,N
        COLSCA(I) = RNOR(I)
        ROWSCA(I) = RNOR(I)
 110   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF DIAGONAL SCALING'
      RETURN
      END SUBROUTINE SMUMPS_238
      SUBROUTINE SMUMPS_240(NSCA,N,NZ,IRN,ICN,VAL,
     *    RNOR,ROWSCA,MPRINT)
      INTEGER N, NZ, NSCA
      INTEGER IRN(NZ), ICN(NZ)
      REAL VAL(NZ),RNOR(N), ROWSCA(N)
      REAL VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO,ONE
      PARAMETER (ZERO=0.0E0, ONE=1.0E0)
      DO 50 J=1,N
       RNOR(J)   = REAL(ZERO)
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     *        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.ABS(RNOR(I))) THEN
              RNOR(I) =  VDIAG
            ENDIF
 100   CONTINUE
      DO 130 J=1,N
       IF (REAL(RNOR(J)).LE.REAL(ZERO)) THEN
         RNOR(J)   = REAL(ONE)
       ELSE
         RNOR(J)   = REAL(ONE)/RNOR(J)
       ENDIF
 130  CONTINUE
      DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I)* RNOR(I)
 110  CONTINUE
      IF ( (NSCA.EQ.4) .OR. (NSCA.EQ.6) ) THEN
        DO 150 K=1,NZ
          I   = IRN(K)
          J   = ICN(K)
          IF (MIN(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 150
          VAL(K) = VAL(K) * RNOR(I)
 150    CONTINUE
      ENDIF
      IF (MPRINT.GT.0)
     *  WRITE(MPRINT,'(A)') '  END OF ROW SCALING'
      RETURN
      END
      SUBROUTINE SMUMPS_216(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      REAL A(NE)
      INTEGER IRN(NE),ICN(NE)
      REAL R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL
      INTRINSIC LOG,ABS,MIN
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      REAL ONE
      REAL SMIN,ZERO
      PARAMETER (ONE=1.0E0,SMIN=0.1E0,ZERO=0.0E0)
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      REAL E,E1,EM,Q,Q1,QM,S,S1,SM,U,V
      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF
      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3
      DO 10 I = 1,M
         R(I) = REAL(ZERO)
         W(I1+I) = REAL(ZERO)
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = REAL(ZERO)
         W(I2+J) = REAL(ZERO)
         W(I3+J) = REAL(ZERO)
         W(I4+J) = REAL(ZERO)
   20 CONTINUE
      DO 30 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*NE
      DO 60 K = 1,NE
         IF (ABS(A(K)).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
      E = REAL(ZERO)
      Q = ONE
      S = REAL(ZERO)
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (ABS(S).LE.ABS(SM)) GO TO 160
      DO 150 ITER = 1,MAXIT
         DO 80 K = 1,NE
            IF (ABS(A(K)).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = REAL(ZERO)
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
         IF (ABS(S).LE.ABS(SM)) E = REAL(ZERO)
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (ABS(S).LE.ABS(SM)) GO TO 180
         EM = E*E1
         DO 110 K = 1,NE
            IF (ABS(A(K)).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = REAL(ZERO)
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
         IF (ABS(S).LE.ABS(SM)) Q = ONE
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
         IF (ABS(S).LE.ABS(SM)) GO TO 160
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
  180 DO 190 K = 1,NE
         IF (ABS(A(K)).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN
  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC29AS ****',' IFAIL =',IFAIL
      END SUBROUTINE SMUMPS_216
      SUBROUTINE SMUMPS_89(N,KEEP28, IW,LIW,A,LA,NIRBDU,
     *       LRLU,IPTRLU,IWPOS,
     *       IWPOSCB,PTRIST,PTRAST,STEP, PIMASTER,PAMASTER)
      IMPLICIT NONE
      INTEGER N,LIW,LA,NIRBDU,LRLU,KEEP28,
     &        IPTRLU,IWPOS,IWPOSCB
      INTEGER IW(LIW),PTRIST(KEEP28),PTRAST(KEEP28),
     &        STEP(N), 
     * PIMASTER(KEEP28),
     * PAMASTER(KEEP28)
      REAL A(LA)
      LOGICAL FREE, SAME_PROC
      INTEGER IPTIW,IPTA,NCOL,NROW,NSLSON,HF,
     &        SIZFI,SIZFR,LONGI,LONGR
      INTEGER I
      IPTIW = IWPOSCB
      IPTA  = IPTRLU
      LONGI = 0
      LONGR = 0
10    FREE = .FALSE.
      SAME_PROC = IW(IPTIW+1) .LT.0
      IF (SAME_PROC.AND.(IW(IPTIW+1) .LT. -N)) THEN
       FREE = .TRUE.
       NROW = -(IW(IPTIW+1) + N)
       NCOL = IW(IPTIW+2)
       SIZFR   = NROW*NCOL
       SIZFI   = 2
      ENDIF
      IF (.NOT.SAME_PROC) THEN
      IF (IW(IPTIW+5).EQ.0) THEN
        FREE = .TRUE.
        NCOL    = IW(IPTIW+1) + IW(IPTIW + 4)
        NROW    = IW(IPTIW+3)
        SIZFR   = NROW * NCOL
        NSLSON  = IW(IPTIW + 6)
        HF      = NSLSON + 6
        SIZFI   = NROW + NCOL + HF
        IF (IW(IPTIW + 4).LT.0) THEN 
         SIZFI   = NROW + IW(IPTIW+1) + HF
        ENDIF
      END IF
      ENDIF
      IF (FREE) THEN
        IF (LONGI.NE.0) THEN
         DO I=0,LONGI-1
           IW(IPTIW + SIZFI - I) = IW (IPTIW - I )
         ENDDO
         DO I=0,LONGR-1
           A(IPTA + SIZFR - I)   = A(IPTA - I )
         ENDDO
          DO I=1,KEEP28
           IF  ((PTRIST(I).LE.IPTIW).AND.
     &       (PTRIST(I).GT.IWPOSCB) ) THEN
            PTRIST(I) = PTRIST(I) + SIZFI
            IF (IW(PTRIST(I)).LT.0) THEN
             PAMASTER(I) = PAMASTER(I) + SIZFR
            ELSE
             PTRAST(I) = PTRAST(I) + SIZFR
            ENDIF
           ENDIF 
           IF ((PIMASTER(I).LE.IPTIW).AND.
     &       (PIMASTER(I).GT.IWPOSCB) ) THEN
            PIMASTER(I) = PIMASTER(I) + SIZFI
            PAMASTER(I) = PAMASTER(I) + SIZFR
           ENDIF
          ENDDO 
        ENDIF
        IWPOSCB = IWPOSCB + SIZFI
        IPTRLU = IPTRLU + SIZFR
        LRLU   = LRLU + SIZFR
      ELSE
        IF (SAME_PROC) THEN
         NROW  = -IW(IPTIW+1)
         NCOL  = IW(IPTIW+2)
         SIZFR = NROW*NCOL
         SIZFI = 2
        ELSE
         NCOL    = IW(IPTIW+1) + IW(IPTIW + 4)
         NROW    = IW(IPTIW+3)
         SIZFR   = NROW * NCOL
         NSLSON  = IW(IPTIW + 6)
         HF      = NSLSON + 6
         SIZFI   = NROW + NCOL + HF
         IF (IW(IPTIW + 4).LT.0) THEN 
          SIZFI   = NROW + IW(IPTIW+1) + HF
         ENDIF
        ENDIF
        LONGI = LONGI + SIZFI
        LONGR = LONGR + SIZFR
      ENDIF
      IPTA  = IPTA  + SIZFR
      IPTIW = IPTIW + SIZFI
      IF (IPTIW.NE.NIRBDU) GOTO 10
       RETURN
       END SUBROUTINE SMUMPS_89
      SUBROUTINE SMUMPS_221(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU, SEUIL)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV,NOFFW
      REAL A(LA) 
      REAL UU, SEUIL
      INTEGER IW(LIW) 
      INTEGER  IOLDPS, POSELT
       include 'mpif.h'
      REAL SWOP
      INTEGER APOS
      REAL AMROW
      REAL ZERO,RMAX
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER ISAMAX
      INTRINSIC MAX
      DATA ZERO /0.0E0/
        INOPV   = 0
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
        NASSW   = IABS(IW(IOLDPS+3))
          DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + NFRONT*(IPIV-1) + NPIV
            JMAX = 1
            IF (UU.GT.ZERO) GO TO 340
            IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
            GO TO 380
  340       AMROW = ZERO
            J1 = APOS
            J2 = APOS - NPIV + NASS - 1
             J3    = NASS -NPIV
             JMAX  = ISAMAX(J3,A(J1),1)
             JJ    = JMAX + J1 - 1
             AMROW = ABS(A(JJ))
            RMAX = AMROW
            J1 = J2 + 1
            J2 = APOS - NPIV + NFRONT - 1
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
  360       CONTINUE
  370       IDIAG = APOS + IPIV - NPIVP1
            IF (ABS(A(IDIAG)).GT. MAX(UU*RMAX,SEUIL)) THEN
              JMAX = IPIV - NPIV
              GO TO 380
            ENDIF
            IF (AMROW.LE. MAX(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV*NFRONT
            J2 = J1 + NFRONT - 1
            J3 = POSELT + (IPIV-1)*NFRONT
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1
  390       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIVP1
            ISWPS2 = IOLDPS + 5 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            J1 = POSELT + NPIV
            J2 = POSELT + NPIV + JMAX - 1
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT
              J2 = J2 + NFRONT
  410       CONTINUE
            ISWPS1 = IOLDPS + 5 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + 5 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
      WRITE(*,*) 'Detected a null pivot, INODE/NPIV=',INODE,NPIV
  420 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_221
      SUBROUTINE SMUMPS_220(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *   INOPV,NOFFW,IOLDPS,POSELT,UU)
      INTEGER NFRONT,NASS,N,LIW,LA,INODE,INOPV
      REAL UU
      REAL A(LA)
      INTEGER IW(LIW)
      REAL AMROW
      REAL ZERO,RMAX
      REAL  SWOP
      INTEGER APOS, POSELT, IOLDPS
      INTEGER NOFFW,NPIV,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER ISAMAX
      INTRINSIC MAX
      DATA ZERO /0.0E0/
        INOPV   = 0
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
          DO 460 IPIV=NPIVP1,NASS
            APOS = POSELT + NFRONT*NPIV + (IPIV-1)
            JMAX = 1
            AMROW = ZERO
            J1 = APOS
            J3    = NASS -NPIV
            JMAX  = ISAMAX(J3,A(J1),NFRONT)
            JJ    = J1 + (JMAX-1)*NFRONT
            AMROW = ABS(A(JJ))
            RMAX = AMROW
            J1 = APOS +  (NASS-NPIV) * NFRONT
            J3 = NFRONT - NASS
            IF (J3.EQ.0) GOTO 370
            DO 360 JJ=1,J3
              RMAX = MAX(ABS(A(J1)),RMAX)
              J1 = J1 + NFRONT
  360       CONTINUE
  370       IF (RMAX.EQ.ZERO) GO TO 460
            IDIAG = APOS + (IPIV - NPIVP1)*NFRONT
            IF (ABS(A(IDIAG)).GE.UU*RMAX) JMAX = IPIV - NPIV
            IF (ABS(A(IDIAG)).GE.UU*RMAX) GO TO 380
            IF (AMROW.LT.UU*RMAX) GO TO 460
            NOFFW = NOFFW + 1
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV
            J3 = POSELT + (IPIV-1)
            DO 390 JJ= 1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J3)
              A(J3) = SWOP
              J1 = J1 + NFRONT
              J3 = J3 + NFRONT
  390       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIVP1 + NFRONT
            ISWPS2 = IOLDPS + 5 + IPIV + NFRONT
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            J1 = POSELT + NPIV*NFRONT
            J2 = POSELT + (NPIV + JMAX - 1)*NFRONT
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + 1
              J2 = J2 + 1
  410       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIV + 1
            ISWPS2 = IOLDPS + 5 + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
       INOPV = 1
  420 RETURN
      END SUBROUTINE SMUMPS_220
      SUBROUTINE SMUMPS_225(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *     IOLDPS,POSELT,IFINB,LKJIB,LKJIT)
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFINB,LKJIB
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    VALPIV
      INTEGER APOS, UUPOS, IOLDPS, POSELT
      INTEGER LKJIT
      REAL ONE, ALPHA
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,KROW,LPOS,NEL
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NPIV   = IW(IOLDPS+1)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
        IF (IW(IOLDPS+3).LE.0) THEN
          IF (NASS.LT.LKJIT) THEN
           IW(IOLDPS+3) = NASS
          ELSE
           IW(IOLDPS+3) = MIN0(NASS,LKJIB)
          ENDIF
        ENDIF
        JROW2 = IW(IOLDPS+3)
        NEL2   = JROW2 - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
          IW(IOLDPS+3) = MIN0(JROW2+LKJIB,NASS)
          IW(IOLDPS+4) = NPIVP1+1
         ENDIF
        ELSE
         APOS   = POSELT + NPIV*(NFRONT + 1)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541     CONTINUE
         LPOS   = APOS + NFRONT
         UUPOS  = APOS+1
         CALL SGER(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     *              A(LPOS+1),NFRONT)
        ENDIF
        RETURN
        END SUBROUTINE SMUMPS_225
      SUBROUTINE SMUMPS_229(NFRONT,N,INODE,IW,LIW,A,LA,IOLDPS,
     *          POSELT)
      INTEGER NFRONT,N,INODE,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA,VALPIV
      INTEGER APOS, POSELT, UUPOS
      INTEGER IOLDPS,NPIV,NEL
      INTEGER LPOS,JROW,IRWPOS
      REAL ONE
      DATA ONE /1.0E0/
        NPIV   = IW(IOLDPS+1)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + (NPIV)*NFRONT + NPIV
        IF (NEL.EQ.0) GO TO 650
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        DO 340 JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
  340   CONTINUE
        LPOS   = APOS + NFRONT
        UUPOS  = APOS+1
        DO 440 JROW = 1,NEL
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL SAXPY(NEL,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
  650   RETURN
        END SUBROUTINE SMUMPS_229
      SUBROUTINE SMUMPS_228(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *       IOLDPS,POSELT,IFINB)
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFINB
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA,VALPIV
      INTEGER APOS, POSELT,UUPOS
      INTEGER IOLDPS,NPIV,KROW
      INTEGER NEL,LPOS,ICOL,NEL2,IRWPOS
      INTEGER NPIVP1
      REAL ONE
      DATA ONE /1.0E0/
        NPIV   = IW(IOLDPS+1)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        NEL2   = NASS - NPIVP1
        IFINB  = 0
        IF (NPIVP1.EQ.NASS) IFINB = 1
        APOS   = POSELT + NPIV*(NFRONT + 1)
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        DO 541 KROW = 1,NEL
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541    CONTINUE
        LPOS   = APOS + NFRONT
        UUPOS  = APOS+1
        DO 440 ICOL = 1,NEL
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL SAXPY(NEL2,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
        RETURN
        END SUBROUTINE SMUMPS_228
      SUBROUTINE SMUMPS_231(A,LA,NFRONT,
     *       NPIV,NASS,POSELT)
      INTEGER LA,POSELT
      REAL    A(LA)
      INTEGER NFRONT, NPIV, NASS
      INTEGER NEL1,NEL11,LPOS2,LPOS1,LPOS
      REAL ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT + NASS*NFRONT
        CALL STRSM('L','L','N','N',NPIV,NEL1,ONE,A(POSELT),NFRONT,
     *              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + NPIV
        LPOS1  = POSELT + NPIV
        CALL SGEMM('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        END SUBROUTINE SMUMPS_231
      SUBROUTINE SMUMPS_232(A,LA,NFRONT,NPIV,NASS,POSELT,LKJIB)
      INTEGER LA, NFRONT, NPIV, NASS, LKJIB
      REAL    A(LA)
      INTEGER POSELT
      INTEGER POSELT_LOCAL
      INTEGER NEL1, NEL11, NPBEG, LPOS, LPOS1, LPOS2
      REAL ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        POSELT_LOCAL = POSELT
        NEL1   = NASS - NPIV
        NPBEG  = NPIV - LKJIB + 1
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT_LOCAL + NPIV*NFRONT + NPBEG - 1
        POSELT_LOCAL = POSELT_LOCAL + (NPBEG-1)*NFRONT + NPBEG - 1
        CALL STRSM('L','L','N','N',LKJIB,NEL1,ONE,A(POSELT_LOCAL),
     *               NFRONT,A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIB
        LPOS1  = POSELT_LOCAL + LKJIB
        CALL SGEMM('N','N',NEL11,NEL1,LKJIB,ALPHA,A(LPOS1),
     *       NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        END SUBROUTINE SMUMPS_232
      SUBROUTINE SMUMPS_233(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *    IOLDPS,POSELT,LKJIB,LKJIT )
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE
      INTEGER POSELT
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1, NEL11
      INTEGER LBP, IPOS, KPOS, LPOS2, HF
      INTEGER LPOS1,LPOS,LBPT,I1,K1,II,ISWOP,LBP1
      INTEGER LKJIT, POSLOCAL
      REAL ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NPIV   = IW(IOLDPS+1)
        JROW2  = IABS(IW(IOLDPS+3))
        NPBEG  = IW(IOLDPS+4)
        HF     = 6 + IW(IOLDPS+5)
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
         LKJIB       = LKJIB + NONEL
         IW(IOLDPS+3)= MIN0(NPIV+LKJIB,NASS)
        ELSE
          IW(IOLDPS+3) = NASS
        ENDIF
        IW(IOLDPS+4) = NPIV + 1
        NEL1   = NASS - JROW2
        LKJIW  = NPIV - NPBEG + 1
        NEL11  = NFRONT - NPIV
        IF ((NEL1.EQ.0).OR.(LKJIW.EQ.0)) GO TO 500
        LPOS2  = POSELT + JROW2*NFRONT + NPBEG - 1
         POSLOCAL = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
         CALL STRSM('L','L','N','N',LKJIW,NEL1,ONE,
     *               A(POSLOCAL),NFRONT,
     *               A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIW
        LPOS1  = POSLOCAL + LKJIW
        CALL SGEMM('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   LBP = JROW2 - NPIV
         LBPT  = LBP + LBP
        IF  ((NEL1.GE.LBPT).AND.(NEL1.GE.LKJIB)) THEN
         I1 = IOLDPS + HF + NPIV
         K1 = IOLDPS + HF + NASS - LBP
         DO 10 II=1,LBP
          ISWOP  = IW(I1)
          IW(I1) = IW(K1)
          IW(K1) = ISWOP
          I1     = I1 +1
          K1     = K1 + 1
  10     CONTINUE
         IPOS = POSELT + NPIV*NFRONT
         KPOS = POSELT + (NASS-LBP)*NFRONT
         LBP1 = LBP * NFRONT
         CALL SSWAP(LBP1,A(IPOS),1,A(KPOS),1)
        ENDIF
        END SUBROUTINE SMUMPS_233
      SUBROUTINE SMUMPS_236(A,LA,NPIVB,NFRONT,
     *                             NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER NPIVB,NASS,LA
      REAL    A(LA)
      INTEGER APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER LPOS, LPOS1, LPOS2, NEL1, NEL11, NPIVE
      REAL    ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + NPIVB*NFRONT + NPIVB
        LPOS2  = APOS + NASSL
        CALL STRSM('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     *              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + NFRONT*NPIVE
        LPOS1  = APOS + NFRONT*NPIVE
        CALL SGEMM('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     *          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
        END SUBROUTINE SMUMPS_236
      SUBROUTINE SMUMPS_94(N,KEEP28,IW,LIW,A,LA,NIRBDU,
     *       LRLU,IPTRLU,IWPOS,
     *       IWPOSCB,PTRIST,PTRAST,STEP,PIMASTER,PAMASTER, ITLOC)
      IMPLICIT NONE
      INTEGER N,LIW,LA,NIRBDU,LRLU,KEEP28,
     &        IPTRLU,IWPOS,IWPOSCB
      INTEGER IW(LIW),PTRIST(KEEP28),PTRAST(KEEP28),
     &        STEP(N), 
     * PIMASTER(KEEP28),
     * PAMASTER(KEEP28), ITLOC(N)
      REAL A(LA)
      INTEGER IPTIW,IPTA, SIZFI,SIZFR, NBLOCKS, SIZE_PTR, I
      LOGICAL FREE, SAME_PROC, REDO
      SIZE_PTR = 0
      DO I = 1, KEEP28
        IF ( PTRIST( I ) .GT. IWPOSCB .and.
     *         PTRIST( I ) .lt. NIRBDU )  THEN
          SIZE_PTR = SIZE_PTR + 1
        END IF
        IF ( PIMASTER( I ) .GT. IWPOSCB .and.
     *         PIMASTER( I ) .lt. NIRBDU ) THEN
          SIZE_PTR = SIZE_PTR + 1
        END IF
      END DO
      NBLOCKS = ( N - 2 * SIZE_PTR ) / 4
      IF ( NBLOCKS .LT. 2 ) THEN
        WRITE(*,*) ' WARNING : INEFFICIENT COMPRESS CALLED'
        CALL SMUMPS_89( N, KEEP28, IW, LIW, A, LA, NIRBDU,
     *       LRLU, IPTRLU, IWPOS, IWPOSCB, PTRIST, PTRAST,
     *       STEP, PIMASTER, PAMASTER )
      ELSE
        REDO = .TRUE.
        DO WHILE ( REDO )
          CALL SMUMPS_91( N, KEEP28, IW, LIW, A, LA,
     *       NIRBDU, LRLU,
     *       IPTRLU, IWPOS, IWPOSCB, PTRIST, PTRAST, STEP,
     *       PIMASTER, PAMASTER, NBLOCKS, SIZE_PTR, ITLOC( 1 ),
     *       ITLOC( 1 + 2 * SIZE_PTR ),
     *       ITLOC( 1 + 2 *SIZE_PTR + NBLOCKS ),
     *       ITLOC( 1 + 2 *SIZE_PTR + 2 * NBLOCKS ),
     *       ITLOC( 1 + 2 *SIZE_PTR + 3 * NBLOCKS ), REDO )
        END DO
        ITLOC( 1: 2 * SIZE_PTR + 4 * NBLOCKS ) = 0
      END IF
      RETURN
      END SUBROUTINE SMUMPS_94
      SUBROUTINE SMUMPS_90
     *( N, IPTIW, LIW, IW, SIZFI, SIZFR, FREE )
      IMPLICIT NONE
      INTEGER LIW, N
      INTEGER IW(LIW) 
      INTEGER SIZFI, SIZFR, IPTIW
      LOGICAL FREE
      LOGICAL SAME_PROC
      INTEGER NROW, NCOL, HF
      SAME_PROC = IW( IPTIW + 1 ) .LT. 0
      IF ( SAME_PROC ) THEN
        IF ( IW( IPTIW + 1 ) .LT. -N ) THEN
          FREE = .TRUE.
          NROW = -( IW( IPTIW + 1 ) + N )
        ELSE
          FREE = .FALSE.
          NROW = -IW( IPTIW + 1 )
        END IF
        NCOL = IW( IPTIW + 2 )
        SIZFI = 2
        SIZFR = NROW * NCOL
      ELSE
        NCOL  = IW( IPTIW + 1 ) + IW( IPTIW + 4 )
        NROW  = IW( IPTIW + 3 )
        SIZFR = NROW * NCOL
        HF = IW( IPTIW + 6 ) + 6
        IF ( IW( IPTIW + 4 ) .LT. 0 ) THEN
          SIZFI = NROW + IW( IPTIW + 1 ) + HF
        ELSE
          SIZFI = NROW + NCOL + HF
        END IF
        FREE = IW( IPTIW + 5 ) .eq. 0
      END IF
      RETURN
      END SUBROUTINE SMUMPS_90
      SUBROUTINE SMUMPS_91( N, KEEP28,
     *       IW, LIW, A, LA, NIRBDU, LRLU,
     *       IPTRLU, IWPOS, IWPOSCB, PTRIST, PTRAST, STEP, PIMASTER,
     *       PAMASTER, NBLOCKS_MAX, SIZE_PTR, PTRI,
     *       IBLOCBEG, IBLOCEND, RBLOCBEG, RBLOCEND, UNFINISHED )
      IMPLICIT NONE
      INTEGER N, LIW, LA, NIRBDU, LRLU, IPTRLU, IWPOS, IWPOSCB,
     *        NBLOCKS_MAX, SIZE_PTR, KEEP28
      INTEGER IW(LIW), PTRIST(KEEP28), PTRAST(KEEP28), STEP(N),
     * PIMASTER(KEEP28),
     *        PAMASTER(KEEP28), IBLOCBEG( NBLOCKS_MAX ),
     *        IBLOCEND( NBLOCKS_MAX ), RBLOCBEG( NBLOCKS_MAX ),
     *        RBLOCEND( NBLOCKS_MAX ), PTRI( SIZE_PTR, 2 )
      REAL A( LA )
      LOGICAL UNFINISHED
      INTEGER SIZFI, SIZFR
      LOGICAL FREE, LAST_FREE, TOP
      INTEGER I, J, K, ISIZEHOLE, RSIZEHOLE, IPTIW, IPTA
      J = 0
      DO I = 1, KEEP28
        IF ( PTRIST(I) .GT. IWPOSCB .and.
     *       PTRIST(I) .LT. NIRBDU )
     *    THEN
            J = J + 1
            PTRI( J, 1 ) = I
            PTRI( J, 2 ) = PTRIST(I)
        END IF
        IF ( PIMASTER( I ) .GT. IWPOSCB .and.
     *       PIMASTER( I ) .LT. NIRBDU )
     *    THEN
            J = J + 1
            PTRI( J, 1 ) = -I
            PTRI( J, 2 ) = PIMASTER( I )
        END IF
      END DO
      IF ( J .NE. SIZE_PTR ) THEN
        WRITE(*,*) 'Error 1 in COMPRESS, ', J, SIZE_PTR
        CALL SMUMPS_ABORT()
      END IF
      CALL SMUMPS_92( SIZE_PTR, PTRI( 1, 1 ), PTRI( 1, 2 ) )
      UNFINISHED = .FALSE.
      LAST_FREE  = .TRUE.
      TOP        = .TRUE.
      IPTIW      = IWPOSCB
      IPTA       = IPTRLU
      K          = 0
      DO WHILE ( IPTIW .NE. NIRBDU )
        CALL SMUMPS_90
     *  ( N, IPTIW, LIW, IW, SIZFI, SIZFR, FREE )
        IF ( TOP .and. FREE ) THEN
          IWPOSCB = IWPOSCB + SIZFI
          IPTRLU  = IPTRLU  + SIZFR
        ELSE
          IF ( .not. FREE ) THEN
            IF ( LAST_FREE ) THEN
              IF ( K .eq. 0 ) TOP = .FALSE.
              IF ( K .LT. NBLOCKS_MAX ) THEN
                K = K + 1
                IBLOCBEG( K ) = IPTIW + 1
                RBLOCBEG( K ) = IPTA  + 1
              ELSE 
                UNFINISHED = .TRUE.
              END IF
            END IF
            IBLOCEND( K ) = IPTIW + SIZFI
            RBLOCEND( K ) = IPTA  + SIZFR
          END IF
          LAST_FREE = FREE
        END IF
        IPTIW = IPTIW + SIZFI
        IPTA  = IPTA  + SIZFR
      END DO
      J = 1
      IF ( IBLOCEND( K ) .eq. NIRBDU ) THEN
        IPTIW = IBLOCBEG( K ) - 1
        IPTA  = RBLOCBEG( K ) - 1
 5      CONTINUE
        IF ( J .LE. SIZE_PTR ) THEN
          IF ( PTRI( J, 2 ) .GE. IBLOCBEG( K ) ) THEN
            J = J + 1
            GOTO 5
          END IF
        END IF
        K = K - 1
      ELSE
        IPTIW = NIRBDU
        IPTA  = LA
      END IF
      IF ( K .eq. 0 ) THEN
        RSIZEHOLE = 0
        ISIZEHOLE = 0
      END IF
      DO WHILE ( K .GE. 1 )
        ISIZEHOLE = IPTIW - IBLOCEND( K )
        RSIZEHOLE = IPTA  - RBLOCEND( K )
        DO I = IBLOCEND( K ), IBLOCBEG( K ), -1
          IW( IPTIW - IBLOCEND( K ) + I ) = IW( I )
        END DO
        IF ( RSIZEHOLE .NE. 0 ) THEN
          DO I = RBLOCEND( K ), RBLOCBEG( K ), -1
            A( IPTA - RBLOCEND( K ) + I ) = A( I )
          END DO
        END IF
 10     CONTINUE
        IF ( J .LE. SIZE_PTR ) THEN
          IF ( PTRI( J, 2 ) .GE. IBLOCBEG( K ) ) THEN
            I = PTRI( J, 1 )
            IF ( I > 0 ) THEN
              PTRIST(I) = PTRIST(I) + ISIZEHOLE
              IF ( PTRIST(I) .LT. NIRBDU ) THEN
                IF (IW(PTRIST(I)).LT.0) THEN
                  PAMASTER(I)=PAMASTER(I) + RSIZEHOLE
                ELSE
                  PTRAST(I)  =PTRAST(I) + RSIZEHOLE
                ENDIF
              END IF
            ELSE
              I = -I
              PIMASTER(I) = PIMASTER(I) + ISIZEHOLE
              PAMASTER(I) = PAMASTER(I) + RSIZEHOLE
            END IF
            J = J + 1
            GOTO 10
          END IF
        END IF
        IPTIW = IBLOCBEG( K ) + ISIZEHOLE - 1
        IPTA  = RBLOCBEG( K ) + RSIZEHOLE - 1
        K = K - 1
      END DO
      LRLU    = LRLU    + RSIZEHOLE
      IWPOSCB = IWPOSCB + ISIZEHOLE
      IPTRLU  = IPTRLU  + RSIZEHOLE
      RETURN
      END SUBROUTINE SMUMPS_91
      SUBROUTINE SMUMPS_92( N, I1, I2 )
      IMPLICIT NONE
      INTEGER N
      INTEGER I1( N ), I2( N )
      LOGICAL DONE
      INTEGER I, ISWAP
      DONE = .FALSE.
      DO WHILE ( .not. DONE )
        DONE = .TRUE.
        DO I =1 , N - 1
          IF ( I2( I ) .LT. I2( I + 1 ) ) THEN
            ISWAP = I2( I )
            I2( I ) = I2( I + 1 )
            I2( I + 1 ) = ISWAP
            ISWAP = I1( I )
            I1( I ) = I1( I + 1 )
            I1( I + 1 ) = ISWAP
            DONE = .FALSE.
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_92
      SUBROUTINE SMUMPS_27( id,  ANORMINF, LSCAL )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      INCLUDE 'smumps_struc.h'
      TYPE(SMUMPS_STRUC), TARGET :: id
      REAL, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      INTEGER, DIMENSION (:), POINTER :: KEEP, INFO
      LOGICAL :: I_AM_SLAVE
      REAL DUMMY(1)
      REAL ZERO
      PARAMETER( ZERO = 0.0E0)
      REAL, ALLOCATABLE :: SUMR(:), SUMR_LOC(:)
      INTEGER :: ALLOCOK, MTYPE, I
      INFO =>id%INFO
      KEEP =>id%KEEP
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               KEEP(46) .eq. 1 ) )
      IF (id%MYID .EQ. MASTER) THEN
       ALLOCATE( SUMR( id%N ), stat =allocok )
       IF (allocok .GT.0 ) THEN
        id%INFO(1)=-13
        id%INFO(2)=id%N
        RETURN
       ENDIF
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
             IF (.NOT.LSCAL) THEN
              CALL SMUMPS_207(id%A(1),
     *          id%NZ, id%N,
     *          id%IRN(1), id%JCN(1),
     *          SUMR, KEEP )
             ELSE
              CALL SMUMPS_289(id%A(1),
     *          id%NZ, id%N,
     *          id%IRN(1), id%JCN(1), 
     *          SUMR, KEEP, 
     *          id%COLSCA(1))
             ENDIF
            ELSE
             MTYPE = 1
             IF (.NOT.LSCAL) THEN
              CALL SMUMPS_119(MTYPE, id%N,
     *           id%NELT, id%ELTPTR,
     *           id%LELTVAR, id%ELTVAR,
     *           id%NA_ELT, id%A_ELT(1),
     *           SUMR, KEEP )
             ELSE
              CALL SMUMPS_135(MTYPE, id%N,
     *           id%NELT, id%ELTPTR(1),
     *           id%LELTVAR, id%ELTVAR(1),
     *           id%NA_ELT, id%A_ELT(1),
     *           SUMR, KEEP, id%COLSCA(1))
             ENDIF
            ENDIF
          ENDIF
      ELSE
          LSCAL = .FALSE.
          ALLOCATE( SUMR_LOC( id%N ), stat =allocok )
          IF (allocok .GT.0 ) THEN
             id%INFO(1)=-13
             id%INFO(2)=id%N
             RETURN
          ENDIF
          IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
           IF (.NOT.LSCAL) THEN
              CALL SMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc, 
     *          SUMR_LOC, id%KEEP )
           ELSE
              CALL SMUMPS_289(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc, 
     *          SUMR_LOC, id%KEEP, 
     *          id%COLSCA)
           ENDIF
          ELSE
           SUMR_LOC = REAL(ZERO)
          ENDIF
          IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SUMR_LOC, SUMR,
     *        id%N, MPI_REAL,
     *        MPI_SUM,MASTER,id%COMM, IERR)
          ELSE
              CALL MPI_REDUCE( SUMR_LOC, DUMMY,
     *        id%N, MPI_REAL,
     *        MPI_SUM,MASTER,id%COMM, IERR)
          END IF
        DEALLOCATE (SUMR_LOC)
      ENDIF
      IF ( id%MYID .eq. MASTER ) THEN
       ANORMINF = REAL(ZERO)
        IF (LSCAL) THEN
         DO I = 1, id%N
          ANORMINF = MAX(ABS(id%ROWSCA(I) * SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ELSE
         DO I = 1, id%N
          ANORMINF = MAX(ABS(SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ENDIF
      ENDIF
      CALL MPI_BCAST(ANORMINF, 1,
     *              MPI_REAL, MASTER,
     *              id%COMM, IERR )
      IF (id%MYID .eq. MASTER) DEALLOCATE (SUMR)
      RETURN
      END SUBROUTINE SMUMPS_27
      SUBROUTINE SMUMPS_248(N, A, LA, IW, LIW, WCB, LWCB, RHS,
     *    PTRICB, IWCB, LIWCB, W1, NE_STEPS, NA, LNA, STEP,
     *    FRERE, FILS,
     *    NSTK_S, IPOOL, LPOOL, PTRIST, MYLEAF, INFO, KEEP,
     *    PROCNODE_STEPS,
     *    SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,
     *    RHS_ROOT, MTYPE, 
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER N, LA, LIW, LWCB, LPOOL, LIWCB, LNA
      INTEGER SLAVEF, MYLEAF, COMM, MYID
      INTEGER INFO( 40 ), KEEP(150)
      INTEGER PROCNODE_STEPS( KEEP(28) )
      REAL A( LA ), RHS( N ), WCB( LWCB )
      REAL W1( N ), RHS_ROOT( * )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NA( LNA ), NE_STEPS( KEEP(28) )
      INTEGER STEP( N ), FRERE( KEEP(28) ), FILS( N )
      INTEGER NSTK_S(KEEP(28)), IPOOL( LPOOL )
      INTEGER PTRIST(KEEP(28))
      INTEGER PTRICB( KEEP(28) )
      INTEGER IW( LIW ), IWCB( LIWCB )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER  MSGTAG, MSGSOU, DUMMY(1)
      LOGICAL FLAG
      INTEGER NBFIN, MYROOT
      INTEGER POSIWCB,POSWCB,PLEFTWCB
      INTEGER INODE
      INTEGER I
      INTEGER III, NBROOT,NBLEAF,LEAF
      LOGICAL BLOQ
      EXTERNAL SMUMPS_275
      INTEGER SMUMPS_275
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      POSIWCB = LIWCB
      POSWCB  = LWCB
      PLEFTWCB= 1
      DO I = 1, KEEP(28)
        NSTK_S(I)   = NE_STEPS(I)
      ENDDO
      PTRICB = 0
      DO 10 I = 1, N
        W1(I)     = REAL(ZERO)
   10 CONTINUE
      CALL SMUMPS_362(N, LEAF, NBROOT, MYROOT, MYID,
     *     SLAVEF, NA, LNA, KEEP, STEP,
     *     PROCNODE_STEPS, IPOOL, LPOOL)
      NBFIN = SLAVEF
      IF ( MYROOT .EQ. 0 ) THEN
        NBFIN = NBFIN - 1
        DUMMY(1) = 1
        CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *       RACINE_SOLVE, SLAVEF)
      END IF
      MYLEAF = LEAF - 1
      III    = 1
   50 CONTINUE
      IF (SLAVEF .EQ. 1) THEN
        LEAF  = LEAF - 1
        INODE = IPOOL( LEAF )
        GOTO 60
      ENDIF
      BLOQ = ( III .EQ. LEAF )
      CALL SMUMPS_303( BLOQ, FLAG,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           MYID, SLAVEF, COMM,
     *           N, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *           IWCB, LIWCB,
     *           WCB, LWCB, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB, INFO, KEEP, STEP,
     *           PROCNODE_STEPS,W1)
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      IF (.not. FLAG) THEN
        IF (III .NE. LEAF) THEN
          LEAF = LEAF - 1
          INODE = IPOOL( LEAF )
          GOTO 60
        ENDIF
      ENDIF
      GOTO 50
 60   CONTINUE
      CALL SMUMPS_302( INODE, BUFR, LBUFR, LBUFR_BYTES,
     *        MSGTAG, MSGSOU, MYID, SLAVEF, COMM,  N,
     *        IPOOL, LPOOL, III, LEAF, NBFIN, NSTK_S,
     *        IWCB, LIWCB, WCB, LWCB, A, LA,
     *        IW, LIW, RHS, W1,
     *        POSWCB, PLEFTWCB, POSIWCB,
     *        PTRICB, PTRIST, PROCNODE_STEPS,
     *        FILS, STEP, FRERE,
     *        MYROOT, INFO, KEEP, RHS_ROOT, MTYPE, 
     *
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *      )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      GOTO 50
  260 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_248
      SUBROUTINE SMUMPS_272(BUFR,LBUFR,
     *     LBUFR_BYTES,
     *     root, N, IW, LIW, A, LA,
     *     NIRBDU, NBPROCFILS, LRLU, IPTRLU, IWPOS, IWPOSCB,
     *     PTRIST, PTLUST_S, PTRAST, STEP, PIMASTER, PAMASTER,
     *     COMP, LRLUS, IPOOL, LPOOL, LEAF,
     *     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *     KEEP, IFLAG, IERROR, COMM, COMM_LOAD, ITLOC,
     *     ND,PROCNODE_STEPS,SLAVEF )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC ) :: ROOT
      INTEGER KEEP( 150 )
      INTEGER LBUFR, LBUFR_BYTES, N, LIW, LA, NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, COMP, COMM, COMM_LOAD, IFLAG,
     *        IERROR, LRLUS
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LEAF )
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), PTLUST_S( KEEP(28) )
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)), ITLOC( N )
      INTEGER BUFR( LBUFR_BYTES ), NBPROCFILS( KEEP(28) )
      INTEGER IW( LIW )
      INTEGER ND(KEEP(28)), PROCNODE_STEPS(KEEP(28)),SLAVEF
      REAL A( LA )
      INTEGER   MYID
      INTEGER FILS( N ), PTRAIW(N), PTRARW( N )
      INTEGER INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
        INCLUDE 'mpif.h'
        INTEGER IERR
        INTEGER POSITION, LOCAL_M, LOCAL_N, POS_ROOT, LREQI, LREQA
        INTEGER NROW_SON, NCOL_SON, IROOT, ISON
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   ISON, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NROW_SON, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NCOL_SON, 1, MPI_INTEGER, COMM, IERR )
        IROOT = KEEP( 38 )
        IF ( PTRIST( STEP(IROOT) ) .NE. 0 .OR.
     *       PTLUST_S( STEP(IROOT)) .NE. 0 ) THEN
            NBPROCFILS(STEP(IROOT)) = NBPROCFILS(STEP(IROOT))-1
            IF ( NBPROCFILS( STEP(IROOT) ) .eq. 0 ) THEN
              CALL SMUMPS_507( N, IPOOL, LPOOL,
     *             PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     *             STEP, IROOT + N)
              IF (KEEP(47) .GE. 3) THEN
                 CALL SMUMPS_500(
     $                IPOOL, LPOOL, 
     *                PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *                MYID, STEP, N, ND, FILS )
              ENDIF
            END IF
        ELSE
            NBPROCFILS(STEP( IROOT ) ) = -1
            CALL SMUMPS_284( root, IROOT, N,
     *                     IW, LIW, A, LA,
     *                     NIRBDU,
     *                     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *                     LRLU, IPTRLU,
     *                     IWPOS, IWPOSCB, PTRIST, PTRAST,
     *                     STEP, PIMASTER, PAMASTER, ITLOC,
     *                     COMP, LRLUS, IFLAG, KEEP, IERROR )
            IF ( IFLAG .LT. 0 ) RETURN
        END IF
        LREQI = NROW_SON + NCOL_SON
        LREQA = NROW_SON * NCOL_SON
        IF ( (LREQA.NE.0) .AND.
     *       (PTRIST(STEP(IROOT)).LT.0)) THEN
         WRITE(*,*) ' Error in SMUMPS_272 '
         CALL SMUMPS_ABORT()
        ENDIF
        IF (LREQA.NE.0) THEN
          CALL SMUMPS_22(.FALSE.,
     *     MYID,N,KEEP,IW,LIW,A, LA, NIRBDU,
     *     LRLU, IPTRLU, IWPOS, IWPOSCB, PTRIST,
     *     PTRAST, STEP, PIMASTER, PAMASTER, ITLOC,
     *     LREQI, LREQA,
     *     COMP, LRLUS, IFLAG, IERROR )
          IF ( IFLAG .LT. 0 ) RETURN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   IW( IWPOSCB + 1 ), LREQI,
     *                   MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   A( IPTRLU + 1 ), LREQA,
     *                   MPI_REAL, COMM, IERR )
          IF ( PTRIST(STEP(IROOT)) .NE. 0 ) THEN
                 LOCAL_N  = -IW( PTRIST(STEP( IROOT ))     )
                 LOCAL_M  =  IW( PTRIST(STEP( IROOT )) + 1 )
                 POS_ROOT = PAMASTER(STEP( IROOT ))
          ELSE
                 LOCAL_N = IW( PTLUST_S(STEP( IROOT ) ) + 1 )
                 LOCAL_M = IW( PTLUST_S(STEP( IROOT ) ) + 2 )
                 POS_ROOT = IW( PTLUST_S(STEP( IROOT )) + 4 )
          END IF
          CALL SMUMPS_38( NROW_SON, NCOL_SON, IW( IWPOSCB + 1 ),
     *                     IW( IWPOSCB + NROW_SON + 1 ),
     *                     A( IPTRLU + 1 ),
     *                     A( POS_ROOT ), LOCAL_M, LOCAL_N )
          IWPOSCB = IWPOSCB + LREQI
          IPTRLU  = IPTRLU  + LREQA
          LRLU    = LRLU    + LREQA
          LRLUS   = LRLUS   + LREQA
          CALL SMUMPS_471(.FALSE.,LA-LRLUS,0,-LREQA,KEEP)
        ENDIF
      RETURN
      END SUBROUTINE SMUMPS_272
      SUBROUTINE SMUMPS_224(NFRONT,NASS,IBEGKJI, LPIV, TIPIV,
     *    N,INODE,IW,LIW,A,LA,
     *    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU,SEUIL)
      IMPLICIT NONE
      INTEGER IBEGKJI, LPIV 
      INTEGER TIPIV(LPIV)
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV,NOFFW
      REAL A(LA) 
      REAL UU, SEUIL
      INTEGER IW(LIW) 
      INTEGER IOLDPS, POSELT
      REAL SWOP
      INTEGER APOS,ILOC
      REAL ZERO, RMAX, AMROW
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW, HF
      INTEGER ISAMAX
      INTRINSIC MAX
      DATA ZERO /0.0E0/
        INOPV   = 0
        NPIV    = IW(IOLDPS+1)
        HF = 6 + IW(IOLDPS+5)
        NPIVP1  = NPIV + 1
        ILOC    = NPIVP1 - IBEGKJI + 1
        TIPIV(ILOC) = ILOC
        NASSW   = IABS(IW(IOLDPS+3))
          DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + NFRONT*(IPIV-1) + NPIV
            JMAX = 1
            IF (UU.GT.ZERO) GO TO 340
            IF (A(APOS).EQ.ZERO) GO TO 630
            GO TO 380
  340       AMROW = ZERO
            J1 = APOS
            J2 = APOS - NPIV + NASS - 1
             J3    = NASS -NPIV
             JMAX  = ISAMAX(J3,A(J1),1)
             JJ    = JMAX + J1 - 1
             AMROW = ABS(A(JJ))
            RMAX = AMROW
            J1 = J2 + 1
            J2 = APOS - NPIV + NFRONT - 1
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
  360       CONTINUE
  370       IDIAG = APOS + IPIV - NPIVP1
            IF (ABS(A(IDIAG)).GT.MAX(UU*RMAX,SEUIL)) THEN
               JMAX = IPIV - NPIV
               GO TO 380
            ENDIF
            IF (AMROW.LE.MAX(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV*NFRONT
            J2 = J1 + NFRONT - 1
            J3 = POSELT + (IPIV-1)*NFRONT
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1
  390       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NPIVP1
            ISWPS2 = IOLDPS + HF - 1 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            TIPIV(ILOC) = ILOC + JMAX - 1
            J1 = POSELT + NPIV
            J2 = POSELT + NPIV + JMAX - 1
            DO 410 KSW=1,NASS
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT
              J2 = J2 + NFRONT
  410       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + HF - 1 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
      WRITE(*,*) 'NIV2:Detected 0 pivot, INODE,NPIV=',INODE,NPIV
  420 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_224
      SUBROUTINE  SMUMPS_294( COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE,
     *             IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDA_FS, 
     *             IBEGKJI, IEND, TIPIV, LPIV, LASTBL, NB_BLOC_FAC,
     *
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL, 
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, STEP, PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, IBEGKJI, IEND, LPIV, 
     *        IOLDPS, POSELT, LA, LDA_FS, NB_BLOC_FAC
      INTEGER IW(LIW), TIPIV(LPIV)
      LOGICAL LASTBL
      REAL A(LA)
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT ) 
      INTEGER KEEP(150)
      INTEGER NBFIN, IFLAG, IERROR, LEAF, LPOOL,
     *        SLAVEF, ICNTL(40)
      INTEGER POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, 
     *        NIRBDU, COMP
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     *        ITLOC(N), FILS(N),
     *        PTRARW(LPTRAR), PTRAIW(LPTRAR), 
     *        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), PTLUST_S(KEEP(28)),
     *        STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)), NSTK_S(KEEP(28)),
     *        NBPROCFILS(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      REAL OPASSW, OPELIW
      REAL DBLARR(MAX(1,KEEP(13)))
      EXTERNAL  SMUMPS_329
      INTEGER NPIV, NCOL, APOS, PDEST, NSLAVES
      INTEGER IERR, IERR_MPI, LREQA, LREQI
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      NSLAVES= IW(IOLDPS+5)
          IF (NSLAVES.EQ.0) THEN
           WRITE(6,*) ' ERROR 1 in SMUMPS_294 '
           CALL SMUMPS_ABORT()
          ENDIF
      NPIV   = IEND - IBEGKJI + 1
      NCOL   = LDA_FS - IBEGKJI + 1
      APOS   = POSELT + LDA_FS*(IBEGKJI-1) + IBEGKJI - 1
      IF ((NPIV.GT.0) .OR. 
     *    ((NPIV.EQ.0).AND.(LASTBL)) ) THEN
        PDEST  = IOLDPS + 6
        IERR = -1
        IF ( NPIV .NE. 0 ) THEN
          NB_BLOC_FAC = NB_BLOC_FAC + 1
        END IF
        DO WHILE (IERR .EQ.-1)
          CALL SMUMPS_65( INODE, LDA_FS, NCOL, 
     *               NPIV, FPERE, LASTBL, TIPIV, A(APOS),
     *               IW(PDEST), NSLAVES, KEEP(50), NB_BLOC_FAC,
     *               COMM, IERR )
        IF (IERR.EQ.-1) THEN
           BLOCKING  = .FALSE.
           SET_IRECV = .TRUE.
           MESSAGE_RECEIVED = .FALSE.
           CALL SMUMPS_329( COMM_LOAD, ASS_IRECV, 
     *      BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *      MPI_ANY_SOURCE, MPI_ANY_TAG,
     *      STATUS, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *      IERROR, COMM,
     *      NBPROCFILS,
     *      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *      root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *      INTARR, DBLARR, ICNTL, KEEP, ND, FRERE,
     *      LPTRAR, NELT, FRTPTR, FRTELT, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE )
           IF (MESSAGE_RECEIVED) POSELT = PTRAST(STEP(INODE))
           IF ( IFLAG .LT. 0 ) GOTO 500
         ENDIF
        ENDDO
        IF (IERR .EQ. -2) GOTO 300
      ENDIF
      GOTO 500
  300 CONTINUE
      IFLAG = -17
      LREQA = NCOL*NPIV  
      LREQI = NPIV + 6 + 2*NSLAVES
      IERROR =  LREQI  * KEEP( 34 )+ LREQA * KEEP( 35 )
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
 500  RETURN
      END SUBROUTINE  SMUMPS_294
      SUBROUTINE SMUMPS_249(N, A, LA, IW, LIW, W, LWC, RHS, RHSSOL, 
     *    PTRICB, PTRACB, IWCB, LIWW, POOLSS, W2, 
     *    NE_STEPS, NA, LNA, STEP,
     *    FRERE, FILS, IPOOL, LPOOL, PTRIST, MYLEAF, INFO, 
     *    PROCNODE_STEPS,
     *    SLAVEF, COMM,MYID, BUFR, LBUFR, LBUFR_BYTES,
     *    KEEP, RHS_ROOT, MTYPE, 
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER N,LA,LIW,LIWW,LWC,LPOOL,LNA
      INTEGER SLAVEF,MYLEAF,COMM,MYID
      INTEGER KEEP( 150 )
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER NA(LNA),NE_STEPS(KEEP(28))
      INTEGER IPOOL(LPOOL),POOLSS(N)
      INTEGER INFO(40)
      INTEGER PTRIST(KEEP(28)),PTRICB(KEEP(28)),PTRACB(KEEP(28))
      REAL A(LA), RHS(N), W(LWC),RHSSOL(N) 
      REAL W2(N)
      INTEGER IW(LIW),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)),FILS(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      REAL RHS_ROOT( * )
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      LOGICAL FLAG
      INTEGER POSIWCB,POSWCB
      INTEGER APOS,APOSCB,NPIV
      INTEGER IPOS,IPOSCB,LIELL,NELIM,IFR,JJ,I
      INTEGER J1,J2,J,IST,NCB,NBFINF
      INTEGER NBLEAF,INODE,NBROOT,NROOT,NBFILS
      INTEGER IN,IF,LONG,IFRR
      INTEGER III,IIPOOL,MYLEAFE
      INTEGER IISSAR
      INTEGER NSLAVES
      REAL ALPHA,ONE,ZERO
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0, ZERO=0.0E0)
      LOGICAL BLOQ,DEBUT
      INTEGER PROCDEST, DEST
      INTEGER SSII,POSII, POSINDICES
      INTEGER DUMMY(1)
      INTEGER PLEFTW
      INTEGER Offset, EffectiveSize, ISLAVE, FirstIndex
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      LOGICAL SMUMPS_283
      INTEGER SMUMPS_330
      EXTERNAL SGEMV, STRSV, SMUMPS_283, SMUMPS_330
      PLEFTW = 1
      POSIWCB = LIWW
      POSWCB = LWC
      NBFINF = SLAVEF
      MYLEAFE = MYLEAF
      IF (MYLEAFE .EQ. 0) THEN
        CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, FEUILLE,
     *                  SLAVEF)
        NBFINF = NBFINF - 1
        IF (NBFINF .EQ. 0) GOTO 340
      ENDIF
      NROOT = 0
      NBLEAF = NA(1)
      NBROOT = NA(2)
      DO I = 1, NBROOT
        INODE = NA(NBLEAF+I+2)
        IF (SMUMPS_275(STEP(INODE),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
          NROOT = NROOT + 1
          IPOOL(NROOT) = INODE
        ENDIF
      END DO
      III = 1
      IIPOOL = NROOT + 1
      IISSAR = 1
   50 CONTINUE
      BLOQ = (III .EQ. IIPOOL .AND. IISSAR .LE. 1 )
      CALL SMUMPS_41( BLOQ, FLAG, BUFR, LBUFR,
     *     LBUFR_BYTES, MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, III, IISSAR,
     *     NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *     IPOOL, LPOOL, STEP,  FRERE, FILS, PROCNODE_STEPS,
     *     PLEFTW, KEEP,
     *     PTRIST, IW, LIW, A, LA, W2, MYLEAFE, 
     *     RHSSOL, RHS )
      IF ( INFO(1) .LT. 0 ) RETURN
      IF ( .NOT. FLAG ) THEN
        IF (III .NE. IIPOOL) THEN
          INODE = IPOOL(III)
          III = III + 1
          GO TO 60
        ELSE IF ( IISSAR .GT. 1 ) THEN
          GO TO 200
        ENDIF
      END IF
      IF ( NBFINF .eq. 0 ) GOTO 340
      GOTO 50
   60 CONTINUE
      IF ( INODE .EQ. KEEP( 38 ) .OR. INODE .EQ. KEEP( 20 ) ) THEN
        IPOS = PTRIST(STEP(INODE))
        LIELL = IW(IPOS+3)
        NPIV  = LIELL
        IPOS =  PTRIST(STEP(INODE)) + 5
        IF ( MTYPE .EQ. 1 ) THEN
          J1   = IPOS + LIELL + 1
          J2   = IPOS + LIELL + NPIV
        ELSE
          J1   = IPOS + 1
          J2   = IPOS + NPIV
        END IF
        IFR  = 0
        DO JJ = J1, J2
          J  = IW( JJ )
          IFR = IFR + 1
          RHSSOL(J) = RHS_ROOT(IFR)
        END DO 
        IN = INODE
  270   IN = FILS(IN)
        IF (IN .GT. 0) GOTO 270
        IF (IN .EQ. 0) THEN
          MYLEAFE = MYLEAFE - 1
          IF (MYLEAFE .EQ. 0) THEN
            CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                       FEUILLE, SLAVEF )
            NBFINF = NBFINF - 1
            IF (NBFINF .EQ. 0) GOTO 340
          ENDIF
          GOTO 50
        ENDIF
        IF   = -IN
        LONG = NPIV
        NBFILS = NE_STEPS(STEP(INODE))
        DEBUT = .TRUE.
        DO I = 0, SLAVEF - 1
          DEJA_SEND( I ) = .FALSE.
        END DO
        DO I = 1, NBFILS
          IF (SMUMPS_275(STEP(IF),PROCNODE_STEPS,SLAVEF)
     *       .EQ. MYID) THEN
            IPOOL(IIPOOL) = IF
            IIPOOL = IIPOOL + 1
          ELSE
            PROCDEST = SMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *                 SLAVEF)
            IF (.NOT. DEJA_SEND( PROCDEST ))  THEN
 600          CALL SMUMPS_78( IF, 0, 0, LONG, IW( J1 ),
     *           RHS_ROOT( 1 ), PROCDEST,
     *           NOEUD, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL SMUMPS_41(
     *          .FALSE., FLAG,
     *          BUFR, LBUFR, LBUFR_BYTES,
     *          MYID, SLAVEF, COMM,
     *          N, IWCB, LIWW, POSIWCB,
     *          W, LWC, POSWCB,
     *          IIPOOL, III, IISSAR,
     *          NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *          IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *          PLEFTW, KEEP,
     *          PTRIST, IW, LIW, A, LA, W2, MYLEAFE, 
     *          RHSSOL, RHS )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 600
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = LONG * KEEP(35) +
     *                      ( LONG + 2 ) * KEEP(34)
                GOTO 330
              END IF
              DEJA_SEND( PROCDEST ) = .TRUE.
            END IF
            IF ( IERR .NE. 0 ) CALL SMUMPS_ABORT()
          ENDIF
          IF = FRERE(STEP(IF))
        ENDDO
        GOTO 50
      END IF
      IF ( SMUMPS_283( STEP(INODE), PROCNODE_STEPS,
     *     SLAVEF ) ) THEN
        POOLSS(IISSAR) = INODE
        IISSAR = IISSAR + 1
        GO TO 50
      END IF
      IPOS = PTRIST(STEP(INODE)) + 2
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      NELIM = IW(IPOS-1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      IPOS = IPOS + 1
      APOS = IW(IPOS)
      NSLAVES = IW( PTRIST(STEP(INODE)) + 5 )
      IPOS = IPOS + 1 + NSLAVES
      IFRR = 0
      LONG = 0
      IF ( MTYPE .eq. 1 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ELSE
        J1 = IPOS + LIELL + 1
        J2 = IPOS + NPIV + LIELL
      END IF
      IF ( POSIWCB - LIELL - 2 .LT. 0 .or.
     *     POSWCB - LIELL .LT. PLEFTW - 1 ) THEN
        CALL SMUMPS_95( N, KEEP(28), IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
        IF ( POSWCB - LIELL .LT. PLEFTW - 1 ) THEN
          INFO( 1 ) = -11
          INFO( 2 ) = LIELL - POSWCB - PLEFTW + 1
          GOTO 330
        END IF
        IF ( POSIWCB - LIELL - 2 .LT. 0 ) THEN
          INFO( 1 ) = -14
          INFO( 2 ) = LIELL + 2 - POSIWCB
          GO TO 330
        END IF
      END IF
      POSIWCB = POSIWCB - LIELL - 2
      POSWCB  = POSWCB - LIELL
      PTRICB(STEP( INODE )) = POSIWCB + 1
      PTRACB(STEP( INODE )) = POSWCB  + 1
      IWCB( PTRICB(STEP( INODE ))     ) = LIELL
      IWCB( PTRICB(STEP( INODE )) + 1 ) = 1
      IF ( MTYPE .eq. 1 ) THEN
         POSINDICES = IPOS + LIELL + 1
      ELSE
         POSINDICES = IPOS + 1
      END IF
      IWCB( PTRICB(STEP( INODE )) + 2 :
     * PTRICB(STEP( INODE )) + 1 + LIELL )
     *   = IW( POSINDICES : POSINDICES + LIELL - 1 )
      IFR = PTRACB(STEP( INODE )) - 1
      DO JJ = J1, J2
        J = IW(JJ)
        IFR = IFR + 1
        W(IFR) = RHS(J)
      END DO
      IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          J1 = IPOS + LIELL + NPIV + 1
          J2 = IPOS + 2 * LIELL
        ELSE
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
        END IF
        DO 130 JJ = J1, J2
          J = IW(JJ)
          IFR = IFR + 1
          W(IFR) = RHSSOL(J)
  130   CONTINUE
        NCB = LIELL - NPIV
        IF (NPIV .EQ. 0) GOTO 160
        IF ( MTYPE .eq. 1 ) THEN
          IST = APOS + NPIV
          CALL SGEMV( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     *                W(NPIV + PTRACB(STEP(INODE))), 1, ONE,
     *                W(PTRACB(STEP(INODE))), 1 )
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
            IST = APOS + NPIV * LIELL
          ELSE
            IST = APOS + NPIV * NPIV
          END IF
          IF ( SMUMPS_330(STEP(INODE),PROCNODE_STEPS,
     *         SLAVEF).eq.2 ) THEN
            IF ( NELIM .ne. 0 )
     *      CALL SGEMV( 'N', NPIV, NELIM, ALPHA, A( IST ), NPIV,
     *                W( NPIV + PTRACB(STEP(INODE)) ), 1, ONE,
     *                W(PTRACB(STEP(INODE))), 1 )
            Offset        = NPIV + NELIM
            IWCB( PTRICB(STEP( INODE )) + 1 ) =
     *      IWCB( PTRICB(STEP( INODE )) + 1 ) + NSLAVES
            DO ISLAVE = 1, NSLAVES
              CALL SMUMPS_49( 
     &                KEEP, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB - NELIM, 
     &                NSLAVES, 
     &                EffectiveSize,
     &                FirstIndex )
 500         DEST = IW( PTRIST(STEP(INODE))+5+ISLAVE)
             CALL SMUMPS_63(INODE,
     *             W(Offset+PTRACB(STEP(INODE))), EffectiveSize,
     *             DEST,
     *             BACKSLV_MASTER2SLAVE,
     *             COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL SMUMPS_41(
     *          .FALSE., FLAG,
     *          BUFR, LBUFR, LBUFR_BYTES,
     *          MYID, SLAVEF, COMM,
     *          N, IWCB, LIWW, POSIWCB,
     *          W, LWC, POSWCB,
     *          IIPOOL, III, IISSAR,
     *          NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *          IPOOL, LPOOL, STEP, FRERE, FILS,
     *          PROCNODE_STEPS, PLEFTW, KEEP,
     *          PTRIST, IW, LIW, A, LA, W2, MYLEAFE, 
     *          RHSSOL, RHS )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 500
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = EffectiveSize * KEEP(35) +
     *                            2 * KEEP(34)
                GOTO 330
              END IF
              Offset = Offset + EffectiveSize
            END DO
            GOTO 50
          ELSE
            CALL SGEMV( 'N', NPIV, NCB, ALPHA, A( IST ), NPIV,
     *                W( NPIV + PTRACB(STEP(INODE)) ), 1, ONE,
     *                W(PTRACB(STEP(INODE))), 1 )
          END IF ! node level =  2
        END IF ! MTYPE == 1
      ENDIF  ! LIELL > NPIV
      IF ( MTYPE .eq. 1 ) THEN
        CALL STRSV('L', 'T', 'N', NPIV, A(APOS), LIELL,
     *              W(PTRACB(STEP(INODE))), 1)
      ELSE
        IF ( KEEP(50) .EQ. 0 ) THEN
          CALL STRSV('U','N','U', NPIV, A(APOS), LIELL,
     *              W(PTRACB(STEP(INODE))), 1)
        ELSE
          CALL STRSV('U','N','U', NPIV, A(APOS), NPIV,
     *              W(PTRACB(STEP(INODE))), 1)
        END IF
      END IF
      IF ( MTYPE .eq. 1 ) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      DO 150 I = 1, NPIV
        JJ = IW(J1 + I - 1)
        RHSSOL(JJ) = W(PTRACB(STEP(INODE))+I-1)
  150 CONTINUE
  160 IN = INODE
  170 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 170
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (MYLEAFE .EQ. 0) THEN
          CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                     FEUILLE, SLAVEF )
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        GOTO 50
      ENDIF
      IF = -IN
      LONG = LONG + NPIV
      NBFILS = NE_STEPS(STEP(INODE))
      DEBUT = .TRUE.
      DO i = 0, SLAVEF - 1
        DEJA_SEND( i ) = .FALSE.
      END DO
      IF ( MTYPE .eq. 1 ) THEN
         POSINDICES = IPOS + LIELL + 1
      ELSE
         POSINDICES = IPOS + 1
      END IF
      DO 190 I = 1, NBFILS
        IF (SMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
          IPOOL(IIPOOL) = IF
          IIPOOL = IIPOOL + 1
          IF = FRERE(STEP(IF))
        ELSE
          PROCDEST = SMUMPS_275(STEP(IF),PROCNODE_STEPS,SLAVEF)
          IF (.not. DEJA_SEND( PROCDEST ))  THEN
 400        CONTINUE
            CALL SMUMPS_78( IF, 0, 0, LIELL,
     *         IWCB( PTRICB(STEP( INODE )) + 2),
     *         W   ( PTRACB(STEP( INODE ))), PROCDEST,
     *         NOEUD, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
              CALL SMUMPS_41(
     *        .FALSE., FLAG,
     *        BUFR, LBUFR, LBUFR_BYTES,
     *        MYID, SLAVEF, COMM,
     *        N, IWCB, LIWW, POSIWCB,
     *        W, LWC, POSWCB,
     *        IIPOOL, III, IISSAR,
     *        NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *        IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *        PLEFTW, KEEP,
     *        PTRIST, IW, LIW, A, LA, W2, MYLEAFE, 
     *        RHSSOL, RHS )
              IF ( INFO( 1 ) .LT. 0 ) GOTO 340
              GOTO 400
            ELSE IF ( IERR .EQ. -2 ) THEN
              INFO( 1 ) = -17
              INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
              GOTO 330
            END IF
            DEJA_SEND( PROCDEST ) = .TRUE.
          END IF
          IF = FRERE(STEP(IF))
        ENDIF
  190 CONTINUE
      IWCB(PTRICB(STEP(INODE))+1) = IWCB(PTRICB(STEP(INODE))+1)-1
      CALL SMUMPS_151(N, KEEP(28), IWCB, LIWW, W, LWC,
     *     POSWCB,POSIWCB,PTRICB,PTRACB)
      GOTO 50
  200 CONTINUE
      IISSAR = IISSAR - 1
      INODE  = POOLSS( IISSAR )
      SSII  = IISSAR
      POSII = IISSAR + 1
  230 CONTINUE
      POSII = POSII - 1
      INODE = POOLSS( POSII )
      IPOS = PTRIST(STEP(INODE)) + 2
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      IF (NPIV .EQ. 0) GOTO 280
      IPOS = IPOS + 1
      APOS = IW(IPOS)
      NSLAVES = IW( PTRIST(STEP(INODE)) +5 )
      IPOS = IPOS + 1 + NSLAVES
      IF ( MTYPE .eq. 1 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ELSE
        J1 = IPOS + 1 + LIELL
        J2 = IPOS + NPIV + LIELL
      END IF
      IFR = 0
      IF ( POSWCB .LT. LIELL ) THEN
          CALL SMUMPS_95( N, KEEP(28), IWCB, LIWW, W, LWC,
     *                 POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB .LT. LIELL ) THEN
            INFO(1) = -11
            INFO(2) = LIELL - POSWCB
            GOTO 330
          END IF
      END IF
      DO 240 JJ = J1, J2
        J = IW(JJ)
        IFR = IFR + 1
        W(IFR) = RHS(J)
  240 CONTINUE
      IF (LIELL .GT. NPIV) THEN
        IF ( MTYPE .eq. 1 ) THEN
          J1 = IPOS + LIELL + NPIV + 1
          J2 = IPOS + 2 * LIELL
        ELSE
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
        END IF
        DO 250 JJ = J1, J2
          J = IW(JJ)
          IFR = IFR + 1
          W(IFR) = RHSSOL(J)
  250   CONTINUE
        NCB = LIELL - NPIV
        IF ( MTYPE .EQ. 1 ) THEN
          IST = APOS + NPIV
          CALL SGEMV( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     *                W(NPIV + 1), 1, ONE, W(1), 1 )
        ELSE
          IF ( KEEP(50) .eq.0 ) THEN
            IST = APOS + NPIV * LIELL
          ELSE
            IST = APOS + NPIV * NPIV
          END IF
          CALL SGEMV( 'N', NPIV, NCB, ALPHA, A(IST), NPIV,
     *                W(NPIV + 1), 1, ONE, W(1), 1 )
        END IF
      ENDIF
      IF ( MTYPE .eq. 1 ) THEN
        CALL STRSV('L', 'T', 'N', NPIV, A(APOS), LIELL, W(1), 1)
      ELSE
        IF ( KEEP(50) .eq. 0 ) THEN
        CALL STRSV('U', 'N', 'U', NPIV, A(APOS), LIELL, W(1), 1)
        ELSE
        CALL STRSV('U', 'N', 'U', NPIV, A(APOS), NPIV, W(1), 1)
        END IF
      END IF
      IFR = 0
      IF ( MTYPE .eq. 1 ) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      J2 = J1 + NPIV - 1
      DO I = 1, NPIV
        JJ = IW(J1 + I - 1)
        RHSSOL(JJ) = W(I)
      ENDDO
  280 IN = INODE
  290 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 290
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (MYLEAFE .EQ. 0) THEN
          CALL SMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID,
     *        COMM, FEUILLE, SLAVEF)
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        IF (SSII .EQ. POSII) GOTO 50
        GOTO 230
      ENDIF
      IF = -IN
      NBFILS = NE_STEPS(STEP(INODE))
      DO 300 I = 1, NBFILS
        POOLSS(POSII) = IF
        POSII = POSII + 1
        IF = FRERE(STEP(IF))
  300 CONTINUE
      GOTO 230
  330 CONTINUE
      CALL SMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, TERREUR,
     * SLAVEF)
  340 RETURN
      END SUBROUTINE SMUMPS_249
      SUBROUTINE SMUMPS_273( ROOT, 
     *    INODE, NELIM, NSLAVES, ROW_LIST,
     *    COL_LIST, SLAVE_LIST, 
     *
     *    PROCNODE_STEPS, IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST, PTLUST_S,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, ITLOC, COMP,
     *    IFLAG, IERROR, 
     *    IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,
     *    COMM,COMM_LOAD,FILS,ND )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: ROOT
      INTEGER INODE, NELIM, NSLAVES 
      INTEGER KEEP( 150 )
      INTEGER ROW_LIST(*), COL_LIST(*), 
     *        SLAVE_LIST(*)
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST( KEEP(28) ), PTLUST_S(KEEP(28)), PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), ITLOC( N ), PROCNODE_STEPS( KEEP(28) )
      INTEGER IFLAG, IERROR
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF
      INTEGER COMM,COMM_LOAD,ND(KEEP(28)),FILS(N)
      INTEGER IROOT, TYPE_INODE, DEB_ROW, DEB_COL,
     *        NOINT, NOREAL
      INCLUDE 'mumps_tags.h'
      INTEGER SMUMPS_330
      EXTERNAL SMUMPS_330
      IROOT        = KEEP(38)
      NSTK_S(STEP(IROOT))= NSTK_S(STEP(IROOT)) - 1
      KEEP(42) = KEEP(42) + NELIM
      TYPE_INODE= SMUMPS_330(STEP(INODE), PROCNODE_STEPS, SLAVEF)
      IF (TYPE_INODE.EQ.1) THEN 
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + 1
        ELSE 
         KEEP(41) = KEEP(41) + 3
        ENDIF
      ELSE
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + NSLAVES
        ELSE 
         KEEP(41) = KEEP(41) + 2*NSLAVES + 1
        ENDIF
      ENDIF
      IF  (NELIM.EQ.0) THEN
        PIMASTER(STEP(INODE)) = 0 
      ELSE
       NOINT = 6 + NSLAVES + NELIM  + NELIM
       NOREAL= 0
       CALL SMUMPS_22(.FALSE.,
     *   MYID,N,KEEP,IW,LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER, PAMASTER, ITLOC,
     *   NOINT, NOREAL,
     *   COMP, LRLUS, IFLAG, IERROR )
       IF ( IFLAG .LT. 0 ) THEN
         WRITE(*,*) ' Failure in int space allocation in CB area ',
     *    ' during assembly of root : SMUMPS_273',
     *    ' size required was :', NOINT,
     *    'INODE=',INODE,' NELIM=',NELIM, ' NSLAVES=', NSLAVES
         RETURN
        ENDIF
        PIMASTER(STEP( INODE )) = IWPOSCB + 1
        PAMASTER(STEP( INODE )) = IPTRLU  + 1
        IW( IWPOSCB + 1 ) = 2*NELIM
        IW( IWPOSCB + 2 ) = NELIM
        IW( IWPOSCB + 3 ) = 0
        IW( IWPOSCB + 4 ) = 0
        IW( IWPOSCB + 5 ) = 1
        IW( IWPOSCB + 6 ) = NSLAVES
        IF (NSLAVES.GT.0) THEN
         IW( IWPOSCB+7: IWPOSCB+7+NSLAVES-1) = 
     *                   SLAVE_LIST(1:NSLAVES)
        ENDIF
        DEB_ROW = IWPOSCB+7+NSLAVES
        IW(DEB_ROW : DEB_ROW+NELIM -1) = ROW_LIST(1:NELIM)
        DEB_COL = DEB_ROW + NELIM
        IW(DEB_COL : DEB_COL+NELIM -1) = COL_LIST(1:NELIM)
      ENDIF
      IF (NSTK_S(STEP(IROOT)) .EQ. 0 ) THEN
          CALL SMUMPS_507(N, IPOOL, LPOOL, PROCNODE_STEPS,
     *         SLAVEF, KEEP(28), KEEP(76), STEP, IROOT )
          IF (KEEP(47) .GE. 3) THEN
             CALL SMUMPS_500(
     $            IPOOL, LPOOL, 
     *            PROCNODE_STEPS, KEEP, SLAVEF, COMM_LOAD,
     *            MYID, STEP, N, ND, FILS )
          ENDIF
      END IF
      RETURN
      END SUBROUTINE SMUMPS_273
