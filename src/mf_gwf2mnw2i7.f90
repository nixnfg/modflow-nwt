MODULE GWFMNW2IMODULE
   INTEGER,SAVE,POINTER  ::Wel1flag,QSUMflag,BYNDflag,MNWOBS
   CHARACTER(LEN=20),SAVE, DIMENSION(:),   POINTER     ::MNWIID
   DOUBLE PRECISION, SAVE, DIMENSION(:,:), POINTER     ::MNWILST
   TYPE GWFMNWITYPE
      INTEGER,POINTER  ::Wel1flag,QSUMflag,BYNDflag,MNWOBS
      CHARACTER(LEN=20),DIMENSION(:),   POINTER     ::MNWIID
      DOUBLE PRECISION, DIMENSION(:,:), POINTER     ::MNWILST
   END TYPE
   TYPE(GWFMNWITYPE), SAVE:: GWFMNWIDAT(10)
END MODULE GWFMNW2IMODULE
   !   GZH  20080208
   !   LFK  March 2011  Revision in GWF2MNW2I7OT to fix WEL1 output file for inactive wells.
   !   LFK  Nov. 2012   Additional revisions--see read_me file.
   !
   ! GWF2MNW2I7AR READ INIT DATA AND ALLOCATE SPACE FOR MNW WELLS DESIGNATED FOR OBSERVATION
   !
   !     ******************************************************************
   !
SUBROUTINE GWF2MNW2I7AR(INMNWI,INMNW2,IGRID)
   !
   !     ******************************************************************
   !
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,        ONLY:IOUT
   USE GWFMNW2IMODULE, ONLY:Wel1flag,QSUMflag,BYNDflag,MNWOBS,&
   &MNWIID,MNWILST
   !
   IF(INMNWI.GT.0.AND.INMNW2.LE.0) THEN
      WRITE(IOUT,*) '***ERROR*** : MNWI PACKAGE CAN ONLY BE&
      &USED IF MNW2 PACKAGE IS ACTIVE'
      STOP 'MNWI ERROR'
   END IF
   !1------Allocate scalar variables, which makes it possible for multiple
   !1------grids to be defined.
   ALLOCATE(Wel1flag,QSUMflag,BYNDflag,MNWOBS)
   !
   IF(INMNWI.EQ.0) THEN
      LCMNIO=1
   ELSE
   !     if transport on, read concflag
      READ(INMNWI,*) Wel1flag,QSUMflag,BYNDflag
      WRITE(IOUT,*) 'MNWI Package input:'
      WRITE(IOUT,*) 'Wel1flag = ',Wel1flag
      WRITE(IOUT,*) 'QSUMflag = ',QSUMflag
      WRITE(IOUT,*) 'BYNDflag = ',BYNDflag
      WRITE(IOUT,*)
   !
      READ(INMNWI,*) MNWOBS
      IF(MNWOBS.LT.0) THEN
         WRITE(IOUT,*) 'MNWOBS MUST BE > 0'
         STOP
      END IF
   !
   !5------ALLOCATE SPACE FOR MNWILST ARRAY.
   !5------FOR EACH OBS WELL, THERE ARE 6 DATA VALUES
      NMNWIVL=6
      ALLOCATE (MNWILST(NMNWIVL,MNWOBS))
   !5------ALLOCATE SPACE FOR MNWIID ARRAY.
      ALLOCATE (MNWIID(MNWOBS+1))
   END IF
   !7------SAVE POINTERS TO DATA AND RETURN.
   CALL SGWF2MNW2IPSV(IGRID)
   !
   RETURN
END
   !
   !_________________________________________________________________________________
   !
   !
   !  GWF2MNW2I7RP READ INPUT FILE FOR MNW2 WELLS DESIGNATED FOR OBSERVATION
   !
   !     ******************************************************************
   !
SUBROUTINE GWF2MNW2I7RP(INMNWI,GWTUNIT,IGRID)
   !
   !     ******************************************************************
   !
   !     READ LOCATIONS OF MNW2 WELLS DESIGNATED FOR OBSERVATION
   !     ******************************************************************
   !
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:IOUT
   USE GWFMNW2MODULE, ONLY:MNWMAX,MNW2,WELLID
   USE GWFMNW2IMODULE, ONLY:Wel1flag,QSUMflag,BYNDflag,MNWOBS,&
   &MNWILST,MNWIID
   !     ------------------------------------------------------------------
   INTEGER GWTUNIT
   CHARACTER*20 SITE,MSITE
   !
   !     ******************************************************************
   !
   CALL SGWF2MNW2IPNT(IGRID)
   !
   IF(MNWOBS.EQ.0) THEN
      RETURN
   ENDIF
   IF(MNWOBS.EQ.1) THEN
      WRITE (IOUT,120) MNWOBS
   ELSEIF(MNWOBS.GT.1) THEN
      WRITE (IOUT,140) MNWOBS
   ELSEIF(MNWOBS.LT.1) THEN
      RETURN
   END IF
   WRITE (IOUT,150)
   IF(MNWOBS.GT.MNWMAX) then
      write(iout,*) '***ERROR*** MNWOBS > MNWMAX'
      STOP 'MNWI ERROR'
   end if
   !
   !  Initialize data array
   MNWILST=0.0
   ! READ THE FIRST RECORD
   IOB=1
   IS_SITE=0
   !
   ! MNWILST(1,IOB) is Well # in MNW list
   ! MNWILST(2,IOB) is net volume in/out well
   ! MNWILST(3,IOB) is unit number for output
   ! MNWILST(4,IOB4) is QNDflag
   ! MNWILST(5,IOB) is QBHflag
   ! MNWILST(6,IOB) is CONCflag
   !
   if(GWTUNIT.GT.0) then
      READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),&
      &MNWILST(5,IOB),MNWILST(6,IOB)
   else
      READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),&
      &MNWILST(5,IOB)
      MNWILST(6,IOB)=0
   end if
   SITE=MNWIID(IOB)
   call UPCASE(SITE)
   ! check site vs list of site names in MNWSITE
   ! Loop over all MNW locations
   !   Loop over all wells
   do iw=1,MNWMAX
      MSITE=WELLID(iw)
      call UPCASE(MSITE)
      IF(SITE.EQ.MSITE) THEN
         IS_SITE=1
         MNWILST(1,IOB)=iw
      END IF
   end do
   !
   WRITE(IOUT,'(I8,3X,A12,3I8)') IOB,SITE,INT(MNWILST(3,IOB)),&
   &INT(MNWILST(4,IOB)),INT(MNWILST(5,IOB))
   IF(IS_SITE.EQ.0) THEN
      WRITE(IOUT,*) '***ERROR***   SITE FOR MNWI ',&
      &'WELL DESIGNATED FOR OBSERVATION NOT FOUND'
      STOP 'MNWI ERROR'
   ENDIF
   ! CYCLE THROUGH THE REMAINING RECORDS
   DO IOB=2,MNWOBS
      IS_SITE=0
      if(GWTUNIT.GT.0) then
         READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),&
         &MNWILST(5,IOB),MNWILST(6,IOB)
      else
         READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),&
         &MNWILST(5,IOB)
         MNWILST(6,IOB)=0
      end if
   ! check site vs list of site names in WELLID
      SITE=MNWIID(IOB)
      call UPCASE(SITE)
   !   Loop over all wells
      do iw=1,MNWMAX
         MSITE=WELLID(iw)
         call UPCASE(MSITE)
         IF(SITE.EQ.MSITE) THEN
            IS_SITE=1
            MNWILST(1,IOB)=iw
         END IF
      end do
   !
      WRITE(IOUT,'(I8,3X,A12,3I8)') IOB,SITE,INT(MNWILST(3,IOB)),&
      &INT(MNWILST(4,IOB)),INT(MNWILST(5,IOB))
      IF(IS_SITE.EQ.0) THEN
         WRITE(IOUT,*) '***ERROR***   SITE FOR MNWI ',&
         &'WELL DESIGNATED FOR OBSERVATION NOT FOUND'
         STOP 'MNWI ERROR'
      ENDIF
   !
   END DO
   WRITE(IOUT,'(140A)') 'DATA FOR MNW WELLS DESIGNATED FOR&
   & OBSERVATION WILL BE WRITTEN ON UNIT NUMBERS LISTED ABOVE'
   WRITE(IOUT,'(/)')
120 FORMAT(///'SITE ID FOR',I4,&
   &' MNW2 WELL DESIGNATED FOR OBSERVATION:')
140 FORMAT(///'SITE IDS FOR',I4,&
   &' MNW2 WELLS DESIGNATED FOR OBSERVATION:')
150 FORMAT(/'  WELL #   SITE ID         UNIT  QNDflag QBHflag')
   RETURN
END
   !
   !_________________________________________________________________________________
   !
SUBROUTINE GWF2MNW2I7OT(nstp,kkstp,kkper,IGRID)
   !     VERSION 20070923 GZH
   !
   !     ******************************************************************
   !    Sort well output into useful tables
   !     ******************************************************************
   !
   !        specifications:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:IOUT,ncol,nrow,nlay,hnew
   USE GWFBASMODULE, ONLY:HDRY,DELT,TOTIM
   USE GWFMNW2MODULE, ONLY:MNWMAX,NMNWVL,MNWAUX,MNW2,WELLID,NODTOT,&
   &NTOTNOD,MNWNOD,SMALL
   USE GWFMNW2IMODULE, ONLY:Wel1flag,QSUMflag,BYNDflag,MNWOBS,&
   &MNWILST,MNWIID
   ALLOCATABLE QBH(:)
   INTEGER firstnode,lastnode,QNDflag,QBHflag,QCONCflag,&
   &iaux,naux
   DOUBLE PRECISION q,hwell,qin,qout,qnet,hcell,&
   &QBH
   !--LFK  Nov. 2012
   DOUBLE PRECISION COND,SEEPFLG
   CHARACTER*20 obssite
   CHARACTER*50 LFRMAT
   !--lfk
   CHARACTER*50 dtext
   !
   !------------------------------------------------------------------
   !-lfk  10/10/2012
   dtext = 'Hwell < value shown because of seepage face calc. '
   !
   CALL SGWF2MNW2PNT(IGRID)
   CALL SGWF2MNW2IPNT(IGRID)
   !
   ALLOCATE(QBH(NODTOT),STAT=ISTAT)
   IF (ISTAT.NE.0) THEN
      WRITE(*,1700)ISTAT
1700  FORMAT(1X,'ALLOCATION OF MNW ROUTING ARRAY FAILED,',&
      &' RETURNED ERROR MESSAGE NUMBER: ',I6)
      CALL USTOP(' ')
   ENDIF
   !
   ! Print WEL1 file
   if(Wel1flag.gt.0) then
   ! print max number of wells, set IWELCB=0
   !-LFK     CHECK IF AUXILIARY VARIABLES PRESENT, & WRITE THEM IF PRESENT
      NAUX=NMNWVL-30
      IF (NAUX.LE.0) THEN
         if(kkper.eq.1.and.kkstp.eq.1)&
         &write(Wel1flag,'(2i10)') NODTOT,0
      ELSE
         if(kkper.eq.1.and.kkstp.eq.1)&
         &write(Wel1flag,1750) NODTOT,0,(mnwaux(iaux),iaux=1,naux)
1750     FORMAT(2I10,2x,5(:(' AUX ',16A,2X)))
      END IF
   ! only write at end of stress period (nstp=kkstp)
      if(nstp.eq.kkstp) then
   ! write number of wells
         write(Wel1flag,'(i10)') ntotnod
   !   Loop over all wells
         do iw=1,MNWMAX
   !   Loop over nodes in well
            firstnode=MNW2(4,iw)
            lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
            do INODE=firstnode,lastnode
               il=MNWNOD(1,INODE)
               ir=MNWNOD(2,INODE)
               ic=MNWNOD(3,INODE)
               q=MNWNOD(4,INODE)
   ! LFK  check if well is inactive; if yes, set q at all nodes = 0 to assure consistent WEL1 file.
               if(MNW2(1,iw).eq.0) q=0.0
   !   Write L,R,C,q (& OPTIONAL AUX values)
               IF (NAUX.LE.0) THEN
   !  LFK 3/26/12                write(Wel1flag,'(i9,2i10,1x,g15.6)') il,ir,ic,q
                  write(Wel1flag,'(i10,2i10,1x,g14.6)') il,ir,ic,q
               ELSE
                  write(Wel1flag,1760)il,ir,ic,q,(MNW2(30+IAUX,IW),IAUX=1,NAUX)
   ! LFK 1760 FORMAT (i9,2i10,1x,g15.6,2X,5(:(f4.0,4x)))
1760              FORMAT (i10,2i10,1x,g14.6,2X,5(:(f4.0,4x)))
               END IF
            end do
         end do
      end if
   end if
   !  Print QSUM file
   if(QSUMflag.gt.0) then
   !   Write header
      if(kkstp.eq.1.and.kkper.eq.1) then
         write(QSUMflag,'(200A)') 'WELLID                    Totim&
         &        Qin           Qout           Qnet          hwell'
      end if
   !   Loop over all wells
      do iw=1,MNWMAX
         qin=0.0D0
         qout=0.0D0
         qnet=0.0D0
   !   Only operate on active wells (MNW2(1,iw)=1)
         if (MNW2(1,iw).EQ.1) then
   !   Loop over nodes in well
   !--lfk
            NNDSIW=ABS(MNW2(2,iw))
   !-lfk
            firstnode=MNW2(4,iw)
   !            lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
            lastnode=MNW2(4,iw)+NNDSIW-1
            hwell=MNW2(17,iw)
            do INODE=firstnode,lastnode
               il=MNWNOD(1,INODE)
               ir=MNWNOD(2,INODE)
               ic=MNWNOD(3,INODE)
               q=MNWNOD(4,INODE)
   !--LFK
               HCELL=hnew(ic,ir,il)
               SEEPFLG=MNWNOD(15,INODE)
               COND=MNWNOD(14,INODE)
   !-LFK
               if(q.lt.0.0D0) then
                  qin=qin+q
               else
                  qout=qout+q
               end if
               qnet=qnet+q
            end do
   !            if(qnet.lt.(small*qin)) qnet=0.d0
   !--lfk  Nov. 2012
            if(SEEPFLG.EQ.hwell.or.&
            &SEEPFLG.eq.Hdry) then
               write(QSUMflag,'(A20,5(1x,1Pg14.6))')&
               &WELLID(iw),totim,qin,qout,qnet,hwell
            else
   !   If seepage face in cell, MNWNOD(15) will hold the bottom elev of the cell
   !--LFK  update Hwell to realistic value (not Hnew)for single-node MNW wells
               if (NNDSIW.eq.1) then
                  hwell=HCELL+(q/COND)
                  write(QSUMflag,'(A20,5(1x,1Pg14.6),3x,A50)')&
                  &WELLID(iw),totim,qin,qout,qnet,hwell,dtext
               else
                  write(QSUMflag,'(A20,5(1x,1Pg14.6))')&
                  &WELLID(iw),totim,qin,qout,qnet,hwell
               end if
            end if
   !--LFK
         end if
      end do
   end if
   !   Print BYND (ByNode) file
   if(BYNDflag.gt.0) then
   !   Write header
      if(kkstp.eq.1.and.kkper.eq.1) then
         write(BYNDflag,'(101A)') 'WELLID                NODE   Lay&
         &  Row   Col        Totim        Q-node         hwell         hcell&
         &   Seepage_elev.'
      end if
   !   Loop over all wells
      do iw=1,MNWMAX
   !   Only operate on active wells (MNW2(1,iw)=1)
         if (MNW2(1,iw).EQ.1) then
   !   Loop over nodes in well
            firstnode=MNW2(4,iw)
            lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
            hwell=MNW2(17,iw)
            do INODE=firstnode,lastnode
               il=MNWNOD(1,INODE)
               ir=MNWNOD(2,INODE)
               ic=MNWNOD(3,INODE)
               q=MNWNOD(4,INODE)
               hcell=hnew(ic,ir,il)
               nd=INODE-firstnode+1
   !   If no seepage face in cell, don't print seepage elev.
               if(MNWNOD(15,INODE).EQ.hwell.or.&
   !--lfk  Nov. 2012
               &MNWNOD(15,INODE).eq.Hdry) then
   !--lfk     &           MNWNOD(15,INODE).eq.Hdry.OR.MNW2(2,iw).eq.1) then
                  write(BYNDflag,'(A20,4i6,1x,1P4e14.6)')&
                  &WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell
               else
   !   If seepage face in cell, MNWNOD(15) will hold the bottom elev of the cell,
   !   which is used with hcell to get the gradient used to calculate the Q for the
   !   seepage face.
   !--LFK  update Hwell to realistic value (not Hnew)for single-node MNW wells
                  if (firstnode.eq.lastnode) then
                     hwell=hcell+(q/MNWNOD(14,INODE))
                     write(BYNDflag,'(A20,4i6,1x,1P5e14.6,3x,A50)')&
                     &WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,&
                     &MNWNOD(15,INODE),dtext
                  else
                     write(BYNDflag,'(A20,4i6,1x,1P5e14.6)')&
                     &WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,&
                     &MNWNOD(15,INODE)
                  end if
   !--LFK
               end if
            end do
         end if
      end do
   end if
   !
   !  Print MNWOBS files
   do iwobs=1,MNWOBS
      qnet=0.d0
      obssite=MNWIID(iwobs)
      iw=MNWILST(1,iwobs)
   !   Loop over nodes in well
   !--lfk  Nov. 2012
      NNDSIW=ABS(MNW2(2,iw))
   !-lfk
      firstnode=MNW2(4,iw)
   !            lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
      lastnode=MNW2(4,iw)+NNDSIW-1
      hwell=MNW2(17,iw)
      qin=0.D0
      qout=0.D0
      do INODE=firstnode,lastnode
   !--LFK
         if (MNW2(1,iw).EQ.1) then
            il=MNWNOD(1,INODE)
            ir=MNWNOD(2,INODE)
            ic=MNWNOD(3,INODE)
            HCELL=hnew(ic,ir,il)
            SEEPFLG=MNWNOD(15,INODE)
            COND=MNWNOD(14,INODE)
         end if
   !-LFK
         q=MNWNOD(4,INODE)
         if(q.lt.0.0D0) then
            qin=qin+q
         else
            qout=qout+q
         end if
         qnet=qnet+q
      end do
   !--lfk  Nov. 2012
      if (MNW2(1,iw).EQ.1) then
         sftest=0.0
         if(SEEPFLG.ne.hwell.and.SEEPFLG.ne.Hdry) then
   !   If seepage face in cell, MNWNOD(15) [=seepflg] will hold the bottom elev of the cell
   !--LFK  update Hwell to realistic value (not Hnew)for single-node MNW wells
            if (NNDSIW.eq.1) then
               hwell=HCELL+(q/COND)
               sftest=1.0
            end if
         end if
      end if
   !--LFK
   !   Cumulative volume for this well
      MNWILST(2,iwobs)=MNWILST(2,iwobs)+qnet*DELT

   ! get NNODES
      NNODES=INT(ABS(MNW2(2,iw)))
   !
   !  Print according to flags
      QNDflag=INT(MNWILST(4,iwobs))
      QBHflag=INT(MNWILST(5,iwobs))
      QCONCflag=INT(MNWILST(6,iwobs))
   !--LFK(Nov.2016)--make sure QNDflag and QBHflag = 0 when NNODES=1
   !         (as described on p. 53 of model documentation).
      IF(NNODES.EQ.1) THEN
         IF(QNDflag.GT.0) then
            QNDflag=0
            MNWILST(4,iwobs)=0
            write(iout,30) WELLID(iwobs)
30          FORMAT('**note** NNODES=1 IN WELL: ',A20,'; QNDflag reset to 0.')
         END IF
         IF(QBHflag.GT.0) then
            QBHflag=0
            MNWILST(5,iwobs)=0
            write(iout,32) WELLID(iwobs)
32          FORMAT('**note** NNODES=1 IN WELL: ',A20,'; QBHflag reset to 0.')
         END IF
   !
      END IF
   !
      if(QBHflag.gt.0)&
      &call GWF2MNW27BH(iw,IGRID)
   !
   !  Create format for header--single node
   !  Ignore ND and BH flags
   !  Create format for header--single node
   !--lfk/Nov.2012:  eliminate separate processing for single-node MNW wells
   !          if(NNODES.eq.1) then
   !            if(kkstp.eq.1.and.kkper.eq.1) then
   !          Write (INT(MNWILST(3,iwobs)),'(A)') 'WELLID
   !     &       TOTIM           Qin
   !     &          Qout          Qnet         QCumu         hwell'
   !            end if
   !            write(INT(MNWILST(3,iwobs)),'(A20,1x,1P6e14.6)')
   !     &             WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell
   !          else
   !  all multi-node well output below
      if(QNDflag.eq.0) then
         if(QBHflag.eq.0) then
            if(QCONCflag.eq.0) then
   !  QNDflag=0, QBHflag=0, QCONCflag=0
   ! write header
               if(kkstp.eq.1.and.kkper.eq.1) then
                  Write (INT(MNWILST(3,iwobs)),'(120A)') 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet      Cum.Vol.         hwell&
                  &   '
               end if
   !   Only operate on active wells (MNW2(1,iw)=1)
               if (MNW2(1,iw).EQ.1) then
   !--lfk
                  if (sftest.lt.1.0) then
                     write(INT(MNWILST(3,iwobs)),'(A20,1x,1P6e14.6)')&
                     &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell
                  else
                     write(INT(MNWILST(3,iwobs)),'(A20,1x,1P6e14.6,3x,A50)')&
                     &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell&
                     &,dtext
                  end if
               else
   !-lfk-Nov 2012
                  write(INT(MNWILST(3,iwobs)),&
                  &'(A20,1x,1Pe14.6,"   (Well is inactive)")')&
                  &WELLID(iw),totim
               end if
            else
   !  QNDflag=0, QBHflag=0, QCONCflag=1
   ! write header
               if(kkstp.eq.1.and.kkper.eq.1) then
                  Write (INT(MNWILST(3,iwobs)),'(120A)') 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet         QCumu         hwell&
                  &   CONC'
               end if
   !gzh debug  concflag not coded yet...separate routine for when GWT active?
               if (MNW2(1,iw).EQ.1) then
                  write(INT(MNWILST(3,iwobs)),'(A20,1x,1P6e14.6,3A)')&
                  &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,&
                  &'N/A'
   !
               end if
            end if
         else
            if(QCONCflag.eq.0) then
   !  QNDflag=0, QBHflag=1, QCONCflag=0
   ! write "smart" header
               if(kkstp.eq.1.and.kkper.eq.1) then
   !  Create format for header
                  WRITE(LFRMAT,14) NNODES-1
14                FORMAT('(A,',I4,'I12)')
   !  WRITE FORMAT FOR HEADER LINE
                  Write (INT(MNWILST(3,iwobs)),LFRMAT) 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet         QCumu         hwell&
                  & QBH_seg-->1',(ibh,ibh=2,NNODES)
               end if
               if (MNW2(1,iw).EQ.1) then

   ! Create format for write
                  WRITE (LFRMAT,355) NNODES
355               FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))')
                  write(INT(MNWILST(3,iwobs)),LFRMAT)&
                  &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,&
                  &(MNWNOD(27,i),i=firstnode,lastnode)
               end if
            else
   !  QNDflag=0, QBHflag=1, QCONCflag=1
   ! write header
               if(kkstp.eq.1.and.kkper.eq.1) then
   !  Create format for header
                  WRITE(LFRMAT,15) NNODES-1
15                FORMAT('(A,',I4,'I12,A)')
   !  WRITE FORMAT FOR HEADER LINE
                  Write (INT(MNWILST(3,iwobs)),LFRMAT) 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet         QCumu         hwell&
                  & QBH_seg-->1',(ibh,ibh=2,NNODES),'     CONC'
               end if
               if (MNW2(1,iw).EQ.1) then
   ! Create format for write
                  WRITE (LFRMAT,255) NNODES
255               FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))',3A)
                  write(INT(MNWILST(3,iwobs)),LFRMAT)&
                  &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,&
                  &(MNWNOD(27,i),i=firstnode,lastnode),'N/A'
               end if
            end if
         end if
      else
         if(QBHflag.eq.0) then
            if(QCONCflag.eq.0) then
   !  QNDflag=1, QBHflag=0, QCONCflag=0
   ! write "smart" header
               if(kkstp.eq.1.and.kkper.eq.1) then
   !  Create format for header
                  WRITE(LFRMAT,14) NNODES-1
                  Write (INT(MNWILST(3,iwobs)),LFRMAT) 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet         QCumu         hwell&
                  & Flow@Nd-->1',(ind,ind=2,NNODES)
               end if
               if (MNW2(1,iw).EQ.1) then
   ! Create format for write
                  WRITE (LFRMAT,155) NNODES
155               FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))')

                  write(INT(MNWILST(3,iwobs)),LFRMAT)&
                  &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,&
                  &(MNWNOD(4,i),i=firstnode,lastnode)
               end if
            else
   !  QNDflag=1, QBHflag=0, QCONCflag=1
   !
   ! write "smart" header
               if(kkstp.eq.1.and.kkper.eq.1) then
   !  Create format for header
                  WRITE(LFRMAT,15) NNODES-1
                  Write (INT(MNWILST(3,iwobs)),LFRMAT) 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet         QCumu         hwell&
                  & Flow@Nd-->1',(ind,ind=2,NNODES),'   CONC'
               end if
               if (MNW2(1,iw).EQ.1) then
   ! Create format for write
                  WRITE (LFRMAT,455) NNODES
455               FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4)),3A')

                  write(INT(MNWILST(3,iwobs)),LFRMAT)&
                  &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,&
                  &(MNWNOD(4,i),i=firstnode,lastnode),'N/A'
   !
               end if
            end if
         else
            if(QCONCflag.eq.0) then
   !  QNDflag=1, QBHflag=1, QCONCflag=0
   ! Create format for write

   ! write "smart" header
               if(kkstp.eq.1.and.kkper.eq.1) then
   !  Create format for header
                  WRITE(LFRMAT,16) NNODES-1,NNODES-1
16                FORMAT('(A,',I4,'I12,A,',I4,'I12)')
                  Write (INT(MNWILST(3,iwobs)),LFRMAT) 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet         QCumu         hwell&
                  & Flow@Nd-->1',(ind,ind=2,NNODES),&
                  &' QBH_seg-->1',(ibh,ibh=2,NNODES)
               end if
               if (MNW2(1,iw).EQ.1) then
                  numvals=NNODES+NNODES
                  WRITE (LFRMAT,156) numvals
156               FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))')
                  write(INT(MNWILST(3,iwobs)),LFRMAT)&
                  &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,&
                  &(MNWNOD(4,i),i=firstnode,lastnode),&
                  &(MNWNOD(27,i),i=firstnode,lastnode)

               end if
            else
   !  QNDflag=1, QBHflag=1, QCONCflag=1
               if(kkstp.eq.1.and.kkper.eq.1) then
                  WRITE(LFRMAT,17) NNODES-1,NNODES-1
17                FORMAT('(A,',I4,'I12,A,',I4,'I12,A)')
                  Write (INT(MNWILST(3,iwobs)),LFRMAT) 'WELLID&
                  &       TOTIM           Qin&
                  &          Qout          Qnet         QCumu         hwell&
                  & Flow@Nd-->1',(ind,ind=2,NNODES),&
                  &' QBH_seg-->1',(ibh,ibh=2,NNODES),'   CONC'
               end if
               if (MNW2(1,iw).EQ.1) then
                  WRITE (LFRMAT,456) NNODES,&
                  &NNODES
456               FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))',I4,'(1pE12.4))')
                  write(INT(MNWILST(3,iwobs)),LFRMAT)&
                  &WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,&
                  &(MNWNOD(4,i),i=firstnode,lastnode),&
                  &(MNWNOD(27,i),i=firstnode,lastnode),'N/A'
               end if
            end if
         end if
      end if
   !-lfk-Nov. 2012          end if
   end do
   !
   return
end
   !
SUBROUTINE GWF2MNW2I7DA(IGRID)
   !  Deallocate MNW MEMORY
   USE GWFMNW2IMODULE
   !
   CALL SGWF2MNW2IPNT(IGRID)
   DEALLOCATE(Wel1flag)
   DEALLOCATE(QSUMflag)
   DEALLOCATE(BYNDflag)
   DEALLOCATE(MNWOBS)
   DEALLOCATE(MNWIID)
   DEALLOCATE(MNWILST)
   !
   RETURN
END
SUBROUTINE SGWF2MNW2IPNT(IGRID)
   !  Change MNW data to a different grid.
   USE GWFMNW2IMODULE
   !
   Wel1flag=>GWFMNWIDAT(IGRID)%Wel1flag
   QSUMflag=>GWFMNWIDAT(IGRID)%QSUMflag
   BYNDflag=>GWFMNWIDAT(IGRID)%BYNDflag
   MNWOBS=>GWFMNWIDAT(IGRID)%MNWOBS
   MNWIID=>GWFMNWIDAT(IGRID)%MNWIID
   MNWILST=>GWFMNWIDAT(IGRID)%MNWILST
   !
   RETURN
END
SUBROUTINE SGWF2MNW2IPSV(IGRID)
   !  Save MNW2 data for a grid.
   USE GWFMNW2IMODULE
   !
   GWFMNWIDAT(IGRID)%Wel1flag=>Wel1flag
   GWFMNWIDAT(IGRID)%QSUMflag=>QSUMflag
   GWFMNWIDAT(IGRID)%BYNDflag=>BYNDflag
   GWFMNWIDAT(IGRID)%MNWOBS=>MNWOBS
   GWFMNWIDAT(IGRID)%MNWIID=>MNWIID
   GWFMNWIDAT(IGRID)%MNWILST=>MNWILST
   !
   RETURN
END
