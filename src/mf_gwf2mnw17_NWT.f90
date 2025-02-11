   !                  KJH  20030327      -- Patched Hyd.K term in LPF option -- cel2wel function
   !                  KJH  20030717      -- Patched budget output switch -- subroutine GWF1MNW1bd
   !                                        Cleaned output so outrageous pointers are not printed
   !                  GZH  20050405      -- Converted calculations to use double precision
   !                  KJH  20050419      -- Array WELL2 dimensioned to 18 to store well id
   !                  AWH  20080411      -- Retrieve HDRY from GWFBASMODULE rather than from
   !                                        LPF, BCF, or HUF
   !                  RGN  20111001      -- Modified to support UPW; removed BC switching is UPW is used.
   !
MODULE GWFMNW1MODULE
   DOUBLE PRECISION, PARAMETER :: TWOPI=2.0D0*3.1415926535897932D0
   DOUBLE PRECISION, PARAMETER :: ZERO25=1.0D-25, ZERO20=1.0D-20
   DOUBLE PRECISION, PARAMETER :: ZERO8=1.0D-8, BIG=1.0D30
   CHARACTER(LEN=200),SAVE,POINTER:: MNWNAME
   INTEGER,          SAVE,POINTER :: NWELL2, MXWEL2, IWL2CB, KSPREF
   INTEGER,          SAVE,POINTER :: IWELPT, NOMOITER
   DOUBLE PRECISION, SAVE,POINTER :: PLOSS
   DOUBLE PRECISION, SAVE,POINTER :: SMALL, HMAX
   CHARACTER(LEN=32),SAVE,DIMENSION(:),    POINTER :: MNWSITE
   INTEGER,          SAVE,DIMENSION(:),    POINTER :: IOWELL2
   DOUBLE PRECISION, SAVE,DIMENSION(:,:),  POINTER :: WELL2
   DOUBLE PRECISION, SAVE,DIMENSION(:,:,:),POINTER :: HREF
   TYPE GWFMNWTYPE
      CHARACTER(LEN=200),    POINTER :: MNWNAME
      INTEGER,               POINTER :: NWELL2, MXWEL2, IWL2CB, KSPREF
      INTEGER,               POINTER :: IWELPT, NOMOITER
      DOUBLE PRECISION,      POINTER :: PLOSS
      DOUBLE PRECISION,      POINTER :: SMALL, HMAX
      CHARACTER(LEN=32),     DIMENSION(:),    POINTER :: MNWSITE
      INTEGER,               DIMENSION(:),    POINTER :: IOWELL2
      DOUBLE PRECISION,      DIMENSION(:,:),  POINTER :: WELL2
      DOUBLE PRECISION,      DIMENSION(:,:,:),POINTER :: HREF
   END TYPE
   TYPE(GWFMNWTYPE), SAVE:: GWFMNWDAT(10)
END MODULE GWFMNW1MODULE
   !
   !-------------------------------------------------------------------------
   !
SUBROUTINE GWF2MNW17AR(In, Iusip, Iude4, Iunwt, Iusor, Iupcg,&
&Iulmg, Iugmg, Fname, Igrid)
   !rgn------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
   !rgn------NEW VERSION NUMBER 1.3.0, 7/01/2022
   !
   !----- MNW by K.J. Halford        1/31/98
   !     ******************************************************************
   !     allocate array storage for well package
   !     ******************************************************************
   !
   !        specifications:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY: IOUT,NCOL,NROW,NLAY
   USE GWFMNW1MODULE
   USE SIPMODULE,ONLY:HCLOSE
   USE DE4MODULE,ONLY:HCLOSEDE4
   USE PCGMODULE,ONLY:HCLOSEPCG
   USE GMGMODULE,ONLY:HCLOSEGMG
   USE GWFNWTMODULE,ONLY:Tol
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   INTRINSIC ABS
   INTEGER, EXTERNAL :: IFRL
   EXTERNAL NCREAD, UPCASE, QREAD, USTOP
   !     ------------------------------------------------------------------
   !     Arguments
   !     ------------------------------------------------------------------
   INTEGER :: In, Iusip, Iude4, Iusor, Iupcg, Iulmg, Iugmg,&
   &Iunwt, Igrid
   CHARACTER(LEN=200) :: Fname                 !!08/19/02KJH-MODIFIED
   !     ------------------------------------------------------------------
   !     Local Variables
   !     ------------------------------------------------------------------
   REAL :: bs
   INTEGER :: ierr, io, iok, jf, ke, kf, ki, kio
   DOUBLE PRECISION :: rn(25)
   CHARACTER(LEN=256) :: txt, tx2
   !     ------------------------------------------------------------------
   !     Static Variables
   !     ------------------------------------------------------------------
   CHARACTER(LEN=6) :: ftag(3)
   INTEGER :: icf(3)
   DATA ftag/'WEL1  ', 'BYNODE', 'QSUM  '/
   DATA icf/4, 6, 4/
   !     ------------------------------------------------------------------
   ALLOCATE (MNWNAME, NWELL2, MXWEL2, IWL2CB, NOMOITER, KSPREF,&
   &IWELPT)
   ALLOCATE (PLOSS, SMALL, HMAX, IOWELL2(3),&
   &HREF(NCOL,NROW,NLAY))
   !
   IOWELL2(1) = 0
   IOWELL2(2) = 0
   IOWELL2(3) = 0
   !
   !1------identify package and initialize nwell2
   WRITE (IOUT, 9001) In
9001 FORMAT (/, ' MNW1 -- MULTI-NODE WELL 1 PACKAGE, VERSION 7,',&
   &' 04/05/2012.', /, '    INPUT READ FROM UNIT', i4)
   NWELL2 = 0
   !
   !2------read max number of wells and
   !2------unit or flag for cell-by-cell flow terms.
   CALL NCREAD(In, txt, ierr)
   CALL UPCASE(txt)
   !
   ki = INDEX(txt, 'REF')
   IF ( ki.GT.0 ) THEN
      tx2 = txt(ki:256)
      CALL QREAD(rn, 1, tx2, ierr)
      IF ( ierr.EQ.0 ) KSPREF = IFRL(rn(1))
      txt(ki:256) = '                                '
   ELSE
      KSPREF = 1
   ENDIF
   !
   CALL QREAD(rn, 4, txt, ierr)
   MXWEL2 = IFRL(rn(1))
   IWL2CB = 0
   IF ( ierr.LE.2 ) IWL2CB = IFRL(rn(2))
   IWELPT = 0
   IF ( ierr.EQ.1 ) IWELPT = IFRL(rn(3))
   NOMOITER = 9999
   IF ( ierr.EQ.0 ) NOMOITER = IFRL(rn(4))
   !
   WRITE (IOUT, 9002) MXWEL2
   IF ( IWL2CB.GT.0 ) WRITE (IOUT, 9003) IWL2CB
   IF ( IWL2CB.LT.0 ) WRITE (IOUT, 9004)
   WRITE (IOUT, 9005) KSPREF
   WRITE (IOUT, 9006) NOMOITER
9002 FORMAT (' MAXIMUM OF', i7, ' WELLS')
9003 FORMAT (' CELL-BY-CELL FLOWS WILL BE RECORDED ON UNIT', i3)
9004 FORMAT (' CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
9005 FORMAT ('  The heads at the beginning of SP:', i4,&
   &' will be the default reference elevations.', /)
9006 FORMAT (' Flow rates will not be estimated after the', i4,&
   &'th iteration')
   !
   !   Define well model to be used
   !
   CALL NCREAD(In, txt, ierr)
   CALL UPCASE(txt)
   PLOSS = 0.0D0   !!  Default use of Skin so linear loss varies with T
   IF ( INDEX(txt, 'LINEAR').GT.0 ) THEN
      PLOSS = 1.0D0 !!  ADD THIS LINE to make sure that the power term is 1 for the linear model
      ki = INDEX(txt, ':') + 1
      tx2 = txt(ki:256)
      CALL QREAD(rn, 1, tx2, ierr)
      IF ( ierr.EQ.0 ) PLOSS = rn(1)
   !   Add error checking to shut down MODFLOW
      bs = 3.6           !!   Maximum limit on power term
      IF ( PLOSS.GT.bs ) THEN
         WRITE (*, *) 'Power term of', PLOSS, ' exceeds maximum of', bs
         WRITE (IOUT, *) 'Power term of', PLOSS, ' exceeds maximum of',&
         &bs
   !
   !         When compiling MNW with Modflow-96, comment out the call to
   !         USTOP and uncomment the STOP statement
         CALL USTOP(' ')
   !         STOP
   !
      ENDIF
   !
   ENDIF
   !
   !   Test for a specified PREFIX NAME  for time series output from MNW7OT
   !
   CALL NCREAD(In, txt, ierr)
   tx2 = txt
   CALL UPCASE(tx2)
   kf = INDEX(tx2, 'PREFIX:')
   IF ( kf.GT.0 ) THEN
      MNWNAME = txt(kf+7:256)
      ke = INDEX(MNWNAME, ' ')
      MNWNAME(ke:200) = '               '
      tx2 = MNWNAME
      CALL UPCASE(tx2)
      IF ( INDEX(tx2, 'FILEPREFIX').GT.0 ) THEN
         MNWNAME = Fname
         ke = INDEX(MNWNAME, '.')
         MNWNAME(ke:200) = '               '
      ENDIF
   ELSE
      MNWNAME = 'OUTput_MNW'
      BACKSPACE (In)
   ENDIF
   !
   !     Test for creation of a WEL1 package and auxillary output files
   !
   iok = 1
   DO WHILE ( iok.EQ.1 )
      CALL NCREAD(In, txt, ierr)
      tx2 = txt
      CALL UPCASE(tx2)
      kf = INDEX(tx2, 'FILE:')
      IF ( kf.GT.0 ) THEN
         kio = 0
         jf = 0
         DO WHILE ( kio.EQ.0 .AND. jf.LT.3 )
            jf = jf + 1
            kio = INDEX(tx2, ftag(jf)(1:icf(jf)))
            IF ( kio.GT.0 ) THEN
               tx2 = txt(kio+1+icf(jf):256)
               CALL QREAD(rn, 1, tx2, ierr)
               IF ( ierr.EQ.0 ) THEN
                  IOWELL2(jf) = IFRL(rn(1))
   !            OC over ride is ALLTIME
                  IF ( INDEX(tx2, 'ALLTIME').GT.0 ) IOWELL2(jf)&
                  &= -IOWELL2(jf)
   !            Find and use file name
                  tx2 = txt(kf+5:256)
                  kf = INDEX(tx2, ' ') - 1
                  CLOSE (ABS(IOWELL2(jf)))
                  OPEN (ABS(IOWELL2(jf)), FILE=tx2(1:kf))
                  WRITE (tx2(253:256), '(i4)') ABS(IOWELL2(jf))
                  txt = ' A '//ftag(jf)&
                  &//' data input file will be written'//' to '//&
                  &tx2(1:kf)//' on unit '//tx2(253:256)
                  WRITE (IOUT, '(/1x,a79)') txt
                  IF ( jf.EQ.1 )&
                  &WRITE (ABS(IOWELL2(jf)),'(3i10)') MXWEL2, IWL2CB, IWELPT
               ENDIF
            ENDIF
         ENDDO
      ELSE
         BACKSPACE (In)
         iok = 0
      ENDIF
   ENDDO
   !
   !  Write header in Auxillary BYNODE file if KPER=1 & IO>0
   !
   IF ( IOWELL2(2).NE.0 ) THEN
      io = ABS(IOWELL2(2))
      WRITE (io, 9008)
   ENDIF
   !
   !  Write header in Auxillary QSUM file if KPER=1 & IO>0
   !
   IF ( IOWELL2(3).NE.0 ) THEN
      io = ABS(IOWELL2(3))
      WRITE (io, 9009)
   ENDIF
   !
9008 FORMAT ('SiteID', 27x, 'Entry  NODE', 5x, 'Total_Time', 8x, 'Q',&
   &5x, 'H-Well', 5x, 'H-Cell', 5x, 'QW-Avg')
9009 FORMAT ('SiteID', 31x, 'Entry', 5x, 'Total_Time', 10x, 'Qin',&
   &10x, 'Qout', 10x, 'Qsum', 5x, 'H-Well', 5x, 'QW-Avg')
   !
   !  4/18/2005 - KJH:  Explicit well tracking addition changed 1st WELL2
   !                    dimension from 17 to 18
   ALLOCATE (WELL2(18, MXWEL2+1), MNWSITE(MXWEL2))
   !
   !-------SET SMALL DEPENDING ON CLOSURE CRITERIA OF THE SOLVER
   SMALL = 0.0D0
   IF ( Iusip.NE.0 ) SMALL = HCLOSE
   IF ( Iude4.NE.0 ) SMALL = HCLOSEDE4
   !     IF ( Iusor.NE.0 ) SMALL = HCLOSESOR
   IF ( Iupcg.NE.0 ) SMALL = HCLOSEPCG
   IF ( Iulmg.NE.0 ) SMALL = 0.0D0  !LMG SETS HCLOSE TO ZERO
   IF ( Iugmg.NE.0 ) SMALL = HCLOSEGMG
   IF ( Iunwt.NE.0 ) SMALL = TOL
   !
   !-----SAVE POINTERS FOR GRID AND RETURN
   CALL SGWF2MNW1PSV(Igrid)
   !
   !7------return
END SUBROUTINE GWF2MNW17AR
   !
   !_________________________________________________________________________________
   !
SUBROUTINE GWF2MNW17RP(In, Iubcf, Iulpf, Iuhuf, Iuupw, Kper,&
&Igrid)
   !     VERSION 20020819 KJH
   !
   !----- MNW by K.J. Halford        1/31/98
   !     ******************************************************************
   !     read new well locations, stress rates, conc, well char., and limits
   !     ******************************************************************
   !
   !        specifications:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NODES,NCOL,NROW,NLAY,IBOUND,HOLD,HNEW,IOUT
   USE GWFBASMODULE,ONLY:TOTIM,HDRY
   USE GWFMNW1MODULE, ONLY:NWELL2,MXWEL2,IWELPT,PLOSS,HMAX,&
   &MNWSITE,IOWELL2,WELL2,HREF,KSPREF,&
   &BIG,ZERO25
   IMPLICIT NONE
   INTRINSIC ABS, MAX, MOD, INT
   INTEGER, EXTERNAL :: IFRL, IDIRECT
   DOUBLE PRECISION, EXTERNAL :: CEL2WELBCF, CEL2WELLPF, CEL2WELHUF,&
   &CEL2WELUPW
   EXTERNAL NCREAD, UPCASE, QREAD, USTOP
   !     ------------------------------------------------------------------
   !     Arguments
   !     ------------------------------------------------------------------
   INTEGER, INTENT(IN) :: Iubcf, Iulpf, Iuhuf, Kper, Igrid, Iuupw
   INTEGER, INTENT(INOUT) :: In
   !     ------------------------------------------------------------------
   !     Local Variables
   !     ------------------------------------------------------------------
   DOUBLE PRECISION :: qfrcmn, qfrcmx, qreject, drytest, ipole, hlim
   DOUBLE PRECISION :: hrfw, qsum, rw, cond, qact, sk, cf, q, rn(25)
   INTEGER :: i, icmn, ierr, igrp, ii, iin, io, iok, ip, ipt, irmx
   INTEGER :: itmp, j, k, kblk, kcp, kfini, ki, kpc, kqc, ksiteid
   INTEGER :: ktab, m, mstep, n, n1, nb, ne, ngrp, nl, nn, node
   INTEGER :: nqreject, nstart
   INTEGER :: idwell,mm
   CHARACTER(LEN=1) :: tab
   CHARACTER(LEN=32) :: tempsite
   CHARACTER(LEN=256) :: txt, tx2, txtraw
   !     ------------------------------------------------------------------
   !-----SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2MNW1PNT(Igrid)
   !swm: SET POINTERS FOR FLOW PACKAGE TO GET K's FOR CEL2WEL
   IF ( Iubcf.NE.0 ) CALL SGWF2BCF7PNT(Igrid)
   IF ( Iulpf.NE.0 ) CALL SGWF2LPF7PNT(Igrid)
   IF ( Iuhuf.NE.0 ) CALL SGWF2HUF7PNT(Igrid)
   IF ( Iuupw.NE.0 ) CALL SGWF2UPW1PNT(Igrid)
   !
   icmn = 1
   kfini = 1
   tab = CHAR(9)
   qfrcmn = ZERO25
   qfrcmx = ZERO25
   qreject = 0.0D0
   nqreject = 0
   nl = 0
   IF ( PLOSS.GT.1.001D0 ) nl = 1 !!  Read NL loss Coefficient after Skin
   !
   !  Check for setting the HREFerence array
   !ERB     IN FIRST STRESS PERIOD, HOLD IS UNDEFINED, SO USE HNEW INSTEAD
   IF ( Kper.EQ.1 ) THEN
      HMAX = ABS(HNEW(1,1,1))
      DO k = 1, NLAY
         DO i = 1, NROW
            DO j = 1, NCOL
               HREF(j,i,k) = HNEW(j,i,k)
               HMAX = MAX(ABS(HREF(j,i,k)),HMAX)
            ENDDO
         ENDDO
      ENDDO
   ELSE IF ( Kper.LE.KSPREF ) THEN
      HMAX = ABS(HOLD(1,1,1))
      DO k = 1, NLAY
         DO i = 1, NROW
            DO j = 1, NCOL
               HREF(j,i,k) = HOLD(j,i,k)
               HMAX = MAX(ABS(HREF(j,i,k)),HMAX)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   !
   !------------------------------------------------------------------
   !     The 18 rows of the well array store:
   !      Row #  = Description
   !------------------------------------------------------------------
   !         1   = Well node locator
   !         2   = Desired flow rate
   !         3   = Actual flow rate used
   !         4   = Water Quality attribute to be averaged by flow
   !         5   = Radius of wellbore
   !         6   = Skin associated with well completion
   !         7   = Minimum/Maximum head or drawdown
   !         8   = Elevation of reference head for computing lift costs
   !         9   = Water Quality Group identifier
   !        10   = Water level in wellbore
   !        11   = HCOF value / QWaverage
   !        12   = RHS  value
   !        13   = Minimum flow rate -  to turn off
   !        14   = Minimum flow rate -- to turn on
   !        15   = Reserve Desired flow rate
   !        16   = Non-linear loss term
   !        17   = Actual flow rate to individual nodes of a multi-node well
   !               kept for transport or other purposes !!7/13/2003 - CZ
   !        18   = Explicit well identifier -- Same value for all nodes in a well
   !------------------------------------------------------------------
   !
   !1------read itmp(number of wells or flag saying reuse well data)
   CALL NCREAD(In, txtraw, ierr)
   txt = txtraw
   CALL UPCASE(txt)
   CALL QREAD(rn, 1, txt, ierr)
   itmp = rn(1)
   !
   IF ( itmp.LT.0 ) THEN
   !        if itmp less than zero reuse data. print message and return.
      WRITE (IOUT, 9001)
9001  FORMAT (1X,/1X,'REUSING MNW7  FROM LAST STRESS PERIOD')
      RETURN
   ELSE
   !  If itmp > 0,  Test if wells are to replace old ones or be added.
   !
      IF ( INDEX(txt, 'ADD').EQ.0 ) NWELL2 = 0
   !
   !   return if there are no wells to read ........
      IF ( itmp.EQ.0 ) RETURN
   !
   !  Redundant well information is allowed in MNW
   !
   !   Read additional well info
      nstart = NWELL2
      DO m = 1, itmp
         CALL NCREAD(In, txtraw, ierr)
         txt = txtraw
         CALL UPCASE(txt)
   !   Attempt read with QREAD first
         CALL QREAD(rn, 4, txt, ierr)
         IF ( ierr.EQ.0 .AND. rn(5).LT.0.5D0 ) THEN
            k = IFRL(rn(1))
            j = IFRL(rn(2))
            i = IFRL(rn(3))
            q = rn(4)
            irmx = IFRL(rn(6)) + 1
         ELSE
   !  Use fixed form reader if errors were detected
            READ (txt(1:40), '(3i10,f10.0)') k, j, i, q
            irmx = 41
         ENDIF
         node = (k-1)*NCOL*NROW + (j-1)*NCOL + i
   !    Test for if well is in active grid ......
         iok = 1
         IF ( i.GT.NCOL .OR. j.GT.NROW .OR. node.GT.NODES ) iok = 0
         drytest = HNEW(i,j,k) - HDRY
         IF (iok.GT.0 .AND. ABS(drytest).GT.ZERO25) iok = IBOUND(i,j,k)
   !
   !  Should MNW wells be allowed in specified-head cells?
         IF ( iok.NE.0 ) THEN      !! Allow SH now, "gt" for no SH
   !    Test for redundant info ......
            ipt = 0
   !    The commented statements prevent having multiple MNW sites in the same cells
   !            nt  = 0
   !            do while (nt.lt.nwell2 .and. ipt.eq.0 )
   !              nt = nt + 1
   !              if( well2(1,nt).eq.node ) ipt = nt
   !            enddo
            IF ( ipt.EQ.0 ) THEN
               NWELL2 = NWELL2 + 1
               ipt = NWELL2
            ENDIF
   !
   !    Assign data now that the pointer is set
            WELL2(1, ipt) = node
            WELL2(2, ipt) = q
            ipole = 0.0D0
            IF ( ABS(q).GT.ZERO25 ) ipole = q/ABS(q)
            WELL2(3, ipt) = WELL2(2, ipt)
            WELL2(13, ipt) = qfrcmn       ! default lower limit
            WELL2(14, ipt) = qfrcmx
   !
   !    Look for limit modifications
            kqc = INDEX(txt, 'QCUT')
            kpc = INDEX(txt, '%CUT')
            IF ( kqc+kpc.GT.0 .AND. ABS(q).GT.ZERO25 ) THEN
               tx2 = txt(kqc+kpc+5:256)
               CALL QREAD(rn, 2, tx2, ierr)
               IF ( kqc.GT.0 ) THEN         !!  Absolute value was provided
                  rn(1) = 100.0D0*rn(1)/q        !!  Convert to percentage
                  rn(2) = 100.0D0*rn(2)/q
               ENDIF
               IF ( ierr.GE.1 ) rn(2) = rn(1)
               WELL2(13, ipt) = rn(1)*0.01D0    !! convert percentages to fractions
               WELL2(14, ipt) = rn(2)*0.01D0
               IF ( INDEX(tx2, 'DEFAULT').GT.0 ) THEN
                  qfrcmn = rn(1)*0.01D0          !!  New default lower limit
                  qfrcmx = rn(2)*0.01D0          !!  New default upper limit
               ENDIF
            ENDIF
   !
   !    Look for NonLinear coefficient
            WELL2(16, ipt) = 0.0D0              !!  NonLinear Loss Coefficient
            kcp = INDEX(txt, 'CP:')
            IF ( kcp.GT.0 .AND. nl.GT.0 ) THEN
               tx2 = txt(kcp+3:256)
               CALL QREAD(rn, 1, tx2, ierr)
               IF ( ierr.EQ.0 ) THEN
                  WELL2(16, ipt) = rn(1)
   !         Could reset default C-term here to a non-zero value
               ENDIF
            ENDIF
   !
   !   Look for Site Identifier   -- Set to NO-PRINT  if not present.
            ksiteid = INDEX(txt, 'SITE')
            IF ( ksiteid.GT.0 ) THEN
               MNWSITE(ipt) = txtraw(ksiteid+5:256)
               kblk = INDEX(MNWSITE(ipt), ' ')
               ktab = INDEX(MNWSITE(ipt), tab)
               IF ( kblk.GT.0 ) kfini = kblk
               IF ( ktab.GT.0 .AND. ktab.LT.kblk ) kfini = ktab
               IF ( kfini.LE.32 ) THEN
                  MNWSITE(ipt)(kfini:32) = '                 '
               ELSE
                  kfini = 32
               ENDIF
               txt(ksiteid:ksiteid+kfini+4) = '                        '
            ELSE
               MNWSITE(ipt) = 'NO-PRINT                     '
            ENDIF
   !
   !    Read remaining info from card to set MNW specific parameters
            tx2 = txt(irmx:256)
            ki = INDEX(tx2, 'ZONE')
            IF ( ki.GT.0 ) tx2(ki:256) = '                         '
            CALL QREAD(rn, 6, tx2, ierr)
   !
   !   Move from well data from temp to permanent locations
            DO ip = 1, 6 - ierr
               WELL2(ip+3, ipt) = rn(ip)
            ENDDO
            IF ( ierr.GE.1 ) WELL2(9, ipt) = ipt
            IF ( ierr.GE.2 .OR. ABS(WELL2(8,ipt)).GT.HMAX )&
            &WELL2(8, ipt) = HREF(i,j,k)
   !  Compute HLIM relative to reference elevation if HLIM read was a DrawDown (DD)
            IF ( INDEX(txt, 'DD').GT.0 )&
            &WELL2(7, ipt) = ipole*WELL2(7, ipt) + WELL2(8, ipt)
            IF ( ierr.GE.3 ) WELL2(7, ipt) = ipole*1.0D+26
            IF ( ierr.GE.4 ) WELL2(6, ipt) = 0.0D0
            IF ( ierr.GE.5 ) WELL2(5, ipt) = 0.0D0
            IF ( ierr.GE.6 ) WELL2(4, ipt) = -1.0D0
   !  Flag as 2-point definition of a multi-node well if MULTI is detected.
            IF ( INDEX(tx2, 'MULTI').GT.0 .AND.&
            &ABS(WELL2(5,ipt)).GT.ZERO25 ) THEN
   !  Define direction and # of points in well
               WELL2(2, ipt-1) = WELL2(2, ipt) + WELL2(2, ipt-1)
               n1 = IFRL(WELL2(1, ipt-1))
               mstep = IDIRECT(n1, node, NCOL, NROW)
               DO nn = n1 + mstep, node, mstep
                  ipt = ipt + 1
                  NWELL2 = NWELL2 + 1
                  WELL2(1, ipt) = nn
                  WELL2(2, ipt) = 0.0D0
                  WELL2(3, ipt) = WELL2(2, ipt)
                  WELL2(4, ipt) = WELL2(4, ipt-1)
                  WELL2(5, ipt) = WELL2(5, ipt-1)
                  WELL2(6, ipt) = WELL2(6, ipt-1)
                  WELL2(16, ipt) = WELL2(16, ipt-1)  !!  NonLinear Loss Coefficient
                  WELL2(9, ipt) = WELL2(9, ipt-1)
                  WELL2(8, ipt) = -1.0D31
                  WELL2(13, ipt) = 0.0D0
                  WELL2(14, ipt) = 0.0D0
                  icmn = icmn + 1
                  WELL2(7, ipt) = icmn
               ENDDO
   !  Flag as part of a multi-node well if MN is detected.
            ELSE IF ( INDEX(tx2, 'MN').GT.0 .AND.&
            &ABS(WELL2(5,ipt)).GT.ZERO25 ) THEN
   !  Set to very large -value to flag MN status
               WELL2(8, ipt) = -1.0D31
               icmn = icmn + 1
               WELL2(7, ipt) = icmn
            ELSE
               icmn = 1
            ENDIF
         ELSE
   !   Sum details on rejected wells
            qreject = qreject + q
            nqreject = nqreject + 1
         ENDIF   !   IBOUND test statement
      ENDDO     !   end of well entry loop
   !
   !   Process wells that are screened across multiple nodes
   !
   ! Check for extreme contrast in conductance
   !
      WELL2(8, NWELL2+1) = 0.0D0
      IF ( nstart.LT.1 ) nstart = 1
   !swm        IF ( nstart.GT.NWELL2 ) nstart = NWELL2 - itmp + 1
      IF ( nstart.GT.NWELL2 ) nstart = itmp - NWELL2 + 1
      DO i = nstart, NWELL2
         IF ( WELL2(8, i).LT.-1.E30 .AND. WELL2(8, i+1).GT.-1.E30 .OR.&
         &WELL2(8, i).LT.-1.E30 .AND. i.EQ.NWELL2 ) THEN
            ngrp = IFRL(WELL2(7, i))
            ne = i
            nb = ne - ngrp + 1
            hlim = WELL2(7, nb)
            hrfw = WELL2(8, nb)
            qsum = 0.0D0
            tempsite = 'NO-PRINT                     '
            DO iin = nb, ne
               qsum = qsum + WELL2(2, iin)
               IF ( MNWSITE(iin)(1:8).NE.'NO-PRINT' )&
               &tempsite = MNWSITE(iin)
               WELL2(2, iin) = 0.0D0
               WELL2(7, iin) = 1.0D31
               WELL2(8, iin) = 1.0D31
   !   Set to very large +value to flag MN status
            ENDDO
   !   Set All SiteIDs in a multinode well to a common tag
            DO iin = nb, ne
               MNWSITE(iin) = tempsite
            ENDDO
            WELL2(7, nb) = ne
            WELL2(2, ne) = qsum
            WELL2(7, ne) = hlim
            WELL2(8, ne) = hrfw
         ENDIF
      ENDDO  !   end of multi-well pointer setting
   !
   ENDIF
   !
   !  nwell2>mxwel2.  print message. stop.
   IF ( NWELL2.GT.MXWEL2 ) THEN
      WRITE (IOUT, 9002) NWELL2, MXWEL2
9002  FORMAT (1X,/&
      &1X,'nwell2(', i4, ') IS GREATER THAN mxwel2(', i4, ')')
   !
   !       When compiling MNW with Modflow-96, comment out the call to
   !       USTOP and uncomment the STOP statement
      CALL USTOP(' ')
   !        STOP
   !
   ENDIF
   !
   !   Place desired flow rates in a reserved location
   !
   DO m = 1, NWELL2
      WELL2(15, m) = WELL2(2, m)
   ENDDO
   !
   !--assign unique well id for use with MT3DMS link package (cdl: 4/19/05)
   m=0
   IDwell=1
   do while (m.lt.nwell2)
      m=m+1
      if(well2(8,m).gt.1.D30) then
         do mm=m,ifrl(well2(7,m))
            well2(18,mm)=IDwell
         enddo
         m=ifrl(well2(7,m))
      else
         well2(18,m)=IDwell
      endif
      IDwell=IDwell+1
   enddo
   !
   !   Echo input to iout file
   !
   IF ( IWELPT.EQ.0 ) THEN
      IF ( nqreject.GT.0 ) THEN
         txt = ' wells were outside of the model domain.'
         WRITE (IOUT, '(1X,/,5x,i5,a50)') nqreject, txt
         txt = 'The rejected pumpage totaled: '
         WRITE (IOUT, '(1X,/1X,a34,g14.5)') txt, qreject
      ENDIF
   !
      WRITE (IOUT, '(1X,/,10x,i5," MNW WELLS")') NWELL2
      WRITE (IOUT, 9003)
9003  FORMAT ('    No.   Lay   Row   Col    Stress   QW param', 6x,&
      &'Rw       Skin    WL Limit    WL Refer   NonLinear Cp',&
      &'  QW Group  Cell-To-Well  Min-Qoff  Min-Qon',&
      &'  Site Identifier')
   !
      DO m = 1, NWELL2
         n = INT(WELL2(1, m))
         k = (n-1)/(NCOL*NROW) + 1
         j = MOD((n-1), NCOL*NROW)/NCOL + 1  !swm: note i,j, are switched
         i = MOD((n-1), NCOL) + 1            !swm: from usual MF2K style
         igrp = INT(WELL2(9, m))
   !
         rw = WELL2(5, m)
         IF ( rw.LT.-ZERO25 ) THEN
            cond = -rw
         ELSE
            qact = WELL2(3, m)
            sk = WELL2(6, m)
            cf = WELL2(16, m)
            IF ( Iubcf.NE.0 ) cond = CEL2WELBCF(i,j,k,rw,sk,qact,cf)
            IF ( Iulpf.NE.0 ) cond = CEL2WELLPF(i,j,k,rw,sk,qact,cf)
            IF ( Iuhuf.NE.0 ) cond = CEL2WELHUF(i,j,k,rw,sk,qact,cf)
            IF ( Iuupw.NE.0 ) cond = CEL2WELUPW(i,j,k,rw,sk,qact,cf)
            IF ( rw.LT.ZERO25 ) cond = cond*1.0D3
         ENDIF
         WELL2(11, m) = cond
   !
   ! ---------Modified OUTPUT to hide internal pointers that "Look Funny" --KJH-- July 10, 2003
         IF ( WELL2(8, m).LT.BIG ) THEN
            hlim = WELL2(7, m)
            hrfw = WELL2(8, m)
         ELSE IF ( WELL2(7, m).LT.BIG ) THEN
            ne = IFRL(WELL2(7, m))
            hlim = WELL2(7, ne)
            hrfw = WELL2(8, ne)
         ENDIF
         WRITE (IOUT, 9004) m, k, j, i, (WELL2(ii, m), ii=3, 6), hlim,&
         &hrfw, WELL2(16, m), igrp, WELL2(11, m),&
         &(WELL2(ii, m)*100.0D0, ii=13, 14),&
         &MNWSITE(m)
9004     FORMAT (1x, 4I6, 6g11.4, g13.6, i10, g13.6, 2F10.3, 2x, a32)
   !
      ENDDO
   ELSE
      WRITE (IOUT, *) 'WELLS WILL NOT BE PRINTED'
   ENDIF
   !
   !  Write blank fields in Auxillary BYNODE file if KPER=1 & IO>0
   !
   IF ( TOTIM.LT.1E-26 .AND. IOWELL2(2).NE.0 ) THEN
      io = ABS(IOWELL2(2))
      DO m = 1, NWELL2
         n = IFRL(WELL2(1, m))
         WRITE (io, '(a32,1x,2i8)') MNWSITE(m), m, n
      ENDDO
   ENDIF
   !
   !  Write blank fields in Auxillary QSUM file if KPER=1 & IO>0
   !
   IF ( TOTIM.LT.1E-26 .AND. IOWELL2(3).NE.0 ) THEN
      io = ABS(IOWELL2(3))
      m = 0
      DO WHILE ( m.LT.NWELL2 )
         m = m + 1
         IF ( WELL2(8, m).GT.BIG ) THEN
            ne = IFRL(WELL2(7, m))
            WRITE (ABS(IOWELL2(3)), '(a32,1x,i5.5,"-",i5.5)') MNWSITE(m)&
            &, m, ne
            m = ne
         ENDIF
      ENDDO
   ENDIF
   !
END SUBROUTINE GWF2MNW17RP
   !
   !_________________________________________________________________________________
   !
SUBROUTINE GWF2MNW17AD(Iubcf, Iulpf, Iuhuf, Iuupw, Igrid)
   !     VERSION 20020819 KJH
   !
   !----- MNW by K.J. Halford
   !
   !     ******************************************************************
   !     Update Qact for wells that were constrained
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,HNEW
   USE GWFMNW1MODULE,ONLY:NWELL2,SMALL,WELL2,ZERO8,BIG
   IMPLICIT NONE
   INTRINSIC ABS, MOD
   INTEGER, EXTERNAL :: IFRL
   DOUBLE PRECISION, EXTERNAL :: CEL2WELBCF, CEL2WELLPF, CEL2WELHUF,&
   &CEL2WELUPW
   ! Arguments
   INTEGER, INTENT(IN) :: Iubcf, Iuhuf, Iulpf, Igrid, Iuupw
   ! Local Variables
   DOUBLE PRECISION :: rw, cond, qact, sk, cf, qoff, qon, qsmall
   DOUBLE PRECISION :: qdes, csum, chsum, hwell, hlim, href, ipole
   DOUBLE PRECISION :: ddmax, ddsim, qpot, ratio
   INTEGER iin, m, n, ne, k, j, i
   !     ------------------------------------------------------------------
   !-----SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2MNW1PNT(Igrid)
   !
   !
   !1------if number of wells <= 0 then return.
   IF ( NWELL2.LE.0 ) RETURN
   !
   !   Compute cell-to-well conductance for each well node
   !
   DO m = 1, NWELL2
      n = IFRL(WELL2(1, m))
      k = (n-1)/(NCOL*NROW) + 1
      j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
      i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
   !-----if the cell is inactive or specified then bypass processing.
      IF ( IBOUND(i,j,k).NE.0 ) THEN
         rw = WELL2(5, m)
         IF ( rw.LT.-ZERO8 ) THEN
            cond = -rw
         ELSE
            qact = WELL2(3, m)
            sk = WELL2(6, m)
            cf = WELL2(16, m)
            IF ( Iubcf.NE.0 ) cond = CEL2WELBCF(i,j,k,rw,sk,qact,cf)
            IF ( Iulpf.NE.0 ) cond = CEL2WELLPF(i,j,k,rw,sk,qact,cf)
            IF ( Iuhuf.NE.0 ) cond = CEL2WELHUF(i,j,k,rw,sk,qact,cf)
            IF ( Iuupw.NE.0 ) cond = CEL2WELUPW(i,j,k,rw,sk,qact,cf)
            IF ( rw.LT.ZERO8 ) cond = cond*1.0D3
         ENDIF
         WELL2(11, m) = cond
      ENDIF
   ENDDO
   !
   !   Allow constrained wells a new chance with the next time step
   !
   m = 0
   DO WHILE ( m.LT.NWELL2 )
      m = m + 1
      qoff = WELL2(13, m)
      qon = WELL2(14, m)
      qact = WELL2(3, m)
      qsmall = SMALL
   !
   !   A very large # in WL reference array (8,m) triggers multi-node calculation
   !
      IF ( WELL2(8, m).GT.BIG ) THEN
   !     Compute hwell / Qpot for multi-node well
         ne = IFRL(WELL2(7, m))
         qdes = WELL2(15, ne)
         csum = 0.0D0
         chsum = 0.0D0
         qact = 0.0D0
         qsmall = SMALL*ABS(qdes)
         DO iin = m, ne
            n = IFRL(WELL2(1, iin))
            k = (n-1)/(NCOL*NROW) + 1
            j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
            i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
            IF ( IBOUND(i,j,k).NE.0 ) THEN
               csum = csum + WELL2(11, iin)
               chsum = chsum + WELL2(11, iin)*HNEW(i,j,k)
               qact = qact + WELL2(3, iin)
            ELSE
               qact = 0.0D0
            ENDIF
         ENDDO
   !---div0 ---  CSUM could go to zero if the entire well is dry
         IF ( csum.GT.ZERO8 ) THEN
            hwell = (qdes+chsum)/csum
         ELSE
            hwell = HNEW(i,j,k)
         ENDIF
         m = ne
   !   Test DD constraints, Hlim is assumed to be a Max/Min for Injection/Production wells
         ipole = 0.0D0
         IF ( ABS(qdes).GT.ZERO8 ) ipole = qdes/ABS(qdes)
         hlim = WELL2(7, ne)
         href = WELL2(8, ne)
         ddmax = ipole*(hlim-href)
         ddsim = ipole*(hwell-href)
         qpot = hlim*csum - chsum
         IF ( ddsim.GT.ddmax ) THEN
            hwell = hlim
            qpot = hwell*csum - chsum
         ENDIF
         cond = csum
      ELSE       !  End of multi-node conditioning IF statement
   !     Compute hwell / Qpot for single-node well
         n = IFRL(WELL2(1, m))
         k = (n-1)/(NCOL*NROW) + 1
         j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
         i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
         cond = WELL2(11, m)
         hlim = WELL2(7, m)
         qpot = (hlim-HNEW(i,j,k))*cond
         qdes = WELL2(15, m)
      ENDIF
   !
   !  Compute ratio of potential/desired flow rates
      ratio = 1.00D0
      IF ( ABS(qdes).GT.SMALL ) ratio = qpot/qdes
      IF ( ratio.GT.0.9999D0 ) THEN
         ratio = 1.0D0
         qpot = qdes
      ENDIF
   !  Check if potential flow rate is below cutoff
      IF ( ratio.LT.qoff ) THEN
         qact = 0.0D0
         qdes = qact
         WELL2(2, m) = qdes
         WELL2(3, m) = qact
   !  Check if potential flow rate is above restart threshold
      ELSE IF ( ratio.GT.qon .AND. ABS(qact).LT.qsmall ) THEN
         qdes = WELL2(15, m)
         WELL2(2, m) = qdes
         WELL2(3, m) = qpot
   !       ELSE
   !  Otherwise leave the flow rate alone
      ENDIF
   !
   ENDDO  ! End of overall test loop
   !
END SUBROUTINE GWF2MNW17AD
   !
   !_________________________________________________________________________________
   !
SUBROUTINE GWF2MNW17FM(Kiter, Iubcf, Iulpf, Iuhuf, Iuupw, Igrid)
   !     VERSION 20020819 KJH
   !
   !----- MNW by K.J. Halford
   !
   !     ******************************************************************
   !     add well flow to source term
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,HNEW,HCOF,RHS
   USE GWFMNW1MODULE,ONLY:NWELL2,NOMOITER,SMALL,WELL2,ZERO20,BIG
   IMPLICIT NONE
   INTRINSIC ABS, MOD
   INTEGER, EXTERNAL :: IFRL
   DOUBLE PRECISION, EXTERNAL :: CEL2WELBCF, CEL2WELLPF, CEL2WELHUF,&
   &CEL2WELUPW
   ! Arguments
   INTEGER, INTENT(IN) :: Iubcf, Iuhuf, Iulpf, Kiter, Igrid, Iuupw
   ! Local Variables
   INTEGER :: iin, iqslv, m, n, ne, k, j, i
   DOUBLE PRECISION :: rw, cond, qact, sk, cf, qdes, csum, chsum
   DOUBLE PRECISION :: hwell, ipole, hlim, href, ddmax, ddsim, ratio
   DOUBLE PRECISION :: dhc2w
   !     ------------------------------------------------------------------
   !-----SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2MNW1PNT(Igrid)
   !
   !                 CR( i, j, k)    ------>   CR  i + 1/2
   !                 CC( i, j, k)    ------>   CC  j + 1/2
   !                 CV( i, j, k)    ------>   CV  k + 1/2
   !
   !1------if number of wells <= 0 then return.
   IF ( NWELL2.LE.0 ) RETURN
   !
   !   Compute cell-to-well conductance for each well node
   !
   DO m = 1, NWELL2
      n = IFRL(WELL2(1, m))
      k = (n-1)/(NCOL*NROW) + 1
      j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
      i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
   !-----if the cell is inactive or specified then bypass processing.
      IF ( IBOUND(i,j,k).NE.0 ) THEN
         rw = WELL2(5, m)
         IF ( rw.LT.-ZERO20 ) THEN
            cond = -rw
         ELSE
            qact = WELL2(3, m)
            sk = WELL2(6, m)
            cf = WELL2(16, m)
            IF ( Iubcf.NE.0 ) cond = CEL2WELBCF(i,j,k,rw,sk,qact,cf)
            IF ( Iulpf.NE.0 ) cond = CEL2WELLPF(i,j,k,rw,sk,qact,cf)
            IF ( Iuhuf.NE.0 ) cond = CEL2WELHUF(i,j,k,rw,sk,qact,cf)
            IF ( Iuupw.NE.0 ) cond = CEL2WELUPW(i,j,k,rw,sk,qact,cf)
            IF ( rw.LT.ZERO20 ) cond = cond*1.0D3
         ENDIF
         WELL2(11, m) = cond
      ENDIF
   ENDDO
   !
   !   Prepare components and limits of a multi-node well
   m = 0
   DO WHILE ( m.LT.NWELL2 )
      m = m + 1
      WELL2(10, m) = 1.0D31
   !
   !   A very large # in WL reference array (8,m) triggers multi-node calculation
   !
      IF ( WELL2(8, m).GT.BIG ) THEN
         ne = IFRL(WELL2(7, m))
         qdes = WELL2(2, ne)
         qact = qdes
         csum = 0.0D0                    !swm: why not initialize DP
         chsum = 0.0D0                   !swm: why not initialize DP
         DO iin = m, ne
            n = IFRL(WELL2(1, iin))
            k = (n-1)/(NCOL*NROW) + 1
            j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
            i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
            IF ( IBOUND(i,j,k).NE.0 ) THEN
               csum = csum + WELL2(11, iin)
               chsum = chsum + WELL2(11, iin)*HNEW(i,j,k)
            ELSE
               WELL2(3, iin) = 0.0D0      !swm: why not initialize DP
            ENDIF
         ENDDO
   !---div0 ---  CSUM could go to zero if the entire well is dry
         IF ( csum.GT.ZERO20 ) THEN
            hwell = (qact+chsum)/csum
         ELSE
            hwell = HNEW(i,j,k)
         ENDIF
   !
   !   Test DD constraints, Hlim is assumed to be a Max/Min for Injection/Production wells
         ipole = 0.0D0
         IF ( ABS(qdes).GT.ZERO20 ) ipole = qdes/ABS(qdes)
         hlim = WELL2(7, ne)
         href = WELL2(8, ne)
         ddmax = ipole*(hlim-href)
         ddsim = ipole*(hwell-href)
   !
         IF ( ddsim.GT.ddmax ) THEN
            hwell = hlim
            qact = hwell*csum - chsum
   !      DD constraints that stop production are not tested until after the 2nd iteration
            IF ( Kiter.GT.2 ) THEN
               ratio = 1.00D0
               IF ( ABS(qdes).GT.SMALL ) ratio = qact/qdes
               IF ( ratio.LT.0.00001 ) THEN
                  qact = 0.0D0
                  IF ( csum.GT.0.0 ) THEN
                     hwell = chsum/csum
                  ELSE
                     hwell = HNEW(i,j,k)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   !
   !   Assign flow rates and water levels to individual nodes
         DO iin = m, ne
            n = IFRL(WELL2(1, iin))
            k = (n-1)/(NCOL*NROW) + 1
            j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
            i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
            WELL2(10, iin) = hwell
            qact = (hwell-HNEW(i,j,k))*WELL2(11, iin)
            WELL2(3, iin) = qact
         ENDDO
         m = ne
      ENDIF      !  End of multi-node conditioning IF statement
   ENDDO      ! End of overall multi-node test loop
   !
   !2------process each well in the well list.
   m = 0
   DO WHILE ( m.LT.NWELL2 )
      m = m + 1
      n = IFRL(WELL2(1, m))
      k = (n-1)/(NCOL*NROW) + 1
      j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
      i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
      qdes = WELL2(2, m)
   !-----if the cell is inactive then bypass processing.
      IF ( IBOUND(i,j,k).GT.0 ) THEN
         qact = WELL2(3, m)
         cond = WELL2(11, m)
   !
         hlim = WELL2(7, m)
         href = WELL2(8, m)
         IF ( WELL2(10, m).GT.BIG .AND. cond.GT.ZERO20 ) THEN
            dhc2w = qact/cond
   !   Process single-node wells
   !   Test DD constraints, Hlim is assumed to be a Max/Min for Injection/Production wells
            ipole = 0.0D0
            IF ( ABS(qdes).GT.ZERO20 ) ipole = qdes/ABS(qdes)
            hwell = HNEW(i,j,k) + dhc2w
            WELL2(10, m) = hwell
            ddsim = ipole*(hwell-href)
            ddmax = ipole*(hlim-href) - SMALL
            ratio = 1.00D0
            IF ( ABS(qdes).GT.ZERO20 ) ratio = qact/qdes
            IF ( ABS(ratio).GT.1.00 ) qact = qdes
            IF ( ratio.LT.ZERO20 ) qact = 0.0D0
   !    Well will be simulated as a specified rate or GHB
            iqslv = 0
            IF ( ddsim.GT.ddmax .AND. ddmax.GT.ZERO20 ) iqslv = 1
            IF ( (qdes-qact)**2.GT.SMALL ) iqslv = 1
            IF ( ABS(qact).LT.ZERO20 .AND. ddsim.GT.ddmax ) iqslv = 0
            IF ( ABS(qact).LT.ZERO20 .AND. ddsim.LT.ddmax ) iqslv = 1
            IF ( ABS(qdes).LT.ZERO20 .OR. ratio.GT.1.0D0-ZERO20 )&
            &iqslv = 0
   !
         ELSE IF ( cond.LT.ZERO20 ) THEN
            qact = 0.0D0
            iqslv = 0
   ! Process multi-node wells, Constraints were already tested when allocating flow
         ELSE IF ( MOD(Kiter, 2).EQ.0 .AND. ABS(qact).GT.SMALL ) THEN
            hlim = WELL2(10, m)
            iqslv = 1
         ELSE
            qact = WELL2(3, m)
            iqslv = 0
         ENDIF
   ! Added next lines to stop BC switching. RGN 7/20/11
   !
         IF ( Iuupw.GT. 0 ) THEN
            hlim = WELL2(10, m)
            iqslv = 1
         END IF
   !
   !   Modify HCOF and RHS arrays
         IF ( iqslv.NE.0 .AND. Kiter.GT.1 .AND. Kiter.LT.NOMOITER )THEN
            qact = (hlim-HNEW(i,j,k))*cond
            HCOF(i,j,k) = HCOF(i,j,k) - cond
            RHS(i,j,k) = RHS(i,j,k) - cond*hlim
         ELSE
   !  Specify Q and solve for head;  add Q to RHS accumulator.
            RHS(i,j,k) = RHS(i,j,k) - qact
         ENDIF
         WELL2(3, m) = qact
      ENDIF
   ENDDO      !    End of DO WHILE loop
   !
END SUBROUTINE GWF2MNW17FM
   !
   !_________________________________________________________________________________
   !
SUBROUTINE GWF2MNW17BD(Nstp, Kstp, Kper, Igrid)
   !     VERSION 20030710 KJH
   !
   !----- MNW by K.J. Halford        1/31/98
   !     ******************************************************************
   !     calculate volumetric budget for wells
   !     ******************************************************************
   !
   !        specifications:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,BUFF,HNEW,IOUT
   USE GWFBASMODULE,ONLY:DELT,PERTIM,TOTIM,ICBCFL,VBVL,VBNM,MSUM,HDRY
   USE GWFMNW1MODULE,ONLY:NWELL2,PLOSS,MNWSITE,IWL2CB,IOWELL2,&
   &WELL2,ZERO25,BIG
   IMPLICIT NONE
   INTRINSIC ABS, MOD
   INTEGER, EXTERNAL :: IFRL
   EXTERNAL UBDSV4, UBDSVB, UBUDSV
   ! Arguments
   INTEGER, INTENT(IN) :: Nstp, Kstp, Kper, Igrid
   ! Local Variables
   DOUBLE PRECISION :: ratin, ratout, qwsum, qsum, qwbar, drytest
   DOUBLE PRECISION :: qd, hlim, href, hwell, dd, s, ipole, sNL, sL
   DOUBLE PRECISION :: qin, qout, qwfsum, q
   REAL :: dummy(5), qsing
   INTEGER :: ibd, igrp1, igrp2, iin, imult, iobynd, ioc
   INTEGER :: ioch, ioqsum, m, m2, n, naux, ne, nwelvl, k, i, j
   CHARACTER(LEN=16) :: text, auxtxt(20)
   !     ------------------------------------------------------------------
   !-----SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2MNW1PNT(Igrid)
   !
   !             ----+----1----+-
   text = '             MNW'
   !     ------------------------------------------------------------------
   !  clear ratin and ratout accumulators.
   ratin = 0.0D0                         !swm: why not initialized DP
   ratout = 0.0D0                        !swm: why not initialized DP
   ibd = 0
   IF ( IWL2CB.GT.0 ) ibd = ICBCFL
   !
   !2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
   IF ( ibd.EQ.2 ) THEN
      auxtxt = ' ' !rsr set auxtxt to blank strings
      naux = 0    !!   Set to zero -- Change order to dump
   !         IF ( IAUXSV.EQ.0 ) NAUX=0
      CALL UBDSV4(Kstp, Kper, text, naux, auxtxt, IWL2CB, NCOL, NROW,&
      &NLAY, NWELL2, IOUT, DELT, PERTIM, TOTIM, IBOUND)
   ENDIF
   !  clear the buffer.
   DO k = 1, NLAY
      DO i = 1, NROW
         DO j = 1, NCOL
            Buff(j,i,k) = 0.0
         ENDDO
      ENDDO
   ENDDO
   ! -----print the header for individual rates if requested(IWL2CB<0).
   IF ( IWL2CB.LT.0 .AND. ICBCFL.NE.0 ) THEN
      WRITE (IOUT, '(/,1x,a16," PERIOD =",i5,8h  STEP =,i5)') text,&
      &Kper, Kstp
      WRITE (IOUT, 9001)
   ENDIF
9001 FORMAT ('  Entry LAY ROW COL', 9x, 'Q', 6x, 'H-Well', 7x,&
   &'H-Cell', 7x, 'DD    ', 7x, 'QW-Avg', 6x, 's-LINEAR', 3x,&
   &'s-NonLINEAR')
   !
   !  Create WEL1 file if iowell2(1) > 0
   IF ( IOWELL2(1).GT.0 .AND. Nstp.EQ.Kstp )&
   &WRITE (IOWELL2(1), '(1i10)') NWELL2
   !
   !2------if there are no wells do not accumulate flow
   IF ( NWELL2.GT.0 ) THEN
   !
   !     Compute flow weighted QW values and store in well2(11,m)
   !
      DO m = 1, NWELL2
         WELL2(11, m) = WELL2(3, m)*WELL2(4, m)
         WELL2(12, m) = 0.0D0          !swm: why not initialized DP
         IF ( WELL2(4, m).LT.0.00 .OR. WELL2(3, m).GT.0.00 ) THEN
            WELL2(11, m) = -1
            WELL2(12, m) = 1
         ENDIF
      ENDDO
   !
      DO m = 1, NWELL2
         igrp1 = IFRL(WELL2(9, m))
         IF ( WELL2(12, m).LT.0.5 ) THEN
            qwsum = 0.0D0
            qsum = 0.0D0
            DO m2 = m, NWELL2
               igrp2 = IFRL(WELL2(9, m2))
               IF ( igrp1.EQ.igrp2 .AND. WELL2(12, m2).LT.0.5 ) THEN
                  qwsum = qwsum + WELL2(11, m2)
                  qsum = qsum + WELL2(3, m2)
                  WELL2(12, m2) = 1
               ENDIF
            ENDDO
   !
            qwbar = qwsum
            IF ( qsum**2.GT.ZERO25 ) qwbar = qwsum/qsum
            DO m2 = m, NWELL2
               igrp2 = IFRL(WELL2(9, m2))
               IF ( igrp1.EQ.igrp2 .AND. WELL2(4, m2).GE.0.0 )&
               &WELL2(11, m2) = qwbar
            ENDDO
         ENDIF
      ENDDO
   !
      imult = 0
   !rsr added next line to be sure ioch has a value
      ioch = 0
      DO m = 1, NWELL2
         n = IFRL(WELL2(1, m))
         k = (n-1)/(NCOL*NROW) + 1
         j = MOD((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
         i = MOD((n-1),NCOL) + 1           !swm: from usual MF2K style
         drytest = HNEW(i,j,k) - HDRY
         IF ( ABS(drytest).LT.ZERO25 ) WELL2(3, m) = 0.0D0  !swm: why not initialized DP
         q = WELL2(3, m)
         WELL2(17, m) = q  !!7/13/2003 - CZ: preserve q
   !
   !    Report all wells with production less than the desired rate......
         IF ( IBOUND(i,j,k).NE.0 .OR. ABS(drytest).LT.ZERO25 ) THEN
            qd = WELL2(2, m)
            hlim = WELL2(7, m)
   ! -----Modified OUTPUT to hide internal pointers that "Look Funny" in DD column--KJH-- July 10, 2003
            IF ( WELL2(8, m).GT.BIG ) THEN
               imult = 1
               IF ( hlim.LT.BIG ) THEN
                  ne = IFRL(hlim)
                  href = WELL2(8, ne)
   !             ELSE
               ENDIF
            ELSE
               href = WELL2(8, m)
            ENDIF
            hwell = WELL2(10, m)
            qwbar = WELL2(11, m)
            dd = hwell - href
   !
            ioch = 0
            IF ( IWL2CB.LT.0 .AND. ICBCFL.NE.0 ) ioch = 1
   ! -----print the individual rates if requested(IWL2CB<0).
            IF ( ioch.EQ.1 ) THEN
               s = HNEW(i,j,k) - hwell
               ipole = 0.0D0
               IF ( ABS(s).GT.ZERO25 ) ipole = s/ABS(s)
               sNL = ipole*WELL2(16, m)*ABS(q)**PLOSS
               sL = s - sNL
               WRITE (IOUT, '(i7,3i4,9g13.6)') m, k, j, i, q, hwell,&
               &HNEW(i,j,k), dd, qwbar, sL, sNL
            ENDIF
   !
   ! -----print the individual rates to auxillary file if requested(IWL2CB<0).
            iobynd = ABS(IOWELL2(2))
            IF ( iobynd.GT.0 ) THEN
               IF ( ioch.EQ.1 .OR. IOWELL2(2).LT.0 )&
               &WRITE (iobynd, '(a32,1x,2i8,6g15.8)')&
               &MNWSITE(m), m, n, TOTIM, q, hwell, HNEW(i,j,k), qwbar
            ENDIF
   !  Create WEL1 file if iowell2(1) > 0
            IF ( IOWELL2(1).GT.0 .AND. Nstp.EQ.Kstp )&
            &WRITE(IOWELL2(1),'(i9,2i10,1X,g11.4,1X,i10,2x,6g11.4)')&
            &k, j, i, q, 0, qd, hwell, HNEW(i,j,k), dd, href, qwbar
   !
            Buff(i,j,k) = Buff(i,j,k) + q
            IF ( q.GE.0.0 ) THEN
   ! -----pumping rate is positive(recharge). add it to ratin.
               ratin = ratin + q
            ELSE
   ! -----pumping rate is negative(discharge). add it to ratout.
               ratout = ratout - q
            ENDIF
         ENDIF
      ENDDO
   !
   !   Sum components of  multi-node wells
   !
   ! -----print the header for multi-node rates if requested(IWL2CB<0).
      IF ( ioch.EQ.1 .AND. imult.EQ.1 ) THEN
         WRITE (IOUT, '(/,5x," Multi-Node Rates & Average QW")')
         WRITE (IOUT, 9002)
9002     FORMAT (' Site Identifier ', 5x, 'ENTRY: Begin - End', 2x,&
         &'Q-Total', 7x, 'H-Well', 7x, 'DD    ', 7x, 'QW-Avg')
      ENDIF
   !
      m = 0
      DO WHILE ( m.LT.NWELL2 )
         m = m + 1
         IF ( WELL2(8, m).GT.BIG ) THEN
            ne = IFRL(WELL2(7, m))
            qwsum = 0.0D0
            qwfsum = 0.0D0
            qsum = 0.0D0
            qin = 0.0D0
            qout = 0.0D0
            DO iin = m, ne
               n = IFRL(WELL2(1, iin))
               k = (n-1)/(NCOL*NROW) + 1
               j = mod((n-1),NCOL*NROW)/NCOL + 1 !swm: note i,j are switched
               i = mod((n-1),NCOL) + 1           !swm: from usual MF2K style
               IF ( IBOUND(i,j,k).EQ.0 ) WELL2(3, iin) = 0.0D0
               IF ( WELL2(4, iin).GE.0.0 .AND. WELL2(3, iin).LE.0.0 )THEN
                  qwfsum = qwfsum + WELL2(3, iin)
                  qwsum = qwsum + WELL2(3, iin)*WELL2(4, iin)
               ENDIF
               IF ( WELL2(3, iin).LE.0.0D0 ) THEN
                  qin = qin + WELL2(3, iin)
               ELSE
                  qout = qout + WELL2(3, iin)
               ENDIF
               qsum = qsum + WELL2(3, iin)
               WELL2(3, iin) = 0.0D0
            ENDDO
            WELL2(3, ne) = qsum
   ! -----print the summed rates if requested(IWL2CB<0).
            qwbar = WELL2(4, ne)
            IF ( qwfsum**2.GT.ZERO25 ) qwbar = qwsum/qwfsum
            href = WELL2(8, ne)
            hwell = WELL2(10, ne)
            dd = hwell - href
            IF ( ioch.EQ.1 ) WRITE (IOUT, '(A26,1x,2i6,6g13.6)')&
            &MNWSITE(m), m, ne, qsum, hwell, dd,&
            &qwbar
   ! -----print the summed rates to auxillary file if requested .
            ioqsum = ABS(IOWELL2(3))
            IF ( ioqsum.GT.0 ) THEN
               IF ( ioch.EQ.1 .OR. IOWELL2(3).LT.0 ) THEN
                  WRITE (ioqsum, 9003) MNWSITE(m), m, ne, TOTIM, qin,&
                  &qout, qsum, hwell, qwbar
9003              FORMAT (A32, i6.5, '-', i5.5, 12g15.8)
               ENDIF
            ENDIF
            m = ne
         ENDIF
      ENDDO
      IF ( ioch.EQ.1 .AND. imult.EQ.1 ) WRITE (IOUT, *)
   !
   !  ----- END  MULTI-NODE reporting section -------------
   !
   !6------if cell-by-cell flows will be saved call ubudsv to record them
      IF ( ABS(IWL2CB).GT.0 .AND. ICBCFL.NE.0 ) THEN          !! BooBoo Fix--July 10,2003  KJH
         ioc = ABS(IWL2CB)
         IF ( ibd.EQ.2 ) THEN  !!  Write COMPACT budget
            nwelvl = 1   !!  Dummy value
            DO m = 1, NWELL2
               n = IFRL(WELL2(1, m))
               qsing = WELL2(17, m)
               CALL UBDSVB(ioc, NCOL, NROW, n, 1, 1, qsing, dummy,&
               &nwelvl, naux, 5, IBOUND, NLAY)
            ENDDO
         ELSE                  !!  Write full 3D array
            CALL UBUDSV(Kstp, Kper, text, ioc, Buff, NCOL, NROW, NLAY,&
            &IOUT)
         ENDIF
      ENDIF
   ENDIF
   !
   !7------move rates into vbvl for printing by module bas1ot.
   VBVL(3, MSUM) = ratin
   VBVL(4, MSUM) = ratout
   !
   !8------move rates times time step length into vbvl accumulators.
   VBVL(1, MSUM) = VBVL(1, MSUM) + ratin*DELT
   VBVL(2, MSUM) = VBVL(2, MSUM) + ratout*DELT
   !
   !9------move budget term labels into vbnm for printing.
   VBNM(MSUM) = text
   !
   !10-----increment budget term counter(msum).
   MSUM = MSUM + 1
   !
   !11-----return
END SUBROUTINE GWF2MNW17BD
   !
   !_________________________________________________________________________________
   !
SUBROUTINE GWF2MNW17OT(Igrid)
   !     VERSION 20020819 KJH
   !
   !     ******************************************************************
   !     Sort well output into useful tables
   !     ******************************************************************
   !
   !        specifications:
   !     ------------------------------------------------------------------
   USE GWFMNW1MODULE,ONLY:NWELL2,MXWEL2,IOWELL2,WELL2
   IMPLICIT NONE
   INTRINSIC CHAR, ABS
   INTEGER, EXTERNAL :: IFRL
   EXTERNAL IOWELLOUT
   ! Arguments
   INTEGER, INTENT(IN) :: Igrid
   ! Local Variables
   DOUBLE PRECISION :: hwell, conc, qt, qin, qout, q, timein, hcell
   DOUBLE PRECISION :: qsum, qwbar, timmult,timeinlast
   INTEGER :: i, icnt, io, iobynd, iopt, ioqsum, iostart, iot, istop
   INTEGER :: me, nb, ne, node
   CHARACTER(LEN=1) :: tab
   CHARACTER(LEN=32) :: temptag, tt, lasttag, eoftag
   !     ------------------------------------------------------------------
   CALL SGWF2MNW1PNT(IGRID)
   !
   tab = CHAR(9)
   iostart = 1000
   eoftag = 'EndOfFile__EndOfFile__EndOfFile_'
   !   Set Flag for printing header info once
   DO i = 1, MXWEL2
      WELL2(16, i) = 0   !! 16 = Header Flag
   ENDDO
   !
   !   Site file names are constructed to be OUTname_MNWSITE.txt
   !   All info for a well are dumped as tab delimited output
   !
   !------------------------------------------------------------------
   !     The 16 rows of the well array store:
   !      Row #  = Description
   !------------------------------------------------------------------
   !         1   = Well node locator
   !         2   = Net Discharge
   !         3   = Water level in wellbore
   !         4   = Water Quality
   !         5   = Qin  ---------------MNW ONLY--------
   !         6   = Qout
   !         7   = Q-node   / Net Discharge
   !         8   = MN-Flag --- Number of nodes / well
   !         9   = I/O unit for well output
   !        16   = Header Flag   Print= 0 / NoPrint = 1
   !------------------------------------------------------------------
   !
   !   Test auxillary output files for cleaning data sets
   iobynd = ABS(IOWELL2(2))
   IF ( iobynd.GT.0 ) THEN
      WRITE (iobynd, '(a32)') eoftag
      REWIND (iobynd)
      READ (iobynd, '(a)')
   ELSE
      RETURN
   ENDIF
   !
   ioqsum = ABS(IOWELL2(3))
   IF ( ioqsum.GT.0 ) THEN
      WRITE (ioqsum, '(a32)') eoftag
      REWIND (ioqsum)
      READ (ioqsum, '(a)')
   ELSE
      RETURN
   ENDIF
   !
   NWELL2 = 0
   icnt = 0
   istop = 0
   lasttag = 'NO-PRINT'
   timeinlast=-1.
   !
   DO WHILE ( istop.EQ.0 )
      READ (iobynd, '(a32,1x,2i8,6g15.8)') temptag, me, node, timein,&
      &q, hwell, hcell, conc
   !
   !   Test for output before accumulating INFO
      IF ( lasttag(1:8).NE.'NO-PRINT' .AND.&
      &(temptag.NE.lasttag .or.timein.ne.timeinlast)) THEN  !! Output
         iot = IFRL(WELL2(9, 1))
   !   Write a Header ?????
         iopt = iot - iostart
         IF ( IFRL(WELL2(16,iopt)).EQ.0 ) THEN
            WELL2(16, iopt) = 1
            IF ( icnt.GT.1 ) THEN
               WRITE (iot,&
               &'(a4,a1,a9,a1,a6,a1,a13,a1,a10,a1,a11,99(a1,a4,i7.7))'&
               &) 'TIME', tab, 'Discharge', tab, 'H-Well', tab,&
               &'Concentration', tab, 'Net-Inflow', tab,&
               &'Net-Outflow',&
               &(tab, 'Node', IFRL(WELL2(1,i)), i=1, icnt)
            ELSE
               WRITE (iot, '(a4,a1,a9,a1,a6,a1,a13)') 'TIME', tab,&
               &'Discharge', tab, 'H-Well', tab, 'Concentration'
            ENDIF
         ENDIF
   !   END of "Write a Header"  Section
         hwell = WELL2(3, icnt)
         conc = WELL2(4, icnt)
         IF ( icnt.GT.1 ) THEN
            qt = WELL2(7, 1)
            qin = WELL2(5, 1)
            qout = WELL2(6, 1)
            WRITE (iot, '(99(g15.8,a1))') WELL2(10, 1), tab, qt, tab,&
            &hwell, tab, conc, tab, qin, tab, qout,&
            &(tab, WELL2(2, i), i=1, icnt)
         ELSE
            qt = WELL2(2, 1)
            WRITE (iot, '(99(g15.8,a1))') WELL2(10, 1), tab, qt, tab,&
            &hwell, tab, conc
   !
         ENDIF
         icnt = 0     !! RESET node counter
      ENDIF

   !   Is this the EOF?
      IF ( temptag.EQ.eoftag ) istop = 1
   !
   !   Skip if this is a NO-PRINT node
      IF ( temptag(1:8).EQ.'NO-PRINT' .OR. temptag.EQ.eoftag ) THEN
         icnt = 0
      ELSE
         icnt = icnt + 1
   !   Identify pointer
         CALL IOWELLOUT(temptag, iostart, io)
         WELL2(1, icnt) = node      !!  1 = Well node locator
         WELL2(2, icnt) = q         !!  2 = Net Discharge
         WELL2(3, icnt) = hwell     !!  3 = Water level in wellbore
         WELL2(4, icnt) = conc      !!  4 = Water Quality
         WELL2(8, 1) = icnt         !!  8 = MN-Flag --- Number of nodes / well
         WELL2(9, 1) = io + iostart !!  9 = IO output
         WELL2(10, 1) = timein      !! 10 = Time
   !
   !     Read MN well output if available
         IF ( temptag.EQ.lasttag ) THEN  !! MN well
            IF ( icnt.EQ.2 .AND. ioqsum.GT.0 ) THEN
               READ (ioqsum, 9001) tt, nb, ne, timmult, qin, qout, qsum,&
               &hwell, qwbar
9001           FORMAT (A32, i6.5, i6.5, 12g15.8)
               WELL2(4, 1) = qwbar     !!  4 = Water Quality
               WELL2(5, 1) = qin       !!  5 = Qin  ---------------MNW ONLY--------
               WELL2(6, 1) = qout      !!  6 = Qout
               WELL2(7, 1) = qsum      !!  7 = Qnet
            ENDIF
         ENDIF
      ENDIF
   !
   !   Save Value of TEMPTAG for comparison
      lasttag = temptag
      timeinlast=timein
   ENDDO
   !
   !   Add IO close routine here if needed
   !
END SUBROUTINE GWF2MNW17OT
   !
   !     ******************************************************************
   !     ******************************************************************
SUBROUTINE IOWELLOUT(Temptag, Iostart, Io)
   USE GWFMNW1MODULE, ONLY : NWELL2,MNWSITE,MNWNAME
   IMPLICIT NONE
   ! Arguments
   INTEGER, INTENT(IN) :: Iostart
   INTEGER, INTENT(OUT) :: Io
   CHARACTER(LEN=32), INTENT(IN) :: Temptag
   !     ------------------------------------------------------------------
   ! Local Variables
   INTEGER :: i, k1, k2
   CHARACTER(LEN=256) :: txt
   !     ------------------------------------------------------------------
   i = 0
   Io = 0
   DO WHILE ( i.LT.NWELL2 )
      i = i + 1
      IF ( MNWSITE(i).EQ.Temptag ) THEN
         Io = i
         i = NWELL2
      ENDIF
   ENDDO
   !
   IF ( Io.EQ.0 ) THEN
      NWELL2 = NWELL2 + 1
      MNWSITE(NWELL2) = Temptag
      Io = NWELL2
   !   Build output file name and Open file
      k1 = INDEX(MNWNAME, ' ') - 1
      k2 = INDEX(Temptag, ' ') - 1
      IF ( k2.EQ.0 .AND. Temptag(32:32).NE.' ' ) k2 = 32
      txt = MNWNAME(1:k1)//'.'//Temptag(1:k2)//'.txt'
      OPEN (Io+Iostart, FILE=txt)
   ENDIF
   !
END SUBROUTINE IOWELLOUT
   !
   !_________________________________________________________________________________
   !
DOUBLE PRECISION FUNCTION CEL2WELBCF(Ix, Iy, Iz, Rw, Skin, Q, Cf)
   !     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
   !
   !----- MNW by K.J. Halford
   !
   !     ******************************************************************
   !     Compute conductance term to define head loss from cell to wellbore
   !      Methodology is described in full by Peaceman (1983)
   !swm: NOTE: MODIFIED FOR USE WITH THE BCF PACKAGE
   !     ******************************************************************
   !     Note: BCF, when LAYCON=0 or 2, does not save cell-by-cell
   !     Transmissivity (T) values.  Instead, it converts the cell-by-cell
   !     T values to branch conductances CR and CC, using harmonic
   !     averaging.  When BCF is used, the method used in this routine to
   !     generate cell-specific values of tx and ty is an approximation
   !     based on CR and CC.  When LPF or HUF is used, cell-by-cell
   !     hydraulic-conductivity values are stored, this approximation is
   !     not needed, and the values generated for tx and ty are exact --
   !     ERB 1/29/01.
   !
   !        SPECIFICATIONS::
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,LAYHDT,DELR,DELC,BOTM,LBOTM,CR,CC,&
   &HNEW
   USE GWFBASMODULE,ONLY:HDRY
   USE GWFBCFMODULE,ONLY:HY,TRPY,LAYCON
   USE GWFMNW1MODULE,ONLY:PLOSS,SMALL,TWOPI,ZERO25
   !
   IMPLICIT NONE
   INTRINSIC LOG, ABS, SQRT
   ! Arguments
   INTEGER, INTENT(IN) :: Ix, Iy, Iz
   DOUBLE PRECISION, INTENT(IN) :: Rw, Skin, Q, Cf
   ! Local Variables
   DOUBLE PRECISION :: ro, ah, tempKX, dx, dy, top, bot, dxp, dxm
   DOUBLE PRECISION :: txm, txp, dyp, dym, typ, tym, div, txx, tyy
   DOUBLE PRECISION :: thick, upper, yx4, xy4, tpi2, a, b, c, cel2wel
   !     ------------------------------------------------------------------
   !1000 FORMAT(/1X,
   !    &'***ERROR: MNW PACKAGE DOES NOT SUPPORT HEAD-DEPENDENT',/,
   !    &' THICKNESS OPTION OF SELECTED FLOW PACKAGE',/,
   !    &' (MNW DOES FULLY SUPPORT BCF, LPF, AND HUF PACKAGES)',/,
   !    &' -- STOP EXECUTION (CEL2WEL)')
   !
   dx = DELR(Ix)
   dy = DELC(Iy)
   top = BOTM(Ix, Iy, LBOTM(Iz)-1)
   bot = BOTM(Ix, Iy, LBOTM(Iz))
   !
   !     FIND HORIZONTAL ANISOTROPY, THE RATIO Ky/Kx
   ah = TRPY(Iz)
   !
   IF ( LAYHDT(Iz).EQ.0 ) THEN
   !       THICKNESS IS NOT HEAD-DEPENDENT
      dxp = dx
      txp = 0.0D0
      IF ( Ix.LT.NCOL ) THEN
         dxp = DELR(Ix+1)
         txp = CR(Ix, Iy, Iz)*(dx+dxp) / 2.0D0
      ENDIF
      dxm = dx
      txm = txp
      IF ( Ix.GT.1  ) THEN
         dxm = DELR(Ix-1)
         txm = CR(Ix-1, Iy, Iz)*(dx+dxm) / 2.0D0
      ENDIF
      IF ( txp.LT.SMALL ) txp = txm
      IF ( txm.LT.SMALL ) txm = txp
   !
      dyp = dy
      typ = 0.0D0
      IF ( Iy.LT.NROW ) THEN
         dyp = DELC(Iy+1)
         typ = CC(Ix, Iy, Iz)*(dy+dyp) / 2.0D0
      ENDIF
      dym = dy
      tym = typ
      IF ( Iy.GT.1 ) THEN
         dym = DELC(Iy-1)
         tym = CC(Ix, Iy-1, Iz)*(dy+dym) / 2.0D0
      ENDIF
      IF ( typ.LT.SMALL ) typ = tym
      IF ( tym.LT.SMALL ) tym = typ
      txp = txp / dy
      txm = txm / dy
      typ = typ / dx
      tym = tym / dx
   !
   !  Eliminate zero values .....
   !
      IF ( typ.LT.SMALL .OR. NROW.LT.2 ) THEN
         typ = txp
         tym = txm
      ENDIF
   !
      IF ( txp.LT.SMALL .OR. NCOL.LT.2 ) THEN
         txp = typ
         txm = tym
      ENDIF
   !
   !   Assuming expansion of grid is slight, if present, & that txx and tyy of the adjacent
   !  cells are about the same value.
      txx = 0.0D0
      div = txp + txm
      IF ( div.GT.SMALL ) txx = 2.0D0*txp*txm / div
      tyy = 0.0D0
      div = typ + tym
      IF ( div.GT.SMALL ) tyy = 2.0D0*typ*tym / div
      IF ( txx.GT.SMALL .AND. tyy.LT.SMALL ) tyy = txx
      IF ( tyy.GT.SMALL .AND. txx.LT.SMALL ) txx = tyy
   ELSE
   !       THICKNESS IS HEAD-DEPENDENT
   !  Estimate T to well in an unconfined system
   !
      upper = HNEW(Ix, Iy, Iz)
      tempKX = HY(Ix, Iy, Iz)        !!  BCF Hydraulic Conductivity array
      IF ( LAYCON(Iz).EQ.3 ) THEN
         IF ( upper.GT.top ) upper = top
      ENDIF
      thick = upper - bot
   !   set thickness / conductance to 0 if cell is dry
      IF ( (HNEW(Ix,Iy,Iz)-HDRY)**2.0D0.LT.ZERO25 ) thick = 0.0D0
      txx = tempKX*thick
      IF ( txx.LT.ZERO25 ) txx = 0.0D0
      tyy = txx*ah
   ENDIF
   !
   IF ( Rw.LT.ZERO25 .OR. txx.LT.ZERO25 .OR. tyy.LT.ZERO25 ) THEN
      cel2wel = SQRT( txx*tyy )
   ELSE
      yx4 = SQRT(SQRT(tyy/txx))
      xy4 = SQRT(SQRT(txx/tyy))
      ro = 0.28D0*SQRT((yx4*dx)**2.0D0 + (xy4*dy)**2.0D0) / (yx4+xy4)
   !
      tpi2 = TWOPI*SQRT(txx*tyy)
      a = LOG(ro/Rw) / tpi2
      IF ( PLOSS.GT.0.99D0 ) THEN
         b = Skin
         c = Cf*ABS(Q)**(PLOSS-1.0D0)
      ELSE
         b = Skin / tpi2
         c = 0.0D0
      ENDIF
      cel2wel = 1.0D0 / (a + b + c)
   ENDIF
   !
   CEL2WELBCF = cel2wel
END FUNCTION CEL2WELBCF
   !
   !_________________________________________________________________________________
   !
DOUBLE PRECISION FUNCTION CEL2WELLPF(Ix, Iy, Iz, Rw, Skin, Q, Cf)
   !     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
   !
   !----- MNW by K.J. Halford
   !
   !     ******************************************************************
   !     Compute conductance term to define head loss from cell to wellbore
   !      Methodology is described in full by Peaceman (1983)
   !swm: NOTE: MODIFIED FOR USE WITH THE LPF PACKAGE
   !     ******************************************************************
   !     Note: BCF, when LAYCON=0 or 2, does not save cell-by-cell
   !     Transmissivity (T) values.  Instead, it converts the cell-by-cell
   !     T values to branch conductances CR and CC, using harmonic
   !     averaging.  When BCF is used, the method used in this routine to
   !     generate cell-specific values of tx and ty is an approximation
   !     based on CR and CC.  When LPF or HUF is used, cell-by-cell
   !     hydraulic-conductivity values are stored, this approximation is
   !     not needed, and the values generated for tx and ty are exact --
   !     ERB 1/29/01.
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:LAYHDT,DELR,DELC,BOTM,LBOTM,HNEW
   USE GWFBASMODULE,ONLY:HDRY
   USE GWFLPFMODULE,ONLY:HK,CHANI,HANI
   USE GWFMNW1MODULE,ONLY:PLOSS,TWOPI,ZERO25
   IMPLICIT NONE
   INTRINSIC LOG, ABS, SQRT
   ! Arguments
   INTEGER, INTENT(IN) :: Ix, Iy, Iz
   DOUBLE PRECISION, INTENT(IN) :: Rw, Skin, Q, Cf
   ! Local Variables
   DOUBLE PRECISION :: ro, ah, tempKX, dx, dy, top, bot, txx, tyy
   DOUBLE PRECISION :: thick, upper, yx4, xy4, tpi2, a, b, c, cel2wel
   !     ------------------------------------------------------------------
   !1000 FORMAT(/1X,
   !    &'***ERROR: MNW PACKAGE DOES NOT SUPPORT HEAD-DEPENDENT',/,
   !    &' THICKNESS OPTION OF SELECTED FLOW PACKAGE',/,
   !    &' (MNW DOES FULLY SUPPORT BCF, LPF, AND HUF PACKAGES)',/,
   !    &' -- STOP EXECUTION (CEL2WEL)')
   !
   dx = DELR(Ix)
   dy = DELC(Iy)
   top = BOTM(Ix, Iy, LBOTM(Iz)-1)
   bot = BOTM(Ix, Iy, LBOTM(Iz))
   !
   !     FIND HORIZONTAL ANISOTROPY, THE RATIO Ky/Kx
   IF ( CHANI(Iz).GT.0.0 ) THEN
      ah = CHANI(Iz)
   ELSE
      ah = HANI(Ix, Iy, Iz)
   ENDIF
   !
   IF ( LAYHDT(Iz).EQ.0 ) THEN
   !       THICKNESS IS NOT HEAD-DEPENDENT
      thick = top - bot
      txx = HK(Ix, Iy, Iz)*thick
      tyy = txx*ah
   ELSE
   !       THICKNESS IS HEAD-DEPENDENT
   !  Estimate T to well in an unconfined system
   !
      upper = HNEW(Ix, Iy, Iz)
      IF ( upper.GT.top ) upper = top
      tempKX = HK(Ix, Iy, Iz)        !!  LPF Hydraulic Conductivity array
      thick = upper - bot
   !   set thickness / conductance to 0 if cell is dry
      IF ( (HNEW(Ix,Iy,Iz)-HDRY)**2.0D0.LT.ZERO25 ) thick = 0.0D0
      txx = tempKX*thick
      IF ( txx.LT.ZERO25 ) txx = 0.0D0
      tyy = txx*ah
   ENDIF
   !
   IF ( rw.LT.ZERO25 .OR. txx.LT.ZERO25 .OR. tyy.LT.ZERO25 ) THEN
      cel2wel = SQRT(txx*tyy)
   ELSE
      yx4 = SQRT(SQRT(tyy/txx))
      xy4 = SQRT(SQRT(txx/tyy))
      ro = 0.28D0*SQRT((yx4*dx)**2.0D0 +(xy4*dy)**2.0D0) / (yx4+xy4)
   !
      tpi2 = TWOPI*SQRT(txx*tyy)
      a = LOG(ro/rw) / tpi2
      IF ( PLOSS.GT.0.99D0 ) THEN
         b = Skin
         c = Cf*ABS(Q)**(PLOSS-1.0D0)
      ELSE
         b = Skin / tpi2
         c = 0.0D0
      ENDIF
      cel2wel = 1.0D0 / (a + b + c)
   ENDIF
   !
   CEL2WELLPF = cel2wel
END FUNCTION CEL2WELLPF
   !
   !_________________________________________________________________________________
   !
DOUBLE PRECISION FUNCTION CEL2WELHUF(Ix, Iy, Iz, Rw, Skin, Q, Cf)
   !     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
   !
   !----- MNW by K.J. Halford
   !
   !     ******************************************************************
   !     Compute conductance term to define head loss from cell to wellbore
   !      Methodology is described in full by Peaceman (1983)
   !swm: NOTE: MODIFIED FOR USE WITH THE HUF PACKAGE
   !     ******************************************************************
   !     Note: BCF, when LAYCON=0 or 2, does not save cell-by-cell
   !     Transmissivity (T) values.  Instead, it converts the cell-by-cell
   !     T values to branch conductances CR and CC, using harmonic
   !     averaging.  When BCF is used, the method used in this routine to
   !     generate cell-specific values of tx and ty is an approximation
   !     based on CR and CC.  When LPF or HUF is used, cell-by-cell
   !     hydraulic-conductivity values are stored, this approximation is
   !     not needed, and the values generated for tx and ty are exact --
   !     ERB 1/29/01.
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:LAYHDT,DELR,DELC,BOTM,LBOTM,HNEW
   USE GWFBASMODULE,ONLY:HDRY
   USE GWFHUFMODULE,ONLY:HK,HKCC
   USE GWFMNW1MODULE,ONLY:PLOSS,TWOPI,ZERO25
   IMPLICIT NONE
   INTRINSIC LOG, ABS, SQRT
   ! Arguments
   INTEGER, INTENT(IN) :: Ix, Iy, Iz
   DOUBLE PRECISION, INTENT(IN) :: Rw, Skin, Q, Cf
   ! Local Variables
   DOUBLE PRECISION :: ro, ah, tempKX, dx, dy, top, bot, txx, tyy, ky
   DOUBLE PRECISION :: thick, upper, yx4, xy4, tpi2, a, b, c, cel2wel
   !     ------------------------------------------------------------------
   !1000 FORMAT(/1X,
   !    &'***ERROR: MNW PACKAGE DOES NOT SUPPORT HEAD-DEPENDENT',/,
   !    &' THICKNESS OPTION OF SELECTED FLOW PACKAGE',/,
   !    &' (MNW DOES FULLY SUPPORT BCF, LPF, AND HUF PACKAGES)',/,
   !    &' -- STOP EXECUTION (CEL2WEL)')
   !
   dx = DELR(Ix)
   dy = DELC(Iy)
   top = BOTM(Ix, Iy, LBOTM(Iz)-1)
   bot = BOTM(Ix, Iy, LBOTM(Iz))
   !
   !     FIND HORIZONTAL ANISOTROPY, THE RATIO Ky/Kx
   tempKX = HK(Ix, Iy, Iz)
   ky = HKCC(Ix, Iy, Iz)
   ah = ky/tempKX
   !
   IF ( LAYHDT(Iz).EQ.0 ) THEN
   !       THICKNESS IS NOT HEAD-DEPENDENT
      thick = top - bot
      txx = HK(Ix, Iy, Iz)*thick
      tyy = txx*ah
   ELSE
   !       THICKNESS IS HEAD-DEPENDENT
   !  Estimate T to well in an unconfined system
   !
      upper = hnew(Ix, Iy, Iz)
      IF ( upper.GT.top ) upper = top
      tempKX = hk(Ix, Iy, Iz)        !!  HUF Hydraulic Conductivity array
      thick = upper - bot
   !   set thickness / conductance to 0 if cell is dry
      IF ( (hnew(Ix, Iy, Iz)-HDRY )**2.0D0.LT.ZERO25 ) thick = 0.0D0
      txx = tempKX*thick
      IF ( txx.LT.ZERO25 ) txx = 0.0D0
      tyy = txx*ah
   ENDIF
   !
   IF ( rw.LT.ZERO25 .OR. txx.LT.ZERO25 .OR. tyy.LT.ZERO25 ) THEN
      cel2wel = SQRT( txx*tyy )
   ELSE
      yx4 = SQRT(SQRT(tyy/txx))
      xy4 = SQRT(SQRT(txx/tyy))
      ro = 0.28D0*SQRT((yx4*dx)**2.0D0 + (xy4*dy)**2.0D0) / (yx4+xy4)
   !
      tpi2 = TWOPI*SQRT(txx*tyy)
      a = LOG(ro/rw) / tpi2
      IF ( PLOSS.GT.0.99D0 ) THEN
         b = Skin
         c = Cf*ABS(Q)**(PLOSS-1.0D0)
      ELSE
         b = Skin / tpi2
         c = 0.0D0
      ENDIF
      cel2wel = 1.0D0 / (a + b + c)
   ENDIF
   !
   CEL2WELHUF = cel2wel
END FUNCTION CEL2WELHUF
   !
   !_________________________________________________________________________________
   !
DOUBLE PRECISION FUNCTION CEL2WELUPW(Ix, Iy, Iz, Rw, Skin, Q, Cf)
   !     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
   !
   !----- MNW by K.J. Halford
   !
   !     ******************************************************************
   !     Compute conductance term to define head loss from cell to wellbore
   !      Methodology is described in full by Peaceman (1983)
   !swm: NOTE: MODIFIED FOR USE WITH THE LPF PACKAGE
   !     ******************************************************************
   !     Note: BCF, when LAYCON=0 or 2, does not save cell-by-cell
   !     Transmissivity (T) values.  Instead, it converts the cell-by-cell
   !     T values to branch conductances CR and CC, using harmonic
   !     averaging.  When BCF is used, the method used in this routine to
   !     generate cell-specific values of tx and ty is an approximation
   !     based on CR and CC.  When LPF or HUF is used, cell-by-cell
   !     hydraulic-conductivity values are stored, this approximation is
   !     not needed, and the values generated for tx and ty are exact --
   !     ERB 1/29/01.
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:LAYHDT,DELR,DELC,BOTM,LBOTM,HNEW
   USE GWFBASMODULE,ONLY:HDRY
   USE GWFUPWMODULE,ONLY:HKUPW,CHANI,HANI
   USE GWFMNW1MODULE,ONLY:PLOSS,TWOPI,ZERO25
   IMPLICIT NONE
   INTRINSIC LOG, ABS, SQRT
   ! Arguments
   INTEGER, INTENT(IN) :: Ix, Iy, Iz
   DOUBLE PRECISION, INTENT(IN) :: Rw, Skin, Q, Cf
   ! Local Variables
   DOUBLE PRECISION :: ro, ah, tempKX, dx, dy, top, bot, txx, tyy
   DOUBLE PRECISION :: thick, upper, yx4, xy4, tpi2, a, b, c, cel2wel
   !     ------------------------------------------------------------------
   !1000 FORMAT(/1X,
   !    &'***ERROR: MNW PACKAGE DOES NOT SUPPORT HEAD-DEPENDENT',/,
   !    &' THICKNESS OPTION OF SELECTED FLOW PACKAGE',/,
   !    &' (MNW DOES FULLY SUPPORT BCF, LPF, AND HUF PACKAGES)',/,
   !    &' -- STOP EXECUTION (CEL2WEL)')
   !
   dx = DELR(Ix)
   dy = DELC(Iy)
   top = BOTM(Ix, Iy, LBOTM(Iz)-1)
   bot = BOTM(Ix, Iy, LBOTM(Iz))
   !
   !     FIND HORIZONTAL ANISOTROPY, THE RATIO Ky/Kx
   IF ( CHANI(Iz).GT.0.0 ) THEN
      ah = CHANI(Iz)
   ELSE
      ah = HANI(Ix, Iy, Iz)
   ENDIF
   !
   IF ( LAYHDT(Iz).EQ.0 ) THEN
   !       THICKNESS IS NOT HEAD-DEPENDENT
      thick = top - bot
      txx = HKUPW(Ix, Iy, Iz)*thick
      tyy = txx*ah
   ELSE
   !       THICKNESS IS HEAD-DEPENDENT
   !  Estimate T to well in an unconfined system
   !
      upper = HNEW(Ix, Iy, Iz)
      IF ( upper.GT.top ) upper = top
      tempKX = HKUPW(Ix, Iy, Iz)        !!  UPW Hydraulic Conductivity array
      thick = upper - bot
      IF ( thick.LT.0.0 ) thick = 0.0   ! RGN NWT can have heads below bottom
   !   set thickness / conductance to 0 if cell is dry
      IF ( (HNEW(Ix,Iy,Iz)-HDRY)**2.0D0.LT.ZERO25 ) thick = 0.0D0
      txx = tempKX*thick
      IF ( txx.LT.ZERO25 ) txx = 0.0D0
      tyy = txx*ah
   ENDIF
   !
   IF ( rw.LT.ZERO25 .OR. txx.LT.ZERO25 .OR. tyy.LT.ZERO25 ) THEN
      cel2wel = SQRT(txx*tyy)
   ELSE
      yx4 = SQRT(SQRT(tyy/txx))
      xy4 = SQRT(SQRT(txx/tyy))
      ro = 0.28D0*SQRT((yx4*dx)**2.0D0 +(xy4*dy)**2.0D0) / (yx4+xy4)
   !
      tpi2 = TWOPI*SQRT(txx*tyy)
      a = LOG(ro/rw) / tpi2
      IF ( PLOSS.GT.0.99D0 ) THEN
         b = Skin
         c = Cf*ABS(Q)**(PLOSS-1.0D0)
      ELSE
         b = Skin / tpi2
         c = 0.0D0
      ENDIF
      cel2wel = 1.0D0 / (a + b + c)
   ENDIF
   !
   CEL2WELUPW = cel2wel
END FUNCTION CEL2WELUPW
   !
   !     ******************************************************************
   !     Define direction of pointer along a row, column, or layer
   !     ******************************************************************
INTEGER FUNCTION IDIRECT(N1, N2, Ncol, Nrow)
   IMPLICIT NONE
   INTRINSIC ABS
   ! Arguments
   INTEGER, INTENT(IN) :: N1, N2, Ncol, Nrow
   !     ------------------------------------------------------------------
   IDIRECT = Ncol
   IF ( ABS(N2-N1).GT.Ncol*Nrow ) IDIRECT = Ncol*Nrow
   IF ( ABS(N2-N1).LT.Ncol ) IDIRECT = 1
   IF ( N2.LT.N1 ) IDIRECT = -IDIRECT
   !
END FUNCTION IDIRECT
   !
   !     ******************************************************************
   !     ******************************************************************
   !     INTEGER FUNCTION ITXEND(Txt)
   !     IMPLICIT NONE
   ! Arguments
   !     CHARACTER(LEN=256), INTENT(IN) :: Txt
   ! Local Variables
   !     INTEGER :: k
   !     ------------------------------------------------------------------
   !     k = 256
   !     DO WHILE ( Txt(k:k).EQ.' ' .AND. k.GT.1 )
   !       k = k - 1
   !     ENDDO
   !     ITXEND = k
   !
   !     END FUNCTION ITXEND
   !
   !     ******************************************************************
   !     ******************************************************************
INTEGER FUNCTION IFRL(R)
   IMPLICIT NONE
   INTRINSIC ABS
   ! Arguments
   DOUBLE PRECISION, INTENT(IN) :: R
   ! Local Variables
   INTEGER :: ip
   !     ------------------------------------------------------------------
   ip = ABS(R) + 0.5D0
   IF ( R.LT.0.0D0 ) ip = -ip
   IFRL = ip
END FUNCTION IFRL
   !
   !     ******************************************************************
   !     NCREAD: reads lines of input and ignores lines that begin with a "#" sign.
   !          All information after a ! is wiped from the input card.
   !     ******************************************************************
SUBROUTINE NCREAD(Io, Txt, Ierr)
   IMPLICIT NONE
   EXTERNAL UPCASE, USTOP
   ! Arguments
   INTEGER, INTENT(INOUT) :: Io
   INTEGER, INTENT(OUT) :: Ierr
   CHARACTER(LEN=256), INTENT(OUT) :: Txt
   ! Local Variables
   INTEGER :: ioalt, ioflip, iohold, ki
   CHARACTER(LEN=128) :: afile
   CHARACTER(LEN=256) :: tx2
   DATA ioflip, ioalt/69, 69/
   !     ------------------------------------------------------------------
   Ierr = 0
5  READ (Io, '(a)', END=10) Txt
   IF ( Txt(1:1).EQ.'#' ) GOTO 5
   !
   ki = INDEX(Txt, '!')
   IF ( ki.GT.0 )&
   &Txt(ki:256) = '                                                '
   !
   tx2 = Txt
   CALL UPCASE(tx2)
   !
   !    Test for switching control to an auxillary input file
   !
   ki = INDEX(Txt, ':')
   IF ( INDEX(tx2, 'REDIRECT').GT.0 .AND. ki.GT.0 ) THEN
      afile = Txt(ki+1:256)
      ki = INDEX(afile, '  ') - 1
      iohold = Io
      Io = ioflip
      ioflip = iohold
      OPEN (Io, FILE=afile(1:ki), STATUS='OLD', ERR=20)
      GOTO 5
   ENDIF
   !
   !    Test for returning io control from auxillary input to master input file
   !
   IF ( INDEX(tx2, 'RETURN').GT.0 .AND.&
   &INDEX(tx2, 'CONTROL').GT.0 ) GOTO 10
   !
   ki = INDEX(tx2, '<END>')
   IF ( ki.GT.0 ) THEN
      Ierr = 1
      Txt(ki+5:256) = '                                           '
   ENDIF
   !
   IF ( INDEX(tx2, '<STOP>').GT.0 ) Ierr = 2
   RETURN
   !
   !    Report error in opening auxillary input file and stop
   !
20 WRITE (*, 25) afile
25 FORMAT (/, '  ERROR opening auxillary input file', //,&
   &'   The file:  ', a40, ' does not exist', /)
   !
   !     When compiling MNW with Modflow-96, comment out the call to
   !     USTOP and uncomment the STOP statement
   CALL USTOP(' ')
   !      STOP
   !
10 Txt(1:3) = 'EOF'
   IF ( Io.EQ.ioalt ) THEN
      CLOSE (Io)
      iohold = Io
      Io = ioflip
      ioflip = iohold
      GOTO 5
   ELSE
      Ierr = -1
   ENDIF
   !
END SUBROUTINE NCREAD
   !
   !     ******************************************************************
   !     ******************************************************************
SUBROUTINE QREAD(R, Ni, Ain, Ierr)
   IMPLICIT NONE
   INTRINSIC CHAR, INDEX
   INTEGER, PARAMETER :: MRNV=25
   ! Arguments
   DOUBLE PRECISION, INTENT(OUT), DIMENSION(MRNV) :: R
   INTEGER, INTENT(IN) :: Ni
   INTEGER, INTENT(OUT) :: Ierr
   CHARACTER(LEN=256), INTENT(IN) :: Ain
   ! Local Variables
   INTEGER :: i, istat, ki, n, nd
   CHARACTER(LEN=1) :: tab
   CHARACTER(LEN=8) :: rdfmt
   CHARACTER(LEN=256) :: a256
   !     ------------------------------------------------------------------
   Ierr = 0
   tab = CHAR(9)           ! sets tab delimiter
   !
   !   r(ni+1) records the number of non-numeric entries that were attempted to be read as a number
   !   r(ni+2) records the last column that was read from the card
   !
   R(Ni+1) = -1.0D0
   a256 = Ain
   DO i = 1, 256
      IF ( a256(i:i).EQ.tab ) a256(i:i) = ' '
      IF ( a256(i:i).EQ.',' ) a256(i:i) = ' '
      IF ( a256(i:i).EQ.':' ) a256(i:i) = ' '
      IF ( a256(i:i).EQ.'=' ) a256(i:i) = ' '
   ENDDO
   n = 1
   i = 0
11 R(Ni+1) = R(Ni+1) + 1.0D0
10 i = i + 1
   IF ( i.GE.256 ) GOTO 15
   IF ( a256(i:i).EQ.' ' ) THEN
      a256(i:i) = '?'
      GOTO 10
   ENDIF
   !
   ki = INDEX(a256, ' ') - 1
   nd = ki - i + 1
   rdfmt = '(F??.0) '
   WRITE (rdfmt(3:4), '(i2.2)') nd
   !ERB  Fix for bug that caused i to be incremented by only 1 position
   !ERB  each time the read statement returns an error.  This bug also
   !ERB  incremented r(ni+1) unnecessarily.  With Lahey-compiled code, the
   !ERB  buggy version would read a final E in a word (without returning an
   !ERB  error) as a zero.
   !ERB      read (a256(i:ki),rdfmt,err=11,end=10) r(n)
   READ (a256(i:ki), rdfmt, ERR=13, IOSTAT=istat) R(n)
13 CONTINUE
   i = ki
   IF ( istat.GT.0 ) GOTO 11 ! PART OF BUG FIX -- ERB
   n = n + 1
   IF ( n.LE.Ni .AND. i.LT.256 ) GOTO 10
   !
15 n = n - 1
   Ierr = Ni - n
   R(Ni+2) = i
   !
END SUBROUTINE QREAD
   !***********************************************************************
SUBROUTINE GWF2MNW17DA(Igrid)
   !     ******************************************************************
   !     DEALLOCATE MNW DATA
   !     ******************************************************************
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFMNW1MODULE
   !     ------------------------------------------------------------------
   ! Arguments
   INTEGER :: Igrid, IDUM
   !     ------------------------------------------------------------------
   IDUM=IGRID
   DEALLOCATE (GWFMNWDAT(Igrid)%NWELL2)
   DEALLOCATE (GWFMNWDAT(Igrid)%MXWEL2)
   DEALLOCATE (GWFMNWDAT(Igrid)%IWL2CB)
   DEALLOCATE (GWFMNWDAT(Igrid)%NOMOITER)
   DEALLOCATE (GWFMNWDAT(Igrid)%KSPREF)
   DEALLOCATE (GWFMNWDAT(Igrid)%IWELPT)
   DEALLOCATE (GWFMNWDAT(Igrid)%PLOSS)
   DEALLOCATE (GWFMNWDAT(Igrid)%SMALL)
   DEALLOCATE (GWFMNWDAT(Igrid)%HMAX)
   DEALLOCATE (GWFMNWDAT(Igrid)%MNWNAME)
   DEALLOCATE (GWFMNWDAT(Igrid)%IOWELL2)
   DEALLOCATE (GWFMNWDAT(Igrid)%HREF)
   DEALLOCATE (GWFMNWDAT(Igrid)%WELL2)
   DEALLOCATE (GWFMNWDAT(Igrid)%MNWSITE)
   !
END SUBROUTINE GWF2MNW17DA
   !***********************************************************************
SUBROUTINE SGWF2MNW1PNT(Igrid)
   !     ******************************************************************
   !     SET MNW POINTER DATA TO CURRENT GRID
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFMNW1MODULE
   !     ------------------------------------------------------------------
   ! Arguments
   INTEGER :: Igrid
   !     ------------------------------------------------------------------
   NWELL2=>GWFMNWDAT(Igrid)%NWELL2
   MXWEL2=>GWFMNWDAT(Igrid)%MXWEL2
   IWL2CB=>GWFMNWDAT(Igrid)%IWL2CB
   NOMOITER=>GWFMNWDAT(Igrid)%NOMOITER
   KSPREF=>GWFMNWDAT(Igrid)%KSPREF
   IWELPT=>GWFMNWDAT(Igrid)%IWELPT
   PLOSS=>GWFMNWDAT(Igrid)%PLOSS
   SMALL=>GWFMNWDAT(Igrid)%SMALL
   HMAX=>GWFMNWDAT(Igrid)%HMAX
   MNWNAME=>GWFMNWDAT(Igrid)%MNWNAME
   IOWELL2=>GWFMNWDAT(Igrid)%IOWELL2
   HREF=>GWFMNWDAT(Igrid)%HREF
   WELL2=>GWFMNWDAT(Igrid)%WELL2
   MNWSITE=>GWFMNWDAT(Igrid)%MNWSITE
   !
END SUBROUTINE SGWF2MNW1PNT
   !***********************************************************************
SUBROUTINE SGWF2MNW1PSV(Igrid)
   !     ******************************************************************
   !     SAVE MNW POINTER DATA
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFMNW1MODULE
   !     ------------------------------------------------------------------
   ! Arguments
   INTEGER :: Igrid
   !     ------------------------------------------------------------------
   GWFMNWDAT(Igrid)%NWELL2=>NWELL2
   GWFMNWDAT(Igrid)%MXWEL2=>MXWEL2
   GWFMNWDAT(Igrid)%IWL2CB=>IWL2CB
   GWFMNWDAT(Igrid)%NOMOITER=>NOMOITER
   GWFMNWDAT(Igrid)%KSPREF=>KSPREF
   GWFMNWDAT(Igrid)%IWELPT=>IWELPT
   GWFMNWDAT(Igrid)%PLOSS=>PLOSS
   GWFMNWDAT(Igrid)%SMALL=>SMALL
   GWFMNWDAT(Igrid)%HMAX=>HMAX
   GWFMNWDAT(Igrid)%MNWNAME=>MNWNAME
   GWFMNWDAT(Igrid)%IOWELL2=>IOWELL2
   GWFMNWDAT(Igrid)%HREF=>HREF
   GWFMNWDAT(Igrid)%WELL2=>WELL2
   GWFMNWDAT(Igrid)%MNWSITE=>MNWSITE
   !
END SUBROUTINE SGWF2MNW1PSV
