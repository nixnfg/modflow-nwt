   ! Time of File Save by ERB: 3/23/2004 1:57PM
   !rgn&dep Revised method of computing lake depth March-August 2009.
   !rgn&dep Previously depth computed using surface area of lake at
   !rgn&dep beginning of time step but the method neglected the change
   !rgn&dep in lake area caused by a change over the time step.
   !rgn&dep This could cause substantial lake mass errors. Added explicit
   !rgn&dep lake stage calculations to determine lake seepage.
   !rgn  Made EVAP, PRECIP, SEEP, and SEEP3 double precision Nov. 6, 2006
   !dep  Converted from MODFLOW-2000 to MODFLOW-2005 May 2006, RGN and DEP
   !dep  Lake Package modified by DEP and RGN May 20 through June 29, 2006
   !dep  to compute lake outflow as a function of lake stage inside the
   !dep  FORMULATE MODULE. Lake outflow had been previously computed in the
   !dep  Streamflow-Routing Package. The Streamflow-Routing (sfr7) Package
   !dep  was also modified to remain compatible with the modifications in
   !dep  the Lake Package.
   !     Modifications made February and March 21, 2004; DEP
   !     Last change:  MLM & LFK  10 Oct 2003;  LFK 21 Jan 2004
   !     Previous change:  ERB  13 Sep 2002    9:22 am
   !
   !
SUBROUTINE GWF2LAK7AR(IN,IUNITSFR,IUNITGWT,IUNITUZF,NSOL,IGRID)
   !
   !------OLD USGS VERSION 7.1; JUNE 2006 GWF2LAK7AR;
   !------UPDATED FOR MF-2005, FEBRUARY 6, 2012
   !rgn------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
   !rgn------NEW VERSION NUMBER FOR NWT 1.3.0, 7/01/2022
   !     ******************************************************************
   !     INITIALIZE POINTER VARIABLES USED BY SFR1 TO SUPPORT LAKE3 AND
   !     GAGE PACKAGES AND THE GWT PROCESS
   !     ******************************************************************
   !
   USE GWFLAKMODULE
   USE GLOBAL,       ONLY: IOUT, NCOL, NROW, NLAY, IFREFM, ITRSS,&
   &NODES, IUNIT, ITMUNI       !EDM
   USE GWFSFRMODULE, ONLY: NSS
   !
   !      ******************************************************************
   !      ALLOCATE ARRAY STORAGE FOR LAKES
   !      ******************************************************************
   !
   !      ------------------------------------------------------------------
   !      SPECIFICATIONS:
   CHARACTER (LEN=40):: CARD
   CHARACTER*200 line
   !      ------------------------------------------------------------------
   !rsr  Allocate lake variables used by SFR even if lakes not active so that
   !       argument lists are defined
   ALLOCATE (NLAKES, NLAKESAR,THETA,LAKUNIT,NSFRLAK,NLKFLWTYP)       !EDM
   ALLOCATE (IGSFLOWLAK)
   IGSFLOWLAK = 0
   ALLOCATE(LKFLOWTYPE(6)) ! POSITION 1: STORAGE; 2: DELVOL; 3: PRECIP; 4: EVAP; 5: RUNOFF; 6: WITHDRAWL
   !
   !--REINITIALIZE LKFLOWTYPE WITH EACH STRESS PERIOD
   NLKFLWTYP=0
   IF(IUNIT(49).NE.0) THEN
      LKFLOWTYPE(1)='NA'
      LKFLOWTYPE(2)='NA'
      LKFLOWTYPE(3)='NA'
      LKFLOWTYPE(4)='NA'
      LKFLOWTYPE(5)='NA'
      LKFLOWTYPE(6)='NA'
   ENDIF
   !
   NLAKES = 0
   LAKUNIT = IN
   NLAKESAR = 1
   THETA = 0.0
   IF (IN.GT.0) THEN
   !dep added SURFDEPTH 3/3/2009
      ALLOCATE (ILKCB, NSSITR, SSCNCR, SURFDEPTH,RAMP,SMALLTOL)
      ALLOCATE (MXLKND, LKNODE, ICMX, NCLS, LWRT, NDV, NTRB)
      ALLOCATE (IRDTAB,ISTARTLAK)
   !
   !1------IDENTIFY PACKAGE AND INITIALIZE LKNODE.
      WRITE(IOUT,1) IN
      LKNODE=0
   !dep  initialize number of iterations and closure criteria to zero.
      DUM = 0.0
      NSSITR = 0
      SSCNCR = 0.0
      SURFDEPTH = 0.0
      RAMP = 1.0
      SMALLTOL = 1.0d-4
      IF ( ITMUNI == 1 ) THEN
         SMALLTOL = SMALLTOL/86400.
      ELSE IF ( ITMUNI == 2) THEN
         SMALLTOL = SMALLTOL/1440.
      ELSE IF ( ITMUNI == 3) THEN
         SMALLTOL = SMALLTOL/24.
      ELSE IF ( ITMUNI == 5 ) THEN
         SMALLTOL = SMALLTOL*365.
      END IF
   !
      lloc = 1
      IRDTAB = 0
      ISTARTLAK = 0
      NPP = 0
      MXVL = 0
      CALL URDCOM(In, IOUT, line)
   ! Check for alternate option to specifiy stage/vol/area tables.
      CALL UPARLSTAL(IN,IOUT,LINE,NPP,MXVL)
      lloc = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'TABLEINPUT') THEN
         IRDTAB = 1
         WRITE(IOUT,32)
32       FORMAT(1X,I10,' Stage, volume and area relationship specified ',&
         &'based on an external tabular input file')
      ELSE
         BACKSPACE IN
         WRITE(IOUT,'(A)') ' Model grid will be used to develop ',&
         &' volume and area relationship. '
      END IF
   !
   !2------READ NLAKES, ILKCB.
   !
   !dep  Revised input statement to read THETA,NSSITR,SSCNCR for
   !dep  transient simulations when THETA is negative.
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(2I10)')NLAKES,ILKCB
         IF (ITRSS.LE.0) THEN
            READ(IN,'(F10.2,I10,F10.2)') THETA,NSSITR,SSCNCR
            IF (THETA.LT.0.0) BACKSPACE IN
         ELSE
            READ(IN,'(F10.2)') THETA
            IF (THETA.LT.0.0) BACKSPACE IN
         END IF
      ELSE
         READ(IN,*) NLAKES,ILKCB
         IF (ITRSS.LE.0) THEN
            READ(IN,*) THETA,NSSITR,SSCNCR
            IF(THETA.LT.0.0) BACKSPACE IN
         ELSE
            READ(IN,*) THETA
            IF(THETA.LT.0.0) BACKSPACE IN
         END IF
      END IF

   !dep    Set default values for number of iterations and closure criteria
   !dep     for transient simulations when using original version of
   !dep     LAKE Package.
      IF(THETA.GE.0.0.AND.NSSITR.EQ.0) THEN
         NSSITR=100
         SSCNCR=1.0E-05
      ELSE IF(THETA.LT.0.0)THEN
         THETA=ABS(THETA)
         IF(IFREFM.EQ.0) THEN
   !dep fixed format can't read in exponent notation
   !rsr, old data sets may not have SURFDEPTH, may need to trap this for some compilers
            READ (IN, '(A)') CARD
            NUMCHAR = LEN(TRIM(CARD))
            IF ( NUMCHAR>30 ) THEN
               READ(CARD,'(F10.2,I10,2F10.5)') DUM,NSSITR,SSCNCR,&
               &SURFDEPTH
            ELSE
               READ(CARD,'(F10.2,I10,F10.5)') DUM,NSSITR,SSCNCR
            ENDIF
         ELSE
            READ(IN,*,IOSTAT=IOS) DUM,NSSITR,SSCNCR,SURFDEPTH
            IF ( IOS.NE.0 ) SURFDEPTH = 0.0
         END IF
      END IF
   !dep   Add check to reset THETA when > 1 or < 0.5.
      IF(THETA.GT.1.0) THEN
         THETA = 1.0
      ELSE IF(THETA.LT.0.5)THEN
         THETA = 0.0
      END IF
   END IF
   !
   !
   !  SET NLAKES ARRAY VARIABLE TO NLAKES IF NLAKES GREATER THAN 0.
   IF (NLAKES.GT.0) NLAKESAR = NLAKES
   ALLOCATE (VOL(NLAKESAR), STGOLD(NLAKESAR), STGNEW(NLAKESAR))
   ALLOCATE (RUNF(NLAKESAR), RUNOFF(NLAKESAR))
   ALLOCATE(STGOLD2(NLAKESAR))
   ALLOCATE (VOLOLDD(NLAKESAR))
   !     ALLOCATE (VOLOLDD(NLAKESAR), VOLOLD(NLAKES), VOLINIT(NLAKES))
   ALLOCATE (STGITER(NLAKESAR))
   ALLOCATE (LAKSEEP(NCOL,NROW),DEADPOOLVOL(NLAKESAR),&
   &RELEASABLE_STOR(NLAKESAR), MXLKVOLF(NLAKESAR))
   STGNEW = 0.0D0
   STGOLD = 0.0D0
   STGOLD2 = 0.0D0
   STGITER = 0.0D0
   VOLOLDD = 0.0D0
   LAKSEEP = 0.0
   RUNF = 0.0D0
   RUNOFF = 0.0D0
   DEADPOOLVOL = 0.0
   MXLKVOLF = 0.0
   RELEASABLE_STOR = 0.0
   !dep initialized VOLOLD and VOLINIT  6/4/2009 (VOLOLD is single precision)
   !     VOLOLD = 0.0
   !     VOLINIT = 0.0
   VOL = 0.0
   CALL SGWF2LAK7PSV1(IGRID)
   IF (IN.LT.1) RETURN
   !
   ! Lakes are active
   ALLOCATE (STAGES(NLAKESAR), CLAKE(NLAKESAR,NSOL))
   STAGES = 0.0
   CLAKE = 0.0
   ! Budget variables for GSFLOW
   ALLOCATE (TOTGWIN_LAK,TOTGWOT_LAK,TOTDELSTOR_LAK,TOTSTOR_LAK)
   ALLOCATE (TOTEVAP_LAK,TOTPPT_LAK,TOTRUNF_LAK,TOTWTHDRW_LAK)
   ALLOCATE (TOTSURFIN_LAK,TOTSURFOT_LAK)
   TOTGWIN_LAK = 0.0
   TOTGWOT_LAK = 0.0
   TOTDELSTOR_LAK = 0.0
   TOTSTOR_LAK = 0.0
   TOTEVAP_LAK = 0.0
   TOTPPT_LAK = 0.0
   TOTRUNF_LAK = 0.0
   TOTWTHDRW_LAK = 0.0
   TOTSURFIN_LAK = 0.0
   TOTSURFOT_LAK = 0.0
   !
   !  VALUE OF MXLKND (NUMBER OF LAKE-AQUIFER INTERFACES) IS AN ESTIMATE.
   !    TO SAVE MEMORY, REDUCE ITS SIZE IF APPROPRIATE.
   !    IF MXLKND TOO SMALL, ERROR MESSAGE WILL BE PRINTED.
   MXLKND=NCOL*NROW*NLAY/2
   IF (NLAKES.LT.1) THEN
      WRITE(IOUT,2)
      IN=0
      NLAKES = 0
   ELSE
      WRITE(IOUT,5) MXLKND,NLAKES
      IF (ILKCB.GT.0) WRITE(IOUT,7) ILKCB
      IF (ILKCB.LE.0) WRITE(IOUT,9)
   !dep   Write THETA, NSSITR, SSCNCR
      IF (ITRSS.GT.0) THEN
         WRITE(IOUT,22) THETA
         WRITE(IOUT,10) NSSITR, SSCNCR
      ELSE
         WRITE(IOUT,11) THETA, NSSITR, SSCNCR
      END IF
   !dep   Changed default values for NSSITR and SSCNCR and revised
   !dep     print statements using format statement 10.
   !dep      IF(ITRSS.LE.0.AND.NSSITR.EQ.0) NSSITR = 50
   !dep      IF(ITRSS.LE.0.AND.SSCNCR.EQ.0.0) SSCNCR = 0.01
   !dep      IF(ITRSS.EQ.0) WRITE(IOUT,23) NSSITR, SSCNCR
   !dep      IF(ITRSS.LT.0) WRITE(IOUT,24) NSSITR, SSCNCR
   !-lfk  1     FORMAT(/1X,'LAK7 -- LAKE PACKAGE, VERSION 7, 2/06/2012',
1     FORMAT(/1X,'LAK7 -- LAKE PACKAGE, VERSION 7, 7/01/2022',&
      &' INPUT READ FROM UNIT',I3)
2     FORMAT(1X,' NUMBER OF LAKES=0, ',&
      &' SO LAKE PACKAGE IS BEING TURNED OFF')
5     FORMAT(1X,'SPACE ALLOCATION FOR',I7,' GRID CELL FACES ADJACENT TO&
      &LAKES'/1X,'MAXIMUM NUMBER OF LAKES IS',I3, ' FOR THIS SIMULATION')
7     FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE RECORDED ON UNIT',I5)
9     FORMAT(1X,'CELL-BY-CELL SEEPAGES WILL NOT BE PRINTED OR SAVED')
   !dep added format statement when starting with transient simulation
10    FORMAT(//1X,'LAKE PACKAGE HAS BEEN MODIFIED TO ITERATIVELY ',&
      &'SOLVE FOR LAKE STAGE DURING TRANSIENT STRESS PERIODS:',/1X,&
      &'MAXIMUM NUMBER OF ITERATIONS (NSSITR) = ',I5,/1X,&
      &'CLOSURE CRITERIA FOR LAKE STAGE (SSCNCR) = ',1PE13.6,/1X,& !gsf
      &'DEFAULT VALUES FOR TRANSIENT ONLY SIMULATIONS ARE: ',&
      &'NSSITR = 100 AND SSCNCR = 0.0001',/1X,'VALUES OTHER THAN ',&
      &'DEFAULT CAN BE READ BY SPECIFYING A THETA LESS THAN ZERO ',&
      &'THEN ADDING NSSITR AND SSCNCR PER ORIGINAL INSTRUCTIONS.',/1X,&
      &'NEGATIVE THETA MUST BE LESS THAN ZERO BUT NOT MORE THAN ',&
      &'ONE. THETA IS CONVERTED TO A POSITIVE VALUE.',/1X,&
      &'MINIMUM AND MAXIMUM LAKE STAGES FOR TRANSIENT ',&
      &'SIMULATIONS ARE SET TO BOTTOM AND TOP ELEVATIONS USED TO ',&
      &'COMPUTE LAKE VOLUME, RESPECTIVELY.',//)
   !dep added format statement for steady state only simulations.
11    FORMAT(//1X,'NEWTON ITERATION METHOD FOR COMPUTING LAKE STAGE ',&
      &'DURING STEADY-STATE STRESS PERIODS HAS BEEN MODIFIED:',/1X,&
      &'SPECIFIED THETA OF ',F6.3,' WILL BE AUTOMATICALLY CHANGED TO ',&
      &'1.0 FOR ALL STEADY STATE STRESS PERIODS.',/1X,&
      &'MAXIMUM NUMBER OF STEADY-STATE ITERATIONS (NSSITR) = ',I5,/1X,&
      &'CLOSURE CRITERIA FOR STEADY-STATE LAKE STAGE (SSCNCR) = ',&
      &1PE13.6,//) !gsf
   !dep revised print statement to note that time weighting of theta can
   !dep  vary only between 0.5 and 1 for transient simulations
   !dep   22 FORMAT(/1X,'THETA = ',F10.2,'  METHOD FOR UPDATING LAKE STAGES IN
   !dep     1ITERATIONS OF THE SOLUTION FOR AQUIFER HEADS.'/20X,'0.0 IS EXPLICI
   !dep     2T, 0.5 IS CENTERED, AND 1.0 IS FULLY IMPLICIT.')
22    FORMAT(/1X,'THETA = ',F6.3,/1X,'THETA IS THE TIME WEIGHTING ',&
      &'FACTOR FOR COMPUTING LAKE STAGE DURING TRANSIENT MODFLOW ',&
      &'TIME STEPS AND ITS DEFINITION HAS BEEN MODIFIED.',/1X,'A THETA ',&
      &'OF LESS THEN 0.5 IS AUTOMATICALLY SET TO 0 AND LAKE STAGE IS ',&
      &'EQUAL TO THE STAGE AT THE END OF THE PREVIOUS TIME STEP. ',/1X,&
      &'TRANSIENT SIMULATIONS OF LAKE STAGE WITH THE CURRENT TIME STEP ',&
      &'REQUIRES A THETA BETWEEN 0.5 AND 1.0. ',/1X,'VALUES GREATER ',&
      &'THAN 1.0 ARE AUTOMATICALLY RESET TO  1.0 AND VALUES LESS ',&
      &'THAN 0.5 ARE RESET TO 0.0.',/1X,'A THETA OF 0.5 REPRESENTS THE ',&
      &'AVERAGE LAKE STAGE DURING A TIME STEP.',/1X,'A THETA OF 1.0 ',&
      &'REPRESENTS THE LAKE STAGE AT THE END OF THE TIME STEP.',//)
   !dep   23 FORMAT(/1X,'STEADY-STATE SOLUTION FOR LAKES.'
   !dep     2/1X,'MAXIMUM NUMBER OF ITERATIONS = ',I4,3X,
   !dep     1'CONVERGENCE CRITERION = ',1PE9.2)
   !dep   24 FORMAT(/1X,'COMBINED STEADY-STATE/TRANSIENT SOLUTION FOR LAKES.'
   !dep     2/1X,'MAXIMUM NUMBER OF ITERATIONS = ',I4,3X,
   !dep     1'CONVERGENCE CRITERION = ',1PE9.2)

      ALLOCATE (ILAKE(5,MXLKND), BEDLAK(MXLKND), CNDFCT(MXLKND))
      ALLOCATE (PRCPLK(NLAKES), EVAPLK(NLAKES), WTHDRW(NLAKES))
      ALLOCATE (RNF(NLAKES), CRNF(NLAKES,NSOL), CUMRNF(NLAKES))
      ALLOCATE (CUMUZF(NLAKES))
      ALLOCATE (ISUB(NLAKES,NLAKES), SILLVT(NLAKES,NLAKES))
      ALLOCATE (IRK(2,NLAKES))
      ALLOCATE (CUMPPT(NLAKES), CUMEVP(NLAKES), CUMGWI(NLAKES))
      ALLOCATE (CUMGWO(NLAKES), CUMSWI(NLAKES), CUMSWO(NLAKES))
      ALLOCATE (CUMWDR(NLAKES), CUMFLX(NLAKES))
      ALLOCATE (CAUG(NLAKES,NSOL), CPPT(NLAKES,NSOL))
      ALLOCATE (CLAKINIT(NLAKESAR,NSOL))
      ALLOCATE (ICS(NLAKES),BOTTMS(NLAKES), BGAREA(NLAKES))
      ALLOCATE (SSMN(NLAKES), SSMX(NLAKES))
      ALLOCATE (LKARR1(NCOL,NROW,NLAY), BDLKN1(NCOL,NROW,NLAY))
      ALLOCATE (EVAP(NLAKES), PRECIP(NLAKES), SEEP(NLAKES),&
      &SEEP3(NLAKES),EVAP3(NLAKES), PRECIP3(NLAKES))
      ALLOCATE (SEEPUZ(NLAKES))
      ALLOCATE (FLWITER(NLAKES),FLWITER3(NLAKES))
      ALLOCATE (SURFA(NLAKES), SURFIN(NLAKES), SURFOT(NLAKES))
      ALLOCATE (SUMCNN(NLAKES), SUMCHN(NLAKES))
      ALLOCATE (NCNCVR(NLAKES), LIMERR(NLAKES), DSRFOT(NLAKES))
   !dep  Allocate arrays that track lake budgets for dry lakes
      ALLOCATE (WITHDRW(NLAKES),FLWIN(NLAKES))
      ALLOCATE (GWRATELIM(NLAKES))
      WITHDRW = 0.0D0
      FLWIN = 0.0
      FLWITER = 0.0D0
      FLWITER3 = 0.0D0
      EVAP = 0.0D0
      PRECIP = 0.0D0
      EVAP3 = 0.0D0
      PRECIP3 = 0.0D0
      IF ( IRDTAB.GT.0 ) THEN
         ALLOCATE(LAKTAB(NLAKES))
      ELSE
         ALLOCATE(LAKTAB(1))
      END IF
      LAKTAB = 0
   !rsr    GWRATLIM= 0.0
   !dep  Allocate space for three arrays used in GAGE Package
   !       when Solute Transport is active
      ALLOCATE (XLAKES(NLAKES,1), XLAKINIT(NLAKES,1))
      ALLOCATE (XLKOLD(NLAKES,1))
   !rsr  Allocate arrays for BD subroutine
      ALLOCATE (LDRY(NODES), FLXINL(NLAKES))
      ALLOCATE (NCNT(NLAKES), NCNST(NLAKES))
      ALLOCATE (SVT(NLAKES), KSUB(NLAKES), STGADJ(NLAKES))
      ALLOCATE (MSUB(NLAKES,NLAKES), MSUB1(NLAKES))
      ALLOCATE (GWIN(NLAKES), GWOUT(NLAKES))
      ALLOCATE (DELH(NLAKES), TDELH(NLAKES))
   !dep   Allocate lake budget error arrays for BD subroutine 6/9/2009
      ALLOCATE (CUMVOL(NLAKES), CMLAKERR(NLAKES))
      ALLOCATE (CUMLKIN(NLAKES), CUMLKOUT(NLAKES))
      ALLOCATE (DELVOL(NLAKES), TSLAKERR(NLAKES))
   !dep initialized VOLOLD and VOLINIT  6/4/2009 (VOLOLD is single precision)
      ALLOCATE (VOLOLD(NLAKES), VOLINIT(NLAKES))
      VOLOLD = 0.0
      VOLINIT = 0.0
   END IF
   !dep   ALLOCATE SPACE FOR CONNECTION WITH STREAMS
   IF (IUNITSFR.LE.0) THEN
      NSSAR = 1
   ELSE
      NSSAR = NSS
   END IF
   !dep   ALLOCATE SPACE FOR FLOB ARRAY WHEN TRANSPORT ACTIVE.
   IF (IUNITGWT.LE.0.AND.IUNIT(49).LE.0) THEN
      MXLKAR = 1
   ELSE
      MXLKAR = MXLKND
   END IF
   !dep    ALLOCATE SPACE FOR OVERLAND FLOW WHEN UNSATURATED FLOW ACTIVE.
   ! RGN Allocate NUZFAR to nlakes for all cases because of the GAG package 5/28/09
   !      IF (IUNITUZF.LE.0) THEN
   !       NUZFAR = 1
   !      ELSE
   NUZFAR = NLAKESAR
   !      END IF

   !rsr, what if NLAKES < 1, sanity check
   IF (NLAKES<1 ) THEN
      print *, 'nlakes dimension problem in lak7', nlakes
      stop
   ENDIF

   ALLOCATE (ITRB(NLAKES,NSSAR), IDIV(NLAKES,NSSAR))
   ALLOCATE (FLOB(MXLKAR))
   ALLOCATE (OVRLNDRNF(NUZFAR), CUMLNDRNF(NUZFAR))
   !dep    ALLOCATE SPACE FOR DEPTHTABLE, AREATABLE, AND VOLUMETABLE
   ALLOCATE (DEPTHTABLE(151,NLAKES), AREATABLE(151,NLAKES))
   ALLOCATE (VOLUMETABLE(151,NLAKES))
   ! Tributary inflow to lakes for LMT
   ALLOCATE (LAKSFR(NSSAR),ILKSEG(NSSAR),ILKRCH(NSSAR),SWLAK(NSSAR),&
   &DELVOLLAK(NLAKES))
   ITRB = 0
   IDIV = 0
   FLOB = 0.0
   OVRLNDRNF = 0.0
   CUMLNDRNF = 0.0
   CUMUZF = 0.0
   DEPTHTABLE = 0.0D0
   AREATABLE = 0.0D0
   VOLUMETABLE = 0.0D0
   !dep initialized lake budget error arrays  6/9/2009
   CUMVOL = 0.0
   DELVOL = 0.0
   CMLAKERR = 0.0
   TSLAKERR = 0.0
   CUMLKOUT = 0.0
   CUMLKIN = 0.0
   !-----SAVE POINTERS FOR GRID AND RETURN
   CALL SGWF2LAK7PSV(IGRID)
   !
   !11-----RETURN.
   RETURN
END
   !
SUBROUTINE GWF2LAK7RP(IN,IUNITBCF,IUNITGWT,IUNITLPF,IUNITHUF,&
&IUNITSFR,IUNITUZF,IUNITUPW,KKPER,NSOL,&
&IOUTS,IGRID)
   !
   !------OLD USGS VERSION 7.1;  JUNE 2006 GWF2LAK7RP
   !        REVISED FEBRUARY 6, 2012
   !------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
   !------NEW VERSION NUMBER 1.3.0, 7/01/2022
   !     ******************************************************************
   !       READ INPUT DATA FOR THE LAKE PACKAGE.
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFLAKMODULE
   USE GLOBAL,       ONLY: IOUT, NCOL, NROW, NLAY, IFREFM, IBOUND,&
   &LBOTM, BOTM, DELR, DELC, ISSFLG,IUNIT
   !
   IMPLICIT NONE
   !
   !     USE GWFSFRMODULE, ONLY: NSS
   !     ------------------------------------------------------------------
   !     FUNCTIONS
   !     ------------------------------------------------------------------
   DOUBLE PRECISION VOLTERP
   EXTERNAL VOLTERP
   !     ------------------------------------------------------------------
   INTEGER IGRID,ISS,LM,IUNITGWT,IN,ISOL,L1,I,J,K,KK,LK,ITMP,ITMP1,&
   &INC,NSOL,N,I2,K1,K2,K3,K4,I1,IC,IS,JK,NSLMS,&
   &IUNITSFR,LAKEFLG,LAKE,NTYP,J2,KKPER,IOUTS,IUNITNUM,L,&
   &IL,IR,ITYPE,IUNITBCF,IUNITLPF,IUNITHUF,IUNITUPW,&
   &IUNITUZF,IC1,II,IFACE,LID,M
   REAL BOTIJ,TBNC,TBELV,TOPMST,GTSDPH,EVOL,BOTLK
   CHARACTER*24 ANAME(2)
   !     CHARACTER*30 LFRMAT
   !dep  added STGINIT as double precision
   DOUBLE PRECISION STGINIT
   DATA ANAME(1)/'           LAKE ID ARRAY'/
   DATA ANAME(2)/'  LAKEBED LEAKANCE ARRAY'/
   !
   !     ------------------------------------------------------------------
   !------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2LAK7PNT(IGRID)
   !
   !1A-----IF MXLKND IS LESS THAN 1, THEN LAKE IS INACTIVE. RETURN.
   IF(MXLKND.LT.1) RETURN
   !
   !1A1----READ INITIAL CONDITIONS FOR ALL LAKES (ONLY READ ONCE)
   ISS = ISSFLG(KKPER)
   IF (KKPER.EQ.1) THEN
      WRITE (IOUT,19)
      IF(ISS.NE.0) WRITE (IOUT,20)
      IF(ISS.EQ.0) WRITE (IOUT,820)
      IF (IUNITGWT.EQ.0) THEN
         DO 30 LM=1,NLAKES
            IF (IFREFM.EQ.0) THEN
               IF ( IRDTAB.GT.0 ) THEN
                  IF(ISS.NE.0) READ (IN,'(3F10.4,I5)') STAGES(LM),&
                  &SSMN(LM),SSMX(LM),LAKTAB(LM)
                  IF(ISS.EQ.0) READ (IN,'(F10.4,I5)') STAGES(LM),&
                  &LAKTAB(LM)
               ELSE
                  IF(ISS.NE.0) READ (IN,'(3F10.4)') STAGES(LM),&
                  &SSMN(LM),SSMX(LM)
                  IF(ISS.EQ.0) READ (IN,'(F10.4)') STAGES(LM)
               END IF
            ELSE
               IF ( IRDTAB.GT.0 ) THEN
                  IF(ISS.NE.0) READ (IN,*)STAGES(LM),SSMN(LM),SSMX(LM),&
                  &LAKTAB(LM)
                  IF(ISS.EQ.0) READ (IN,*) STAGES(LM),LAKTAB(LM)
               ELSE
                  IF(ISS.NE.0) READ (IN,*) STAGES(LM),SSMN(LM),SSMX(LM)
                  IF(ISS.EQ.0) READ (IN,*) STAGES(LM)
               END IF
            END IF
            IF(ISS.NE.0) WRITE (IOUT,22) LM,STAGES(LM),SSMN(LM),SSMX(LM)
            IF(ISS.EQ.0) WRITE (IOUT,22) LM,STAGES(LM)
30       CONTINUE
      ELSE
         WRITE (IOUTS,21) NSOL
   !            WRITE (LFRMAT,23) NSOL  !LFRMAT is not set
         DO 35 LM=1,NLAKES
            IF (IFREFM.EQ.0) THEN
               IF ( IRDTAB.GT.0 ) THEN
                  IF(ISS.NE.0) READ(IN,'(100F10.4)') STAGES(LM),&
                  &SSMN(LM),SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL),&
                  &LAKTAB(LM)
                  IF(ISS.EQ.0) READ (IN,'(100F10.4)') STAGES(LM),&
                  &(CLAKE(LM,ISOL),ISOL=1,NSOL),LAKTAB(LM)
               ELSE
                  IF(ISS.NE.0) READ(IN,'(100F10.4)') STAGES(LM),&
                  &SSMN(LM),SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL)
                  IF(ISS.EQ.0) READ (IN,'(100F10.4)') STAGES(LM),&
                  &(CLAKE(LM,ISOL),ISOL=1,NSOL)
               END IF
            ELSE
               IF ( IRDTAB.GT.0 ) THEN
                  IF(ISS.NE.0) READ (IN,*) STAGES(LM),SSMN(LM),&
                  &SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL),&
                  &LAKTAB(LM)
                  IF(ISS.EQ.0) READ (IN,*) STAGES(LM),&
                  &(CLAKE(LM,ISOL),ISOL=1,NSOL),LAKTAB(LM)
               ELSE
                  IF(ISS.NE.0) READ (IN,*) STAGES(LM),SSMN(LM),&
                  &SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL)
                  IF(ISS.EQ.0) READ (IN,*) STAGES(LM),&
                  &(CLAKE(LM,ISOL),ISOL=1,NSOL)
               END IF
            END IF
            IF(ISS.NE.0) WRITE (IOUT,22) LM,STAGES(LM),SSMN(LM),SSMX(LM)
            IF(ISS.EQ.0) WRITE (IOUT,22) LM,STAGES(LM)
35       WRITE (IOUTS,*) LM,(CLAKE(LM,ISOL),ISOL=1,NSOL)
   !gage
   !            CLAKINIT=CLAKE
      END IF
   END IF
   !
   WRITE (IOUT,'(/)')
   WRITE(IOUT,822)
19 FORMAT(//1X,'LAKE PACKAGE ACTIVE:  CALCULATED LAKE STAGE FOR EACH&
   &TIME STEP WILL BE STORED IN HNEW ARRAY.')
20 FORMAT(///1X,'INITIAL LAKE STAGE:  LAKE    STAGE    SS MIN    SS M&
   &AX'/)
21 FORMAT (//1X,'INITIAL LAKE CONCENTRATIONS:  LAKE   CONCENTRATION (&
   &NSOL =',I3,')'/)
22 FORMAT (22X,I3,3F10.3)
23 FORMAT ('(31X,I3,3X,1P',I3,'(E12.3))')
820 FORMAT (/1X,'INITIAL LAKE STAGE:  LAKE    STAGE'/)
822 FORMAT(//1X,'If any subsequent steady-state stress periods, min. a&
   &nd max. stages for each lake will be read in Record 9a.'//)
   !
   ! RGN 9/25/12 moved this to read lake bathymetry before stress period information.
   IF ( KKPER==1 .AND. IRDTAB.GT.0 ) THEN
      DO L1=1,NLAKES
         WRITE(IOUT,1399) L1
         iunitnum = LAKTAB(L1)
1399     FORMAT(//1X,'STAGE/VOLUME RELATION FOR LAKE',I3//6X,'STAGE',&
         &8X,'VOLUME',8X,'AREA'/)
         DO  INC=1,151
            READ(iunitnum,*) DEPTHTABLE(INC,L1), VOLUMETABLE(INC,L1),&
            &AREATABLE(INC,L1)
            WRITE(IOUT,1315) DEPTHTABLE(INC,L1), VOLUMETABLE(INC,L1),&
            &AREATABLE(INC,L1)
         END DO
      END DO
   END IF
   !1B-----READ ITMP (FLAG TO REUSE LAKE-GEOMETRY DATA).
   IF(IFREFM.EQ.0) THEN
      READ(IN,'(3I10)') ITMP, ITMP1, LWRT
   ELSE
      READ(IN,*) ITMP, ITMP1, LWRT
   END IF
   !
   !2A-----IF ITMP < 0 THEN REUSE LAKE CONFIGURATION DATA FROM LAST STRESS
   !       PERIOD.
   IF(ITMP.GE.0) GO TO 50
   WRITE (IOUT,'(/)')
   WRITE(IOUT,2)
2  FORMAT(1H ,'REUSING LAKE CONFIGURATION DATA FROM LAST STRESS PERIO&
   &D'/)
   GO TO 800
   !
   !4------IF THERE ARE NO LAKE NODES THEN RETURN.
50 LKNODE = 0
   IF(ITMP.EQ.0) GOTO 900
   !
   !   INITIALIZE BGAREA
   DO 60 LK=1,NLAKES
      BGAREA(LK)=0.0
60 CONTINUE
   !
   !5------READ INTEGER ARRAYS THAT DEFINE THE POSITIONS OF ALL LAKES IN
   !5A     EACH MODEL GRID LAYER.  THEN READ ARRAYS OF LAKEBED CONDUCTANCES
   !5B     IN EACH LAYER.
   !
   !   READ ARRAY OF LAKE ID'S, LAYER BY LAYER
   !   REVISED 11/30/2005 DEP
   DO 125 K=1,NLAY
      KK = K
      CALL U2DINT(LKARR1(:,:,KK),ANAME(1),NROW,NCOL,KK,IN,IOUT)
125 CONTINUE
   !
   !   CHECK THAT ALL ENTRIES ARE VALID LAKE ID NUMBERS OR ZERO
   !
   DO 130 K=1,NLAY
      DO 130 I=1,NCOL
         DO 130 J=1,NROW
            IF(LKARR1(I,J,K).GT.0.AND.LKARR1(I,J,K).LE.NLAKES) GO TO 130
            LKARR1(I,J,K)=0
130 CONTINUE
   !
   !   CHECK IF LAKE CELLS HAVE VALUES OF IBOUND=0; WARN IF INCONSISTENT
   !
   WRITE (IOUT,'(/)')
   DO 132 K=1,NLAY
      DO 132 I=1,NCOL
         DO 132 J=1,NROW
            IF(LKARR1(I,J,K).GT.0.AND.IBOUND(I,J,K).NE.0) THEN
               WRITE (IOUT,232) IBOUND(I,J,K),LKARR1(I,J,K),I,J,K
232            FORMAT (7X,'*** WARNING: IBOUND = ',I2,&
               &' & LKARR = ',I2,' at CELL I=',I3,&
               &', J=',I3,', K=',I3,' ***')
            END IF
132 CONTINUE
   !
   !   READ ARRAY OF BED LEAKANCES, LAYER BY LAYER
   !dep    REVISED 11/30/2005
   WRITE (IOUT,'(/)')
   DO 135 K=1,NLAY
      KK = K
      CALL U2DREL(BDLKN1(:,:,KK),ANAME(2),NROW,NCOL,KK,IN,IOUT)
135 CONTINUE
   !
   WRITE(IOUT,36)
   WRITE(IOUT,4)
36 FORMAT(/7X,'LOCATIONS, LAKE #, INTERFACE TYPE FOR GRID CELLS',&
   &' ADJACENT TO LAKES:',5X,/&
   &5X,71('-'))
4  FORMAT(5X,'LAYER #',4X,'ROW #',4X,'COLUMN #',3X,'LAKE #',&
   &2X,'INTERFACE TYPE',2X,'LAKEBED LEAKANCE')
   !
   !   IDENTIFY LAKE BORDER CELLS, ASSIGN CELL TYPE ID'S, COMPUTE AND
   !     ASSIGN LAKE-AQUIFER INTERFACE CONDUCTANCES.
   !
   M = 0
   DO 180 I=1,NCOL
      DO 180 J=1,NROW
         K = 1
         IF(LKARR1(I,J,K).EQ.0) GO TO 150
         IF(NLAY.EQ.1) GO TO 145
   !   Keep searching in vertical direction until non-lake cell is found,
   !     and define interface there ("K" for interface is layer below
   !     bottom of lake)
         DO 140 K=2,NLAY
            IF(LKARR1(I,J,K).EQ.0) GO TO 145
140      CONTINUE
   !   Make sure that K=NLAY if lake extends to bottom cell of grid:
         K=NLAY
   !      GO TO 145
   !
   !   VERTICAL LAKEBED INTERFACE (TYPE 0) DETECTED
   !
145      M = M + 1
         IF(M.LE.MXLKND) GO TO 147
         WRITE(IOUT,149) I,J,K
149      FORMAT(/1X,'MAXIMUM NUMBER OF GRID CELLS ADJACENT TO LAKES HAS BEE&
         &N EXCEEDED WITH CELL ',3I5,'  REDEFINE VARIABLE MXLKND TO A LARGER&
         & VALUE IN MODULE GWF2LAK7AR')
         CALL USTOP(' ')
147      ILAKE(1,M) = K
         ILAKE(2,M) = J
         ILAKE(3,M) = I
   !dep  changed if statement August 24, 2009
   !dep      IF(K.GT.1.AND.LKARR1(I,J,K).EQ.0) LID = LKARR1(I,J,K-1)
   !dep      IF(LKARR1(I,J,K).NE.0) LID = LKARR1(I,J,K)
         IF(K.GT.1) THEN
            IF(LKARR1(I,J,K).EQ.0) THEN
               LID = LKARR1(I,J,K-1)
            ELSE
               LID = LKARR1(I,J,K)
            END IF
         ELSE IF (K.EQ.1) THEN
            IF(LKARR1(I,J,K).EQ.0) THEN
               LID = 0
            ELSE
               LID = LKARR1(I,J,K)
            END IF
         END IF
         ILAKE(4,M) = LID
         ILAKE(5,M) = 6
         IF ( K.GT.1 ) THEN             !RGN 5/21/12 added IF test
            BEDLAK(M) = BDLKN1(I,J,K-1)
         ELSE                           !RGN
            BEDLAK(M) = BDLKN1(I,J,K)    !RGN
         END IF                         !RGN
         IF(K.EQ.NLAY.AND.LKARR1(I,J,K).NE.0) BEDLAK(M) = 0.0
         BGAREA(LID) = BGAREA(LID) + DELC(J)*DELR(I)
   !-LFK-JAN. 2013
   !-lfk        WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
         WRITE(IOUT,6) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
5        FORMAT(5I10,10X,F10.5)
6        FORMAT(5I10,12X,1PE10.3)
   !-LFK
         IF(LKARR1(I,J,K).NE.0) GO TO 180
   !
   !   SEARCH FOR CELL(S) ADJACENT TO LAKE
   !
150      K2 = K
         DO 175 K1=K2,NLAY
   !gzh fix for 2D-problems
            IF(NCOL.EQ.1) GO TO 165
            IF(I.NE.1) GO TO 1151
            IF(LKARR1(I+1,J,K1).EQ.0) GO TO 165
            GO TO 1153
1151        IF(I.NE.NCOL) GO TO 1152
            IF(LKARR1(I-1,J,K1).EQ.0) GO TO 165
            GO TO 1153
1152        IF(LKARR1(I+1,J,K1).EQ.0.AND.LKARR1(I-1,J,K1).EQ.0) GO TO 165
   !
   !   CELL(S) LATERALLY ADJACENT TO LAKE IN X-DIRECTION (TYPE 1) DETECTED
   !
1153        DO 160 N=1,2
               IF(N.EQ.2) GO TO 155
               IF(I.EQ.1) GO TO 160
               IF(LKARR1(I-1,J,K1).EQ.0) GO TO 160
               I2 = I-1
               IFACE=1
               GO TO 157
155            IF(I.EQ.NCOL) GO TO 160
               IF(LKARR1(I+1,J,K1).EQ.0) GO TO 160
               I2 = I + 1
               IFACE=2
157            M = M + 1
               IF(M.LE.MXLKND) GO TO 158
               WRITE(IOUT,149) I,J,K1
               CALL USTOP(' ')
158            ILAKE(1,M) = K1
               ILAKE(2,M) = J
               ILAKE(3,M) = I
               ILAKE(4,M) = LKARR1(I2,J,K1)
               ILAKE(5,M) = IFACE
               BEDLAK(M) = BDLKN1(I,J,K1)
               K4 = K1 - 1
               DO 3158 K3=1,K4
                  IF(LKARR1(I,J,K3).EQ.0) GO TO 3158
                  GO TO 3162
3158           CONTINUE
               BEDLAK(M) = BDLKN1(I,J,1)
3162           CONTINUE
   !-lfk      WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
   !-LFK-JAN. 2013
   !-lfk        WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
               WRITE(IOUT,6) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
   !-LFK
160         CONTINUE
   !gzh fix for 2D-problems
165         IF(NROW.EQ.1) GO TO 175
            IF(J.NE.1) GO TO 1161
            IF(LKARR1(I,J+1,K1).EQ.0) GO TO 175
            GO TO 1163
1161        IF(J.NE.NROW) GO TO 1162
            IF(LKARR1(I,J-1,K1).EQ.0) GO TO 175
            GO TO 1163
1162        IF(LKARR1(I,J+1,K1).EQ.0.AND.LKARR1(I,J-1,K1).EQ.0) GO TO 175
   !
   !   CELL(S) LATERALLY ADJACENT TO LAKE IN Y-DIRECTION (TYPE 2) DETECTED
   !
1163        DO 170 N=1,2
               IF(N.EQ.2) GO TO 172
               IF(J.EQ.1) GO TO 170
               IF(LKARR1(I,J-1,K1).EQ.0) GO TO 170
               J2 = J - 1
               IFACE=4
               GO TO 174
172            IF(J.EQ.NROW) GO TO 170
               IF(LKARR1(I,J+1,K1).EQ.0) GO TO 170
               J2 = J + 1
               IFACE=3
174            M = M + 1
               IF(M.LE.MXLKND) GO TO 176
               WRITE(IOUT,149) I,J,K1
               CALL USTOP(' ')
176            ILAKE(1,M) = K1
               ILAKE(2,M) = J
               ILAKE(3,M) = I
               ILAKE(4,M) = LKARR1(I,J2,K1)
               ILAKE(5,M) = IFACE
               BEDLAK(M) = BDLKN1(I,J,K1)
               K4 = K1 - 1
               DO 4158 K3=1,K4
                  IF(LKARR1(I,J,K3).EQ.0) GO TO 4158
                  GO TO 4162
4158           CONTINUE
               BEDLAK(M) = BDLKN1(I,J,1)
4162           CONTINUE
   !-lfk      WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
   !-LFK-JAN. 2013
   !-lfk        WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
               WRITE(IOUT,6) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
   !-LFK
170         CONTINUE
175      CONTINUE
180 CONTINUE
   WRITE(IOUT,195) M
195 FORMAT(/5X,'NUMBER OF LAKE-AQUIFER CELL INTERFACES = ',I5)
   LKNODE = M
   !
   !   SET LAKE BOTTOM ELEVATIONS
   DO 295 LK=1,NLAKES
295 BOTTMS(LK) = 999999
   !
   DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
   !  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2,
   !    6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      IF(NTYP.EQ.0) THEN
         LAKE = ILAKE(4,II)
   !dep  changed if statement August 24, 2009
   !dep        IF(K.GT.1) BOTLK = BOTM(I,J,LBOTM(K-1))
   !dep        IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) BOTLK = BOTM(I,J,LBOTM(K))
         IF(K.EQ.1.OR.K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) THEN
            BOTLK = BOTM(I,J,LBOTM(K))
         ELSE IF (K.EQ.0) THEN
            BOTLK = BOTM(I,J,LBOTM(1))
         ELSE
            BOTLK = BOTM(I,J,LBOTM(K-1))
         END IF
         IF(BOTLK.LT.BOTTMS(LAKE)) BOTTMS(LAKE) = BOTLK
      END IF
350 CONTINUE
   !
   !-- COMPUTE AND PRINT STAGE/VOLUME TABLES WHEN MORE THAN ONE LAYER
   !dep  revised print statement to include stage/area tables
   !
   IF ( IRDTAB.EQ.0 ) THEN
   !      IF(NLAY.EQ.1) GO TO 1331       !RGN 5/21/12
      DO 1330 L1=1,NLAKES
         WRITE(IOUT,1306) L1
   !dep  revised print statement to include area
1306     FORMAT(//1X,'STAGE/VOLUME RELATION FOR LAKE',I3//6X,'STAGE',&
         &8X,'VOLUME',8X,'AREA'/)
         DO  INC=1,151
            AREATABLE(INC,L1) = 0.D0
         END DO
         EVOL = 0.0
         GTSDPH = 40.0
         TOPMST = BOTTMS(L1)+GTSDPH
         TBELV = BOTTMS(L1)
         DO 1340 I=1,NCOL
            DO 1340 J=1,NROW
               IF(LKARR1(I,J,1).NE.L1) GO TO 1340
   !dep Revised estimate of DTHK to be thickness of top most
   !     layer 6/09/2009
               IF(BOTM(I,J,0).GT.TOPMST) TOPMST = BOTM(I,J,0)
   !      DTHK = BOTM(I,J,0) - BOTM(I,J,1)   RGN this was causing problems 7/8/11
   !      IF (DTHK.LE.GTSDPH) THEN
   !        TOPMST = BOTM(I,J,1)+DTHK
   !      ELSE
   !        TOPMST = BOTM(I,J,1)+GTSDPH
   !      END IF
1340     CONTINUE
         TBNC = (TOPMST-BOTTMS(L1))/150.0
   !dep Revised looping for computing lake stage, volume,
   !dep   and area Apr 2009.
   !dep   WRITE(IOUT,1315) TBELV, EVOL
         DO  INC=1,151
            IF (INC.GT.1) THEN
               VOLUMETABLE(INC,L1)=VOLUMETABLE(INC-1,L1)
            END IF
            DO I=1,NCOL
               DO J=1,NROW
                  LAKEFLG = 0
                  K = 1
                  MOSTBOT: DO WHILE (LAKEFLG.EQ.0)
                     IF(LKARR1(I,J,K).EQ.L1) THEN
                        LAKEFLG = K
                     END IF
                     IF(K.EQ.NLAY)EXIT MOSTBOT
                     K = K + 1
                  END DO MOSTBOT
                  IF(LAKEFLG.GT.0) THEN
                     K=LAKEFLG
                     FINDBOT: DO WHILE(LKARR1(I,J,K).GT.0)
                        K=K+1
                        IF(K.EQ.NLAY+1) EXIT
                     END DO FINDBOT
                     BOTIJ = BOTM(I,J,LBOTM(K-1))
                     IF(INC.EQ.1) THEN
                        IF(TBELV+1.0E-03.GT.BOTIJ) THEN
                           AREATABLE(INC,L1)=AREATABLE(INC,L1)+DELC(J)*DELR(I)
                           DEPTHTABLE(INC,L1)=TBELV
                        END IF
                     ELSE
                        IF (TBELV-BOTIJ.GT.1.0E-03) THEN
                           AREATABLE(INC,L1)=AREATABLE(INC,L1)+DELC(J)*DELR(I)
                           DEPTHTABLE(INC,L1)=TBELV
                           IF(ABS(TBELV-BOTIJ).GT.1.0E-03) THEN
                              VOLUMETABLE(INC,L1)=VOLUMETABLE(INC,L1)+&
                              &(DELC(J)*DELR(I))*TBNC
                           END IF
                        END IF
                     END IF
                  END IF
               END DO
            END DO
   !dep PRINT TABLE OF ELEVATION, VOLUME, AND AREA
            WRITE(IOUT,1315) DEPTHTABLE(INC,L1), VOLUMETABLE(INC,L1),&
            &AREATABLE(INC,L1)
            TBELV = TBELV + TBNC
         END DO
1315     FORMAT(3(1X,1PE13.5))
         WRITE(IOUT,1326)
1326     FORMAT(120X)
   !dep  set minimum and maximum lake stages for transient simulations
         IF(ISS.EQ.0) THEN
            SSMN(L1)=BOTTMS(L1)
            SSMX(L1)=TBELV
         END IF
1330  CONTINUE
1331  CONTINUE
   END IF
   IF(IUNITSFR.LE.0) THEN
      NDV=0
      NTRB=0
   END IF
   !
   !
   !--  READ LINKAGE PARAMETERS FOR COALESCING LAKES
   !
   !    FOR EACH CONNECTED LAKE SYSTEM, READ LAKE NUMBERS OF CENTER LAKES
   !    AND ADJOINING LAKES AND SILL ELEVATIONS.  ENTER CARD IMAGES
   !    FOR SUBLAKE SYSTEMS EVEN IF LINKED TO MAIN LAKE SYSTEM.  SYSTEMS
   !    MUST BE ORDERED HIERARCHICALLY.
   !
   ICMX = 0
   NCLS=0
   IF(IFREFM.EQ.0) THEN
      READ(IN,'(I5)') NSLMS
   ELSE
      READ(IN,*) NSLMS
   END IF
   WRITE(IOUT,680) NSLMS
680 FORMAT(/1X,'NUMBER OF CONNECTED LAKE SYSTEMS IN SIMULATION IS ',I3&
   &)
   IF(NSLMS.LE.0) GO TO 760
   DO 700 IS=1,NSLMS
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(16I5)',END=750) IC,(ISUB(IS,I),I=1,IC)
      ELSE
         READ(IN,*,END=750) IC,(ISUB(IS,I),I=1,IC)
      END IF
      IF(IC.LE.0) GO TO 750
      IF(IC.GT.ICMX) ICMX=IC
      ICS(IS)=IC
      IC1 = IC - 1
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(100F10.2)') (SILLVT(IS,I),I=1,IC1)
      ELSE
         READ(IN,*) (SILLVT(IS,I),I=1,IC1)
      END IF
      WRITE(IOUT,18) IS, ICS(IS), ISUB(IS,1)
18    FORMAT(/10X,'SYSTEM',I3//2X,'NUMBER OF LAKES IN SYSTEM',I5,&
      &'  CENTER LAKE NUMBER',I5//1X,'SUBLAKE NUMBER',3X,&
      &'SILL ELEVATION'/)
      DO 715 JK=2,IC
715   WRITE(IOUT,717) ISUB(IS,JK), SILLVT(IS,JK-1)
717   FORMAT(8X,I2,8X,F10.2)
700 CONTINUE
750 CONTINUE
   NCLS=IS-1
   WRITE(IOUT,751) NCLS
751 FORMAT(/1X,'READ DATA FOR',I5,' LAKE SYSTEMS'/)
760 CONTINUE
   !
   !----- READ LAKE PRECIPITATION, EVAPORATION, RUNOFF, AND WITHDRAWAL RATES.
   !      IF ITMP1 LT 0, SPECIFICATIONS FROM LAST STRESS PERIOD ARE USED.
   !
800 IF(ITMP1.GE.0) GO TO 801
   WRITE(IOUT,802)
802 FORMAT(1H0,'REUSING RECH,ET,WITHDRAWAL RATES FROM LAST STRESS PERI&
   &OD'/)
   GOTO 900
801 IF(ISS.NE.0.AND.KKPER.GT.1) WRITE(IOUT,7)
7  FORMAT(/1X,'LAKE',7X,'PRECIP',5X,'EVAP',5X,'RUNOFF',&
   &3X,'WITHDRAW',3X,'BOTTOM',5X,'AREA',5X,'SS MIN',3X,'SS MAX'&
   &/90('-'))
   IF(ISS.EQ.0.OR.KKPER.EQ.1) WRITE(IOUT,77)
77 FORMAT(/1X,'LAKE',7X,'PRECIP',5X,'EVAP',5X,'RUNOFF',&
   &3X,'WITHDRAW',3X,'BOTTOM',5X,'AREA',5X,/70('-'))
   IF (IUNITGWT.GT.0) WRITE (IOUTS,8)
8  FORMAT (//1X,'LAKE',4X,'SOLUTE',6X,'CPPT',6X,'CRNF',6X,'CAUG'/)
   !
   DO 300 LM=1,NLAKES
      IF(IFREFM.EQ.0) THEN
         IF(ISS.NE.0.AND.KKPER.GT.1) READ(IN,'(6F10.4)') PRCPLK(LM),&
         &EVAPLK(LM),RNF(LM),WTHDRW(LM),SSMN(LM),SSMX(LM)
         IF(ISS.EQ.0.OR.KKPER.EQ.1) READ(IN,'(6F10.4)') PRCPLK(LM),&
         &EVAPLK(LM),RNF(LM),WTHDRW(LM)
      ELSE
         IF(ISS.NE.0.AND.KKPER.GT.1) READ(IN,*) PRCPLK(LM),EVAPLK(LM),&
         &RNF(LM),WTHDRW(LM),SSMN(LM),SSMX(LM)
         IF(ISS.EQ.0.OR.KKPER.EQ.1) READ(IN,*) PRCPLK(LM),EVAPLK(LM),&
         &RNF(LM),WTHDRW(LM)
      END IF
   !
   !--EDM: SET FOLLOWING VALUES FOR LMT
      IF(IUNIT(49).NE.0) THEN
         IF(PRCPLK(LM).NE.0.AND.LKFLOWTYPE(3).EQ.'NA') THEN
            LKFLOWTYPE(3)='PRECIP'
            NLKFLWTYP = NLKFLWTYP + 1
         ENDIF
         IF(EVAPLK(LM).NE.0.AND.LKFLOWTYPE(4).EQ.'NA') THEN
            LKFLOWTYPE(4)='EVAP'
            NLKFLWTYP = NLKFLWTYP + 1
         ENDIF
         IF(RNF(LM).NE.0.AND.LKFLOWTYPE(5).EQ.'NA') THEN
            LKFLOWTYPE(5)='RUNOFF'
            NLKFLWTYP = NLKFLWTYP + 1
         ENDIF
         IF(WTHDRW(LM).NE.0.AND.LKFLOWTYPE(6).EQ.'NA') THEN
            LKFLOWTYPE(6)='WITHDRAW'
            NLKFLWTYP = NLKFLWTYP + 1
         ENDIF
      ENDIF
   !
      IF(ISS.NE.0.AND.KKPER.GT.1)WRITE(IOUT,9)LM,PRCPLK(LM),EVAPLK(LM)&
      &,RNF(LM),WTHDRW(LM),BOTTMS(LM),BGAREA(LM),SSMN(LM),SSMX(LM)
9     FORMAT(1X,I3,4X,1P,3E10.3,1X,5E10.3)
      IF(ISS.EQ.0.OR.KKPER.EQ.1)WRITE(IOUT,9)LM,PRCPLK(LM),EVAPLK(LM),&
      &RNF(LM),WTHDRW(LM),BOTTMS(LM),BGAREA(LM)
      IF(IUNITGWT.LE.0) GO TO 300
      DO 850 ISOL=1,NSOL
         IF(IFREFM.EQ.0) THEN
            IF(WTHDRW(LM).LT.0.0) THEN
               READ(IN,'(3F10.4)')CPPT(LM,ISOL),CRNF(LM,ISOL),&
               &CAUG(LM,ISOL)
            ELSE
               READ(IN,'(2F10.4)')CPPT(LM,ISOL),CRNF(LM,ISOL)
            END IF
         ELSE
            IF(WTHDRW(LM).LT.0.0) THEN
               READ(IN,*) CPPT(LM,ISOL),CRNF(LM,ISOL),CAUG(LM,ISOL)
            ELSE
               READ(IN,*) CPPT(LM,ISOL),CRNF(LM,ISOL)
            END IF
         END IF
         IF(WTHDRW(LM).LT.0.0)WRITE(IOUTS,840) LM,ISOL,&
         &CPPT(LM,ISOL),CRNF(LM,ISOL),CAUG(LM,ISOL)
         IF(WTHDRW(LM).GE.0.0)&
         &WRITE(IOUTS,841) LM,ISOL,CPPT(LM,ISOL),CRNF(LM,ISOL)
840      FORMAT(1X,I3,6X,I3,4X,1P,3E10.2)
841      FORMAT(1X,I3,6X,I3,4X,1P,2E10.2)
850   CONTINUE
   !        WRITE (IOUTS,'(/)')
300 CONTINUE
   WRITE (IOUT,'(/)')
   !
   !------Define Initial Lake Volume & Initialize Cumulative Budget Terms
   IF(KKPER.EQ.1) THEN
   !dep revised calculation of initial lake volume July 2009
      STGINIT=0.0D0
      DO 8400 LK=1,NLAKES
   !dep 8400    VOL(LK)=0.0
         STGINIT=STAGES(LK)
         VOL(LK)=VOLTERP(STGINIT,LK)
         VOLINIT(LK)=VOL(LK)
         VOLOLDD(LK)=VOL(LK)
8400  CONTINUE
      DO 8450 LK=1,NLAKES
         CUMPPT(LK)=0.0
         CUMEVP(LK)=0.0
         CUMRNF(LK)=0.0
         CUMGWI(LK)=0.0
         CUMGWO(LK)=0.0
         CUMSWI(LK)=0.0
         CUMSWO(LK)=0.0
         CUMWDR(LK)=0.0
         CUMFLX(LK)=0.0
8450  CONTINUE
      DO 8900 L=1,LKNODE
         IL=ILAKE(1,L)
         IR=ILAKE(2,L)
         IC=ILAKE(3,L)
         LAKE=ILAKE(4,L)
   !------Convert ILAKE(5,L):  1 and 2 are type 1,  3 and 4 are type 2,
   !        6 is type 0
         ITYPE = (ILAKE(5,L)+1)/2
         IF(ITYPE.EQ.3) ITYPE=0
         IF(ITYPE.NE.0) GO TO 8900
         IF(IL.GT.1) BOTLK = BOTM(IC,IR,LBOTM(IL-1))
         IF(IL.EQ.NLAY.AND.LKARR1(IC,IR,IL).GT.0)&
         &BOTLK = BOTM(IC,IR,LBOTM(IL))
8900  CONTINUE
   ENDIF

900 IF (IUNITBCF.GT.0) THEN  ! rsr, moved if block from main
      CALL SGWF2LAK7BCF7RPS(LWRT)
   ELSE IF (IUNITLPF.GT.0) THEN
      CALL SGWF2LAK7LPF7RPS(LWRT)
   ELSE IF (IUNITHUF.GT.0) THEN
      CALL SGWF2LAK7HUF7RPS(LWRT)
   ELSE IF (IUNITUPW.GT.0) THEN
      CALL SGWF2LAK7UPW1RPS(LWRT)
   ELSE
      WRITE (IOUT, *) 'LAK Package requires BCF, LPF, UPW, or HUF'
      CALL USTOP(' ')
   END IF
   IF (IUNITSFR.GT.0) CALL SGWF2LAK7SFR7RPS()
   !
   !7------RETURN
   RETURN
END
   !
SUBROUTINE GWF2LAK7AD(KKPER,KKSTP,IUNITGWT,IGRID)
   !
   !------OLD VERSION 7.1 JUNE 2006 GWF2LAK7AD; REVISED FEBRUARY 6, 2012
   !------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
   !------NEW VERSION NUMBER 1.3.0, 7/01/2022
   !
   !     ******************************************************************
   !     ADVANCE TO NEXT TIME STEP FOR TRANSIENT LAKE SIMULATION, AND COPY
   !             INITIAL LAKE STAGES TO STGOLD FOR STEADY STATE.
   !     ******************************************************************
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFLAKMODULE, ONLY: NLAKES, LKNODE, FLOB, STAGES,&
   &STGNEW, STGOLD, VOLOLDD, VOLOLD, VOLINIT,&
   &BOTTMS, IDIV, STGOLD2, NDV, ISTARTLAK
   USE GWFSFRMODULE, ONLY: DLKSTAGE
   USE GLOBAL,       ONLY: IOUT
   !     ------------------------------------------------------------------
   !     FUNCTIONS
   !     ------------------------------------------------------------------
   DOUBLE PRECISION VOLTERP
   EXTERNAL VOLTERP
   !     ------------------------------------------------------------------
   !
   !------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2LAK7PNT(IGRID)
   !
   !1 --- COPY INITIAL LAKE STAGES TO STGOLD.
   ! RGN COMBINED IF AND ADDED VOLOLDD 4/17/09
   !dep  initialized VOLINIT and VOLOLD to VOLOLDD 6/4/2009
   DO I=1,NLAKES
      IF( ISTARTLAK==0 ) THEN
         IF ( I==NLAKES ) ISTARTLAK = 1
         STGOLD(I)=STAGES(I)
         VOLOLDD(I)=VOLTERP(STGOLD(I),I)
         VOLOLD(I) = VOLOLDD(I)
         VOLINIT(I) = VOLOLDD(I)
         STGNEW(I)=STAGES(I)
      ELSE
         STGOLD2(I)=STGNEW(I)
         STGOLD(I)=STGNEW(I)
         VOLOLDD(I)=VOLTERP(STGOLD(I),I)
         VOLOLD(I)=VOLOLDD(I)
      END IF
   ! Moved this code from 7FM  10/19/10
      DO IDV=1,NDV
         INODE=IDIV(I,IDV)
         IF (INODE.GT.0) THEN
            IF( DLKSTAGE(1,INODE).LT.DBLE(BOTTMS(I))) THEN
               WRITE(IOUT,971)I,BOTTMS(I),&
               &DLKSTAGE(1,INODE),INODE
               CALL USTOP(' ')
            END IF
         END IF
      END DO
      ! To hear.
   END DO
971 FORMAT(' BOTTOM ELEVATION OF LAKE ',I5,' IS ', F10.2,&
   &' AND IS ABOVE OUTLET ELEVATION OF ', F10.2,&
   &' FOR STREAM SEGMENT ',I5,/1X,&
   &' THIS WILL CAUSE PROBLEMS IN COMPUTING LAKE',&
   &' STAGE USING THE NEWTON METHOD. '/1X,&
   &' ELEVATION OF STREAM OUTLET MUST BE GREATER'&
   &' THAN OR EQUAL TO THE LOWEST ELEVATION OF THE',&
   &' LAKE.',/1X,'*****PROGRAM STOPPING'/)
   !2 ----- IF NOT FIRST TIME STEP, OR FIRST STRESS PERIOD, UPDATE
   !           STGOLD BY STGNEW.
   ! RGN MOVED TO ABOVE. STGOLD SHOULD BE UPDATED EVERY TIME STEP! 4/17/09
   !      IF (KKPER.NE.1.OR.KKSTP.NE.1) THEN
   !            DO 30 K=1,NLAKES
   !               STGOLD(K)=STGNEW(K)
   !               VOLOLD(K)=VOLTERP(STGOLD(K),K))
   !30             STGOLD2(K)=STGNEW(K)
   !      ENDIF
   !
   !-----Initialize FLOB array (stores cell by cell flux between lake and
   !                            aquifer)
   IF (IUNITGWT.GT.0) THEN
      DO 50 LK=1,LKNODE
50    FLOB(LK)=0.0
   END IF
   !
   !3------RETURN
   RETURN
END
   !
SUBROUTINE GWF2LAK7ST(NFLG,IGRID)
   !   ********************************************************************
   !   SET IBOUND VALUES SO THAT RECHARGE AND EVAPOTRANSPIRATION (ET) WILL
   !   BE ASSIGNED CORRECTLY UNDERNEATH DRYING LAKES (NFLG = 0), OR RESET
   !   IBOUND AFTER RECHARGE AND ET ARE COMPUTED (NFLG = 1).
   !   ********************************************************************
   !
   !   SPECIFICATIONS:
   !
   !-----------------------------------------------------------------------
   USE GWFLAKMODULE, ONLY: LKNODE, ILAKE, STGOLD
   USE GLOBAL,       ONLY: IBOUND, LBOTM, BOTM
   !-----------------------------------------------------------------------
   !
   !------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2LAK7PNT(IGRID)

   IF(LKNODE.EQ.0) RETURN
   DO 10 L=1,LKNODE
   !  Convert ILAKE(5,L):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      ITYPE = (ILAKE(5,L)+1)/2
      IF(ITYPE.EQ.3) ITYPE=0
   !
   !-------ONLY CHANGE IBOUND FOR VERTICALLY ADJACENT NODE FACES
      IF(ITYPE.NE.0) GO TO 10
      IL = ILAKE(1,L)
      IR = ILAKE(2,L)
      IC = ILAKE(3,L)
   !
   !-------RESET AFTER EXECUTING RECHARGE OR ET ROUTINES
      IF(NFLG.EQ.1) GO TO 8
   !
   !-------RESET BEFORE EXECUTING RECHARGE OR ET ROUTINES
      IF ( IL.GT.1 ) THEN       !RGN 5/21/12 added IF test
         IBOUND(IC,IR,IL-1) = -7
      ELSE                      !RGN
         IBOUND(IC,IR,IL) = -7   !RGN
      END IF                    !RGN
   !
   !-------THIS IS THE CORRECT ASSIGNMENT IF PORTION OF LAKE IN COLUMN
   !       IS WET.
      LAKE = ILAKE(4,L)
      IF(STGOLD(LAKE).GT.BOTM(IC,IR,LBOTM(IL)-1)) GO TO 10
   !
   !-------IF PORTION OF LAKE IN NODE IS DRY, LET RECHARGE AND ET BE
   !       APPLIED TO THE AQUIFER NODE UNDERNEATH THE LAKE BY SETTING
   !       IBOUND EQUAL TO 0.
8     IF ( il.GT.1 ) THEN        !RGN 5/21/12 added IF test
   !     8  IBOUND(IC,IR,IL-1) = 0  !RGN
         IBOUND(IC,IR,IL-1) = 0   !RGN
      ELSE                       !RGN
         IBOUND(IC,IR,IL) = 0     !RGN
      END IF                     !RGN
10 CONTINUE
   !
   !3------RETURN
   RETURN
END
   !
SUBROUTINE GWF2LAK7FM(KKITER,KKPER,KKSTP,IUNITSFR,IUNITUZF,IGRID)
   !
   !------OLD USGS VERSION 7.1; JUNE 2006 GWF2LAK7FM;
   !------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
   !------NEW VERSION NUMBER 1.3.0, 7/01/2022
   !     ******************************************************************
   !     ADD LAKE TERMS TO RHS AND HCOF IF SEEPAGE OCCURS IN MODEL CELLS
   !     ******************************************************************
   !
   USE GWFLAKMODULE
   USE GLOBAL,       ONLY: NLAY, IBOUND, IOUT, ISSFLG,&
   &DELR, DELC, LBOTM, BOTM, HNEW, HCOF, RHS
   USE GWFBASMODULE, ONLY: DELT
   USE GWFSFRMODULE, ONLY: STRIN, STROUT, FXLKOT, SEG
   USE GWFUZFMODULE, ONLY: IUZFBND
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !dep  Added functions for interpolating between areas, derivatives,
   !dep     and outflow rates
   !     ------------------------------------------------------------------
   !     FUNCTIONS
   !     -----------------------------------------------------------------
   IMPLICIT NONE
   DOUBLE PRECISION FINTERP, DERIVTERP, OUTFLWTERP, VOLTERP, SURFTERP
   EXTERNAL FINTERP, DERIVTERP, OUTFLWTERP, VOLTERP, SURFTERP
   DOUBLE PRECISION STGTERP, FXLKOT_TERP
   EXTERNAL STGTERP, FXLKOT_TERP
   !     -----------------------------------------------------------------
   !     ARGUMENTS
   !     -----------------------------------------------------------------
   INTEGER, INTENT(IN) :: KKITER, KKPER, IUNITSFR, IUNITUZF, IGRID,&
   &KKSTP
   !dep  added runoff and flobo3
   !      REAL :: RUNOFF
   !dep  added unsaturated flow beneath lakes flag as a local variable
   INTEGER ISS, LK, ITRIB, INODE, LAKE, MTER, IICNVG, L1, MAXITER
   INTEGER NCNV, LL, II, L, IC, IR, IL, ITYPE
   INTEGER INOFLO, IDV, IL1
   !dep  added SURFDPTH, CONDMX,BOTLKUP,BOTLKDN  3/3/2009
   DOUBLE PRECISION BOTLK,BOTCL,CONDUC,H,FLOBOT,STGON,&
   &FLOBO3,THET1,CLOSEZERO,&
   &SURFDPTH,CONDMX,BOTLKUP,BOTLKDN, FLOTOUZF,&
   &VOL2,WITHDRW3
   !!     3                 VOL2,RAMPGW,RAMPSTGO,RAMPSTGN,
   !!     4                 RAMPSTGON,HTEMP,WITHDRW3
   !dep  added double precision variables
   DOUBLE PRECISION RESID1, RESID2, DERIV, DSTAGE, DLSTG, Botlake,&
   &Splakout, dy, SRFPT, HD
   DOUBLE PRECISION AREA, RAIN, EV, THCK, SSMN1, SSMX1,&
   &OUTFLOW, DSTG, VOLNEW1, VOLNEW2, STAGE2
   !      DOUBLE PRECISION RUNF
   !      PARAMETER(CLOSEZERO = 1.0E-07)
   !     ------------------------------------------------------------------
   !------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2LAK7PNT(IGRID)
   CLOSEZERO = 1.0D-09
   DLSTG = 0.00001D0
   SURFDPTH = DBLE(SURFDEPTH)
   VOLNEW1 = 0.0
   VOLNEW2 = 0.0
   !1------IF LKNODE<=0 THERE ARE NO LAKE NODES. RETURN.
   IF (LKNODE.LE.0) RETURN
   ISS = ISSFLG(KKPER)
   !
   !2------PROCESS EACH CELL IN THE ILAKE LIST.
   !dep   added STGITER, and STGNEW to INITIALIZATION.
   DO LK=1,NLAKES
      IF(KKITER.EQ.1)THEN
         STGITER(LK) = STGOLD(LK)
         STGNEW(LK) = STGOLD(LK)
      END IF
      NCNCVR(LK) = 0
      LIMERR(LK) = 0
      SURFIN(LK)=0.0
      DSRFOT(LK)=0.0
      GWRATELIM(LK) = 0.0
      SURFOT(LK)=0.0
   END DO
   BOTLKUP = 0.0D0
   BOTLKDN = 0.0D0
   CONDMX = 0.0D0
   VOL2 = 0.0D0
   WITHDRW3 = 0.0D0
   !
   !2A --- SUM UP INFLOWS FROM INFLOWING STREAM REACHES.
   IF (IUNITSFR.GT.0) THEN
      DO 200 LK=1,NLAKES
         DO 200 ITRIB=1,NTRB
            INODE=ITRB(LK,ITRIB)
            IF (INODE.LE.0) GO TO 200
            SURFIN(LK)=SURFIN(LK)+STRIN(INODE)
200   CONTINUE
   END IF
   !
   !2B --- SUM UP OVERLAND RUNOFF INTO LAKE.
   DO LAKE = 1,NLAKES
   ! EDM - Add indices to RUNF and RUNOFF because these terms need to be
   !       saved for each lake in the simulation for writing to the FTL file by LMT

      IF(RNF(LAKE).GE.0.0) RUNF(LAKE) = RNF(LAKE)
      IF(RNF(LAKE).LT.0.0) RUNF(LAKE) =-RNF(LAKE)*PRCPLK(LAKE)&
      &*BGAREA(LAKE)
      IF (IUNITUZF.GT.0) THEN
         RUNOFF(LAKE) = OVRLNDRNF(LAKE)
      ELSE
         RUNOFF(LAKE) = 0.0
      END IF
   !
   !2C --- SUM OF BOTH STREAMFLOW IN AND OVERLAND RUNOFF.
   !         (INCLUDES LAKE VOLUME).
      FLWIN(LAKE) = SURFIN(LAKE)+RUNF(LAKE)+RUNOFF(LAKE)+&
      &VOLTERP(STGOLD(LAKE),LAKE)/DELT
   END DO
   !
   !3  --- TIME WEIGHTING FACTOR.

   THET1 = DBLE(THETA)
   IF ( ISS==1 ) THET1 = 1.0
   IF(THET1-0.5.LT.-CLOSEZERO) THET1=0.0D0
   MTER = NSSITR
   IICNVG = 0
   L1 = 0
   MAXITER = MTER
   IF ( THET1.LT.CLOSEZERO ) MAXITER = 0
   !
   !4------MASTER LOOP FOR COMPUTING GROUNDWATER INTERACTION WITH LAKES.
   !dep&rgn   Revised April-August 2009
   CONVERGE: DO WHILE (L1<=MAXITER)
      L1 = L1 + 1
   !
   !4B-----NCNV IS LAKE CONVERGENCE FLAG. IT IS 0 WHEN ALL LAKES HAVE
   !         CONVERGED. NCNCVR(LAKE) IS CONVERGENCE FLAG FOR EACH LAKE.
      NCNV = 0
      DO LAKE=1,NLAKES
         IF(NCNCVR(LAKE).EQ.0) NCNV = 1
      END DO
      IF ( L1.GT.MAXITER ) NCNV = 0
      IF ( THET1.LT.CLOSEZERO) NCNV = 0
      IF ( NCNV.EQ.0 ) IICNVG = 1
   !
   !4C-----INITIALIZE VARIABLES.
      DO LL=1,NLAKES
         SUMCNN(LL) = 0.0
         SUMCHN(LL) = 0.0
         EVAP(LL)=0.0D0
         PRECIP(LL)=0.0D0
         SEEP(LL)=0.0D0
         SEEP3(LL)=0.0D0
         SEEPUZ(LL)=0.0D0
         SURFA(LL)=0.0
         WITHDRW(LL) = WTHDRW(LL)
         FLWITER(LL) = FLWIN(LL)
         FLWITER3(LL) = FLWIN(LL)
         IF ( ISS==1 ) THEN
            FLWITER(LL) = 1.0E10
            FLWITER3(LL) = 1.0E10
         END IF
      END DO
   !
   !5------II LOOP BALANCES INFLOWS AND OUTFLOWS TO/FROM A LAKE.
   !   WHEN II=1, FLOW INTO LAKE FROM ALL SOURCES ARE CALCULATED.
   !   WHEN II=2, SEEPAGE TO AND FROM LAKE AND RESIDUAL TERMS ARE ADDED TO
   !   GROUNDWATER MATRIX.LAKE SEEPAGE TO GROUNDWATER LIMITED TO AVAILABLE
   !   WATER IN LAKE.
      DO II = 1,2
         DO L=1,LKNODE
            IL=ILAKE(1,L)
            IR=ILAKE(2,L)
            IC=ILAKE(3,L)
   !
   !5B------DETERMINE LAKE AND NODAL LAYER,ROW,COLUMN NUMBER.
            LAKE=ILAKE(4,L)
            IF( ISS.EQ.1 ) STGOLD(LAKE)=STGNEW(LAKE)
            STGON = (1.0D0-THET1)*STGOLD(LAKE) + THET1*STGNEW(LAKE)
            ITYPE = (ILAKE(5,L)+1)/2
            IF(ITYPE.EQ.3) ITYPE=0
            AREA = DELR(IC)*DELC(IR)
            IF(IL.GT.1) THEN
               BOTLK = BOTM(IC,IR,LBOTM(IL-1))
            ELSE
               BOTLK = BOTM(IC,IR,LBOTM(IL))
            END IF
            BOTCL = BOTM(IC,IR,LBOTM(IL))
            IF(IL.EQ.NLAY.AND.CNDFCT(L).EQ.0.0) BOTLK = BOTCL
            RAIN = PRCPLK(LAKE)
            EV = EVAPLK(LAKE)
            CONDUC=CNDFCT(L)
            FLOBOT = 0.0D0
            FLOBO3 = 0.0D0
            FLOTOUZF = 0.0D0
            IF ( NCNCVR(LAKE).NE.2 ) THEN
               IL1 = IL
               IF ( ITYPE.EQ.0 ) THEN
                  DO WHILE (IL1 .LE. NLAY)
                     IF( IBOUND(IC,IR,IL1).GT.0 ) THEN
                        EXIT
                     ELSE
                        IL1 = IL1 + 1
                     END IF
                  END DO
                  IF ( IL1.GT.NLAY ) IL1 = NLAY
               END IF
               IF( IBOUND(IC,IR,IL1).LE.0 ) THEN
                  !               IF ( CONDUC.GT.CLOSEZERO ) WRITE(IOUT,506) L,IC,IR,IL
                  CONDUC = 0.0
               END IF
               IF(CONDUC.GT.0.0) THEN
                  H=HNEW(IC,IR,IL1)   !RGN this was IL and not IL1. 1/9/17
                  INOFLO = 0
   !
   !9A------CALCULATE SEEPAGE.
   !
                  CALL GET_FLOBOT(IC, IR, IL1, ITYPE, INOFLO,CONDUC,&
                  &FLOBOT,FLOBO3,FLOTOUZF,DLSTG,CLOSEZERO,H,&
                  &THET1,ISS,LAKE,II,SURFDPTH,AREA,IUNITUZF,&
                  &BOTLK,BOTCL,L1)
   !
   !9B------ADD SEEPAGE RATES AND RESIDUAL TERMS TO GROUNDWATER MATRIX
   !         WHEN ITYPE = 0.
                  IF( NCNV == 0 .AND. II==2 ) THEN
                     IF ( L==LKNODE ) NCNCVR(LAKE) = 2
                     IF ( ITYPE.EQ.0 ) THEN
                        IF ( INOFLO==1 ) THEN
                           RHS(IC,IR,IL1)=RHS(IC,IR,IL1)-FLOBOT
                        ELSE
                           IF (STGON-BOTLK.GT.CLOSEZERO)THEN
                              IF (H.LE.BOTLK) THEN
                                 IF (IUNITUZF.EQ.0) THEN
                                    RHS(IC,IR,IL1)=RHS(IC,IR,IL1)-FLOBOT
                                 ELSE IF (IUZFBND(IC,IR).EQ.0 ) THEN
                                    RHS(IC,IR,IL1)=RHS(IC,IR,IL1)-FLOBOT
                                 END IF
                              ELSE
                                 RHS(IC,IR,IL1)=RHS(IC,IR,IL1) - STGON*CONDUC
                                 HCOF(IC,IR,IL1)=HCOF(IC,IR,IL1) - CONDUC
                              END IF
                           ELSE
                              IF ( H.GT.BOTLK ) THEN
                                 RHS(IC,IR,IL1)=RHS(IC,IR,IL1) - BOTLK*CONDUC
                                 HCOF(IC,IR,IL1)=HCOF(IC,IR,IL1) - CONDUC
                              END IF
                           END IF
                        END IF
   !
   !9BC-----ADD SEEPAGE RATES AND RESIDUAL TERMS TO GROUNDWATER MATRIX
   !         WHEN ITYPE = 1 OR 2.
                     ELSE IF ( ITYPE.EQ.1.OR.ITYPE.EQ.2 ) THEN
                        IF ( INOFLO==1 ) THEN
                           RHS(IC,IR,IL1)=RHS(IC,IR,IL1)-FLOBOT
                        ELSE
                           IF( IBOUND(IC,IR,IL1).GT.0 ) THEN
                              HD = H
                              IF( H.GT.BOTM(IC,IR,LBOTM(IL1)-1) )&
                              &HD = BOTM(IC,IR,LBOTM(IL1)-1)
   !
   !9D------CONDUCTANCE ACROSS VERTICAL CELL FACE DEPENDENT ON
   !          SATURATED THICKNESS.
                              THCK = HD - BOTCL
                              IF( THCK.LE.0.0 ) THCK = 0.0
                           END IF
                           IF (STGON.GT.BOTCL.OR.THCK.GT.0.0) THEN
                              IF(STGON-BOTCL.GT.CLOSEZERO)THEN
                                 IF ( H-BOTCL.GT.CLOSEZERO ) THEN
                                    RHS(IC,IR,IL1)=RHS(IC,IR,IL1) - STGON*CONDUC
                                    HCOF(IC,IR,IL1)=HCOF(IC,IR,IL1) - CONDUC
                                 END IF
                              ELSE IF( H-BOTCL.GT.CLOSEZERO ) THEN
                                 RHS(IC,IR,IL1)=RHS(IC,IR,IL1) - BOTCL*CONDUC
                                 HCOF(IC,IR,IL1)=HCOF(IC,IR,IL1) - CONDUC
                              END IF
                           END IF
                        END IF
                     END IF
   !
   !10------COMPUTE NET SEEPAGE THROUGH LAKEBED SEDIMENTS
   !          WHEN II=2.
                  END IF
                  IF (II==2)THEN
                     SEEP(LAKE)=SEEP(LAKE)-FLOBOT-FLOTOUZF
                     SEEP3(LAKE)=SEEP3(LAKE)-FLOBO3
                  END IF
               END IF
            END IF
         END DO
      END DO
   !
   !11------ONLY COMPUTE LAKE LEVEL AFTER SCANNING THRU ALL NODES OF
   !          A LAKE.
      DO LAKE=1,NLAKES
   !
   !11B-----SET STGITER TO STGNEW WHEN THET1>0 AND TO STGOLD
   !          WHEN THET=0.
         IF( THET1.GT.CLOSEZERO ) THEN
            IF ( L1.LT.MTER) STGITER(LAKE) = STGNEW(LAKE)
         ELSE
            STGITER(LAKE) = STGOLD(LAKE)
         END IF
   !
   !12------COMPUTE EVAPORATION AND PRECIPITATION USING STGOLD AND
   !          ADD PRECIPITATION TO FLWITER AND FLWITER3.
         SURFA(LAKE)=FINTERP(STGNEW(LAKE),LAKE)
   ! The following lines were added to make sure the volume of precip
   ! is a constant from PRMS
         IF ( IGSFLOWLAK == 0 ) THEN
            EVAP(LAKE)=EVAPLK(LAKE)*SURFA(LAKE)
            PRECIP(LAKE)=PRCPLK(LAKE)*SURFA(LAKE)
            SRFPT=FINTERP(STGNEW(LAKE)+DLSTG,LAKE)
            EVAP3(LAKE)=EVAPLK(LAKE)*SRFPT
            PRECIP3(LAKE)=PRCPLK(LAKE)*SRFPT
         ELSE
            EVAP(LAKE) = EVAPLK(LAKE)
            PRECIP(LAKE) = PRCPLK(LAKE)
            EVAP3(LAKE) = EVAPLK(LAKE)
            PRECIP3(LAKE)=PRCPLK(LAKE)
         END IF
         FLWITER(LAKE) = FLWITER(LAKE) + PRECIP(LAKE)
         FLWITER3(LAKE) = FLWITER3(LAKE) + PRECIP3(LAKE)
   !
   !13------LIMIT WITHDRW TO LAKE INFLOW WHEN WITHDRAWALS EXCEED
   !           INFLOW (INCLUDING AVAILABLE LAKE STORAGE).
         IF(WITHDRW(LAKE).GE.FLWITER(LAKE)) THEN
            WITHDRW(LAKE) = FLWITER(LAKE)
            FLWITER(LAKE) = 0.0D0
         ELSE
            FLWITER(LAKE)  = FLWITER(LAKE) - WITHDRW(LAKE)
         END IF
         IF(WITHDRW(LAKE).GE.FLWITER3(LAKE)) THEN
            WITHDRW3 = FLWITER3(LAKE)
            FLWITER3(LAKE) = 0.0
         ELSE
            WITHDRW3 = WITHDRW(LAKE)
            FLWITER3(LAKE)  = FLWITER3(LAKE) - WITHDRW3
         END IF
   !
   !14------LIMIT EVAPORATION TO LAKE INFLOW WHEN EVAPORATION EXCEEDS
   !          INFLOW (INCLUDING AVAILABLE LAKE STORAGE AND WITHDRAWALS).
         IF ( EVAP(LAKE)>=FLWITER(LAKE) ) THEN
            EVAP(LAKE)=FLWITER(LAKE)
            FLWITER(LAKE) = 0.0D0
         ELSE
            FLWITER(LAKE) = FLWITER(LAKE) - EVAP(LAKE)
         END IF
         IF ( EVAP3(LAKE)>=FLWITER3(LAKE) ) THEN
            EVAP3(LAKE)=FLWITER3(LAKE)
            FLWITER3(LAKE) = 0.0D0
         ELSE
            FLWITER3(LAKE) = FLWITER3(LAKE) - EVAP3(LAKE)
         END IF
         SSMN1 = SSMN(LAKE)
         SSMX1 = SSMX(LAKE)
   !
   !15-----SUM UP OUTFLOWS FROM OUTFLOWING STREAM REACHES.
         DSRFOT(LAKE) = 0.0D0
         SURFOT(LAKE) = 0.0D0
         ! 11/10 Outflow is zero when all is lost to ET. ******* NEW(.AND.FLWITER(LAKE).GT.CLOSEZERO)
         IF(IUNITSFR.GT.0.AND.FLWITER(LAKE).GT.CLOSEZERO) THEN
            DO IDV=1,NDV
               INODE=IDIV(LAKE,IDV)
               IF (INODE.GT.0) THEN
                  Splakout = DBLE(SEG(2,INODE))
                  DSTAGE = STGITER(LAKE)
                  IF ( SEG(2,INODE).LE.1.0e-6 ) THEN
                     Botlake = DBLE(BOTTMS(LAKE))
                     DSRFOT(LAKE) = DSRFOT(LAKE) + DERIVTERP(DSTAGE,&
                     &INODE)
                     STROUT(INODE) = OUTFLWTERP(DSTAGE,INODE)
                  ELSE
   ! Set Botlake to elevation of stream channel for specified flow diversion.
                     Botlake = SEG(8,INODE)
                     STROUT(INODE) = 0.0
                     FXLKOT(INODE) = FXLKOT_TERP(RAMP,DSTAGE,Botlake,&
                     &Splakout,dy)
                     DSRFOT(LAKE) = DSRFOT(LAKE) + dy
                  END IF
                  SURFOT(LAKE)= SURFOT(LAKE) + STROUT(INODE)+&
                  &FXLKOT(INODE)
                  IF(SURFOT(LAKE).LT.CLOSEZERO )SURFOT(LAKE)=0.0
               END IF
            END DO
         END IF
         IF ( THET1.GT.CLOSEZERO ) THEN
            IF ( NCNV == 1 ) THEN
   ! Calc overland flow again.
   ! EDM - Add indices to RUNF and RUNOFF because these terms need to be
   !       saved for each lake in the simulation for writing to the FTL file by LMT
               IF(RNF(LAKE).GE.0.0) RUNF(LAKE) = RNF(LAKE)
               IF(RNF(LAKE).LT.0.0) RUNF(LAKE) =-RNF(LAKE)*&
               &PRCPLK(LAKE)*BGAREA(LAKE)
               IF (IUNITUZF.GT.0) THEN
                  RUNOFF(LAKE) = OVRLNDRNF(LAKE)
               ELSE
                  RUNOFF(LAKE) = 0.0
               END IF
   !
   !2C --- SUM OF BOTH STREAMFLOW IN AND OVERLAND RUNOFF.
   !         (INCLUDES LAKE VOLUME).
               FLWIN(LAKE) = SURFIN(LAKE)+RUNF(LAKE)+RUNOFF(LAKE)+&
               &VOLTERP(STGOLD(LAKE),LAKE)/DELT
   !
   !16-----COMPUTE NEW LAKE STAGE USING NEWTON METHOD AND THET1>0.
   !
   !16B----COMPUTE RESIDUALS FOR TRANSIENT SIMULATIONS.
               IF(ISS.EQ.0) THEN
                  VOLNEW1 = VOLTERP(STGNEW(LAKE),LAKE)
                  RESID1 = (PRECIP(LAKE)-EVAP(LAKE)+RUNF(LAKE)+&
                  &RUNOFF(LAKE)-WITHDRW(LAKE)+SURFIN(LAKE)-&
                  &SURFOT(LAKE)+SEEP(LAKE))-&
                  &(VOLNEW1-VOLOLDD(LAKE))/DELT
                  OUTFLOW = SURFOT(LAKE)+ DSRFOT(LAKE)*DLSTG
                  IF(OUTFLOW.LT.0.0)SURFOT(LAKE)=0.0
                  STAGE2 = STGNEW(LAKE)+DLSTG
                  VOLNEW2 = VOLTERP(STAGE2,LAKE)
                  RESID2 = (PRECIP3(LAKE)-EVAP3(LAKE)+RUNF(LAKE)+&
                  &RUNOFF(LAKE)-WITHDRW3+SURFIN(LAKE)-OUTFLOW+&
                  &SEEP3(LAKE))-(VOLNEW2-VOLOLDD(LAKE))/DELT
   !
   !16C----COMPUTE RESIDUALS FOR STEADY STATE SIMULATIONS.
               ELSE
                  RESID1 = (PRECIP(LAKE)-EVAP(LAKE)+RUNF(LAKE)+&
                  &RUNOFF(LAKE)-WITHDRW(LAKE)+SURFIN(LAKE)-&
                  &SURFOT(LAKE)+SEEP(LAKE))
                  OUTFLOW = SURFOT(LAKE)+ DSRFOT(LAKE)*DLSTG
                  RESID2 = (PRECIP3(LAKE)-EVAP3(LAKE)+RUNF(LAKE)+&
                  &RUNOFF(LAKE)-WITHDRW3+SURFIN(LAKE)-&
                  &OUTFLOW+SEEP3(LAKE))
               END IF
   !
   !16D----DETERMINE DERIVATIVE AND COMPUTE NEW LAKE STAGE.
               IF ( DABS(RESID2-RESID1).GT.SMALLTOL ) THEN
                  DERIV = (RESID2-RESID1)/(DLSTG)
                  DSTG = RESID1/DERIV
                  STGNEW(LAKE) = STGITER(LAKE) - DSTG
                  DSTG = ABS(DSTG)
               ELSE
   !16E----LINEAR CASE. SIMPLY CALCULATE STAGE BASED ON VOLUME.
                  VOL2 = RESID1*DELT
                  IF ( VOL2.LT.0.0 ) VOL2 = 0.0
                  STGNEW(LAKE) = STGTERP(VOL2,LAKE)
                  DSTG = ABS(STGNEW(LAKE) - STGITER(LAKE))
                  NCNCVR(LAKE) = 1
               END IF
               !     IF (lake==1 .and. kkper==1)then
               !     write(521,222)PRECIP(LAKE),EVAP(LAKE),RUNF(LAKE),RUNOFF(LAKE),
               !    1                WITHDRW(LAKE),SURFIN(LAKE),SURFOT(LAKE),
               !    2                SEEP(LAKE),VOLNEW1,VOLOLDD(LAKE),STGITER(LAKE),
               !    3                resid1,SURFA(LAKE),deriv,dstg,l1,kkiter,kkstp
               !     write(521,222)PRECIP3(LAKE),EVAP3(LAKE),RUNF(LAKE),RUNOFF(LAKE),
               !    1                WITHDRW3,SURFIN(LAKE),OUTFLOW,
               !    2                SEEP3(LAKE),VOLNEW2,VOLOLDD(LAKE),STGNEW(LAKE),
               !    3                resid2,SRFPT,deriv,dstg,l1,kkiter,kkstp
               !     END IF
               !222  format(15e20.10,3i5)
               IF(STGNEW(LAKE).LT.BOTTMS(LAKE))&
               &STGNEW(LAKE)=BOTTMS(LAKE)
               IF(DSTG.LE.SSCNCR) NCNCVR(LAKE) = 1
            END IF
   !
   !17-----COMPUTE NEW LAKE STAGE EXPLICITLY WITH THET1=0.
         ELSE
            IF(ISS.EQ.0) THEN
   !
   !17B----COMPUTE LAKE VOLUME FOR TRANSIENT SIMULATIONS.
   ! EDM - Add indices to RUNF and RUNOFF because these terms need to be
   !       saved for each lake in the simulation for writing to the FTL file by LMT
               VOL2 = DELT*(PRECIP(LAKE)-EVAP(LAKE)+RUNF(LAKE)+&
               &RUNOFF(LAKE)-WITHDRW(LAKE)+SURFIN(LAKE)-&
               &SURFOT(LAKE)+SEEP(LAKE))+ VOLOLDD(LAKE)
   !
   !17C----COMPUTE LAKE VOLUME FOR STEADY STATE SIMULATIONS.
            ELSE
               VOL2 = (PRECIP(LAKE)-EVAP(LAKE)+RUNF(LAKE)+&
               &RUNOFF(LAKE)-WITHDRW(LAKE)+SURFIN(LAKE)-&
               &SURFOT(LAKE)+SEEP(LAKE))
            END IF
   !
   !17D----NEW LAKE STAGE COMPUTED FROM LAKE VOLUME.
            STGNEW(LAKE) = STGTERP(VOL2,LAKE)
            IF(STGNEW(LAKE).LT.BOTTMS(LAKE))&
            &STGNEW(LAKE)=BOTTMS(LAKE)
         END IF
         VOL(LAKE) = VOLTERP(STGNEW(LAKE),LAKE)
      END DO
      IF ( IICNVG==1 ) EXIT CONVERGE
   END DO CONVERGE
   !
   !18-----PRINT NEW LAKE STAGE AND PREVIOUS LAKE ITERATION STAGE
   !        WHEN LAKE STAGE DOES NOT CONVERGE.
   DO LAKE=1,NLAKES
      IF ( ISS > 0 ) VOLOLDD(LAKE) = VOL(LAKE)   !RGN for modsim 5/15/18
      IF( L1.GE.MTER.AND.NCNCVR(LAKE).EQ.0 ) THEN
         WRITE(IOUT,1004) KKITER,&
         &LAKE,STGNEW(LAKE), STGITER(LAKE)
      END IF
   END DO
   !
   !19------FORMAT STATEMENTS
   !dep  101   FORMAT(4I5,3E20.10)    !format used for debugging
   !dep  202  FORMAT(i5,8(1X,E20.10)) !format used for debugging
506 FORMAT(1X,'ERROR - NO AQUIFER UNDER LAKE CELL ',4I5)
1004 FORMAT(1X,'ITERATION ',I4,2X,'LAKE ',I4,2X,'NEW STAGE ',1PE15.8,& !gsf
   &'  DID NOT CONVERGE-- PREVIOUS INTERNAL ITERATION STAGE  ',&
   &1PE15.8,/) !gsf
END SUBROUTINE GWF2LAK7FM
   !
SUBROUTINE GWF2LAK7BD(KSTP,KPER,IUNITGWT,IUNITGAGE,IUNITSFR,&
&IUNITUZF,NSOL,IGRID)
   !
   !------OLD USGS VERSION 7.1; JUNE 2006 GWF2LAK7BD;
   !------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
   !------NEW VERSION NUMBER 1.3.0, 7/01/2022
   !     ******************************************************************
   !     CALCULATE VOLUMETRIC BUDGET FOR LAKES
   !     ******************************************************************
   !
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFLAKMODULE
   USE GLOBAL,       ONLY: NCOL, NROW, NLAY, NODES, IBOUND, IOUT,&
   &ISSFLG, DELR, DELC, LBOTM, BOTM, HNEW,&
   &BUFF,IUNIT

   USE GWFBASMODULE, ONLY: MSUM, ICBCFL, DELT, PERTIM, TOTIM,&
   &HNOFLO, VBVL, VBNM
   USE GWFSFRMODULE, ONLY: STRIN
   IMPLICIT NONE
   !rsr: argument IUNITSFR not used
   CHARACTER*16 TEXT
   !dep  Added functions for interpolating between areas, derivatives,
   !dep     and outflow rates
   !     ------------------------------------------------------------------
   !     FUNCTIONS
   !     -----------------------------------------------------------------
   DOUBLE PRECISION FINTERP, DERIVTERP, OUTFLWTERP, VOLTERP, STGTERP,&
   &SURFTERP
   EXTERNAL FINTERP, DERIVTERP, OUTFLWTERP, VOLTERP, STGTERP,&
   &SURFTERP
   !     -----------------------------------------------------------------
   !     -----------------------------------------------------------------
   !     ARGUMENTS
   !     -----------------------------------------------------------------
   DOUBLE PRECISION BOTLK,BOTCL,CONDUC,H,FLOBOT,STGON,&
   &RATE,RATIN,RATOUT
   DOUBLE PRECISION THET1,SURFDPTH,CONDMX,BOTLKUP,BOTLKDN,VOL2
   DOUBLE PRECISION FLOTOUZF, SILLELEV, ADJSTAGE
   DOUBLE PRECISION CLOSEZERO, FLOBO2, FLOBO3, DLSTG
   DOUBLE PRECISION RUNFD, AREA, RAIN
   REAL zero, FACE, R, WDRAW, OLDSTAGE, AVHD, TOTARE, SUM
   REAL STGTST, SVT1, TVOLM, STO, FLSUM, TV, PPTIN, EOUT
   REAL SEEPUZF, QIN, QOUT, QSIN, QSOUT, DENOM
   INTEGER L1, IGRID, ISS, KPER, IBD, KCNT, LDR, NAUX, KSTP, IL
   INTEGER IR, IC, LK, ITRIB, INODE, LAKE, IUNITUZF, II, L
   INTEGER IUNITGWT, LL, K, J, I, JCLS, ICL, IC4, ICM4
   INTEGER ITYPE, INOFLO, IC5, ICM, IS1, IS2, ICl1, LK3, LK1
   INTEGER ICM2, IC2, IC3, ICNR, ICM3, ICNT, ICM1, L11, ICNR1
   INTEGER ILB, IRB, ICB, LDR1, NN, IIC, JIC, LIC, IUNITGAGE
   INTEGER NSOL, ILL, IL2, IUNITSFR, IL1
   DIMENSION JCLS(NCLS,ICMX)
   DIMENSION ILB(5),IRB(5),ICB(5)
   CHARACTER*16 LAKAUX(20)
   DIMENSION FACE(1)
   DATA TEXT /'   LAKE  SEEPAGE'/
   DATA LAKAUX(1)/'IFACE'/
   !      PARAMETER (CLOSEZERO=1.0E-7)
   !     ------------------------------------------------------------------
   !
   !------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2LAK7PNT(IGRID)
   ISS = ISSFLG(KPER)
   CLOSEZERO = 1.0D-09
   SURFDPTH = DBLE(SURFDEPTH)
   !
   !1------SET IBD=1 IF BUDGET TERMS SHOULD BE SAVED ON DISK.
   ZERO = 0.0
   IBD=0
   KCNT = 0
   RATIN = 0.
   RATOUT =0.
   LAKSEEP = 0.0
   DLSTG = 0.00001D0
   !1A-----Set Lake budget terms for GSFLOW to zero.
   TOTGWIN_LAK = 0.0
   TOTGWOT_LAK = 0.0
   TOTDELSTOR_LAK = 0.0
   TOTSTOR_LAK = 0.0
   TOTEVAP_LAK = 0.0
   TOTPPT_LAK = 0.0
   TOTRUNF_LAK = 0.0
   TOTWTHDRW_LAK = 0.0
   TOTSURFIN_LAK = 0.0
   TOTSURFOT_LAK = 0.0
   !dep  initialize CONDMX, BOTLKUP, AND BOTLKDN TO ZERO. 3/3/2009
   FLOBO2 = 0.0D0
   CONDMX = 0.0D0
   BOTLKUP = 0.0D0
   BOTLKDN = 0.0D0
   L1 = 0
   !dep  initialize SURFDPTH TO 1.0D0 3/3/2009
   SURFDPTH = DBLE(SURFDEPTH)
   !
   DO 104 LDR = 1,NODES
104 LDRY(LDR) = 0
   LDR = 0
   !
   !1B-----TEST TO SEE IF CELL-BY-CELL TERMS ARE NEEDED.
   IF(ILKCB.GT.0) IBD=ICBCFL
   !1C-----IF COMPACT BUDGET, WRITE LIST HEADER
   IF(IBD.EQ.2) THEN
   !-LFK         NAUX=0
   !-LFK         IF(IAUXSV.NE.0) NAUX=1
      NAUX=1
      CALL UBDSV4(KSTP,KPER,TEXT,NAUX,LAKAUX,ILKCB,NCOL,NROW,NLAY,&
      &LKNODE,IOUT,DELT,PERTIM,TOTIM,IBOUND)
   END IF
   !
   !1D------IF NO LAKE NODES, KEEP ZERO IN ACCUMULATORS.
   IF (LKNODE.EQ.0) GO TO 1200
   !
   !1E-----CLEAR CELL-BY-CELL BUFFER.
   DO 5 IL=1,NLAY
      DO 5 IR=1,NROW
         DO 5 IC=1,NCOL
5  BUFF(IC,IR,IL)=ZERO
   !
   !2------PROCESS EACH CELL IN THE ILAKE LIST.
   DO 100 LK=1,NLAKES
      FLXINL(LK)=ZERO
      SURFIN(LK)=ZERO
   !dep   Initialize flow limiting and 3 dummy arrays to zero
   !dep      for gage package
      GWRATELIM(LK) = ZERO
      XLAKES(lk,1) = ZERO
      XLKOLD(lk,1) = ZERO
      XLAKINIT(lk,1) = ZERO
100 CONTINUE
   !2A------SUM UP INFLOWS FROM INFLOWING STREAM REACHES.
   !dep   removed delta from computation of SURFIN
   ! - EDM Initialize variables needed for LMT
   NSFRLAK = 0
   ILKSEG = 0
   ILKRCH = 0
   LAKSFR = 0
   SWLAK = 0.
   IF (IUNITSFR.GT.0) THEN
      DO 200 LK=1,NLAKES
         DO 200 ITRIB=1,NTRB
            INODE=ITRB(LK,ITRIB)
            IF (INODE.LE.0) GO TO 200
            SURFIN(LK)=SURFIN(LK)+STRIN(INODE)
   !--EDM - SET UP VARIABLES NEEDED BY LMT
            IF (IUNIT(49).NE.0) THEN
               NSFRLAK = NSFRLAK + 1
               ILKSEG(NSFRLAK) = INODE
               LAKSFR(NSFRLAK) = LK
               SWLAK(NSFRLAK) = -STRIN(INODE)
            ENDIF
200   CONTINUE
   END IF
   !
   !2B------SET DOUBLE PRECISION THET1 TO THETA AND IF LESS THAN
   !          0.5, SET THET1 TO 0.0D0 (EXPLICIT LAKE STAGE).
   THET1 = DBLE(THETA)
   IF ( ISS==1 ) THET1 = 1.0
   IF(THET1-0.5.LT.-CLOSEZERO) THET1=0.0D0
   !2C------INITIALIZE SUMMATION PARAMETERS.
   DO LAKE = 1,NLAKES
      SUMCNN(LAKE) = ZERO
      SUMCHN(LAKE) = ZERO
      EVAP(LAKE)=ZERO
      PRECIP(LAKE)=ZERO
      SEEP(LAKE)=ZERO
      SEEPUZ(LAKE)=ZERO
   !            VOL(LAKE)=ZERO
      SURFA(LAKE)=ZERO
      GWIN(LAKE)=ZERO
      GWOUT(LAKE)=ZERO
      WITHDRW(LAKE) = WTHDRW(LAKE)
      IF(RNF(LAKE).GE.0.0) RUNF(LAKE) = RNF(LAKE)
      IF(RNF(LAKE).LT.0.0) RUNF(LAKE) =-RNF(LAKE)*PRCPLK(LAKE)*&
      &BGAREA(LAKE)
      IF (IUNITUZF.GT.0) THEN
         RUNOFF(LAKE) = OVRLNDRNF(LAKE)
      ELSE
         RUNOFF(LAKE) = 0.0
      END IF
      FLWIN(LAKE) = SURFIN(LAKE)+RUNF(LAKE)+RUNOFF(LAKE)+&
      &VOLTERP(STGOLD(LAKE),LAKE)/DELT
      IF ( ISS==1 ) THEN
         FLWIN(LAKE) = 1.0E10
      END IF
      FLWITER(LAKE) = FLWIN(LAKE)   !RGN 5/21/12
   END DO
   !
   !3------MASTER NODE LOOP -- COMPUTE LAKEBED SEEPAGE TERMS AND
   !         BUDGET TERMS.
   IF (ILKCB.LT.0.AND.ICBCFL.NE.0) WRITE (IOUT,'(//)')
   !3B-----II LOOP IS TO BALANCES INFLOWS AND OUTFLOWS TO/FROM A LAKE.
   !         WHEN II=1, FLOW INTO LAKE FROM ALL SOURCES ARE CALCULATED.
   !         WHEN II=2, SEEPAGE TO AND FROM LAKE AND RESIDUAL TERMS ARE
   !         ADDED TO GROUNDWATER MATRIX.LAKE SEEPAGE TO GROUNDWATER
   !         LIMITED TO AVAILABLE WATER IN LAKE.
   DO II = 1,2
      DO L=1,LKNODE
         IL=ILAKE(1,L)
         IR=ILAKE(2,L)
         IC=ILAKE(3,L)
   !
   !4------DETERMINE LAKE AND NODAL LAYER,ROW,COLUMN NUMBER.
         LAKE=ILAKE(4,L)
   !
   !5-------SET STGOLD TO STGNEW FOR STEADY STATE SIMULATIONS.
   !          COMPUTE STGON AS A FRACTION OF STGOLD AND STGNEW.
         IF(ISS.EQ.1 ) STGOLD(LAKE)=STGNEW(LAKE)
         STGON = (1.0D0-THET1)*STGOLD(LAKE) + THET1*STGNEW(LAKE)
         ITYPE = (ILAKE(5,L)+1)/2
         IF(ITYPE.EQ.3) ITYPE=0
         AREA = DELR(IC)*DELC(IR)
         IF(IL.GT.1) THEN
            BOTLK = BOTM(IC,IR,LBOTM(IL-1))
         ELSE
            BOTLK = BOTM(IC,IR,LBOTM(IL))
         END IF
         BOTCL = BOTM(IC,IR,LBOTM(IL))
         IF(IL.EQ.NLAY.AND.CNDFCT(L).EQ.0.0) BOTLK = BOTCL
         RAIN = PRCPLK(LAKE)
   !             EV = EVAPLK(LAKE)
   !
   !5B------CONDUCTANCE FACTOR NEEDED FOR SEEPAGE CALCULATIONS.
   !         FLOB01 USED TO CALCULATE SEEPAGE WITH STGOLD.
   !         FLOBO2 USED TO CALCULATE SEEPAGE WITH STGNEW.
   !         FLOBOT IS A FRACTION OF BOTH FLOBO1 AND FLOBO2 AND
   !           IS DEPENDENT ON VALUE OF THET1.
         CONDUC=CNDFCT(L)
         FLOBOT = 0.0D0
         FLOBO3 = 0.0D0
         FLOTOUZF = 0.0D0
         RATE = 0.0
         IL1 = IL
         IF ( ITYPE.EQ.0 ) THEN       !RGN 2-18-2014 added this if check
            DO WHILE (IL1 .LE. NLAY)
               IF( IBOUND(IC,IR,IL1).GT.0 ) THEN
                  EXIT
               ELSE
                  IL1 = IL1 + 1
               END IF
            END DO
            IF ( IL1.GT.NLAY ) IL1 = NLAY
         END IF
         IF( IBOUND(IC,IR,IL1).LE.0 ) THEN
            !  Commented next line out 12/27/10
            !              IF ( CONDUC.GT.CLOSEZERO )WRITE(IOUT,506) L,IC,IR,IL
            CONDUC = 0.0
            ! 506          FORMAT(1X,'ERROR - NO AQUIFER UNDER LAKE CELL ',4I5)
         END IF
         IF(CONDUC.GT.0.0) THEN
            H=HNEW(IC,IR,IL1)  !RGN set IL to IL1 1/9/16
   !
   !5C------DETERMINE UPPERMOST ACTIVE CELL IF NOT CELL(IL)
   !
            CALL GET_FLOBOT(IC, IR, IL1, ITYPE, INOFLO,CONDUC,&
            &FLOBOT,FLOBO3,FLOTOUZF,DLSTG,CLOSEZERO,H,&
            &THET1,ISS,LAKE,II,SURFDPTH,AREA,IUNITUZF,&
            &BOTLK,BOTCL,L1)
   !
   !9-------SET RATE TO FLOBOT AND SET FLOB(L) TO FLOBOT WHEN
   !          SOLUTE TRANSPORT IS ACTIVE.
            IF ( II==2 ) THEN
               SEEP(LAKE)=SEEP(LAKE)-FLOBOT
   ! Save seepage to UZF for writing Lake-to-UZF in GAG Package.
               SEEPUZ(LAKE)=SEEPUZ(LAKE)+FLOTOUZF
               RATE=FLOBOT
               IF (IUNITGWT.GT.0.OR.IUNIT(49).GT.0) FLOB(L)=FLOBOT !EDM
               IF (ILKCB.LT.0.AND.ICBCFL.NE.0) WRITE(IOUT,880)&
               &TEXT,KPER,KSTP,L,IL,IR,IC,RATE
880            FORMAT(1X,A,'   PERIOD',I3,'   STEP',I6,'   NODE',I4,& !gsf
               &'   LAYER',I3,'   ROW',I4,'   COL',I4,'   RATE',&
               &G15.7)
   !
   !10------ADD RATE TO BUFFER.
               BUFF(IC,IR,IL1)=BUFF(IC,IR,IL1)+RATE
               LAKSEEP(IC,IR) = LAKSEEP(IC,IR)+ RATE + FLOTOUZF  !10/4/2014
   !
   !10B-----CHECK IF RATE IS DISCHARGING FROM AQUIFER (NEGATIVE RATE).
   !dep            IF (RATE) 885,899,890
               IF(RATE.LT.0.0D0) THEN
   !
   !10C-----SUBTRACT RATE FROM RATOUT.
   !dep 885         RATOUT=RATOUT-RATE
                  RATOUT=RATOUT-RATE
                  GWIN(LAKE)=GWIN(LAKE)-RATE
   !dep         GO TO 899
   !
   !10D------CHECK IF RATE IS RECHARGING AQUIFER (POSITIVE RATE).
               ELSE IF (RATE+FLOTOUZF.GT.0.0D0) THEN  !11/3/2014 added flotouzf
   !
   !10E------ADD RATE TO RATIN.
   !dep 890         RATIN=RATIN+RATE
                  RATIN=RATIN+RATE
                  GWOUT(LAKE)=GWOUT(LAKE)+RATE
               END IF
   !11-------IF SAVING COMPACT BUDGET, WRITE FLOW FOR ONE LAKE FACE.
            END IF
         END IF
899      IF(IBD.EQ.2.and.II.EQ.2) THEN
            FACE(1)=ILAKE(5,L)
            R=RATE
            CALL UBDSVB(ILKCB,NCOL,NROW,IC,IR,IL1,R,FACE(1),1,&
            &NAUX,1,IBOUND,NLAY)
         END IF
      END DO
   END DO
   !
   !12------COMPUTE EVAPORATION AND PRECIPITATION USING STGOLD AND
   !          ADD PRECIPITATION TO FLWIN.
   DO LL = 1,NLAKES
      SURFA(LL)=FINTERP(STGNEW(LL),LL)
      IF ( IGSFLOWLAK == 0 ) THEN
         EVAP(LL)=EVAPLK(LL)*SURFA(LL)
         PRECIP(LL)=PRCPLK(LAKE)*SURFA(LL)
      ELSE
         EVAP(LL) = EVAPLK(LL)
         PRECIP(LL) = PRCPLK(LL)
      END IF
      FLWIN(LL) = FLWIN(LL) + PRECIP(LL)
   !
   !13------LIMIT WITHDRW TO LAKE INFLOW WHEN WITHDRAWALS EXCEED
   !           INFLOW (INCLUDING AVAILABLE LAKE STORAGE).
      IF( WITHDRW(LL)-DBLE(FLWIN(LL)).GT.1.0E-04 ) THEN
         WITHDRW(LL) = FLWIN(LL)
         FLWIN(LL) = 0.0
      ELSE
         FLWIN(LL)  = FLWIN(LL) - WITHDRW(LL)
      END IF
   !
   !14------LIMIT EVAPORATION TO LAKE INFLOW WHEN EVAPORATION EXCEEDS
   !          INFLOW (INCLUDING AVAILABLE LAKE STROAGE AND WITHDRAWALS).
      IF ( EVAP(LL)-DBLE(FLWIN(LL)).GT.1.0E-04 ) THEN
         EVAP(LL)=FLWIN(LL)
         FLWIN(LL) = 0.0D0
      ELSE
         FLWIN(LL)=FLWIN(LL)-EVAP(LL)
      END IF
   END DO
   !dep  August 28, 2009  moved to end of do loop
   !15-----COMPUTE CHANGE IN STAGE DURING TIME STEP AND FOR SIMULATION.
   !       SKIP IF STEADY STATE SIMULATION.
   !      IF(ISS.LE.0) GO TO 905
   !      DO 903 LAKE=1,NLAKES
   !            DELH(LAKE)=STGNEW(LAKE)-STGOLD(LAKE)
   !            TDELH(LAKE)=STGNEW(LAKE)-STAGES(LAKE)
   !  903 CONTINUE
   ! RGN commented next line out. 5/4/09.
   !      GO TO 1350
   !
   !16-----COMPUTE STGON AS FRACTION OF STGOLD AND STGNEW.
905 DO 1000 LAKE=1,NLAKES
      VOL2 = 0.0
   !dep       STGON = (1.0-THETA)*STGOLD(LAKE) + THETA*STGNEW(LAKE)
   !dep       Changed THETA to THET1
      STGON = (1.0D0-THET1)*STGOLD(LAKE) + THET1*STGNEW(LAKE)
   !
   !17-----COMPUTE RUNOFF INTO LAKE FROM LAKE PACKAGE AND FROM UZF.
   !dep   Changed WTHDRW(LAKE) TO WITHDRW(LAKE)
      WDRAW=WITHDRW(LAKE)
      IF(RNF(LAKE).GE.ZERO) RUNF(LAKE) = RNF(LAKE)
      IF(RNF(LAKE).LT.ZERO) RUNF(LAKE) =-RNF(LAKE)*PRCPLK(LAKE)&
      &*BGAREA(LAKE)
   !dep  Added runoff from Unsaturated Flow Package
      IF (IUNITUZF.GT.0) THEN
         RUNOFF(LAKE) = OVRLNDRNF(LAKE)
      ELSE
         RUNOFF(LAKE) = 0.0
      END IF
   !dep  Created RUNFD and added to STGNEW
   !-LFK         RUNFD = RUNF(LAKE)+ RUNOFF(LAKE)
      RUNFD = RUNF(LAKE) + RUNOFF(LAKE)
   !
   !18------COMPUTE LAKE VOLUME FROM ALL INFLOWS AND OUTFLOWS FOR
   !          TRANSIENT SIMULATION AND THEN COMPUTE STGNEW FROM
   !          NEW VOLUME.
   !RGN Volume made equal to sum of inflows and outflows plus  ! 6/30/2017 commenting this out as it causes issues for small lakes
   !RGN   plus lake storage from previous time step  4/17/09
      IF(ISS.EQ.0)THEN
         VOL2 = VOLOLDD(LAKE)+(PRECIP(LAKE)-EVAP(LAKE)&
         &-WDRAW+RUNFD+SURFIN(LAKE)-SURFOT(LAKE)+GWIN(LAKE)&    !10/4/2014 added SEEPUZ(LAKE)
         &-GWOUT(LAKE)-SEEPUZ(LAKE))*DELT
         IF(VOL2.LE.0.0) VOL2=0.0
   !          VOL(LAKE) = VOL2
   !          STGNEW(LAKE) = STGTERP(VOL2,LAKE)
   !C
   !C18B-----COMPUTE LAKE VOLUME FROM ALL INFLOWS AND OUTFLOWS FOR
   !C          STEADY STATE SIMULATION.
      ELSE
         VOL2 = VOLTERP(STGNEW(LAKE),LAKE)
         IF(VOL2.LE.0.0D0) VOL2 = 0.0D0
         VOL(LAKE) = VOL2
      END IF
   !18C-----STGON IS FRACTION OF STGOLD AND STGNEW AND SURFACE AREA
   !
   !          IS BASED ON STGOLD.
      STGON = (1.0D0-THET1)*STGOLD(LAKE) + THET1*STGNEW(LAKE)
      SURFA(LAKE)=FINTERP(STGNEW(LAKE),LAKE)
      IF(STGNEW(LAKE)-BOTTMS(LAKE).LT.CLOSEZERO) GO TO 1110
   !
   !19------COMPUTE LAKE BUDGET VOLUMES FOR GSFLOW CSV FILE.
   !dep  EVAP, PRECIP,WDRW AND RUNFD are volumetric rates 4/19/2009
      TOTGWIN_LAK = TOTGWIN_LAK + GWIN(LAKE)*DELT
      TOTGWOT_LAK = TOTGWOT_LAK - GWOUT(LAKE)*DELT
      TOTDELSTOR_LAK = TOTDELSTOR_LAK + VOL(LAKE) - VOLOLDD(LAKE)
      TOTSTOR_LAK = TOTSTOR_LAK + VOL(LAKE)
      TOTEVAP_LAK = TOTEVAP_LAK - EVAP(LAKE)*DELT
      TOTPPT_LAK = TOTPPT_LAK + PRECIP(LAKE)*DELT
      TOTRUNF_LAK = TOTRUNF_LAK + RUNFD*DELT
      TOTWTHDRW_LAK = TOTWTHDRW_LAK - WDRAW*DELT
      TOTSURFIN_LAK = TOTSURFIN_LAK + SURFIN(LAKE)*DELT
      TOTSURFOT_LAK = TOTSURFOT_LAK - SURFOT(LAKE)*DELT
   !
   !20------WRITE WHEN LAKE VOLUME HAS INITIALLY GONE DRY.
      IF(VOL(LAKE).LE.ZERO) WRITE(IOUT,1114) LAKE
1114  FORMAT(1X,'..........LAKE',I3,' JUST GONE DRY..........')
   !dep  August 27, 2009   set delh and tdelh to zero for steady state
      IF (ISS.EQ.1) THEN
         IF(KPER.EQ.1) THEN
            STAGES(LAKE)=STGNEW(LAKE)
         END IF
         DELH(LAKE)= 0.0
         TDELH(LAKE)= 0.0
   !dep  August 27, 2009  set delh and tdelh for transient simulations
      ELSE
         OLDSTAGE= STGOLD(LAKE)
         DELH(LAKE)=STGNEW(LAKE)-OLDSTAGE
         TDELH(LAKE)=STGNEW(LAKE)-STAGES(LAKE)
      END IF
      GO TO 1000
   !
   !21----WRITE WHEN LAKE CONTINUES TO BE DRY.
1110  AVHD = ZERO
   !      BOTARE = ZERO
      WRITE(IOUT,1112) LAKE
1112  FORMAT(1X,'..........LAKE',I3,' IS DRY..........')
      IF(NLAY.EQ.1) GO TO 1000
      DO 1115 L=1,LKNODE
         L1 = ILAKE(4,L)
   !  Convert ILAKE(5,L):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
         ITYPE = (ILAKE(5,L)+1)/2
         IF(ITYPE.EQ.3) ITYPE=0
         IF(L1.NE.LAKE.OR.ITYPE.NE.0) GO TO 1115
         K = ILAKE(1,L)
         J = ILAKE(2,L)
         I = ILAKE(3,L)
         IF(K.EQ.NLAY.AND.IBOUND(I,J,K).EQ.0) GO TO 1000
         IF(K.EQ.1) GO TO 1115
         IF(BOTM(I,J,LBOTM(K-1)).GT.BOTTMS(LAKE)) GO TO 1115
   !dep   lake no longer refills on basis of average head beneath aquifer.
   !      CHECK FOR CONDITION (AVERAGE HEAD IN UNDERLYING AQUIFER
   !      GREATER THAN BOTTOM OF LAKE) FOR REWETTING A DRY LAK
   ! lfk  eliminate DO 1117 loop, ignore LAYHDT value below lake; and only
   !           check layer below lake interface for IBOUND value:
   ! orig lines
   !      DO 1117 K1=K,NLAY
   !      IF(LAYHDT(K1).EQ.0) GO TO 1117
   !      K2 = K1
   !      IF(IBOUND(I,J,K1).LE.0) GO TO 1117
   !      K2 = K
   !      IF(IBOUND(I,J,K).LE.0) GO TO 1117
   !      GO TO 1119
   ! 1117 CONTINUE
   !      AVHD = AVHD + BOTM(I,J,LBOTM(K2))*DELR(I)*DELC(J)
   !      GO TO 1121
   ! 1119 AVHD = AVHD + HNEW(I,J,K1)*DELR(I)*DELC(J)
   !
   ! new lines
   !dep      IF(IBOUND(I,J,K).LE.0) THEN
   !dep        AVHD = AVHD + BOTM(I,J,LBOTM(K))*DELR(I)*DELC(J)
   !dep      ELSE
   !dep        AVHD = AVHD + HNEW(I,J,K)*DELR(I)*DELC(J)
   !dep     END IF
   !
   !dep  1121 BOTARE = BOTARE + DELR(I)*DELC(J)
   !dep     BOTARE = BOTARE + DELR(I)*DELC(J)
1115  CONTINUE
   !dep      IF(BOTARE.LE.ZERO) GO TO 1000
   !dep      AVHD = AVHD/BOTARE
   !dep     IF(AVHD.LT.BOTTMS(LAKE)) GO TO 1125
   !dep      WRITE(IOUT,1122)
   !dep 1122 FORMAT(/1X,'AQUIFER HEAD UNDERNEATH BOTTOM OF LAKE IS HIGHER THAN
   !dep     1LAKE BOTTOM ELEVATION')
   !dep  Revised to not reset lake to average ground-water level beneath lake.
   !dep 1120 STGNEW(LAKE) = BOTTMS(LAKE)
   !dep      SURFA(LAKE) = BOTARE
   !dep      VOL(LAKE) = (STGNEW(LAKE)-BOTTMS(LAKE))*SURFA(LAKE)
   !dep  revised program to allow lakes to rewet while keeping track of lake
   !dep   inflows and outflow
   !dep      WRITE(IOUT,1184) LAKE, STGNEW(LAKE)
   !dep 1184 FORMAT(/1X,'LAKE',I3,' HAS REWET.  SET STAGE OF LAKE TO',F10.2,
   !dep     1  2X,'FT'/)
   !dep     GO TO 1000
   !
   !-----CHECK FOR STREAM OR AUGMENTATION INFLOWS
   !
   !dep 1125 TSTVOL = SURFIN(LAKE)
   !dep     IF(WDRAW.LT.ZERO) TSTVOL = TSTVOL - WDRAW
   !dep     IF(TSTVOL.LE.ZERO) GO TO 1000
   !dep     STGNEW(LAKE) = BOTTMS(LAKE)
   !dep    SURFA(LAKE) = BOTARE
   !dep      AVHD = BOTTMS(LAKE) + TSTVOL/BOTARE
   !dep  revised program to allow lakes to rewet while keeping track of lake
   !dep   inflows and outflow
   !dep   WRITE(IOUT,1123)
   !dep 1123 FORMAT(/1X,'THERE ARE NON-ZERO SURFACE-WATER INFLUXES INTO THE DRY
   !dep     1 LAKE')
   !dep      GO TO 1120
   !
1000 CONTINUE
   !
   !22------ADJUST STAGES OF COALESCENT MULTIPLE-LAKE SYSTEMS.
   !
   KCNT = 0
   IF(NCLS.LE.0) GO TO 1350
   DO 1205 I=1,NCLS
      DO 1205 J=1,ICMX
1205 JCLS(I,J) = 0
   !
   !22B-----CHECK EACH LAKE SYSTEM (ICL) FOR CURRENT CONNECTIONS TO
   !          SUBLAKES.
   DO 1300 ICL=1,NCLS
      DO 1206 K=1,NLAKES
         SVT(K) = ZERO
         NCNST(K) = 0
1206  NCNT(K) = 0
   !
   !22C-----ELIMINATE (JCLS=2) ALL LAKES THAT HAVE ALREADY HAD THEIR STAGES
   !          ADJUSTED AS PART OF A CONNECTED SYSTEM OF LAKES AND SUBLAKES.
      DO 1210 IC4=ICL,NCLS
         ICM4 = ICS(IC4)
         DO 1210 IC5=1,ICM4
            IF(JCLS(IC4,IC5).EQ.1) JCLS(IC4,IC5) = 2
1210  CONTINUE
   !dep 1215 IF(JCLS(ICL,1).GE.2) GO TO 1300
      IF(JCLS(ICL,1).GE.2) GO TO 1300
   !
   !22D-----TAG CENTER LAKE BY SETTING JCLS=1 AND THEN CHECK SUBLAKES FOR
   !   CONNECTIONS.  IF CONNECTED, SET JCLS=1 AND NCNT=1.
      ICM = ICS(ICL)
      IS1 = ISUB(ICL,1)
      JCLS(ICL,1) = 1
      NCNT(IS1) = 1
      SVT(IS1) = -99999.
      DO 1220 J=2,ICM
         IS2 = ISUB(ICL,J)
         IF(IS2.LE.0) GO TO 1225
         IF(STGNEW(IS1).LE.SILLVT(ICL,J-1).AND.STGNEW(IS2).LE.&
         &SILLVT(ICL,J-1)) GO TO 1220
         JCLS(ICL,J) = 1
         NCNT(IS2) = 1
         SVT(IS2) = SILLVT(ICL,J-1)
1220  CONTINUE
1225  IF(ICL.EQ.NCLS) GO TO 1240
   !
   !22E-----CHECK TO SEE IF CENTER LAKES OF REMAINING LAKE SYSTEMS
   !          ARE THE SAMEAS CONNECTED SUBLAKES OF THE PRESENT LAKE SYSTEM.
   !          IF SO, CHECK THEIR SUBLAKES FOR CONNECTIONS TO THE CENTER
   !          LAKES. IF SUBLAKES ARE ADDED TO THE SYSTEM, THEN CHECK
   !          REMAINING CENTER LAKES FOR AN IDENTITY WITH THE NEWLY
   !          ADDED SUBLAKES. ALL CONNECTED LAKES ARE APPROPRIATELY
   !          TAGGED (JCLS=1 AND NCNT=1).
      ICL1 = ICL + 1
      DO 1230 LK=ICL1,NCLS
         IF(JCLS(LK,1).EQ.2) GO TO 1230
         LK3 = LK - 1
         DO 1227 LK1=1,LK3
            ICM2 = ICS(LK1)
            DO 1227 IC2=2,ICM2
               IF(ISUB(LK1,IC2).NE.ISUB(LK,1).OR.JCLS(LK1,IC2).NE.1) GO TO 1227
               JCLS(LK,1) = 1
               IS1 = ISUB(LK,1)
               NCNT(IS1) = 1
               ICM3 = ICS(LK)
               DO 1226 IC3=2,ICM3
                  IS2 = ISUB(LK,IC3)
                  IF(IS2.LE.0) GO TO 1227
                  IF(STGNEW(IS1).LE.SILLVT(LK,IC3-1).AND.STGNEW(IS2).LE.&
                  &SILLVT(LK,IC3-1)) GO TO 1226
                  JCLS(LK,IC3) = 1
                  NCNT(IS2) = 1
                  SVT(IS2) = SILLVT(LK,IC3-1)
1226           CONTINUE
1227     CONTINUE
1230  CONTINUE
   !
   !22F-----COUNT NUMBER OF LAKES IDENTIFIED AS A CONNECTED PART
   !          OF THE PRESENT LAKE SYSTEM, STORE LAKE NUMBERS IN KSUB,
   !          AND CUMULATE TOTAL SURFACE AREA.
1240  ICNT = 0
      KCNT = KCNT + 1
      TOTARE = ZERO
      DO 1245 L=1,NLAKES
         IF(NCNT(L).NE.1) GO TO 1245
         ICNT = ICNT + 1
         TOTARE = TOTARE + SURFA(L)
         KSUB(ICNT) = L
         MSUB(ICNT,KCNT) = L
1245  CONTINUE
      MSUB1(KCNT) = ICNT
      IF(ICNT.LE.1) KCNT = KCNT - 1
      IF(ICNT.LE.1) GO TO 1300
      IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 1251
      WRITE(IOUT,1250) KSTP, ICNT, TOTARE
1250  FORMAT(/1X,80('-')/1X,'TIME STEP',I6,5X,'NUMBER OF CONNECTED LAKE !gsf&
      &S IS',I3,5X,'TOTAL AREA = ',D16.9/)
1251  CONTINUE
   !
   !22G-----COMPUTE STAGE ADJUSTMENTS (STGADJ) REQUIRED FOR CONNECTED
   !           LAKES TO HAVE THE SAME STAGE.
   !dep     changed adjustment on basis of volumes exchanged among lakes.
      DO 1270 I=1,ICNT
         L1 = KSUB(I)
         SUM = ZERO
         DO 1265 J=1,ICNT
            IF(J.EQ.I) GO TO 1265
            L = KSUB(J)
            SUM = SUM + VOLTERP(STGNEW(L),L)-VOLTERP(STGNEW(L1),L)
   !     SUM = SUM + SURFA(L)*(STGNEW(L)-STGNEW(L1))
1265     CONTINUE
         STGADJ(L1) = SUM/TOTARE
1270  CONTINUE
   !
   !22H-----CHECK FOR NEWLY COMPUTED LAKE STAGES (STGNEW) LESS
   !          THAN SILL ELEVATIONS OF THE LAKES (SVT).
1272  ICNR = 0
      ICM3 = ICS(ICL)
      DO 1275 IC3=2,ICM3
         L = ISUB(ICL,IC3)
         IF(L.LE.0) GO TO 1275
         IF(NCNST(L).EQ.1) GO TO 1275
         IF(NCNT(L).EQ.0) GO TO 1275
         IF(STGNEW(L).LT.SVT(L)) GO TO 1275
         STGTST = STGNEW(L) + STGADJ(L)
         IF(STGTST.GE.SVT(L)) GO TO 1275
   !
   !22I-----ADJUST STAGE TO SILL ELEVATION.
   !
         NCNST(L) = 1
   !dep  revised calculation of FLXINL using volumes 6/7/2009
   !dep  created a double precision local variable named SILLELEV
         SILLELEV  = SVT(L)
         FLXINL(L) = VOLTERP(SILLELEV,L)-VOLTERP(STGNEW(L),L)
         STGADJ(L) = SVT(L) - STGNEW(L)
         STGNEW(L) = SVT(L)
   !dep   commented out calculation of FLXINL using surface areas
   !      FLXINL(L) = STGADJ(L)*SURFA(L)
         VOL(L) = VOL(L) + FLXINL(L)
         TOTARE = TOTARE - SURFA(L)
         NCNT(L) = 2
         ICNR = 1
         WRITE(IOUT,2238) L,L,STGADJ(L),STGNEW(L)
1275  CONTINUE
      IF(ICL.EQ.NCLS) GO TO 1277
   !
   !22J-----IF A LAKE STAGE IS ADJUSTED TO THE SILL ELEVATION, CHECK TO SEE
   !          WHETHER THERE ARE SUBLAKES OF THIS LAKE AND ADJUST THEM TO
   !          THE SILL ELEVATION UNLESS THE ORIGINAL STAGE IS ALREADY LOWER,
   !          IN WHICH CASE THEY ARE NO LONGER CONNECTED.
      ICL1 = ICL + 1
      DO 2230 LK=ICL1,NCLS
         IS1 = ISUB(LK,1)
         IF(NCNT(IS1).EQ.0) GO TO 2230
         ICM1 = ICS(LK)
         DO 2225 IC2=2,ICM1
            IS2 = ISUB(LK,IC2)
            IF(NCNST(IS2).EQ.1) GO TO 2225
            IF(NCNT(IS2).EQ.0) GO TO 2225
            IF(STGNEW(IS2).LT.SVT(IS2)) GO TO 2225
            SVT1 = SVT(IS2)
            IF(SVT(IS1).GT.SVT1.AND.NCNT(IS1).EQ.2) SVT1 = SVT(IS1)
            STGTST = STGNEW(IS2) + STGADJ(IS2)
            IF(STGTST.GE.SVT1) GO TO 2225
            ICNR = 1
            NCNST(IS2) = 1
            NCNT(IS2) = 2
            STGADJ(IS2) = SVT1 - STGNEW(IS2)
   !dep  revised calculation of FLXINL using volumes 6/7/2009
   !dep  created a double precision local variable named SILLELEV
            SILLELEV = SVT1
            FLXINL(IS2)=VOLTERP(SILLELEV,IS2)-VOLTERP(STGNEW(IS2),IS2)
            STGNEW(IS2) = SVT1
   !dep commented calculation of FLXINL using surface area
   !      FLXINL(IS2) = STGADJ(IS2)*SURFA(IS2)
            VOL(IS2) = VOL(IS2) + FLXINL(IS2)
            TOTARE = TOTARE - SURFA(IS2)
            IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 2225
            L11 = L
            IF(SVT(IS2).GT.SVT(L)) L11 = IS2
            WRITE(IOUT,2238) IS2,L11,STGADJ(IS2),STGNEW(IS2)
2238        FORMAT(1X,'READJUST STAGE OF LAKE ',I3,' TO LAKE ',I3,&
            &' SILL ELEVATION BY ',F5.2,' TO ',F7.2)
2225     CONTINUE
2230  CONTINUE
1277  IF(ICNR.LE.0) GO TO 1280
   !
   !22K-----RECOMPUTE STAGE ADJUSTMENTS CONSTRAINED NOT TO LOWER LAKES
   !          BELOW SILL ELEVATIONS.
   !
      ICNR1 = 0
      DO 1370 I=1,ICNT
         L1 = KSUB(I)
         IF(NCNT(L1).EQ.0) GO TO 1370
         IF(NCNST(L1).EQ.1) GO TO 1370
         SUM = ZERO
         DO 1365 J=1,ICNT
            IF(J.EQ.I) GO TO 1365
            L = KSUB(J)
            IF(NCNT(L).EQ.0) GO TO 1365
   !dep changed computation of volume change 6/09/2009
            IF(NCNST(L).EQ.0) SUM = SUM +&
            &VOLTERP(STGNEW(L),L)-VOLTERP(STGNEW(L1),L)
   !      IF(NCNST(L).EQ.0) SUM = SUM + SURFA(L)*(STGNEW(L)-STGNEW(L1))
            IF(NCNST(L).EQ.1) SUM = SUM - SURFA(L)*STGADJ(L)
1365     CONTINUE
         STGADJ(L1) = SUM/TOTARE
         STGTST = STGNEW(L1) + STGADJ(L1)
         IF(STGTST.GE.SVT(L1)) GO TO 1370
         ICNR1 = 1
1370  CONTINUE
      IF(ICNR1.NE.0) GO TO 1272
1280  IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 1281
      WRITE(IOUT,1286)
1286  FORMAT(//11X,'SURFACE',7X,'SILL',3X,'WATER BUDGET',2X,'STAGE',2X,&
      &'CORRECTED',3X,'LAKE VOLUME'/2X,'LAKE',7X,'AREA',6X,'ELEVATION',&
      &3X,'STAGE',3X,'CORRECTION',2X,'STAGE',5X,'CORRECTION')
1281  CONTINUE
      TVOLM = ZERO
      DO 1290 I=1,ICNT
         L = KSUB(I)
         STO = STGNEW(L)
         IF(NCNST(L).EQ.1) STO = STGNEW(L) - STGADJ(L)
         IF(NCNST(L).EQ.1) GO TO 1285
   !dep  revised calculation of FLXINL using volumes 6/7/2009
   !dep  created a double precision local variable named ADJSTAGE
         ADJSTAGE = STGNEW(L)+STGADJ(L)
         FLXINL(L) = VOLTERP(ADJSTAGE,L)-VOLTERP(STGNEW(L),L)
         STGNEW(L) = STGNEW(L) + STGADJ(L)
   !dep commented calculation of FLXINL using surface area
   !      FLXINL(L) = STGADJ(L)*SURFA(L)
         VOL(L) = VOL(L) + FLXINL(L)
1327     FORMAT(/10X,'WARNING -- SUM OF INTERLAKE FLUXES ',F10.0,' EXCEEDS&
         &10**6 OF THE TOTAL VOLUME'/)
         WRITE(IOUT,1301)
1301     FORMAT(1X,80('-')/)
1285     TVOLM = TVOLM + VOL(L)
         IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 1290
         WRITE(IOUT,1269) L,SURFA(L),SVT(L),STO,STGADJ(L),STGNEW(L),&
         &FLXINL(L)
1269     FORMAT(1X,I3,1X,G15.5,4F10.2,G15.5)
1290  CONTINUE
   !
   !22L-----RECOMPUTE TIME STEP AND CUMULATIVE STAGE CHANGES FOR
   !          CONNECTED LAKES.
      DO 1295 I=1,ICNT
         L = KSUB(I)
   !dep August 27, 2009 Transient only
         IF (ISS.NE.1) THEN
            OLDSTAGE = STGOLD(L)
            DELH(L) = STGNEW(L) - OLDSTAGE
            TDELH(L) = STGNEW(L) - STAGES(L)
         END IF
1295  CONTINUE
1300 CONTINUE
   !
   !22M-----CHECK ON SUM OF CONNECTED-LAKE INTERCHANGE VOLUMES.
   FLSUM = ZERO
   DO 1325 L=1,NLAKES
      FLSUM = FLSUM + FLXINL(L)
1325 CONTINUE
   IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 1350
   TV = TVOLM/1000000.
   IF(FLSUM.GE.TV) WRITE(IOUT,1327) FLSUM
1350 CONTINUE
   !
   !23------CHECK FOR LAKE CELLS GOING DRY.
   LDR = 0
   DO 950 L=1,LKNODE
      IL=ILAKE(1,L)
      IR=ILAKE(2,L)
      IC=ILAKE(3,L)
      LAKE=ILAKE(4,L)
   !  Convert ILAKE(5,L):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      ITYPE = (ILAKE(5,L)+1)/2
      IF(ITYPE.EQ.3) ITYPE=0
      IF(ITYPE.NE.0) GO TO 950
      IF(IBOUND(IC,IR,IL).GT.0) BOTLK = BOTM(IC,IR,LBOTM(IL-1))
      IF(IBOUND(IC,IR,IL).EQ.0) BOTLK = BOTM(IC,IR,LBOTM(IL))
   !dep  revised to set STGNEW to BOTTOM of LAKE when LESS than BOTLK.
      IF (STGNEW(LAKE).LE.BOTLK) THEN
         LDR = LDR + 1
         LDRY(LDR) = L
   !dep           STGNEW(LAKE)=BOTLK
      END IF
950 CONTINUE
   IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 951
   IF(LDR.EQ.0) WRITE(IOUT,875) LDR
875 FORMAT(//1X,I5,2X,'LAKE CELLS ARE DRY.')
   IF(LDR.EQ.0) GO TO 951
   WRITE(IOUT,874) LDR
874 FORMAT(/5X,'SECTIONS OF THE LAKE BOTTOM HAVE BECOME DRY.  THE DRY&
   &SECTIONS'/5X,'LIE ABOVE THE FOLLOWING',I3,' AQUIFER CELLS (LAYER,R&
   &OW,COLUMN):')
   LDR1 = 0
   DO 952 L=1,LDR
      LDR1 = LDR1 + 1
      L1 = LDRY(L)
      ILB(LDR1) = ILAKE(1,L1)
      IRB(LDR1) = ILAKE(2,L1)
      ICB(LDR1) = ILAKE(3,L1)
      IF(LDR1.LT.5) GO TO 952
      WRITE(IOUT,876) (ILB(I),IRB(I),ICB(I),I=1,5)
876   FORMAT(5X,5('(',I3,',',I3,',',I3,')',2X))
      LDR1 = 0
952 CONTINUE
   IF(LDR1.GT.0) WRITE(IOUT,876) (ILB(I),IRB(I),ICB(I),I=1,LDR1)
951 CONTINUE
   !dep  Added following do loop.
   !24------Set Lake Stage to bottom of lake when lake is dry and set
   !       lake volume to zero.
   DO LAKE=1,NLAKES
      IF(STGNEW(LAKE).LE.BOTTMS(LAKE)) THEN
         STGNEW(LAKE) = BOTTMS(LAKE)
         VOL(LAKE) = 0.0
      END IF
   END DO
   !dep Moved call to GAG7LO to just above comment line 13 6/9/2009
   !dep      IF(IUNITGWT.GT.0) GO TO 1086
   !dep      NDUM=1
   !dep  Buff was passed to GAGE package as a dummy array because
   !dep     transport is inactive and the variables were not defined.
   !dep      IF(IUNITGWT.LE.0 .AND. IUNITGAGE.GT.0)
   !dep     *  CALL SGWF2GAG7LO(IUNITGWT,IUNITUZF,XLAKES,TOTIM,
   !dep     *  GWIN,GWOUT,SEEP,FLXINL,VOLOLD,XLKOLD,XLAKINIT,NSOL)
   !
   !25------WRITE WARNINGS IF LAKE STAGES EXCEED SPECIFIED MINIMUMS
   !          AND MAXIMUMS FOR STEADY STATE SIMULATIONS.
1086 DO LAKE=1,NLAKES
      IF(ISS.GT.0) THEN
         IF (STGNEW(LAKE).LT.SSMN(LAKE)) THEN
            WRITE(IOUT,972) STGNEW(LAKE), SSMN(LAKE), LAKE
972         FORMAT(/1X,'WARNING-- COMPUTED STAGE OF ',F10.2,&
            &' IS LESS THAN SPECIFIED MINIMUM ',F10.2,&
            &' FOR LAKE ',I5)
         ELSE IF (STGNEW(LAKE).GT.SSMX(LAKE)) THEN
            WRITE(IOUT,973) STGNEW(LAKE), SSMX(LAKE), LAKE
973         FORMAT(/1X,'WARNING-- COMPUTED STAGE OF ',F10.2,&
            &' IS GREATER THAN SPECIFIED MAXIMUM ',F10.2,&
            &' FOR LAKE ',I5)
         END IF
      END IF
   END DO
   IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 1061
   !
   !26------WRITE BUDGET SUMMARIES.
   WRITE(IOUT,1025) KPER, KSTP, DELT, PERTIM, TOTIM
1025 FORMAT(/1X,'PERIOD ',I5,5X,'TIME STEP',I6,5X,'TIME STEP LENGTH ',& !gsf
   &1PE11.4/1X,'PERIOD TIME ',E11.4,5X,'TOTAL SIMULATION TIME ',&
   &E11.4)
   WRITE(IOUT,1040)
1040 FORMAT(//17X,'HYDROLOGIC BUDGET SUMMARIES FOR SIMULATED LAKES'&
   &,/,5X,'(ALL FLUID FLUXES ARE VOLUMES ADDED TO THE LAKE DURING'&
   &,' PRESENT TIME STEP)'&
   &,/,5X,'------------------------------------------------------'&
   &,'----------------------------------')
   !dep REVISED PRINT STATEMENT WHEN THERE IS RUNOFF FROM UZF PACKAGE
   IF (IUNITUZF.EQ.0) THEN
      WRITE(IOUT,1045)
1045  FORMAT(1X,'LAKE',4X,'STAGE',9X,'VOLUME',5X,'VOL. CHANGE',6X,&
      &'PRECIP',5X,'EVAPORATION',8X,'RUNOFF')
   ELSE
      WRITE(IOUT,3045)
3045  FORMAT(' LAKE',4x,'STAGE',9X,'VOLUME',5X,'VOL. CHANGE',6X,&
      &'PRECIP',5X,'EVAPORATION',3X,' SPECIFIED RUNOFF',5X,&
      &'COMPUTED RUNOFF',4X, 'TOTAL RUNOFF')
   END IF
   !
   !27-----WRITE LAKE BUDGETS FOR A TIMES STEP (VOLUMES PER TIME STEP).
1061 DO 1100 NN=1,NLAKES
      PPTIN=PRECIP(NN)*DELT
      EOUT=EVAP(NN)*DELT
      SEEPUZF = SEEPUZ(NN)*DELT
      IF(RNF(NN).GE.ZERO) RUNF(NN) = RNF(NN)
      IF(RNF(NN).LT.ZERO) RUNF(NN) =-RNF(NN)*PRCPLK(NN)&
      &*BGAREA(NN)
      RUNFD = RUNF(NN)*DELT
      IF (IUNITUZF.GT.0) THEN
         RUNOFF(NN) = OVRLNDRNF(NN)*DELT
      ELSE
         RUNOFF(NN) = 0.0
      END IF
   !
      CUMPPT(NN)=CUMPPT(NN)+PPTIN
      CUMEVP(NN)=CUMEVP(NN)+EOUT
      CUMRNF(NN)=CUMRNF(NN)+RUNFD
      CUMUZF(NN)=CUMUZF(NN)+SEEPUZF
      IF (IUNITUZF.GT.0) THEN
         CUMLNDRNF(NN) = CUMLNDRNF(NN) + RUNOFF(NN)
      END IF
   !-lfk
      IF(ISS.NE.0) THEN
         DELVOL(NN)=0.0
         VOLINIT(NN)=VOL(NN)
      ELSE
         DELVOL(NN)=VOL(NN)-VOLOLD(NN)
      END IF
   !
      DELVOLLAK(NN)=DELVOL(NN)/DELT
   !-EDM
      IF(IUNIT(49).NE.0 ) THEN
         IF ( LKFLOWTYPE(1).EQ.'NA' ) THEN
            LKFLOWTYPE(1)='VOLUME'
            NLKFLWTYP = NLKFLWTYP + 1
         END IF
      ENDIF
      IF(IUNIT(49).NE.0 ) THEN
         IF ( LKFLOWTYPE(2).EQ.'NA' ) THEN
            LKFLOWTYPE(2)='DELVOL'
            NLKFLWTYP = NLKFLWTYP + 1
         END IF
      ENDIF
   !
      IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 1100
      IF(IUNITUZF.EQ.0) THEN
         IF(ISS.NE.0) THEN
            WRITE(IOUT,1049) NN,STGNEW(NN),VOL(NN),&
            &PPTIN,EOUT,RUNFD
         ELSE
            WRITE(IOUT,1050) NN,STGNEW(NN),VOL(NN),DELVOL(NN),&
            &PPTIN,EOUT,RUNFD
         END IF
      ELSE
         IF(ISS.NE.0) THEN
            WRITE(IOUT,3049) NN,STGNEW(NN),VOL(NN),&
            &PPTIN,EOUT,RUNFD,RUNOFF(NN),RUNFD+&
            &RUNOFF(NN)
         ELSE
            WRITE(IOUT,3050) NN,STGNEW(NN),VOL(NN),DELVOL(NN),&
            &PPTIN,EOUT,RUNFD,RUNOFF(NN),RUNFD+&
            &RUNOFF(NN)
         END IF
      END IF
1100 CONTINUE
1049 FORMAT(1X,I3,2(1X,1PE13.6),'     N/A (SS) ',3(1X,1PE13.6))
1050 FORMAT(1X,I3,6(1X,1PE13.6))
3049 FORMAT(1X,I3,2(1X,1PE13.6),'     N/A (SS) ',2(1X,1PE13.6),&
   &2X,1PE13.6,9X,1PE13.6,5X,1PE13.6)
3050 FORMAT(I4,5(1X,1PE13.6),2X,1PE13.6,9X,1PE13.6,5X,1PE13.6)
   IF(LWRT.LE.0.AND.ICBCFL.GT.0) THEN
      IF(IUNITUZF.EQ.0) THEN
         WRITE(IOUT,1046)
      ELSE
         WRITE(IOUT,1047)
      END IF
   END IF
1046 FORMAT(/12X,'GROUND WATER',12X,'SURFACE WATER',8X,'WATER'/1X,&
   &'LAKE',4X,'INFLOW',5X,'OUTFLOW',6X,'INFLOW',5X,'OUTFLOW'7X,'USE')
1047 FORMAT(/12X,'GROUND WATER',12X,'SURFACE WATER',8X,'WATER',10X,&
   &'UZF INFIL.'/1X,'LAKE',4X,'INFLOW',5X,'OUTFLOW',6X,'INFLOW',5X,&
   &'OUTFLOW',7X,'USE',8X,'FROM LAKE')
   !
   !28-----DETERMINE LAKE BUDGET ERROR FOR A TIME STEP.
   DO 1101 NN=1,NLAKES
      PPTIN=PRECIP(NN)*DELT
      EOUT=EVAP(NN)*DELT
      SEEPUZF = SEEPUZ(NN)*DELT
      IF(RNF(NN).GE.ZERO) RUNF(NN) = RNF(NN)
      IF(RNF(NN).LT.ZERO) RUNF(NN) =-RNF(NN)*PRCPLK(NN)&
      &*BGAREA(NN)
      RUNFD = RUNF(NN)*DELT
      IF (IUNITUZF.GT.0) THEN
         RUNOFF(NN) = OVRLNDRNF(NN)*DELT
      ELSE
         RUNOFF(NN) = 0.0
      END IF
      QIN=GWIN(NN)*DELT
      QOUT=GWOUT(NN)*DELT
      QSIN=SURFIN(NN)*DELT
      QSOUT=SURFOT(NN)*DELT
   !
      CUMGWI(NN)=CUMGWI(NN)+QIN
      CUMGWO(NN)=CUMGWO(NN)+QOUT
      CUMSWI(NN)=CUMSWI(NN)+QSIN
      CUMSWO(NN)=CUMSWO(NN)+QSOUT
   !-LFK
      WDRAW=WITHDRW(NN)*DELT
   !-LFK      Calculate accuracy of lake budget FOR TIME STEP
      CUMLKIN(NN)=PPTIN+RUNFD+RUNOFF(NN)+QIN+QSIN
      CUMLKOUT(NN)=EOUT+WDRAW+QOUT+QSOUT+SEEPUZF
      IF (CUMLKIN(NN).GT.CUMLKOUT(NN)) THEN
         DENOM=CUMLKIN(NN)
      ELSE IF (CUMLKOUT(NN).NE.0.0) THEN
         DENOM=CUMLKOUT(NN)
      ELSE
         DENOM=ABS(DELVOL(NN))
      END IF
      IF (abs(DENOM).GT.Closezero) THEN
         TSLAKERR(NN)=(CUMLKIN(NN)-CUMLKOUT(NN)-DELVOL(NN)&
         &+FLXINL(NN))*100./DENOM
      ELSE
   !            Note: if both ins & outs = 0.0, err set = 1e20
   !                 TSLAKERR(NN)=1.0E20
         TSLAKERR(NN)=0.0
      END IF
      IF(LWRT.LE.0.AND.ICBCFL.GT.0)THEN
         IF(IUNITUZF.EQ.0) THEN
            WRITE(IOUT,1051) NN,QIN,QOUT,QSIN,QSOUT,WDRAW
         ELSE
            WRITE(IOUT,1052) NN,QIN,QOUT,QSIN,QSOUT,WDRAW,SEEPUZF
         END IF
      END IF
1101 END DO
1051 FORMAT(1X,I3,1P,5E12.4)
1052 FORMAT(1X,I3,1P,5E12.4,5X,E12.4)
   IF(LWRT.LE.0.AND.ICBCFL.GT.0) WRITE(IOUT,1035)
   !-LFK
1035 FORMAT(/6X,'CONNECTED LAKE',3X,'TIME-STEP'&
   &,9X,'STAGE-CHANGE',10X,'PERCENT'/1X,'LAKE',4X,'INFLUX',&
   &6X,'SURFACE AREA',3X,'TIME STEP',2X,'CUMULATIVE',4X,&
   &'DISCREPANCY')
   !
   !29-----DETERMINE CUMULATIVE LAKE BUDGET ERRORS.
   DO 1105 NN=1,NLAKES
   !dep  All values are volumes per time (times DELT)
      WDRAW=WITHDRW(NN)*DELT
   !
      CUMWDR(NN)=CUMWDR(NN)+WDRAW
      CUMFLX(NN)=CUMFLX(NN)+FLXINL(NN)
   !-LFK
      IF(ISS.NE.0) THEN
         CUMVOL(NN)=0.0
      ELSE
         CUMVOL(NN)=VOL(NN)-VOLINIT(NN)
      END IF
   !-LFK      Calculate accuracy of CUMULATIVE lake budget
      CUMLKIN(NN)=CUMPPT(NN)+CUMRNF(NN)+CUMGWI(NN)+CUMSWI(NN)&
      &+CUMLNDRNF(NN)
      CUMLKOUT(NN)=CUMEVP(NN)+CUMWDR(NN)+CUMGWO(NN)+CUMSWO(NN)&
      &+CUMUZF(NN)
      IF (CUMLKIN(NN).GT.CUMLKOUT(NN)) THEN
         DENOM=CUMLKIN(NN)
      ELSE IF (CUMLKOUT(NN).NE.0.0) THEN
         DENOM=CUMLKOUT(NN)
      ELSE
         DENOM=ABS(CUMVOL(NN))
      END IF
      IF (DENOM.NE.0.0) THEN
         CMLAKERR(NN)=(CUMLKIN(NN)-CUMLKOUT(NN)-CUMVOL(NN)&
         &+CUMFLX(NN))*100./DENOM
      ELSE
   !     Note: if both ins & outs = 0.0, err set = 1e20
   !                 CMLAKERR(NN)=1.0E20
         CMLAKERR(NN)=0.0
      END IF
      IF(LWRT.LE.0.AND.ICBCFL.GT.0) THEN
         IF(ISS.NE.0) THEN
            WRITE(IOUT,1054) NN,FLXINL(NN),SURFA(NN),&
            &TSLAKERR(NN)
         ELSE
            WRITE(IOUT,1055) NN,FLXINL(NN),SURFA(NN),&
            &DELH(NN),TDELH(NN),TSLAKERR(NN)
         END IF
1054     FORMAT(1X,I3,1P,1E13.4,2X,E13.4,'     N/A (SS) ','  N/A (SS) ',&
         &3X,0P,F9.3)
1055     FORMAT(1X,I3,1P,1E13.4,2X,2E13.4,E12.4,3X,0P,F9.3)
      END IF
1105 END DO
   ! 1055   FORMAT(1X,I3,1PE13.4,2X,2PE13.4,1PE12.4,3X,0P,F9.3)
   IF(LWRT.GT.0.OR.ICBCFL.LE.0) GO TO 1041
   WRITE(IOUT,1301)
   !
   !30------WRITE CUMULATIVE BUDGET SUMMARIES.
   WRITE(IOUT,2040)
2040 FORMAT(//12X,'CUMULATIVE HYDROLOGIC BUDGET SUMMARIES FOR SIMULATED&
   & LAKES'&
   &,/,5X,'(ALL FLUID FLUXES ARE SUMS OF VOLUMES ADDED TO THE'&
   &,' LAKE SINCE INITIAL TIME)'&
   &,/,5X,'------------------------------------------------------'&
   &,'---------------------')
   !dep  added computed runoff from UZF Package to lake budget
   IF(IUNITUZF.LE.0) THEN
      WRITE(IOUT,2045)
2045  FORMAT(1X,'LAKE',7X,'PRECIP',7X,'EVAP',7X,'RUNOFF')
   ELSE
      WRITE(IOUT,4045)
4045  FORMAT(1X,'LAKE',7X,'PRECIP',7X,'EVAP',5X,'SPECIFIED RUNOFF',&
      &3X,'COMPUTED RUNOFF',3X,'TOTAL RUNOFF')
   END IF
   DO 2100 NN=1,NLAKES
      IF(IUNITUZF.LE.0) THEN
         WRITE(IOUT,2050) NN,CUMPPT(NN),CUMEVP(NN),CUMRNF(NN)
      ELSE
         WRITE(IOUT,4050) NN,CUMPPT(NN),CUMEVP(NN),CUMRNF(NN),&
         &CUMLNDRNF(NN),CUMRNF(NN)+CUMLNDRNF(NN)
      END IF
2100 CONTINUE
2050 FORMAT(1X,I3,3X,1P,3E12.4)
4050 FORMAT(1X,I3,3X,1P,2E12.4,3(5X,E12.4))
   IF(IUNITUZF.LE.0) THEN
      WRITE(IOUT,2046)
   ELSE
      WRITE(IOUT,2047)
   END IF
2046 FORMAT(/12X,'GROUND WATER',12X,'SURFACE WATER'/1X,'LAKE',4X,&
   &'INFLOW',5X,'OUTFLOW',6X,'INFLOW',5X,'OUTFLOW')
2047 FORMAT(/12X,'GROUND WATER',12X,'SURFACE WATER',7X,'UZF INFIL.'&
   &/1X,'LAKE',4X,'INFLOW',5X,'OUTFLOW',6X,&
   &'INFLOW',5X,'OUTFLOW',5X,'FROM LAKE')
   DO 2101 NN=1,NLAKES
      IF(IUNITUZF.LE.0) THEN
         WRITE(IOUT,2051) NN,CUMGWI(NN),CUMGWO(NN),CUMSWI(NN),&
         &CUMSWO(NN)
      ELSE
         WRITE(IOUT,2052) NN,CUMGWI(NN),CUMGWO(NN),CUMSWI(NN),&
         &CUMSWO(NN),CUMUZF(NN)
      END IF
2101 END DO
2051 FORMAT(1X,I3,1P,4E12.4)
2052 FORMAT(1X,I3,1P,4E12.4,2X,E12.4)
   WRITE(IOUT,2035)
   !-LFK
2035 FORMAT(/9X,'WATER',4X,'CONNECTED LAKE',3X,'CHANGE',7X,'PERCENT'/&
   &1X,'LAKE',5X,'USE',9X,'INFLUX',7X,'IN VOL.',4X,&
   &'DISCREPANCY')
   DO 2105 NN=1,NLAKES
   !-LFK
      IF(ISS.NE.0) THEN
         WRITE(IOUT,2054) NN,CUMWDR(NN),CUMFLX(NN),CMLAKERR(NN)
      ELSE
         WRITE(IOUT,2055) NN,CUMWDR(NN),CUMFLX(NN),CUMVOL(NN),CMLAKERR(NN)
      END IF
2105 CONTINUE
   !-LFK
2054 FORMAT(1X,I3,1P,2E13.4,'    N/A (SS)   ',0P,F9.3)
2055 FORMAT(1X,I3,1P,3E13.4,2X,0P,F9.3)
   WRITE(IOUT,1301)
   IF(KCNT.LE.0) GO TO 11041
   IF (KCNT.GT.1) THEN
      WRITE(IOUT,11055) KCNT
11055 FORMAT(/1X,I3,' CONNECTED LAKE SETS'/)
      DO 11056 IIC=1,KCNT
         JIC = MSUB1(IIC)
         WRITE(IOUT,11057) JIC, (MSUB(LIC,IIC),LIC=1,JIC)
11057    FORMAT(1X,I3,' LAKES:  ',25I3)
11056 CONTINUE
   ELSE
      WRITE(IOUT,21055) KCNT
21055 FORMAT(/1X,I3,' CONNECTED LAKE SET'/)
   !-LFK
      IIC=1
      JIC = MSUB1(IIC)
      WRITE(IOUT,11057) JIC, (MSUB(LIC,KCNT),LIC=1,JIC)
   END IF
11041 CONTINUE
1041 CONTINUE
   !
   !31-----dep   Moved Call to gage to follow budget summaries 6/9/2009
   !-LFK   don't call GAG5LO from here in Lake Package if GWT is active
   !dep replaced stgold2 and stages with delh and tdelh
   IF(IUNITGWT.LE.0) THEN
      IF (IUNITGAGE.GT.0)THEN
         CALL SGWF2GAG7LO(IUNITGWT,IUNITUZF,XLAKES,TOTIM,&
         &GWIN,GWOUT,SEEP,FLXINL,VOLOLD,XLKOLD,XLAKINIT,NSOL)
      END IF
   END IF
   !
   !32-----IF C-B-C TERMS WILL BE SAVED THEN WRITE TO DISK.
   IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,ILKCB,BUFF,NCOL,NROW,&
   &NLAY,IOUT)
   !
   !32A-----MOVE RATES INTO VBVL FOR PRINTING BY MODULE BAS OT.
1200 VBVL(3,MSUM)=RATIN
   VBVL(4,MSUM)=RATOUT
   !
   !32B-----MOVE PRODUCT OF RATE AND TIME STEP INTO VBVL ACCUMULATORS.
   VBVL(1,MSUM)=VBVL(1,MSUM)+RATIN*DELT
   VBVL(2,MSUM)=VBVL(2,MSUM)+RATOUT*DELT
   !
   !32C-----MOVE BUDGET TERM LABELS INTO VBVM FOR PRINTING BY BAS OT.
   VBNM(MSUM)=TEXT
   !33------INCREASE BUDGET COUNTER.
   MSUM=MSUM+1
   !
   !        Substitute Lake Stage values for HNOFLO values at non-dry lake
   !        cells in HNEW array; loop over all lake nodes.  If Lake Stage
   !        is below bottom of lake cell, set HNEW = HNOFLO.
   !
   DO 1900 L=1,LKNODE
      IL=ILAKE(1,L)
      ILL=IL-1
      IR=ILAKE(2,L)
      IC=ILAKE(3,L)
      LAKE=ILAKE(4,L)
      ITYPE = (ILAKE(5,L)+1)/2
      IF(ITYPE.EQ.3) ITYPE=0
      IF(ITYPE.NE.0) GO TO 1900
      IF(IBOUND(IC,IR,IL).GT.0) THEN
         BOTLK = BOTM(IC,IR,LBOTM(ILL))
         IF (STGNEW(LAKE).GT.BOTLK) HNEW(IC,IR,ILL)=STGNEW(LAKE)
         IF (ILL.GT.1) THEN
            ILL=ILL-1
            DO 1890 IL2=1,ILL
               IF (STGNEW(LAKE).GT.BOTM(IC,IR,LBOTM(IL2))) THEN
                  HNEW(IC,IR,IL2)=STGNEW(LAKE)
               ELSE
                  HNEW(IC,IR,IL2)=HNOFLO
               END IF
1890        CONTINUE
         END IF
      ELSE IF(IBOUND(IC,IR,IL).EQ.0) THEN
         BOTLK = BOTM(IC,IR,LBOTM(IL))
         IF (STGNEW(LAKE).GT.BOTLK) HNEW(IC,IR,IL)=STGNEW(LAKE)
         IF (IL.GT.1) THEN
            ILL=IL-1
            DO 1892 IL2=1,ILL
               IF (STGNEW(LAKE).GT.BOTM(IC,IR,LBOTM(IL2))) THEN
                  HNEW(IC,IR,IL2)=STGNEW(LAKE)
               ELSE
                  HNEW(IC,IR,IL2)=HNOFLO
               END IF
1892        CONTINUE
         END IF
      END IF
   !
1900 CONTINUE
   !
   !34-----RETURN.
   RETURN
END
   !
SUBROUTINE SGWF2LAK7SFR7RPS()
   !
   !    *******************************************************************
   !--  IF STREAMS EXIST, DEFINE CONNECTIONS BETWEEN LAKES AND STREAMS
   !    *******************************************************************
   !
   !    -------------------------------------------------------------------
   !        SPECIFICATIONS:
   !    -------------------------------------------------------------------
   USE GWFLAKMODULE, ONLY: NLAKES, NTRB, NDV, ITRB, IDIV, IRK, RAMP,&
   &DEADPOOLVOL
   USE GLOBAL,       ONLY: IOUT, NODES
   USE GWFSFRMODULE, ONLY: NSS, IDIVAR, IOTSG, SEG, ISEG, CLOSEZERO
   DOUBLE PRECISION :: SPILL,V
   DOUBLE PRECISION VOLTERP
   EXTERNAL VOLTERP
   !
   !-- DOUBLE CHECK SIZE OF IRK (STORED IN BUFF) vs. NLAKES
   !
   IF ((NLAKES*2).GT.NODES) THEN
      WRITE (IOUT,*) '***NLAKES too large for BUFF in Subroutine GWF2&
      &LAK7SFR7RPS***  STOP EXECUTION'
      CALL USTOP(' ')
   END IF
   !
   !-- INITIALIZE ARRAYS
   !
   !     DO 50 I=1,NSS
   !     DO 50 LK=1,NLAKES
   !     ITRB(LK,I) = 0
   !  50 IDIV(LK,I) = 0
   DO 55 LK=1,NLAKES
      IRK(1,LK) = 0
55 IRK(2,LK) = 0
   NTRB = 0
   NDV = 0
   !
   !-- Build arrays to define lake tributary & diversion links ...
   !        based on stream package input data
   !
   !---  Stream Inflow to Lakes
   DO 100 LSEG=1,NSS
      IF(IOTSG(LSEG).LT.0) THEN
         LAKE = -IOTSG(LSEG)
         IRK(1,LAKE) = IRK(1,LAKE) + 1
         K1 = IRK(1,LAKE)
         ITRB(LAKE,K1) = LSEG
         IF(IRK(1,LAKE).GT.NTRB) NTRB = IRK(1,LAKE)
      ENDIF
   !
   !---  Stream Outflow from Lakes
      IF(IDIVAR(1,LSEG).LT.0) THEN
         LAKE = -IDIVAR(1,LSEG)
         IRK(2,LAKE) = IRK(2,LAKE) + 1
         K1 = IRK(2,LAKE)
         IDIV(LAKE,K1) = LSEG
         IF(IRK(2,LAKE).GT.NDV) NDV = IRK(2,LAKE)
   ! CALCULATE DEAD POOL STORAGE
         SPILL = SEG(8,LSEG) + RAMP
         V = VOLTERP(SPILL,LAKE)
   !---   Account for multiple diversions out of a lake
         IF(V.GT.CLOSEZERO) THEN
            IF(V.LT.DEADPOOLVOL(LAKE).OR.&
            &DEADPOOLVOL(LAKE).LT.CLOSEZERO) THEN
               DEADPOOLVOL(LAKE) = V
            ENDIF
         ENDIF
         WRITE(IOUT,*)
         WRITE(IOUT,9008)LAKE,LSEG,V
         WRITE(IOUT,*)
9008     FORMAT(6X,'DEAD POOL STORAGE FOR LAKE ',I4,' BELOW SEGMENT ',&
         &I4,' IS EQUAL TO ',E20.10)
      ENDIF
100 CONTINUE
   !
   !--  PRINT LAKE INFLOW STREAM SEGMENTS.
   WRITE(IOUT,10)
10 FORMAT(6X,'LAKE ',4X,'INFLOWING STREAM SEGMENT')
   DO 520 IK=1,NLAKES
      DO 519 JK=1,NSS
         IF(ITRB(IK,JK).LE.0) GO TO 521
519   CONTINUE
521   JK1 = JK - 1
      IF(JK1.GT.0) WRITE(IOUT,15) IK,(ITRB(IK,JK),JK=1,JK1)
15    FORMAT(5X,I5,14X,100I5)
520 CONTINUE
   WRITE(IOUT,103) NTRB
103 FORMAT(/1X,'MAXIMUM NUMBER OF STREAMS INFLOWING TO A',&
   &' LAKE IS',I5/)
   !
   !--  PRINT LAKE STREAM OUTFLOW SEGMENT (FROM A LAKE) NUMBERS.
   !
   WRITE(IOUT,13)
13 FORMAT(6X,'LAKE ',4X,'OUTFLOWING STREAM',' SEGMENT')
   DO 600 IK=1,NLAKES
      DO 523 JK=1,NSS
         IF(IDIV(IK,JK).LE.0) GO TO 527
523   CONTINUE
527   JK1 = JK - 1
      IF(JK1.GT.0) WRITE(IOUT,15) IK,(IDIV(IK,JK),JK=1,JK1)
600 CONTINUE
   !
   !dep-- PRINT WARNING IF OUTFLOWING STREAM IS ASSIGNED ICALC =0.
   !dep    ADDED OCTOBER 15, 2004; DAVID PRUDIC
   DO ls = 1, NSS
      IF (IDIVAR(1,ls).LT.0) THEN
         lk = -IDIVAR(1,ls)
         IF (ISEG(1,ls).LE.0 .AND. SEG(2,ls).LE.0.0) THEN
            WRITE (IOUT, 9007) ls, lk, ISEG(1,ls), SEG(2,ls)
         END IF
      END IF
   END DO
   WRITE(IOUT,133) NDV
133 FORMAT(/1X,'MAXIMUM NUMBER OF STREAMS OUTFLOWING',&
   &' FROM A LAKE IS',I5/)
9007 FORMAT(/, ' WARNING****  OUTFLOWING STREAM SEGMENT', I6,&
   &' FROM LAKE', I6, ' HAS AN ICALC VALUE OF', I6,&
   &' AND FLOW INTO THE SEGMENT IS', E12.4, /,&
   &' NO OUTFLOW FROM THE LAKE INTO ',&
   &'SEGMENT WILL BE SIMULATED', /,&
   &' SUGGEST CHANGING ICALC TO ANOTHER OPTION')
   !
   !-- RETURN
   RETURN
END
SUBROUTINE SGWF2LAK7BCF7RPS(LWRT)
   !
   !     ******************************************************************
   !     COMPUTE VERTICAL CONDUCTANCES AND HORIZONTAL CONDUCTANCES PER UNIT
   !     THICKNESS FOR LAKES WHEN BCF PACKAGE IS USED
   !     ******************************************************************
   !
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFLAKMODULE, ONLY: LKNODE, BEDLAK, LKARR1, ILAKE, CNDFCT
   USE GLOBAL,       ONLY: NLAY, IOUT, DELR, DELC, LAYHDT
   !!      USE GLOBAL,       ONLY: NLAY, IOUT, DELR, DELC, LAYHDT,NCOL,NROW
   USE GWFBCFMODULE, ONLY: IWDFLG, HY, CVWD, TRPY
   !
   IF ( LWRT <= 0 ) WRITE(IOUT,108)
108 FORMAT(//9X,'C',15X,'INTERFACE CONDUCTANCES BETWEEN LAKE AND ',&
   &'AQUIFER CELLS'/&
   &3X,'L',5X,'O',10X,'(IF TYPE = 6, CONDUCTANCE (L^2/T) IS ',&
   &'BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.)',/&
   &3X,'A',5X,'L',2X,'L',2X,'T',&
   &4X,'(IF TYPE = 1 TO 4, CONDUCTANCES ARE PER UNIT SATURATED ',&
   &'THICKNESS (L/T).)'/&
   &3X,'Y',2X,'R',2X,'U',2X,'A',2X,'Y'/&
   &3X,'E',2X,'O',2X,'M',2X,'K',2X,'P',&
   &24X,'LAKEBED',6X,'C O N D U C T A N C E S'/3X,'R',2X,'W',2X,&
   &'N',2X,'E',&
   &2X,'E',5X,'DELTA Y',3X,'DELTA X',2X,'LEAKANCE',3X,'LAKEBED',3X,&
   &'AQUIFER',2X,'COMBINED'/1X,79('_'))
   !
   IWRN = 0
   IWRN1 = 0
   DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
      CNDFCT(II) = 0.0
   !  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      NTYP = NTYP + 1
      IF(NTYP.EQ.1) THEN
   !
   !  Vertical Conductance
   !    for vertical interface, "K" is layer below bottom of lake
   !
         CNDFC1=0.0
         IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) GO TO 315
         IF(BEDLAK(II).LE.0.0) GO TO 315
         IWRN1 = 1
         CNDFC1 = BEDLAK(II)*DELR(I)*DELC(J)
         IF (IWDFLG.EQ.0) THEN
            CNDFCT(II) = CNDFC1
         ELSE
            IF(CVWD(I,J,K-1).LE.0.0.OR.CNDFC1.LE.0.0) GO TO 315
            CNDFCT(II) = 1.0/(0.5/CVWD(I,J,K-1)+1.0/CNDFC1)
         END IF
315      IF (IWDFLG.EQ.0) THEN
            IF ( LWRT <= 0 ) WRITE(IOUT,7324) (ILAKE(I1,II),I1=1,5),&
            &DELC(J),DELR(I),BEDLAK(II),CNDFC1,CNDFCT(II)
   !-lfk
7324        FORMAT(1X,5I3,2X,1P,4E10.2,10X,E11.3)
   ! 7324     FORMAT(1X,5I3,2X,1P,4E10.2,10X,E10.2)
         ELSE
            CVWD2= 2.0*CVWD(I,J,K-1)
            IF ( LWRT <= 0 ) WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),&
            &DELC(J),DELR(I),BEDLAK(II),CNDFC1,CVWD2,CNDFCT(II)
   !-lfk
7325        FORMAT(1X,5I3,2X,1P,5E10.2,E11.3)
   ! 7325     FORMAT(1X,5I3,2X,1P,6E10.2)
         END IF
      ELSE
   !
   !  Horizontal conductance
   !
   !  HY not read in, thus unavailable.
   !
   !dep  348   IF(LAYHDT(K).EQ.0) THEN
         IF(LAYHDT(K).EQ.0) THEN
            IF(NTYP.EQ.2) CNDFCT(II) = BEDLAK(II)*DELC(J)
            IF(NTYP.EQ.3) CNDFCT(II) = BEDLAK(II)*DELR(I)
            IF ( LWRT <= 0 ) WRITE(IOUT,7324) (ILAKE(I1,II),I1=1,5),&
            &DELC(J),DELR(I),BEDLAK(II),CNDFCT(II),CNDFCT(II)
            IWRN = 1
         ELSE
   !
   !  HY read in, thus available.
   !
            TT = HY(I,J,K)
            IF(NTYP.EQ.2) CNDFC2 = 2.0*TT*DELC(J)/DELR(I)
            IF(NTYP.EQ.3) CNDFC2 = 2.0*TRPY(K)*TT*DELR(I)/DELC(J)
            IF(NTYP.EQ.2) CNDFC1 = BEDLAK(II)*DELC(J)
            IF(NTYP.EQ.3) CNDFC1 = BEDLAK(II)*DELR(I)
            IF (CNDFC1.GT.0.0.AND.CNDFC2.GT.0.0)&
            &CNDFCT(II) = 1.0/(1.0/CNDFC2+1.0/CNDFC1)
            IF ( LWRT <= 0 )WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),&
            &DELR(I),BEDLAK(II),CNDFC1,CNDFC2,CNDFCT(II)
         END IF
      END IF
350 CONTINUE
   !
   !  WRITE WARNINGS ON LAKE/AQUIFER CONDUCTANCES, IF NECESSARY
   IF(IWRN.EQ.1.OR.IWRN1.EQ.1) WRITE(IOUT,345)
345 FORMAT(//5X,'NOTE: INFORMATION ABOUT CALCULATED LAKE/AQUIFER C&
   &ONDUCTANCES WHEN USING BCF PACKAGE FOLLOWS: '/)
   IF(IWRN.EQ.1) WRITE(IOUT,346)
346 FORMAT(1X,'NODE(S) ADJACENT TO LAKE IN CONFINED LAYER:'/&
   &1X,'LAKE/AQUIFER CONDUCTANCES BASED SOLELY ON LAKEBED SPECIFIC&
   &ATION'/)
   IF(IWRN1.EQ.1) WRITE(IOUT,347)
347 FORMAT(1X,'IF WETDRY FLAG NOT TURNED ON, VERTICAL LEAKANCES AR&
   &E NOT SAVED:'/1X,'THEREFORE, LAKE/AQUIFER CONDUCTANCES ARE BASED S&
   &OLELY ON LAKEBED SPECIFICATION'/)
   IF(IWRN.EQ.1.OR.IWRN1.EQ.1) WRITE(IOUT,'(//)')
   !
   RETURN
END
   !
SUBROUTINE SGWF2LAK7LPF7RPS(LWRT)
   !
   !     ******************************************************************
   !     COMPUTE VERTICAL CONDUCTANCES AND HORIZONTAL CONDUCTANCES PER UNIT
   !     THICKNESS FOR LAKES WHEN LPF PACKAGE IS USED
   !     ******************************************************************
   !
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFLAKMODULE, ONLY: LKNODE, BEDLAK, LKARR1, ILAKE, CNDFCT
   USE GLOBAL,       ONLY: NLAY, IOUT, LBOTM, LAYCBD, DELR, DELC,&
   &BOTM
   USE GWFLPFMODULE, ONLY: CHANI, LAYVKA, VKA, VKCB, HANI, HK
   !
   IF ( LWRT <= 0 )WRITE(IOUT,108)
108 FORMAT(//9X,'C',15X,'INTERFACE CONDUCTANCES BETWEEN LAKE AND ',&
   &'AQUIFER CELLS'/&
   &3X,'L',5X,'O',10X,'(IF TYPE = 6, CONDUCTANCE (L^2/T) IS ',&
   &'BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.)',/&
   &3X,'A',5X,'L',2X,'L',2X,'T',&
   &4X,'(IF TYPE = 1 TO 4, CONDUCTANCES ARE PER UNIT SATURATED ',&
   &'THICKNESS (L/T).)'/&
   &3X,'Y',2X,'R',2X,'U',2X,'A',2X,'Y'/&
   &3X,'E',2X,'O',2X,'M',2X,'K',2X,'P',&
   &24X,'LAKEBED',6X,'C O N D U C T A N C E S'/3X,'R',2X,'W',2X,&
   &'N',2X,'E',&
   &2X,'E',5X,'DELTA Y',3X,'DELTA X',2X,'LEAKANCE',3X,'LAKEBED',3X,&
   &'AQUIFER',2X,'COMBINED'/1X,79('_'))
   !
   DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
      CAQ = 0.0
      CNDFCT(II) = 0.0
   !  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      NTYP=NTYP + 1
      IF(NTYP.EQ.1) THEN
   !
   !  Vertical Conductance
   !    for vertical interface, "K" is layer below bottom of lake
         CNDFC1=0.0
         IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) GO TO 315
         IF(BEDLAK(II).LE.0.0) GO TO 315
         CNDFC1 = BEDLAK(II)*DELR(I)*DELC(J)
         VK = 0.0
         IF(LAYVKA(K).EQ.0) THEN
            VK=VKA(I,J,K)
         ELSE IF ( VKA(I,J,K) > 0.0 ) THEN    !RGN 8/21/17 CHECK DIVIDE BY ZERO
            VK=HK(I,J,K)/VKA(I,J,K)
         END IF
   !   skip if zero vk
         IF(VK.LE.0.0) GO TO 350
         BBOT=BOTM(I,J,LBOTM(K))
         TTOP=BOTM(I,J,LBOTM(K)-1)
         CAQ=VK*DELR(I)*DELC(J)/((TTOP-BBOT)*0.5)
         IF(LAYCBD(K-1).GT.0) THEN
   !   skip if zero vkcb
            IF(VKCB(I,J,LAYCBD(K-1)).LE.0.0) GO TO 350
            BBOT=BOTM(I,J,LBOTM(K)-1)
            TTOP=BOTM(I,J,LBOTM(K-1))
            CCB=VKCB(I,J,LAYCBD(K-1))*DELR(I)*DELC(J)/(TTOP-BBOT)
            !include VKCB
            CAQ = 1.0/(1.0/CAQ + 1.0/CCB)
         END IF
         CNDFCT(II) = 1.0/(1.0/CAQ+1.0/CNDFC1)
315      IF ( LWRT <= 0 )WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),&
         &DELC(J),DELR(I),BEDLAK(II),CNDFC1,CAQ,CNDFCT(II)
      ELSE
   !
   !  Horizontal conductance
   !
         TT = HK(I,J,K)
   ! X-DIRECTION
         IF(NTYP.EQ.2) CNDFC2 = 2.0*TT*DELC(J)/DELR(I)
   ! Y-DIRECTION
         IF(NTYP.EQ.3) THEN
            IF(CHANI(K).LE.0) THEN
               KH=-CHANI(K)
               CNDFC2 = 2.0*HANI(I,J,KH)*TT*DELR(I)/DELC(J)
            ELSE
               CNDFC2 = 2.0*CHANI(K)*TT*DELR(I)/DELC(J)
            END IF
         END IF
         IF(NTYP.EQ.2) CNDFC1 = BEDLAK(II)*DELC(J)
         IF(NTYP.EQ.3) CNDFC1 = BEDLAK(II)*DELR(I)
         IF (CNDFC1.GT.0.0.AND.CNDFC2.GT.0.0)&
         &CNDFCT(II) = 1.0/(1.0/CNDFC2+1.0/CNDFC1)
         IF ( LWRT <= 0 )WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),&
         &DELC(J),DELR(I),BEDLAK(II),CNDFC1,CNDFC2,CNDFCT(II)
   !-lfk
7325     FORMAT(1X,5I3,2X,1P,5E10.2,E11.3)
   ! 7325   FORMAT(1X,5I3,2X,1P,6E10.2)
      END IF
350 CONTINUE
   !
   RETURN
END
   !
SUBROUTINE SGWF2LAK7UPW1RPS(LWRT)
   !
   !     ******************************************************************
   !     COMPUTE VERTICAL CONDUCTANCES AND HORIZONTAL CONDUCTANCES PER UNIT
   !     THICKNESS FOR LAKES WHEN UPW PACKAGE IS USED
   !     ******************************************************************
   !
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFLAKMODULE, ONLY: LKNODE, BEDLAK, LKARR1, ILAKE, CNDFCT
   USE GLOBAL,       ONLY: NLAY, IOUT, LBOTM, LAYCBD, DELR, DELC,&
   &BOTM
   USE GWFUPWMODULE, ONLY: CHANI, LAYVKAUPW, VKAUPW, VKCB, HANI,&
   &HKUPW
   !
   IF ( LWRT <= 0 )WRITE(IOUT,108)
108 FORMAT(//9X,'C',15X,'INTERFACE CONDUCTANCES BETWEEN LAKE AND ',&
   &'AQUIFER CELLS'/&
   &3X,'L',5X,'O',10X,'(IF TYPE = 6, CONDUCTANCE (L^2/T) IS ',&
   &'BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.)',/&
   &3X,'A',5X,'L',2X,'L',2X,'T',&
   &4X,'(IF TYPE = 1 TO 4, CONDUCTANCES ARE PER UNIT SATURATED ',&
   &'THICKNESS (L/T).)'/&
   &3X,'Y',2X,'R',2X,'U',2X,'A',2X,'Y'/&
   &3X,'E',2X,'O',2X,'M',2X,'K',2X,'P',&
   &24X,'LAKEBED',6X,'C O N D U C T A N C E S'/3X,'R',2X,'W',2X,&
   &'N',2X,'E',&
   &2X,'E',5X,'DELTA Y',3X,'DELTA X',2X,'LEAKANCE',3X,'LAKEBED',3X,&
   &'AQUIFER',2X,'COMBINED'/1X,79('_'))
   !
   DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
      CAQ = 0.0
      CNDFCT(II) = 0.0
   !  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      NTYP=NTYP + 1
      IF(NTYP.EQ.1) THEN
   !
   !  Vertical Conductance
   !    for vertical interface, "K" is layer below bottom of lake
         CNDFC1=0.0
         IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) GO TO 315
         IF(BEDLAK(II).LE.0.0) GO TO 315
         CNDFC1 = BEDLAK(II)*DELR(I)*DELC(J)
         IF(LAYVKAUPW(K).EQ.0) THEN
            VK=VKAUPW(I,J,K)
         ELSE IF ( VKAUPW(I,J,K) > 0.0 ) THEN    !RGN 8/21/17 CHECK DIVIDE BY ZERO
            VK=HKUPW(I,J,K)/VKAUPW(I,J,K)
         END IF
   !   skip if zero vk
         IF(VK.LE.0.0) GO TO 350
         BBOT=BOTM(I,J,LBOTM(K))
         TTOP=BOTM(I,J,LBOTM(K)-1)
         IF(TTOP-BBOT.LT.1.0e-20) GO TO 350
         CAQ=VK*DELR(I)*DELC(J)/((TTOP-BBOT)*0.5)
         IF(LAYCBD(K-1).GT.0) THEN
   !   skip if zero vkcb
            IF(VKCB(I,J,LAYCBD(K)).LE.0.0) GO TO 350
            BBOT=BOTM(I,J,LBOTM(K)-1)
            TTOP=BOTM(I,J,LBOTM(K-1))
            IF(TTOP-BBOT.LT.1.0e-20) GO TO 350
            CCB=VKCB(I,J,LAYCBD(K-1))*DELR(I)*DELC(J)/(TTOP-BBOT)
            !include VKCB
            CAQ = 1.0/(1.0/CAQ + 1.0/CCB)
         END IF
         CNDFCT(II) = 1.0/(1.0/CAQ+1.0/CNDFC1)
315      IF ( LWRT <= 0 )WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),&
         &DELC(J),DELR(I),BEDLAK(II),CNDFC1,CAQ,CNDFCT(II)
      ELSE
   !
   !  Horizontal conductance
   !
         TT = HKUPW(I,J,K)
   ! X-DIRECTION
         IF(NTYP.EQ.2) CNDFC2 = 2.0*TT*DELC(J)/DELR(I)
   ! Y-DIRECTION
         IF(NTYP.EQ.3) THEN
            IF(CHANI(K).LE.0) THEN
               KH=-CHANI(K)
               CNDFC2 = 2.0*HANI(I,J,KH)*TT*DELR(I)/DELC(J)
            ELSE
               CNDFC2 = 2.0*CHANI(K)*TT*DELR(I)/DELC(J)
            END IF
         END IF
         IF(NTYP.EQ.2) CNDFC1 = BEDLAK(II)*DELC(J)
         IF(NTYP.EQ.3) CNDFC1 = BEDLAK(II)*DELR(I)
         IF (CNDFC1.GT.0.0.AND.CNDFC2.GT.0.0)&
         &CNDFCT(II) = 1.0/(1.0/CNDFC2+1.0/CNDFC1)
         IF ( LWRT <= 0 )WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),&
         &DELC(J),DELR(I),BEDLAK(II),CNDFC1,CNDFC2,CNDFCT(II)
7325     FORMAT(1X,5I3,2X,1P,6E10.2)
      END IF
350 CONTINUE
   !
   RETURN
END
   !
   !
SUBROUTINE SGWF2LAK7HUF7RPS(LWRT)
   !
   !     ******************************************************************
   !     COMPUTE VERTICAL CONDUCTANCES AND HORIZONTAL CONDUCTANCES PER UNIT
   !     THICKNESS FOR LAKES WHEN HUF PACKAGE IS USED
   !     ******************************************************************
   !
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFLAKMODULE, ONLY: LKNODE, BEDLAK, LKARR1, ILAKE, CNDFCT
   USE GLOBAL,       ONLY: NLAY, IOUT, DELR, DELC
   !!      USE GLOBAL,       ONLY: NLAY, IOUT, LBOTM, DELR, DELC, BOTM
   !-lfk      USE GWFLPFMODULE, ONLY: VKA, HK
   USE GWFHUFMODULE, ONLY: VKAH !!,HK,HKCC
   !
   IF ( LWRT <= 0 )WRITE(IOUT,108)
108 FORMAT(//9X,'C',15X,'INTERFACE CONDUCTANCES BETWEEN LAKE AND ',&
   &'AQUIFER CELLS'/&
   &3X,'L',5X,'O',10X,'(IF TYPE = 6, CONDUCTANCE (L^2/T) IS ',&
   &'BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.)',/&
   &3X,'A',5X,'L',2X,'L',2X,'T',&
   &4X,'(IF TYPE = 1 TO 4, CONDUCTANCES ARE PER UNIT SATURATED ',&
   &'THICKNESS (L/T).)'/&
   &3X,'Y',2X,'R',2X,'U',2X,'A',2X,'Y'/&
   &3X,'E',2X,'O',2X,'M',2X,'K',2X,'P',&
   &24X,'LAKEBED',6X,'C O N D U C T A N C E S'/3X,'R',2X,'W',2X,&
   &'N',2X,'E',&
   &2X,'E',5X,'DELTA Y',3X,'DELTA X',2X,'LEAKANCE',3X,'LAKEBED',3X,&
   &'AQUIFER',2X,'COMBINED'/1X,79('_'))
   !
   DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
      CAQ = 0.0
      CNDFCT(II) = 0.0
   !  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      NTYP=NTYP + 1
      IF(NTYP.EQ.1) THEN
   !
   !  Vertical Conductance
   !    for vertical interface, "K" is layer below bottom of lake
         CNDFC1=0.0
         IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) GO TO 315
         IF(BEDLAK(II).LE.0.0) GO TO 315
         CNDFC1 = BEDLAK(II)*DELR(I)*DELC(J)
   !-LFK
         VK=VKAH(I,J,K)
   !-LFK   VK=VKA(I,J,K)
   !   skip if zero vk
   !-LFK        IF(VK.LE.0.0) GO TO 350
         IF(VK.LE.0.0) THEN
            CNDFCT(II) = CNDFC1
         ELSE
   !          BBOT=BOTM(I,J,LBOTM(K))
   !          TTOP=BOTM(I,J,LBOTM(K)-1)
   !          CAQ=VK*DELR(I)*DELC(J)/((TTOP-BBOT)*0.5)
   !          CNDFCT(II) = 1.0/(1.0/CAQ+1.0/CNDFC1)
   !-lfk    When HUF is active, set effective lakebed conductance only on basis of user-specified lakebed leakance value
            CNDFCT(II) = CNDFC1
         END IF
   !-LFK   (ABOVE BLOCK MODIFIED 1/2/2013 WITH IF-THEN-ELSE-ENDIF CONTROLS TO AVOID ZERO LKNCE)
315      IF ( LWRT <= 0 ) WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),&
   !-lfk
         &DELR(I),BEDLAK(II),CNDFC1,CNDFCT(II)
   !    1             BEDLAK(II),CNDFC1,CAQ,CNDFCT(II)
      ELSE
   !
   !  Horizontal conductance
   !
   !        TT = HK(I,J,K)
   !-LFK
   !        TY = HKCC(I,J,K)
   !gsf    TY = HKCC(I,J,K)
   !-LFK        TY = 0.0     !gsf
   ! X-DIRECTION
   !        IF(NTYP.EQ.2) CNDFC2 = 2.0*TT*DELC(J)/DELR(I)
   ! Y-DIRECTION
   !        IF(NTYP.EQ.3) CNDFC2 = 2.0*TY*DELC(J)/DELR(I)
         IF(NTYP.EQ.2) CNDFC1 = BEDLAK(II)*DELC(J)
         IF(NTYP.EQ.3) CNDFC1 = BEDLAK(II)*DELR(I)
   !        IF (CNDFC1.GT.0.0.AND.CNDFC2.GT.0.0)
   !     *         CNDFCT(II) = 1.0/(1.0/CNDFC2+1.0/CNDFC1)
   !-lfk    When HUF is active, set effective lakebed conductance only on basis of user-specified lakebed leakance value
         CNDFCT(II) = CNDFC1
         IF ( LWRT <= 0 )WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),&
   !-lfk
         &DELR(I),BEDLAK(II),CNDFC1,CNDFCT(II)
   !    1    BEDLAK(II),CNDFC1,CNDFC2,CNDFCT(II)
   !-lfk
7325     FORMAT(1X,5I3,2X,1P,4E10.2,10x,E11.3)
   ! 7325   FORMAT(1X,5I3,2X,1P,6E10.2)
      END IF
350 CONTINUE
   !
   !--LFK   INFO FOR HUF-LAKE COMBO:
   !  WRITE WARNINGS ON LAKE/AQUIFER CONDUCTANCES
   WRITE(IOUT,345)
345 FORMAT(//5X,'NOTE: INFORMATION ABOUT CALCULATED LAKE/AQUIFER C&
   &ONDUCTANCES WHEN USING HUF PACKAGE FOLLOWS: '/)
   WRITE(IOUT,346)
346 FORMAT(1X,'FOR ALL NODE(S) ADJACENT TO LAKE: '/&
   &1X,'COMBINED LAKE/AQUIFER CONDUCTANCES BASED SOLELY ON LAKEBED&
   & LEAKANCE SPECIFICATION'/)
   WRITE(IOUT,'(//)')
   !
   RETURN
END
   !dep  Added function statements to compute derivatives for Newton method
   !dep     used in solving lake stage in the FORMULATE SUBROUTINE (LAK7FM).
DOUBLE PRECISION FUNCTION FINTERP (STAGE,LN)
   !dep&rgn  FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
   !          OF LAKE STAGE TO CACULATE LAKE AREA.
   !         ADDED 5/16/2006-- changed 12/2007 from "DOUBLE PRECISION FUNCTION"
   !          to "FUNCTION"
   USE GWFLAKMODULE, ONLY: AREATABLE, DEPTHTABLE
   IMPLICIT NONE
   DOUBLE PRECISION STAGE, AREA, TOLF2, FOLD
   DOUBLE PRECISION a1, a2, d1, d2, SLOPE
   INTEGER LN, I
   TOLF2=1.0E-7
   AREA = 0.0D0
   FINTERP = 0.0D0
   IF (STAGE.GT.DEPTHTABLE(151,LN))THEN
      a1 = AREATABLE(150,LN)
      a2 = AREATABLE(151,LN)
      d1 = DEPTHTABLE(150,LN)
      d2 = DEPTHTABLE(151,LN)
      SLOPE = (a2-a1)/(d2-d1)
      FINTERP =  SLOPE*STAGE+a2-SLOPE*d2
      RETURN
   END IF
   I = 1
   DO
      a1 = AREATABLE(I,LN)
      a2 = AREATABLE(I+1,LN)
      d1 = DEPTHTABLE(I,LN)
      d2 = DEPTHTABLE(I+1,LN)
      SLOPE = (a2-a1)/(d2-d1)
      FOLD= STAGE-d1
      IF (FOLD .LE. 0.0D0) THEN
         AREA=A1
         EXIT
      ELSEIF (STAGE.GE.d1 .AND. STAGE.LE.d2)THEN
         AREA=SLOPE*STAGE+a2-SLOPE*d2
         EXIT
      END IF
      I = I + 1
      IF( I.GT.150 ) THEN
         AREA=SLOPE*STAGE+a2-SLOPE*d2
         EXIT
      END IF
   END DO
   FINTERP = AREA
   RETURN
END FUNCTION FINTERP
   !  RGN Added function statements to compute calculate surface area form volume
DOUBLE PRECISION FUNCTION SURFTERP (VOL,LN)
   !     FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
   !          OF LAKE VOLUME TO CACULATE LAKE AREA.
   USE GWFLAKMODULE, ONLY: AREATABLE, VOLUMETABLE
   IMPLICIT NONE
   DOUBLE PRECISION a1, a2, V1, V2, SLOPE, TOLF2, FOLD, AREA
   INTEGER LN, I
   DOUBLE PRECISION VOL
   TOLF2=1.0E-7
   AREA = 0.0D0
   SURFTERP = 0.0D0
   IF (VOL.GT.VOLUMETABLE(151,LN))THEN
      a1 = AREATABLE(150,LN)
      a2 = AREATABLE(151,LN)
      V1 = VOLUMETABLE(150,LN)
      V2 = VOLUMETABLE(150,LN)
      SLOPE = (a2-a1)/(V2-V1)
      SURFTERP =  SLOPE*VOL+a2-SLOPE*V2
      RETURN
   END IF
   I = 1
   DO
      FOLD=ABS(VOL-VOLUMETABLE(I,LN))
      a1 = AREATABLE(I,LN)
      a2 = AREATABLE(I+1,LN)
      V1 = VOLUMETABLE(I,LN)
      V2 = VOLUMETABLE(I+1,LN)
      SLOPE = (a2-a1)/(V2-V1)
      FOLD = VOL-V1
      IF (FOLD .LE. 0.0D0) THEN
         AREA=A1
         EXIT
      ELSEIF (VOL.GE.V1 .AND. VOL.LE.V2)THEN
         AREA=SLOPE*VOL+a2-SLOPE*V2
         EXIT
      END IF
      I = I + 1
      IF( I.GT.150 ) THEN
         AREA=SLOPE*VOL+a2-SLOPE*V2
         EXIT
      END IF
   END DO
   SURFTERP = AREA
   RETURN
END FUNCTION SURFTERP
   !
   !     Interpolate lake volume as a function of lake stage
   !     used in solving lake stage in the FORMULATE SUBROUTINE (LAK7FM).
DOUBLE PRECISION FUNCTION VOLTERP (STAGE,LN)
   !     FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
   !          OF LAKE STAGE TO CACULATE LAKE VOLUME.
   USE GWFLAKMODULE, ONLY: VOLUMETABLE, DEPTHTABLE
   IMPLICIT NONE
   INTEGER LN, I
   DOUBLE PRECISION STAGE, VOL, TOLF2, FOLD
   DOUBLE PRECISION D1, D2, V1, V2, SLOPE
   TOLF2=1.0E-7
   VOLTERP = 0.0D0
   VOL = 0.0D0
   IF (STAGE.GT.DEPTHTABLE(151,LN))THEN
      V1 = VOLUMETABLE(150,LN)
      V2 = VOLUMETABLE(151,LN)
      D1 = DEPTHTABLE(150,LN)
      D2 = DEPTHTABLE(151,LN)
      SLOPE = (V2-V1)/(D2-D1)
      VOLTERP = SLOPE*STAGE+V2-SLOPE*D2
      RETURN
   END IF
   I = 1
   DO
      V1 = VOLUMETABLE(I,LN)
      V2 = VOLUMETABLE(I+1,LN)
      D1 = DEPTHTABLE(I,LN)
      D2 = DEPTHTABLE(I+1,LN)
      SLOPE = (V2-V1)/(D2-D1)
      FOLD = STAGE-D1
      IF (FOLD .LE. 0.0d0) THEN
         VOL=V1
         EXIT
      ELSEIF (STAGE.GE.D1 .AND. STAGE.LE.D2)THEN
         VOL = SLOPE*STAGE+V2-SLOPE*D2
         EXIT
      END IF
      I = I + 1
      IF( I.GT.150 ) THEN
         VOL = SLOPE*STAGE+V2-SLOPE*D2
         EXIT
      END IF
   END DO
   VOLTERP = VOL
   RETURN
END FUNCTION VOLTERP
   !     Interpolate lake STAGE as a function of lake VOLUME
   !     used in solving lake stage in the FORMULATE SUBROUTINE (LAK7FM).
DOUBLE PRECISION FUNCTION STGTERP (VOL,LN)
   !     FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
   !          OF LAKE VOLUME TO CACULATE LAKE STAGE.
   USE GWFLAKMODULE, ONLY: VOLUMETABLE, DEPTHTABLE
   IMPLICIT NONE
   DOUBLE PRECISION VOL, TOLF2, FOLD
   DOUBLE PRECISION D1, D2, V1, V2, SLOPE
   INTEGER LN, I
   !!      DOUBLE PRECISION VOLUME, STAGE
   TOLF2=1.0E-7
   STGTERP = 0.0D0
   IF (VOL.GT.VOLUMETABLE(151,LN))THEN
      D1 = DEPTHTABLE(150,LN)
      D2 = DEPTHTABLE(151,LN)
      V1 = VOLUMETABLE(150,LN)
      V2 = VOLUMETABLE(151,LN)
      SLOPE = (D2-D1)/(V2-V1)
      STGTERP =  SLOPE*VOL+D2-SLOPE*V2
      RETURN
   END IF
   I = 1
   DO
      D1 = DEPTHTABLE(I,LN)
      D2 = DEPTHTABLE(I+1,LN)
      V1 = VOLUMETABLE(I,LN)
      V2 = VOLUMETABLE(I+1,LN)
      SLOPE = (D2-D1)/(V2-V1)
      FOLD = VOL-V1
      IF (FOLD .LE. 0.0d0) THEN
         STGTERP=D1
         EXIT
      ELSEIF (VOL.GE.V1 .AND. VOL.LE.V2)THEN
         STGTERP =  SLOPE*VOL+D2-SLOPE*V2
         EXIT
      END IF
      I = I + 1
      IF( I.GT.150 ) THEN
         STGTERP =  SLOPE*VOL+D2-SLOPE*V2
         EXIT
      END IF
   END DO
   RETURN
END FUNCTION STGTERP
   !------FUNCTION DERIVTERP FOR INTERPOLATING DERIVATIVE OF LAKE OUTFLOW.
DOUBLE PRECISION FUNCTION DERIVTERP (STAGE,LSEG)
   !dep&rgn  FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
   !          OF LAKE STAGE TO CACULATE LAKE OUTFLOW DERIVATIVE.
   USE GWFSFRMODULE, ONLY: DLKOTFLW, DLKSTAGE
   IMPLICIT NONE
   DOUBLE PRECISION STAGE, FOLD, TOLF2
   DOUBLE PRECISION DS1, DS2, DF1, DF2, SLOPE
   INTEGER LSEG, I
   TOLF2=1.0E-7
   DERIVTERP = 0.0D0
   IF (STAGE.GT.DLKSTAGE(200,LSEG))THEN
      DF1 = DLKOTFLW(199,LSEG)
      DF2 = DLKOTFLW(200,LSEG)
      DS1 = DLKSTAGE(199,LSEG)
      DS2 = DLKSTAGE(200,LSEG)
      SLOPE = (DF2-DF1)/(DS2-DS1)
      DERIVTERP =  SLOPE*STAGE+DF2-SLOPE*DS2
      RETURN
   END IF
   I = 1
   DO
      DF1 = DLKOTFLW(I,LSEG)
      DF2 = DLKOTFLW(I+1,LSEG)
      DS1 = DLKSTAGE(I,LSEG)
      DS2 = DLKSTAGE(I+1,LSEG)
      SLOPE = (DF2-DF1)/(DS2-DS1)
      FOLD= STAGE-DS1
      IF (FOLD .LE. 0.0D0) THEN
         DERIVTERP=0.0D0
         EXIT
      ELSEIF (STAGE.GE.DS1 .AND. STAGE.LE.DS2)THEN
         DERIVTERP=SLOPE*STAGE+DF2-SLOPE*DS2
         EXIT
      END IF
      I = I + 1
      IF( I.GT.199) THEN
         DERIVTERP=SLOPE*STAGE+DF2-SLOPE*DS2
         EXIT
      END IF
   END DO
   RETURN
END FUNCTION DERIVTERP
   !------FUNCTION OUTFLWTERP FOR INTERPOLATING DERIVATIVE OF LAKE OUTFLOW.
DOUBLE PRECISION FUNCTION OUTFLWTERP (STAGE,LSEG)
   !dep&rgn  FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
   !          OF LAKE OUTFLOW STORED IN SLKOTFLW ARRAY.
   !         ADDED 5/16/2006-- changed 12/2007 from "DOUBLE PRECISION FUNCTION"
   !          to "FUNCTION"
   USE GWFSFRMODULE, ONLY: SLKOTFLW, DLKSTAGE
   IMPLICIT NONE
   DOUBLE PRECISION STAGE, FOLD, TOLF2
   DOUBLE PRECISION SL1, SL2, DL1, DL2, SLOPE
   INTEGER LSEG, I
   TOLF2=1.0E-9
   OUTFLWTERP = 0.0D0
   IF (STAGE.GT.DLKSTAGE(200,LSEG))THEN
      SL1 = SLKOTFLW(199,LSEG)
      SL2 = SLKOTFLW(200,LSEG)
      DL1 = DLKSTAGE(199,LSEG)
      DL2 = DLKSTAGE(200,LSEG)
      SLOPE = (SL2-SL1)/(DL2-DL1)
      OUTFLWTERP =  SLOPE*STAGE+SL2-SLOPE*DL2
      RETURN
   END IF
   I = 1
   DO
      SL1 = SLKOTFLW(I,LSEG)
      SL2 = SLKOTFLW(I+1,LSEG)
      DL1 = DLKSTAGE(I,LSEG)
      DL2 = DLKSTAGE(I+1,LSEG)
      SLOPE = (SL2-SL1)/(DL2-DL1)
      FOLD= STAGE-DL1
      IF (FOLD .LE. 0.0d0) THEN
         OUTFLWTERP=0.0d0
         EXIT
      ELSEIF (STAGE.GE.DL1 .AND. STAGE.LE.DL2)THEN
         OUTFLWTERP = SLOPE*STAGE+SL2-SLOPE*DL2
         EXIT
      END IF
      I = I + 1
      IF( I.GT.199) THEN
         OUTFLWTERP = SLOPE*STAGE+SL2-SLOPE*DL2
         EXIT
      END IF
   END DO
   RETURN
END FUNCTION OUTFLWTERP
   !
   !------FUNCTION FXLKOT_TERP FOR SMOOTHING SPECIFIED LAKE OUTFLOWS TO STREAMS.
   !
DOUBLE PRECISION FUNCTION FXLKOT_TERP(RAMP,DSTAGE,Botlake,&
&Splakout,dy)
   IMPLICIT NONE
   DOUBLE PRECISION DSTAGE,Botlake,Splakout, s, aa, ad, b, x, y, dy
   DOUBLE PRECISION RAMP
   FXLKOT_TERP = 0.0D0
   s = RAMP
   x = DSTAGE-Botlake
   aa = -1.0d0/(s**2.0d0)
   ad = -2.0D0/(s**2.0d0)
   b = 2.0d0/s
   y = aa*x**2.0d0 + b*x
   dy = (ad*x + b)*Splakout
   IF ( x.LE.0.0 ) THEN
      y = 0.0D0
      dy = 0.0D0
   ELSE IF ( x-s.GT.-1.0e-14 ) THEN
      y = 1.0D0
      dy = 0.0D0
   END IF
   FXLKOT_TERP = y*Splakout
END FUNCTION FXLKOT_TERP
   !
SUBROUTINE GET_FLOBOT(IC, IR, IL, ITYPE, INOFLO,CONDUC,&
&FLOBOT,FLOBO3,FLOTOUZF,DLSTG,CLOSEZERO,H,&
&THET1,ISS,LAKE,II,SURFDPTH,AREA,IUNITUZF,&
&BOTLK,BOTCL,L1)
   !
   !     ******************************************************************
   !     CALCULATE SEEPAGE BETWEEN LAKE AND GW CELLS
   !     ******************************************************************
   !
   USE GWFLAKMODULE
   USE GLOBAL,       ONLY: IBOUND, LBOTM, BOTM, LAYHDT
   !!      USE GLOBAL,       ONLY: IBOUND, IOUT, LBOTM, BOTM, NLAY,LAYHDT
   USE GWFUZFMODULE, ONLY: IUZFBND,FINF,VKS
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   !     FUNCTIONS
   !     -----------------------------------------------------------------
   !     -----------------------------------------------------------------
   !     ARGUMENTS
   DOUBLE PRECISION FLOBO3,FLOBOT,CONDUC,H,THET1,CLOSEZERO,DLSTG,&
   &SURFDPTH,AREA,BOTLK,BOTCL,HH
   INTEGER ISS, LAKE, II, IC, IR, IL, ITYPE, IUNITUZF, L1
   !     -----------------------------------------------------------------
   INTEGER INOFLO
   !!      INTEGER ICHECK, LI, INOFLO
   DOUBLE PRECISION FLOBO1,FLOBO2,CONDMX,BOTLKUP,&
   &BOTLKDN,FLOTOUZF,RAMPGW,RAMPSTGO,RAMPSTGN,&
   &HTEMP,HD,THCK,RAMPUP
   !!     2                 RAMPSTGON,HTEMP,HD,THCK,RAMPUP
   !
   !5C-----INITIALIZE GROUNDWATER SEEPAGE VARIABLES AND CONDUCTANCE FACTOR.
   FLOBO1 = 0.0D0
   FLOBO2 = 0.0D0
   !
   !6------COMPUTE SEEPAGE INTO OR OUT OF A LAKE BED NODE WHEN ITYPE=0.
   !       HEAD CANNOT FALL BELOW LAKE BOTTOM
   IF (ITYPE.EQ.0) THEN
   !
   !6B------RAMP CONDUCTANCE ACROSS HORIZONTAL CELL FACE WHEN
   !          LAKE STAGE AND GROUNDWATER HEAD NEAR LAKEBED.
      BOTLKUP = BOTLK + SURFDPTH
      BOTLKDN = BOTLK
      CONDMX = CONDUC
      HH = H
      IF ( HH.LT.BOTLKDN ) THEN
         HH = BOTLKDN
         INOFLO = 1
      END IF
      IF(SURFDPTH.GT.CLOSEZERO) THEN
         RAMPGW = CONDMX-(CONDMX/SURFDPTH)*&
         &(BOTLKUP-HH)
         IF ( RAMPGW-CONDMX.GT.0.0D0 ) RAMPGW = CONDMX
         IF ( RAMPGW.LE.0.0D0 ) RAMPGW = 0.0D0
         RAMPSTGO = CONDMX-(CONDMX/SURFDPTH)*&
         &(BOTLKUP-STGOLD(LAKE))
         IF ( RAMPSTGO-CONDMX.GT.0.0D0 ) RAMPSTGO = CONDMX
         IF ( RAMPSTGO.LE.0.0D0 ) RAMPSTGO = 0.0D0
         RAMPSTGN = CONDMX-(CONDMX/SURFDPTH)*&
         &(BOTLKUP-STGNEW(LAKE))
         IF ( RAMPSTGN-CONDMX.GT.0.0D0 ) RAMPSTGN = CONDMX
         IF ( RAMPSTGN.LE.0.0D0 ) RAMPSTGN = 0.0D0
      ELSE
         RAMPGW=CONDMX
         RAMPSTGO=CONDMX
         RAMPSTGN=CONDMX
      END IF
      IF( HH-BOTLKDN.GT.CLOSEZERO ) THEN
         HTEMP = HH
      ELSE
         HTEMP=BOTLKDN
      END IF
   !
   !6C------COMPUTE LAKE SEEPAGE FOR STGOLD USING FLOBO1.
   !        USE UPSTREAM WEIGHTING
      IF ( HH.LT.STGOLD(LAKE) ) THEN
         RAMPUP = RAMPSTGO
      ELSE
         RAMPUP = RAMPGW
      END IF
      CONDUC = RAMPUP
      IF( STGOLD(LAKE)-BOTLKDN.GT.CLOSEZERO ) THEN
         FLOBO1=CONDUC*(STGOLD(LAKE)-HTEMP)
      ELSE
         FLOBO1=CONDUC*(BOTLKDN-HTEMP)
      END IF
      IF ( IUNITUZF.GT.0 ) THEN
         IF ( IUZFBND(IC,IR).GT.0 )THEN
            IF (H-BOTLK.LT.-0.5*SURFDPTH) THEN           !11/3/14 changed to H for UZF
               IF ( VKS(IC,IR)*AREA-FLOBO1.LT.CLOSEZERO )&
               &THEN
                  FLOBO1 = VKS(IC,IR)*AREA
               END IF
            END IF
         END IF
      END IF
   !
   !6D------COMPUTE LAKE SEEPAGE FOR STGNEW USING FLOBO2 AND FLOBO3.
   !        USE UPSTREAM WEIGHTING
      IF ( HH.LT.STGNEW(LAKE) ) THEN
         RAMPUP = RAMPSTGN
      ELSE
         RAMPUP =  RAMPGW
      END IF
      CONDUC = RAMPUP
      IF( STGNEW(LAKE)-BOTLKDN.GT.CLOSEZERO ) THEN
         FLOBO2 = CONDUC*(STGNEW(LAKE)-HTEMP)
         FLOBO3 = CONDUC*(STGNEW(LAKE)+DLSTG-HTEMP)
      ELSE
         FLOBO2 = CONDUC*(BOTLKDN-HTEMP)
         FLOBO3 = CONDUC*(BOTLKDN+DLSTG-HTEMP)
      END IF
      IF ( IUNITUZF.GT.0 ) THEN
         IF ( IUZFBND(IC,IR).GT.0 )THEN
            IF ( H-BOTLK.LT.-0.5*SURFDPTH ) THEN
               IF ( VKS(IC,IR)*AREA-FLOBO2.LT.CLOSEZERO )&
               &THEN
                  FLOBO2 = VKS(IC,IR)*AREA
                  FLOBO3 = VKS(IC,IR)*AREA
               END IF
            END IF
         END IF
      END IF
   !
   !6E------COMPUTE LAKE SEEPAGE (FLOBOT) AS A FRACTION OF FLOBO1 AND
   !          FLOB02 AND FLOBO3 AS A FRACTION OF FLOBO1 AND FLOBO3.
      FLOBOT = THET1*FLOBO2 + (1.0D0-THET1)*FLOBO1
      FLOBO3 = THET1*FLOBO3 + (1.0D0-THET1)*FLOBO1
   !          CONDUC = THET1*RAMPSTGN + (1.0D0-THET1)*RAMPSTGO
      IF ( IUNITUZF.GT.0 ) THEN
         IF ( IUZFBND(IC,IR).GT.0 )THEN
            IF ( H-BOTLK.LT.-0.5*SURFDPTH ) THEN  !11/3/14 changed to H for UZF, was HH
               IF ( FLOBOT/AREA.GT.VKS(IC,IR) ) THEN
                  FLOBOT = VKS(IC,IR)*AREA
                  FLOBO3 = FLOTOUZF
               END IF
               FLOTOUZF = FLOBOT
               FLOBOT = 0.0D0
               IF ( abs((STGNEW(LAKE)-BOTLK)) > closezero ) THEN
                  CONDUC = FLOTOUZF/(STGNEW(LAKE)-BOTLK)
               ELSE
                  CONDUC = 0.0
               END IF
               FINF(IC,IR)=FLOTOUZF/AREA
            END IF
         END IF
      END IF
   !
   !7------COMPUTE SEEPAGE INTO OR OUT OF A LAKE WALL NODE
   !         WHEN ITYPE=1 OR 2.
   ELSE IF ( ITYPE.EQ.1.OR.ITYPE.EQ.2 ) THEN
      IF( IBOUND(IC,IR,IL).GT.0 ) THEN
         HD = H
         IF( H.GT.BOTM(IC,IR,LBOTM(IL)-1) )&
         &HD = BOTM(IC,IR,LBOTM(IL)-1)
   !
   !7B------CONDUCTANCE ACROSS VERTICAL CELL FACE DEPENDENT ON
   !          SATURATED THICKNESS.
         IF ( LAYHDT(il).GT.0 ) THEN
            THCK = HD - BOTCL
         ELSE
            THCK = BOTM(IC,IR,LBOTM(IL)-1) - BOTCL
         END IF
         IF( THCK.LE.0.0 ) THCK = 0.0
         CONDUC = CONDUC*THCK
         IF ( H.LT.BOTM(IC,IR,LBOTM(IL)) )&
         &H = BOTM(IC,IR,LBOTM(IL))
   !
   !7C------COMPUTE LAKE SEEPAGE FOR STGOLD USING FLOBO1.
         IF( STGOLD(LAKE)-BOTCL.GT.CLOSEZERO ) THEN
            FLOBO1 = CONDUC*(STGOLD(LAKE)-H)
         ELSE IF ( H-BOTCL.GT.CLOSEZERO ) THEN
            FLOBO1 = CONDUC*(BOTCL-H)
         END IF
   !
   !7D------COMPUTE LAKE SEEPAGE FOR STGNEW USING FLOBO2 AND FLOBO3.
         IF( STGNEW(LAKE)-BOTCL.GT.CLOSEZERO )THEN
            FLOBO3 = CONDUC*(STGNEW(LAKE)+DLSTG-H)
            FLOBO2 = CONDUC*(STGNEW(LAKE)-H)
         ELSE IF ( H-BOTCL.GT.CLOSEZERO ) THEN
            FLOBO3 = CONDUC*(BOTCL+DLSTG-H)
            FLOBO2 = CONDUC*(BOTCL-H)
         ELSE IF ( STGNEW(LAKE)+DLSTG.GE.BOTCL )THEN
            FLOBO3 = CONDUC*(STGNEW(LAKE)+DLSTG-H)
         END IF
   !
   !7E------COMPUTE LAKE SEEPAGE (FLOBOT) AS A FRACTION OF FLOBO1 AND
   !         FLOB02 AND FLOBO3 AS A FRACTION OF FLOBO1 AND FLOBO3.
         FLOBOT = THET1*FLOBO2 + (1.0D0-THET1)*FLOBO1
         FLOBO3  = THET1*FLOBO3 + (1.0D0-THET1)*FLOBO1
         SUMCNN(LAKE) = SUMCNN(LAKE) + CONDUC
      END IF
   END IF
   !
   !8-------SEEPAGE RATES ADDED TO MATRIX AND RESIDUAL TERMS.
   !8B------COMPUTE FLWITER AND FLWITER3 DURING FIRST LOOP THROUGH
   !          CALCULATIONS. NEGATIVE FLOBOT MEANS INTO LAKE
   IF ( II==1 ) THEN
      IF ( FLOBOT.LT.0.0D0 ) FLWITER(LAKE) =&
      &FLWITER(LAKE) - FLOBOT
      IF ( FLOBO3.LT.0.0D0 ) FLWITER3(LAKE) =&
      &FLWITER3(LAKE) - FLOBO3
   END IF
   !8C------COMPUTE FLWITER AND FLOWITER3 DURING SECOND LOOP THROUGH
   !          CALCULATIONS.
   IF ( II==2 ) THEN
      IF ( FLOBOT>=FLWITER(LAKE) ) THEN
         IF ( FLOBOT.GT.CLOSEZERO ) THEN
   !              FLOBO2=FLWITER(LAKE)
   !              FLOBOT = THET1*FLOBO2 + (1.0D0-THET1)*FLOBO1
            FLOBOT = FLWITER(LAKE)
            FLWITER(LAKE) = 0.0
            INOFLO = 1
         END IF
      ELSE IF ( FLOBOT.GT.CLOSEZERO )THEN
         FLWITER(LAKE) = FLWITER(LAKE) - FLOBOT
      END IF
      IF ( FLOTOUZF>=FLWITER(LAKE) ) THEN
         IF ( FLOTOUZF.GT.CLOSEZERO ) THEN
            FLOTOUZF=FLWITER(LAKE)
            !             FLOTOUZF = THET1*FLOTOUZF + (1.0D0-THET1)*FLOBO1
            FLWITER(LAKE) = 0.0
            INOFLO = 1
         END IF
      ELSE IF ( FLOTOUZF.GT.CLOSEZERO )THEN
         FLWITER(LAKE) = FLWITER(LAKE) - FLOTOUZF
      END IF
      IF ( FLOBO3>=FLWITER3(LAKE) ) THEN
         IF ( FLOBO3.GT.CLOSEZERO ) THEN
            FLOBO3=FLWITER3(LAKE)
            !             FLOBO3  = THET1*FLOBO3 + (1.0D0-THET1)*FLOBO1
            FLWITER3(LAKE) = 0.0
            INOFLO = 1
         END IF
      ELSE IF ( FLOBO3.GT.CLOSEZERO )THEN
         FLWITER3(LAKE) = FLWITER3(LAKE) - FLOBO3
      END IF
   END IF
   !
   !6E------COMPUTE LAKE SEEPAGE (FLOBOT) AS A FRACTION OF FLOBO1 AND
   !          FLOBO2 AND FLOBO3 AS A FRACTION OF FLOBO1 AND FLOBO3.
   RETURN
END SUBROUTINE GET_FLOBOT
   !
   !-------SUBROUTINE LAK2MODSIM
SUBROUTINE LAK2MODSIM(DELTAVOL, LAKEVOL, Diversions, Nsegshold)
   !     *******************************************************************
   !     SET VOLUMES, SFR INFLOWS, AND SFR OUTFLOWS FOR MODSIM
   !--------MARCH 8, 2017
   !     *******************************************************************
   USE GWFLAKMODULE, ONLY: NLAKES, SURFIN, SURFOT, VOLOLDD, VOL,&
   &DEADPOOLVOL, RELEASABLE_STOR, IDIV,&
   &MXLKVOLF
   USE GWFSFRMODULE, ONLY: IDIVAR, FXLKOT
   USE GWFBASMODULE, ONLY: DELT
   IMPLICIT NONE
   !     -------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     -------------------------------------------------------------------
   !     ARGUMENTS
   INTEGER,          INTENT(IN)    :: Nsegshold
   DOUBLE PRECISION, INTENT(INOUT) :: DELTAVOL(NLAKES),&
   &LAKEVOL(NLAKES),&
   &Diversions(Nsegshold)
   !     -------------------------------------------------------------------
   !      INTEGER
   !      DOUBLE PRECISION
   !     -------------------------------------------------------------------
   !     LOCAL VARIABLES
   !     -------------------------------------------------------------------
   INTEGER LAKE, M, LAK_ID, INODE
   DOUBLE PRECISION DELTAQ, start_vol, inflow, outflow, gwsw
   !     -------------------------------------------------------------------
   !
   !
   !1-------SET FLOWS IN AND OUT OF LAKES AND CHANGE IN LAKE VOLUME.
   !
   !    OPEN(222, FILE='LAKE_DEBUG_.TXT')
   DO LAKE=1,NLAKES
      DELTAQ = (SURFIN(LAKE) - SURFOT(LAKE))*DELT
      DELTAVOL(LAKE) = VOL(LAKE) - VOLOLDD(LAKE) - DELTAQ
      LAKEVOL(LAKE) = VOL(LAKE)
      INODE=IDIV(LAKE,1)
      IF(INODE.EQ.0) INODE=1
      !      WRITE(222,333)LAKE,PRECIP(LAKE),EVAP(LAKE),
      !   1     RUNF(LAKE),RUNOFF(LAKE),WITHDRW(LAKE),SURFIN(LAKE),
      !   2     SURFOT(LAKE),SEEP(LAKE),VOL(LAKE),VOLOLDD(LAKE),
      !   3     SEG(8,INODE),BOTTMS(LAKE),STGNEW(LAKE),Diversions(505)    ! SEG(8,)=BotLake
   END DO
   !333 FORMAT(I5,14E15.7)
   !
   !-----UPDATE RELEASABLE STORAGE ARRAY RETURNED TO MODSIM FOR RESETTING POTENTIAL RELEASE AMOUNT
   DO LAKE=1,NLAKES
      RELEASABLE_STOR(LAKE) = VOL(LAKE) - DEADPOOLVOL(LAKE)
      IF (RELEASABLE_STOR(LAKE).LT.0.0) RELEASABLE_STOR(LAKE) = 0.0
   ENDDO
   !
   !-----LOOP OVER REACHES AND OVERRIDE LAKE RELEASES IF WATER LIMITED
   IF (Nsegshold.GT.-1) THEN
      DO M = 1, Nsegshold
         IF (IDIVAR(1,M).LT.0) THEN
            LAK_ID = ABS(IDIVAR(1,M))
            ! Set some variables for readability of ELSEIF statement
            start_vol = VOLOLDD(LAK_ID)
            inflow = SURFIN(LAK_ID)
            outflow = SURFOT(LAK_ID)
            gwsw = DELTAVOL(LAK_ID)
   !
            ! The following bit of code added to handle stress period 54
            IF (VOL(LAK_ID).GT.(MXLKVOLF(LAK_ID)+1.0).AND.&
            &.NOT.MXLKVOLF(LAK_ID).LT.0.0) THEN
               !Diversions(M) = (VOL(LAK_ID) - MXLKVOLF(LAK_ID)) / DELT
            ELSEIF((start_vol + inflow + gwsw).LT.&
            &DEADPOOLVOL(LAK_ID)) THEN
               Diversions(M) = 0.0001
               !INODE = IDIV(LAKE,IDV)
   !           THE FOLLOWING CONDITION WILL BE TRIGGERED WHEN NEAR DEADPOOL STORAGE
   !           CHECKS WHAT'S AVAILABLE AGAINST WHAT MODSIM IS ASKING FOR ("Diversions")
            ELSEIF((start_vol + inflow - outflow + gwsw + 1.0).LT.&
            &DEADPOOLVOL(LAK_ID)) THEN
               Diversions(M)=MAX((RELEASABLE_STOR(LAK_ID)/DELT),FXLKOT(M))
               !ELSEIF((RELEASABLE_STOR(LAK_ID)/DELT).LT.Diversions(M)) THEN
               !Diversions(M)=RELEASABLE_STOR(LAK_ID) / DELT
               !Diversions(M)=MAX((RELEASABLE_STOR(LAK_ID)/DELT),FXLKOT(M))
               ! NEED TODO: Look into "/ DELT", could be an issue in Deschutes where seconds are used.
               ! Need to check this calculation by hand to make sure units are as expected for MODSIM
               !          Diversions(M) = (RELEASABLE_STOR(LAK_ID) / DELT) +
               !  &                       FXLKOT(M)
            ENDIF
         ENDIF
      ENDDO
   ENDIF
   !
   !5------RETURN.
   RETURN
END SUBROUTINE LAK2MODSIM
   !
   !-------SUBROUTINE LAK2MODSIM, But directly callable by MODSIM
SUBROUTINE LAK2MODSIM_InitLakes(DELTAVOL,LAKEVOL, MXLKVOL)&
&BIND(C,NAME="LAK2MODSIM_InitLakes")

   !DEC$ ATTRIBUTES DLLEXPORT :: LAK2MODSIM_InitLakes

   !     *******************************************************************
   !     SET VOLUMES, SFR INFLOWS, AND SFR OUTFLOWS FOR MODSIM
   !--------MARCH 8, 2017
   !     *******************************************************************
   USE GWFLAKMODULE, ONLY: NLAKES, DEADPOOLVOL, MXLKVOLF
   IMPLICIT NONE
   !     -------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     -------------------------------------------------------------------
   !     ARGUMENTS
   DOUBLE PRECISION, INTENT(INOUT) :: DELTAVOL(NLAKES),&
   &LAKEVOL(NLAKES),&
   &MXLKVOL(NLAKES)
   !      INTEGER, INTENT(IN) :: KITER,KSTP,KPER
   !     -------------------------------------------------------------------
   !      INTEGER
   !      DOUBLE PRECISION
   !     -------------------------------------------------------------------
   !     LOCAL VARIABLES
   !     -------------------------------------------------------------------
   INTEGER LAKE
   DOUBLE PRECISION :: dmy(1)
   !     -------------------------------------------------------------------
   !
   !0----FILL A NEW VARIABLE CALLED MXLKVOLF CONTAINING THE MODSIM MAX LAKE STORAGE
   DO LAKE=1, NLAKES
      MXLKVOLF(LAKE) = MXLKVOL(LAKE)
   ENDDO
   !
   !1-------SET FLOWS IN AND OUT OF LAKES AND CHANGE IN LAKE VOLUME.
   !
   dmy(1) = 0.0D0
   CALL LAK2MODSIM(DELTAVOL,LAKEVOL, dmy(1), -1) ! rsr, the 0 needs to be an array, values
   !
   !2------STUFF DELTAVOL WITH DEADPOOL INFORMATION CALCULATED BY MODFLOW
   DO LAKE=1, NLAKES
      DELTAVOL(LAKE) = DEADPOOLVOL(LAKE)
   ENDDO
   !
   !5------RETURN.
   RETURN
END SUBROUTINE LAK2MODSIM_InitLakes
   !
   !
SUBROUTINE GWF2LAK7DA(IUNITLAK, IGRID)
   !dep  End of FUNCTIONS used for Newton method in
   !dep     FORMULATE SUBROUTINE (LAK7FM).
   !  Deallocate LAK data
   USE GWFLAKMODULE
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     ------------------------------------------------------------------
   INTEGER, INTENT(IN) :: IUNITLAK, IGRID
   !
   DEALLOCATE (GWFLAKDAT(IGRID)%NLAKES)
   DEALLOCATE (GWFLAKDAT(IGRID)%NLAKESAR)
   DEALLOCATE (GWFLAKDAT(IGRID)%THETA)
   DEALLOCATE (GWFLAKDAT(IGRID)%STGNEW)
   DEALLOCATE (GWFLAKDAT(IGRID)%STGOLD)
   DEALLOCATE (GWFLAKDAT(IGRID)%STGOLD2)
   DEALLOCATE (GWFLAKDAT(IGRID)%STGITER)
   DEALLOCATE (GWFLAKDAT(IGRID)%VOL)
   DEALLOCATE (GWFLAKDAT(IGRID)%LAKUNIT)
   IF ( IUNITLAK.LT.1 ) RETURN

   DEALLOCATE (GWFLAKDAT(IGRID)%ILKCB)
   DEALLOCATE (GWFLAKDAT(IGRID)%LAKTAB)
   DEALLOCATE (GWFLAKDAT(IGRID)%IRDTAB)
   DEALLOCATE (GWFLAKDAT(IGRID)%ISTARTLAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%NSSITR)
   !dep  deallocate SURFDEPTH 3/3/2009
   DEALLOCATE (GWFLAKDAT(IGRID)%SURFDEPTH)
   DEALLOCATE (GWFLAKDAT(IGRID)%RAMP)
   DEALLOCATE (GWFLAKDAT(IGRID)%SMALLTOL)
   DEALLOCATE (GWFLAKDAT(IGRID)%DEADPOOLVOL)
   DEALLOCATE (GWFLAKDAT(IGRID)%RELEASABLE_STOR)
   DEALLOCATE (GWFLAKDAT(IGRID)%MXLKVOLF)
   DEALLOCATE (GWFLAKDAT(IGRID)%MXLKND)
   DEALLOCATE (GWFLAKDAT(IGRID)%LKNODE)
   DEALLOCATE (GWFLAKDAT(IGRID)%ICMX)
   DEALLOCATE (GWFLAKDAT(IGRID)%NCLS)
   DEALLOCATE (GWFLAKDAT(IGRID)%LWRT)
   DEALLOCATE (GWFLAKDAT(IGRID)%NDV)
   DEALLOCATE (GWFLAKDAT(IGRID)%NTRB)
   DEALLOCATE (GWFLAKDAT(IGRID)%SSCNCR)
   DEALLOCATE (GWFLAKDAT(IGRID)%ICS)
   DEALLOCATE (GWFLAKDAT(IGRID)%NCNCVR)
   DEALLOCATE (GWFLAKDAT(IGRID)%LIMERR)
   DEALLOCATE (GWFLAKDAT(IGRID)%ILAKE)
   DEALLOCATE (GWFLAKDAT(IGRID)%ITRB)
   DEALLOCATE (GWFLAKDAT(IGRID)%IDIV)
   DEALLOCATE (GWFLAKDAT(IGRID)%ISUB)
   DEALLOCATE (GWFLAKDAT(IGRID)%IRK)
   DEALLOCATE (GWFLAKDAT(IGRID)%LKARR1)
   DEALLOCATE (GWFLAKDAT(IGRID)%STAGES)
   DEALLOCATE (GWFLAKDAT(IGRID)%FLOB)
   DEALLOCATE (GWFLAKDAT(IGRID)%DSRFOT)
   DEALLOCATE (GWFLAKDAT(IGRID)%PRCPLK)
   DEALLOCATE (GWFLAKDAT(IGRID)%EVAPLK)
   DEALLOCATE (GWFLAKDAT(IGRID)%BEDLAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%WTHDRW)
   DEALLOCATE (GWFLAKDAT(IGRID)%RNF)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMRNF)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMPPT)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMEVP)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMGWI)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMGWO)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMSWI)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMSWO)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMWDR)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMFLX)
   DEALLOCATE (GWFLAKDAT(IGRID)%CNDFCT)
   DEALLOCATE (GWFLAKDAT(IGRID)%VOLINIT)
   DEALLOCATE (GWFLAKDAT(IGRID)%BOTTMS)
   DEALLOCATE (GWFLAKDAT(IGRID)%BGAREA)
   DEALLOCATE (GWFLAKDAT(IGRID)%SSMN)
   DEALLOCATE (GWFLAKDAT(IGRID)%SSMX)
   DEALLOCATE (GWFLAKDAT(IGRID)%EVAP)
   DEALLOCATE (GWFLAKDAT(IGRID)%PRECIP)
   DEALLOCATE (GWFLAKDAT(IGRID)%EVAP3)
   DEALLOCATE (GWFLAKDAT(IGRID)%PRECIP3)
   DEALLOCATE (GWFLAKDAT(IGRID)%SEEP)
   DEALLOCATE (GWFLAKDAT(IGRID)%SEEP3)
   DEALLOCATE (GWFLAKDAT(IGRID)%SURFA)
   DEALLOCATE (GWFLAKDAT(IGRID)%SURFIN)
   DEALLOCATE (GWFLAKDAT(IGRID)%SURFOT)
   DEALLOCATE (GWFLAKDAT(IGRID)%SUMCNN)
   DEALLOCATE (GWFLAKDAT(IGRID)%SUMCHN)
   DEALLOCATE (GWFLAKDAT(IGRID)%CLAKE)
   DEALLOCATE (GWFLAKDAT(IGRID)%CRNF)
   DEALLOCATE (GWFLAKDAT(IGRID)%SILLVT)
   DEALLOCATE (GWFLAKDAT(IGRID)%CAUG)
   DEALLOCATE (GWFLAKDAT(IGRID)%CPPT)
   DEALLOCATE (GWFLAKDAT(IGRID)%CLAKINIT)
   DEALLOCATE (GWFLAKDAT(IGRID)%BDLKN1)
   !dep  Added arrays that track lake budgets for dry lakes
   DEALLOCATE (GWFLAKDAT(Igrid)%WITHDRW)
   DEALLOCATE (GWFLAKDAT(Igrid)%FLWIN)
   DEALLOCATE (GWFLAKDAT(Igrid)%FLWITER)
   DEALLOCATE (GWFLAKDAT(Igrid)%FLWITER3)
   DEALLOCATE (GWFLAKDAT(Igrid)%GWRATELIM)
   !dep  Deallocate arrays used in conjunction with UZF Package
   DEALLOCATE (GWFLAKDAT(Igrid)%OVRLNDRNF)
   DEALLOCATE (GWFLAKDAT(Igrid)%CUMLNDRNF)
   DEALLOCATE (GWFLAKDAT(Igrid)%CUMUZF)
   !dep  Deallocate arrays for storing depth, and area arrays
   DEALLOCATE (GWFLAKDAT(Igrid)%DEPTHTABLE)
   DEALLOCATE (GWFLAKDAT(Igrid)%AREATABLE)
   DEALLOCATE (GWFLAKDAT(Igrid)%VOLUMETABLE)
   DEALLOCATE (GWFLAKDAT(Igrid)%XLAKES)
   DEALLOCATE (GWFLAKDAT(Igrid)%XLAKINIT)
   DEALLOCATE (GWFLAKDAT(Igrid)%XLKOLD)
   !rsr allocate BD arrays
   DEALLOCATE (GWFLAKDAT(IGRID)%LDRY)
   DEALLOCATE (GWFLAKDAT(IGRID)%NCNT)
   DEALLOCATE (GWFLAKDAT(IGRID)%NCNST)
   DEALLOCATE (GWFLAKDAT(IGRID)%KSUB)
   DEALLOCATE (GWFLAKDAT(IGRID)%MSUB1)
   DEALLOCATE (GWFLAKDAT(IGRID)%MSUB)
   DEALLOCATE (GWFLAKDAT(IGRID)%FLXINL)
   DEALLOCATE (GWFLAKDAT(IGRID)%VOLOLD)
   DEALLOCATE (GWFLAKDAT(IGRID)%GWIN)
   DEALLOCATE (GWFLAKDAT(IGRID)%GWOUT)
   DEALLOCATE (GWFLAKDAT(IGRID)%DELH)
   DEALLOCATE (GWFLAKDAT(IGRID)%TDELH)
   DEALLOCATE (GWFLAKDAT(IGRID)%SVT)
   DEALLOCATE (GWFLAKDAT(IGRID)%STGADJ)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTGWIN_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTGWOT_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTDELSTOR_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTSTOR_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTEVAP_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTPPT_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTRUNF_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTWTHDRW_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTSURFIN_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%TOTSURFOT_LAK)
   DEALLOCATE (GWFLAKDAT(IGRID)%VOLOLDD)
   !dep  Added arrays that calculate lake budgets 6/9/2009
   DEALLOCATE (GWFLAKDAT(IGRID)%DELVOL)
   DEALLOCATE (GWFLAKDAT(IGRID)%TSLAKERR)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMVOL)
   DEALLOCATE (GWFLAKDAT(IGRID)%CMLAKERR)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMLKIN)
   DEALLOCATE (GWFLAKDAT(IGRID)%CUMLKOUT)
   DEALLOCATE (GWFLAKDAT(IGRID)%LAKSEEP)
   DEALLOCATE (GWFLAKDAT(IGRID)%RUNF)      !EDM
   DEALLOCATE (GWFLAKDAT(IGRID)%RUNOFF)    !EDM
   DEALLOCATE (GWFLAKDAT(IGRID)%LKFLOWTYPE)
   DEALLOCATE (GWFLAKDAT(IGRID)%NLKFLWTYP)
   DEALLOCATE (GWFLAKDAT(IGRID)%IGSFLOWLAK)
END SUBROUTINE GWF2LAK7DA

SUBROUTINE SGWF2LAK7PNT(IGRID)
   !  Set pointers to LAK data for grid
   USE GWFLAKMODULE
   !
   LKFLOWTYPE=>GWFLAKDAT(IGRID)%LKFLOWTYPE
   NLKFLWTYP=>GWFLAKDAT(IGRID)%NLKFLWTYP
   NLAKES=>GWFLAKDAT(IGRID)%NLAKES
   NLAKESAR=>GWFLAKDAT(IGRID)%NLAKESAR
   ILKCB=>GWFLAKDAT(IGRID)%ILKCB
   LAKTAB=>GWFLAKDAT(IGRID)%LAKTAB
   IRDTAB=>GWFLAKDAT(IGRID)%IRDTAB
   ISTARTLAK=>GWFLAKDAT(IGRID)%ISTARTLAK
   NSSITR=>GWFLAKDAT(IGRID)%NSSITR
   MXLKND=>GWFLAKDAT(IGRID)%MXLKND
   LKNODE=>GWFLAKDAT(IGRID)%LKNODE
   ICMX=>GWFLAKDAT(IGRID)%ICMX
   NCLS=>GWFLAKDAT(IGRID)%NCLS
   LWRT=>GWFLAKDAT(IGRID)%LWRT
   NDV=>GWFLAKDAT(IGRID)%NDV
   NTRB=>GWFLAKDAT(IGRID)%NTRB
   THETA=>GWFLAKDAT(IGRID)%THETA
   SSCNCR=>GWFLAKDAT(IGRID)%SSCNCR
   !dep  added SURFDEPTH 3/3/2009
   SURFDEPTH=>GWFLAKDAT(IGRID)%SURFDEPTH
   RAMP=>GWFLAKDAT(IGRID)%RAMP
   SMALLTOL=>GWFLAKDAT(IGRID)%SMALLTOL
   DEADPOOLVOL=>GWFLAKDAT(IGRID)%DEADPOOLVOL
   RELEASABLE_STOR=>GWFLAKDAT(IGRID)%RELEASABLE_STOR
   MXLKVOLF=>GWFLAKDAT(IGRID)%MXLKVOLF
   ICS=>GWFLAKDAT(IGRID)%ICS
   NCNCVR=>GWFLAKDAT(IGRID)%NCNCVR
   LIMERR=>GWFLAKDAT(IGRID)%LIMERR
   ILAKE=>GWFLAKDAT(IGRID)%ILAKE
   ITRB=>GWFLAKDAT(IGRID)%ITRB
   IDIV=>GWFLAKDAT(IGRID)%IDIV
   ISUB=>GWFLAKDAT(IGRID)%ISUB
   IRK=>GWFLAKDAT(IGRID)%IRK
   LKARR1=>GWFLAKDAT(IGRID)%LKARR1
   STAGES=>GWFLAKDAT(IGRID)%STAGES
   STGNEW=>GWFLAKDAT(IGRID)%STGNEW
   STGOLD=>GWFLAKDAT(IGRID)%STGOLD
   STGOLD2=>GWFLAKDAT(IGRID)%STGOLD2
   STGITER=>GWFLAKDAT(IGRID)%STGITER
   VOL=>GWFLAKDAT(IGRID)%VOL
   FLOB=>GWFLAKDAT(IGRID)%FLOB
   DSRFOT=>GWFLAKDAT(IGRID)%DSRFOT
   PRCPLK=>GWFLAKDAT(IGRID)%PRCPLK
   EVAPLK=>GWFLAKDAT(IGRID)%EVAPLK
   BEDLAK=>GWFLAKDAT(IGRID)%BEDLAK
   WTHDRW=>GWFLAKDAT(IGRID)%WTHDRW
   RNF=>GWFLAKDAT(IGRID)%RNF
   CUMRNF=>GWFLAKDAT(IGRID)%CUMRNF
   CUMPPT=>GWFLAKDAT(IGRID)%CUMPPT
   CUMEVP=>GWFLAKDAT(IGRID)%CUMEVP
   CUMGWI=>GWFLAKDAT(IGRID)%CUMGWI
   CUMGWO=>GWFLAKDAT(IGRID)%CUMGWO
   CUMSWI=>GWFLAKDAT(IGRID)%CUMSWI
   CUMSWO=>GWFLAKDAT(IGRID)%CUMSWO
   CUMWDR=>GWFLAKDAT(IGRID)%CUMWDR
   CUMFLX=>GWFLAKDAT(IGRID)%CUMFLX
   CNDFCT=>GWFLAKDAT(IGRID)%CNDFCT
   VOLINIT=>GWFLAKDAT(IGRID)%VOLINIT
   BOTTMS=>GWFLAKDAT(IGRID)%BOTTMS
   BGAREA=>GWFLAKDAT(IGRID)%BGAREA
   SSMN=>GWFLAKDAT(IGRID)%SSMN
   SSMX=>GWFLAKDAT(IGRID)%SSMX
   EVAP=>GWFLAKDAT(IGRID)%EVAP
   PRECIP=>GWFLAKDAT(IGRID)%PRECIP
   EVAP3=>GWFLAKDAT(IGRID)%EVAP3
   PRECIP3=>GWFLAKDAT(IGRID)%PRECIP3
   SEEP=>GWFLAKDAT(IGRID)%SEEP
   SEEP3=>GWFLAKDAT(IGRID)%SEEP3
   SURFA=>GWFLAKDAT(IGRID)%SURFA
   SURFIN=>GWFLAKDAT(IGRID)%SURFIN
   SURFOT=>GWFLAKDAT(IGRID)%SURFOT
   SUMCNN=>GWFLAKDAT(IGRID)%SUMCNN
   SUMCHN=>GWFLAKDAT(IGRID)%SUMCHN
   CLAKE=>GWFLAKDAT(IGRID)%CLAKE
   CRNF=>GWFLAKDAT(IGRID)%CRNF
   SILLVT=>GWFLAKDAT(IGRID)%SILLVT
   CAUG=>GWFLAKDAT(IGRID)%CAUG
   CPPT=>GWFLAKDAT(IGRID)%CPPT
   CLAKINIT=>GWFLAKDAT(IGRID)%CLAKINIT
   BDLKN1=>GWFLAKDAT(IGRID)%BDLKN1
   !dep  Added arrays that track lake budgets for dry lakes
   WITHDRW=>GWFLAKDAT(Igrid)%WITHDRW
   FLWIN=>GWFLAKDAT(Igrid)%FLWIN
   FLWITER=>GWFLAKDAT(Igrid)%FLWITER
   FLWITER3=>GWFLAKDAT(Igrid)%FLWITER3
   GWRATELIM=>GWFLAKDAT(Igrid)%GWRATELIM
   !dep  added two variable arrays
   OVRLNDRNF=>GWFLAKDAT(Igrid)%OVRLNDRNF
   CUMLNDRNF=>GWFLAKDAT(Igrid)%CUMLNDRNF
   CUMUZF=>GWFLAKDAT(Igrid)%CUMUZF
   !dep  added three variable arrays for depth,area, and volume
   DEPTHTABLE=>GWFLAKDAT(Igrid)%DEPTHTABLE
   AREATABLE=>GWFLAKDAT(Igrid)%AREATABLE
   VOLUMETABLE=>GWFLAKDAT(Igrid)%VOLUMETABLE
   XLAKES=>GWFLAKDAT(Igrid)%XLAKES
   XLAKINIT=>GWFLAKDAT(Igrid)%XLAKINIT
   XLKOLD=>GWFLAKDAT(Igrid)%XLKOLD
   !rsr allocate BD arrays
   LDRY=>GWFLAKDAT(IGRID)%LDRY
   NCNT=>GWFLAKDAT(IGRID)%NCNT
   NCNST=>GWFLAKDAT(IGRID)%NCNST
   KSUB=>GWFLAKDAT(IGRID)%KSUB
   MSUB1=>GWFLAKDAT(IGRID)%MSUB1
   MSUB=>GWFLAKDAT(IGRID)%MSUB
   FLXINL=>GWFLAKDAT(IGRID)%FLXINL
   VOLOLD=>GWFLAKDAT(IGRID)%VOLOLD
   GWIN=>GWFLAKDAT(IGRID)%GWIN
   GWOUT=>GWFLAKDAT(IGRID)%GWOUT
   DELH=>GWFLAKDAT(IGRID)%DELH
   TDELH=>GWFLAKDAT(IGRID)%TDELH
   SVT=>GWFLAKDAT(IGRID)%SVT
   STGADJ=>GWFLAKDAT(IGRID)%STGADJ
   TOTGWIN_LAK=>GWFLAKDAT(IGRID)%TOTGWIN_LAK
   TOTGWOT_LAK=>GWFLAKDAT(IGRID)%TOTGWOT_LAK
   TOTDELSTOR_LAK=>GWFLAKDAT(IGRID)%TOTDELSTOR_LAK
   TOTSTOR_LAK=>GWFLAKDAT(IGRID)%TOTSTOR_LAK
   TOTEVAP_LAK=>GWFLAKDAT(IGRID)%TOTEVAP_LAK
   TOTPPT_LAK=>GWFLAKDAT(IGRID)%TOTPPT_LAK
   TOTRUNF_LAK=>GWFLAKDAT(IGRID)%TOTRUNF_LAK
   TOTWTHDRW_LAK=>GWFLAKDAT(IGRID)%TOTWTHDRW_LAK
   TOTSURFIN_LAK=>GWFLAKDAT(IGRID)%TOTSURFIN_LAK
   TOTSURFOT_LAK=>GWFLAKDAT(IGRID)%TOTSURFOT_LAK
   LAKUNIT=>GWFLAKDAT(IGRID)%LAKUNIT
   VOLOLDD=>GWFLAKDAT(IGRID)%VOLOLDD
   !dep  Allocate lake budget error arrays 6/9/2009
   DELVOL=>GWFLAKDAT(IGRID)%DELVOL
   TSLAKERR=>GWFLAKDAT(IGRID)%TSLAKERR
   CUMVOL=>GWFLAKDAT(IGRID)%CUMVOL
   CMLAKERR=>GWFLAKDAT(IGRID)%CMLAKERR
   CUMLKOUT=>GWFLAKDAT(IGRID)%CUMLKOUT
   CUMLKIN=>GWFLAKDAT(IGRID)%CUMLKIN
   LAKSEEP=>GWFLAKDAT(IGRID)%LAKSEEP
   RUNF=>GWFLAKDAT(IGRID)%RUNF        !EDM
   RUNOFF=>GWFLAKDAT(IGRID)%RUNOFF    !EDM
   IGSFLOWLAK=>GWFLAKDAT(IGRID)%IGSFLOWLAK
END SUBROUTINE SGWF2LAK7PNT

SUBROUTINE SGWF2LAK7PSV1(IGRID)
   !  Save LAK data for a grid for data shared with SFR
   USE GWFLAKMODULE
   !
   GWFLAKDAT(IGRID)%NLAKES=>NLAKES
   GWFLAKDAT(IGRID)%NLAKESAR=>NLAKESAR
   GWFLAKDAT(IGRID)%THETA=>THETA
   GWFLAKDAT(IGRID)%STGOLD=>STGOLD
   GWFLAKDAT(IGRID)%STGOLD2=>STGOLD2
   GWFLAKDAT(IGRID)%STGNEW=>STGNEW
   GWFLAKDAT(IGRID)%STGITER=>STGITER
   GWFLAKDAT(IGRID)%VOL=>VOL
   GWFLAKDAT(IGRID)%LAKUNIT=>LAKUNIT
END SUBROUTINE SGWF2LAK7PSV1

SUBROUTINE SGWF2LAK7PSV(IGRID)
   !  Save LAK data for a grid
   USE GWFLAKMODULE
   !
   GWFLAKDAT(IGRID)%LKFLOWTYPE=>LKFLOWTYPE
   GWFLAKDAT(IGRID)%NLKFLWTYP=>NLKFLWTYP
   GWFLAKDAT(IGRID)%ILKCB=>ILKCB
   GWFLAKDAT(IGRID)%NSSITR=>NSSITR
   GWFLAKDAT(IGRID)%MXLKND=>MXLKND
   GWFLAKDAT(IGRID)%LKNODE=>LKNODE
   GWFLAKDAT(IGRID)%ICMX=>ICMX
   GWFLAKDAT(IGRID)%LAKTAB=>LAKTAB
   GWFLAKDAT(IGRID)%IRDTAB=>IRDTAB
   GWFLAKDAT(IGRID)%ISTARTLAK=>ISTARTLAK
   GWFLAKDAT(IGRID)%NCLS=>NCLS
   GWFLAKDAT(IGRID)%LWRT=>LWRT
   GWFLAKDAT(IGRID)%NDV=>NDV
   GWFLAKDAT(IGRID)%NTRB=>NTRB
   GWFLAKDAT(IGRID)%SSCNCR=>SSCNCR
   !dep  Added SURFDEPTH 3/3/2009
   GWFLAKDAT(IGRID)%SURFDEPTH=>SURFDEPTH
   GWFLAKDAT(IGRID)%RAMP=>RAMP
   GWFLAKDAT(IGRID)%SMALLTOL=>SMALLTOL
   GWFLAKDAT(IGRID)%DEADPOOLVOL=>DEADPOOLVOL
   GWFLAKDAT(IGRID)%RELEASABLE_STOR=>RELEASABLE_STOR
   GWFLAKDAT(IGRID)%MXLKVOLF=>MXLKVOLF
   GWFLAKDAT(IGRID)%ICS=>ICS
   GWFLAKDAT(IGRID)%NCNCVR=>NCNCVR
   GWFLAKDAT(IGRID)%LIMERR=>LIMERR
   GWFLAKDAT(IGRID)%ILAKE=>ILAKE
   GWFLAKDAT(IGRID)%ITRB=>ITRB
   GWFLAKDAT(IGRID)%IDIV=>IDIV
   GWFLAKDAT(IGRID)%ISUB=>ISUB
   GWFLAKDAT(IGRID)%IRK=>IRK
   GWFLAKDAT(IGRID)%LKARR1=>LKARR1
   GWFLAKDAT(IGRID)%STAGES=>STAGES
   GWFLAKDAT(IGRID)%FLOB=>FLOB
   GWFLAKDAT(IGRID)%DSRFOT=>DSRFOT
   GWFLAKDAT(IGRID)%PRCPLK=>PRCPLK
   GWFLAKDAT(IGRID)%EVAPLK=>EVAPLK
   GWFLAKDAT(IGRID)%BEDLAK=>BEDLAK
   GWFLAKDAT(IGRID)%WTHDRW=>WTHDRW
   GWFLAKDAT(IGRID)%RNF=>RNF
   GWFLAKDAT(IGRID)%CUMRNF=>CUMRNF
   GWFLAKDAT(IGRID)%CUMPPT=>CUMPPT
   GWFLAKDAT(IGRID)%CUMEVP=>CUMEVP
   GWFLAKDAT(IGRID)%CUMGWI=>CUMGWI
   GWFLAKDAT(IGRID)%CUMGWO=>CUMGWO
   GWFLAKDAT(IGRID)%CUMSWI=>CUMSWI
   GWFLAKDAT(IGRID)%CUMSWO=>CUMSWO
   GWFLAKDAT(IGRID)%CUMWDR=>CUMWDR
   GWFLAKDAT(IGRID)%CUMFLX=>CUMFLX
   GWFLAKDAT(IGRID)%CNDFCT=>CNDFCT
   GWFLAKDAT(IGRID)%VOLINIT=>VOLINIT
   GWFLAKDAT(IGRID)%BOTTMS=>BOTTMS
   GWFLAKDAT(IGRID)%BGAREA=>BGAREA
   GWFLAKDAT(IGRID)%SSMN=>SSMN
   GWFLAKDAT(IGRID)%SSMX=>SSMX
   GWFLAKDAT(IGRID)%EVAP=>EVAP
   GWFLAKDAT(IGRID)%PRECIP=>PRECIP
   GWFLAKDAT(IGRID)%EVAP3=>EVAP3
   GWFLAKDAT(IGRID)%PRECIP3=>PRECIP3
   GWFLAKDAT(IGRID)%SEEP=>SEEP
   GWFLAKDAT(IGRID)%SEEP3=>SEEP3
   GWFLAKDAT(IGRID)%SURFA=>SURFA
   GWFLAKDAT(IGRID)%SURFIN=>SURFIN
   GWFLAKDAT(IGRID)%SURFOT=>SURFOT
   GWFLAKDAT(IGRID)%SUMCNN=>SUMCNN
   GWFLAKDAT(IGRID)%SUMCHN=>SUMCHN
   GWFLAKDAT(IGRID)%CLAKE=>CLAKE
   GWFLAKDAT(IGRID)%CRNF=>CRNF
   GWFLAKDAT(IGRID)%SILLVT=>SILLVT
   GWFLAKDAT(IGRID)%CAUG=>CAUG
   GWFLAKDAT(IGRID)%CPPT=>CPPT
   GWFLAKDAT(IGRID)%CLAKINIT=>CLAKINIT
   GWFLAKDAT(IGRID)%BDLKN1=>BDLKN1
   !dep  Added arrays that track lake budgets for dry lakes
   GWFLAKDAT(Igrid)%WITHDRW=>WITHDRW
   GWFLAKDAT(Igrid)%FLWIN=>FLWIN
   GWFLAKDAT(Igrid)%FLWITER=>FLWITER
   GWFLAKDAT(Igrid)%FLWITER3=>FLWITER3
   GWFLAKDAT(Igrid)%GWRATELIM=>GWRATELIM
   !dep  added two variable arrays
   GWFLAKDAT(Igrid)%OVRLNDRNF=>OVRLNDRNF
   GWFLAKDAT(Igrid)%CUMLNDRNF=>CUMLNDRNF
   GWFLAKDAT(Igrid)%CUMUZF=>CUMUZF
   !dep  added three variable arrays for depth, area, and volume
   GWFLAKDAT(Igrid)%DEPTHTABLE=>DEPTHTABLE
   GWFLAKDAT(Igrid)%AREATABLE=>AREATABLE
   GWFLAKDAT(Igrid)%VOLUMETABLE=>VOLUMETABLE
   GWFLAKDAT(Igrid)%XLAKES=>XLAKES
   GWFLAKDAT(Igrid)%XLAKINIT=>XLAKINIT
   GWFLAKDAT(Igrid)%XLKOLD=>XLKOLD
   !rsr allocate BD arrays
   GWFLAKDAT(IGRID)%LDRY=>LDRY
   GWFLAKDAT(IGRID)%NCNT=>NCNT
   GWFLAKDAT(IGRID)%NCNST=>NCNST
   GWFLAKDAT(IGRID)%KSUB=>KSUB
   GWFLAKDAT(IGRID)%MSUB1=>MSUB1
   GWFLAKDAT(IGRID)%MSUB=>MSUB
   GWFLAKDAT(IGRID)%FLXINL=>FLXINL
   GWFLAKDAT(IGRID)%VOLOLD=>VOLOLD
   GWFLAKDAT(IGRID)%GWIN=>GWIN
   GWFLAKDAT(IGRID)%GWOUT=>GWOUT
   GWFLAKDAT(IGRID)%DELH=>DELH
   GWFLAKDAT(IGRID)%TDELH=>TDELH
   GWFLAKDAT(IGRID)%SVT=>SVT
   GWFLAKDAT(IGRID)%STGADJ=>STGADJ
   !dep  Allocate lake budget error arrays 6/9/2009
   GWFLAKDAT(IGRID)%DELVOL=>DELVOL
   GWFLAKDAT(IGRID)%TSLAKERR=>TSLAKERR
   GWFLAKDAT(IGRID)%CUMVOL=>CUMVOL
   GWFLAKDAT(IGRID)%CMLAKERR=>CMLAKERR
   GWFLAKDAT(IGRID)%CUMLKOUT=>CUMLKOUT
   GWFLAKDAT(IGRID)%CUMLKIN=>CUMLKIN
   !rgn Allocate budget arrays for GSFLOW CSV file
   GWFLAKDAT(IGRID)%TOTGWIN_LAK=>TOTGWIN_LAK
   GWFLAKDAT(IGRID)%TOTGWOT_LAK=>TOTGWOT_LAK
   GWFLAKDAT(IGRID)%TOTDELSTOR_LAK=>TOTDELSTOR_LAK
   GWFLAKDAT(IGRID)%TOTSTOR_LAK=>TOTSTOR_LAK
   GWFLAKDAT(IGRID)%TOTEVAP_LAK=>TOTEVAP_LAK
   GWFLAKDAT(IGRID)%TOTPPT_LAK=>TOTPPT_LAK
   GWFLAKDAT(IGRID)%TOTRUNF_LAK=>TOTRUNF_LAK
   GWFLAKDAT(IGRID)%TOTWTHDRW_LAK=>TOTWTHDRW_LAK
   GWFLAKDAT(IGRID)%TOTSURFIN_LAK=>TOTSURFIN_LAK
   GWFLAKDAT(IGRID)%TOTSURFOT_LAK=>TOTSURFOT_LAK
   GWFLAKDAT(IGRID)%VOLOLDD=>VOLOLDD
   GWFLAKDAT(IGRID)%LAKSEEP=>LAKSEEP
   GWFLAKDAT(IGRID)%RUNF=>RUNF        !EDM
   GWFLAKDAT(IGRID)%RUNOFF=>RUNOFF    !EDM
   GWFLAKDAT(IGRID)%IGSFLOWLAK=>IGSFLOWLAK
END SUBROUTINE SGWF2LAK7PSV
