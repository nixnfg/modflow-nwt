   !     ******************************************************************
   !     MAIN CODE FOR U.S. GEOLOGICAL SURVEY MODULAR MODEL -- MODFLOW-NWT
   !rgn------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
   !rgn------NEW VERSION NUMBER 1.3.0, 7/01/2022
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   !1------USE package modules.
USE GLOBAL
USE GWFBASMODULE
USE GWFHUFMODULE, ONLY:IOHUFHDS,IOHUFFLWS
USE GWFEVTMODULE, ONLY:NEVTOP
USE GWFRCHMODULE, ONLY:NRCHOP
USE GWFLAKMODULE, ONLY:NLAKESAR,THETA,STGOLD,STGNEW,VOL
USE GWFNWTMODULE, ONLY:LINMETH,ICNVGFLG,ITREAL
USE PCGMODULE
USE SIPMODULE
USE DE4MODULE
USE GMGMODULE
   !      USE PCGN
use openspec
   !
   !-------ASSIGN VERSION NUMBER AND DATE
CHARACTER*40 VERSION,VERSION2,VERSION3
CHARACTER*10 MFVNAM
PARAMETER (VERSION='1.3.0 07/01/2022')
PARAMETER (VERSION2='1.12.0 02/03/2017')
PARAMETER (VERSION3='1.05.0 03/10/2022')
PARAMETER (MFVNAM='-NWT-SWR1')
   !
CHARACTER*80 HEADNG(2)
CHARACTER*200 FNAME
INTEGER IBDT(8), AGCONVERGE
   !
CHARACTER*4 CUNIT(NIUNIT)
DATA CUNIT/'BCF6', 'WEL ', 'DRN ', 'RIV ', 'EVT ', 'gfd ', 'GHB ',&  !  7
&'RCH ', 'SIP ', 'DE4 ', '    ', 'OC  ', 'PCG ', 'lmg ',&  ! 14
&'gwt ', 'FHB ', 'RES ', 'STR ', 'IBS ', 'CHD ', 'HFB6',&  ! 21
&'LAK ', 'LPF ', 'DIS ', '    ', 'PVAL', '    ', 'HOB ',&  ! 28
&'    ', '    ', 'ZONE', 'MULT', 'DROB', 'RVOB', 'GBOB',&  ! 35
&'STOB', 'HUF2', 'CHOB', 'ETS ', 'DRT ', '    ', 'GMG ',&  ! 42
&'HYD ', 'SFR ', '    ', 'GAGE', 'LVDA', '    ', 'LMT6',&  ! 49
&'MNW2', 'MNWI', 'MNW1', 'KDEP', 'SUB ', 'UZF ', 'gwm ',&  ! 56
&'SWT ', 'cfp ', 'pcgn', '    ', 'fmp ', 'UPW ', 'NWT ',&  ! 63
&'SWR ', 'SWI2', 'AG ', '    ', '    ', '    ', '    ',&  ! 70     - SWR - JDH
&30*'    '/                                               ! 71-100 - SWR - JDH
   !     ------------------------------------------------------------------
   !
   !2------WRITE BANNER TO SCREEN AND DEFINE CONSTANTS.
WRITE (*,1) MFVNAM,VERSION,VERSION2,VERSION3
1 FORMAT (/,34X,'MODFLOW',A,/,&
&4X,'U.S. GEOLOGICAL SURVEY MODULAR FINITE-DIFFERENCE',&
&' GROUNDWATER-FLOW MODEL',/,29X,'WITH NEWTON FORMULATION',&
&/,29X,'Version ',A/,20X,'BASED ON MODFLOW-2005 Version ',A/,&
&/,20X,'SWR1 Version ',A/)
INUNIT = 99
NCVGERR=0
AGCONVERGE=1
   !
   !3------GET THE NAME OF THE NAME FILE
CALL GETNAMFIL(FNAME)
MAXUNIT= INUNIT
   !
   !4------OPEN NAME FILE.
OPEN (UNIT=INUNIT,FILE=FNAME,STATUS='OLD',ACTION=ACTION(1))
NC=INDEX(FNAME,' ')
WRITE(*,490)' Using NAME file: ',FNAME(1:NC)
490 FORMAT(A,A)
   !
   !5------Get current date and time, assign to IBDT, and write to screen
CALL DATE_AND_TIME(VALUES=IBDT)
WRITE(*,2) (IBDT(I),I=1,3),(IBDT(I),I=5,7)
2 FORMAT(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ',&
&I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2,/)
   !
   !6------ALLOCATE AND READ (AR) PROCEDURE
IGRID=1
NSOL=1
CALL GWF2BAS7AR(INUNIT,CUNIT,VERSION,24,31,32,MAXUNIT,IGRID,12,&
&HEADNG,26,MFVNAM)
IF(IUNIT(50).GT.0 .AND. IUNIT(52).GT.0) THEN
   WRITE(IOUT,'(1X,/,1X,A)')&
   &'MNW1 and MNW2 cannot both be active in the same simulation'
   CALL USTOP(' ')
END IF
IF(IUNIT(50).GT.0 .AND. IUNIT(52).GT.0) THEN
   WRITE(IOUT,'(1X,/,1X,A)')&
   &'MNW1 and MNW2 cannot both be active in the same simulation'
   CALL USTOP(' ')
END IF
IF(IUNIT(1).GT.0) CALL GWF2BCF7AR(IUNIT(1),IGRID)
IF(IUNIT(23).GT.0) CALL GWF2LPF7AR(IUNIT(23),IGRID)
IF(IUNIT(37).GT.0) CALL GWF2HUF7AR(IUNIT(37),IUNIT(47),&
&IUNIT(53),IGRID)
   ! Allocate arrays for Newton Solver
IF(IUNIT(63).GT.0) CALL GWF2NWT1AR(IUNIT(63),MXITER,&
&IUNIT(22),IGRID)
IF(IUNIT(62).GT.0) CALL GWF2UPW1AR(IUNIT(62), Igrid)
IF(IUNIT(2).GT.0) CALL GWF2WEL7AR(IUNIT(2),IUNIT(63),IGRID)
IF(IUNIT(3).GT.0) CALL GWF2DRN7AR(IUNIT(3),IGRID)
IF(IUNIT(4).GT.0) CALL GWF2RIV7AR(IUNIT(4),IGRID)
IF(IUNIT(5).GT.0) CALL GWF2EVT7AR(IUNIT(5),IGRID)
IF(IUNIT(7).GT.0) CALL GWF2GHB7AR(IUNIT(7),IGRID)
IF(IUNIT(8).GT.0) CALL GWF2RCH7AR(IUNIT(8),IGRID)
IF(IUNIT(16).GT.0) CALL GWF2FHB7AR(IUNIT(16),IGRID)
IF(IUNIT(17).GT.0) CALL GWF2RES7AR(IUNIT(17),IGRID)
IF(IUNIT(18).GT.0) CALL GWF2STR7AR(IUNIT(18),IGRID)
IF(IUNIT(19).GT.0) CALL GWF2IBS7AR(IUNIT(19),IUNIT(54),IGRID)
IF(IUNIT(20).GT.0) CALL GWF2CHD7AR(IUNIT(20),IGRID)
IF(IUNIT(21).GT.0) CALL GWF2HFB7AR(IUNIT(21),IGRID)
   ! Modify conductance for HFB when using UPW.
IF ( IUNIT(62).GT.0 ) THEN
   IF(IUNIT(21).GT.0) CALL GWF2HFB7UPW(IGRID)
END IF
IF(IUNIT(44).GT.0) CALL GWF2SFR7AR(IUNIT(44),IUNIT(1),IUNIT(23),&
&IUNIT(37),IUNIT(15),NSOL,IOUTS,&
&IUNIT(62),IUNIT(55),IGRID)
IF(IUNIT(55).GT.0) CALL GWF2UZF1AR(IUNIT(55),IUNIT(1),&
&IUNIT(23),IUNIT(37),&
&IUNIT(63),IGRID)
IF(IUNIT(22).GT.0 .OR. IUNIT(44).GT.0) CALL GWF2LAK7AR(&
&IUNIT(22),IUNIT(44),IUNIT(15),IUNIT(55),NSOL,IGRID)
IF(IUNIT(46).GT.0) CALL GWF2GAG7AR(IUNIT(46),IUNIT(44),&
&IUNIT(22),IGRID)
IF(IUNIT(39).GT.0) CALL GWF2ETS7AR(IUNIT(39),IGRID)
IF(IUNIT(40).GT.0) CALL GWF2DRT7AR(IUNIT(40),IGRID)
IF(IUNIT(54).GT.0) CALL GWF2SUB7AR(IUNIT(54),IGRID)
IF(IUNIT(9).GT.0) CALL SIP7AR(IUNIT(9),MXITER,IGRID)
IF(IUNIT(10).GT.0) CALL DE47AR(IUNIT(10),MXITER,IGRID)
IF(IUNIT(13).GT.0) CALL PCG7AR(IUNIT(13),MXITER,IGRID)
   !      IF(IUNIT(14).GT.0) CALL LMG7AR(IUNIT(14),MXITER,IGRID)
IF(IUNIT(42).GT.0) CALL GMG7AR(IUNIT(42),MXITER,IGRID)
   !      IF(IUNIT(59).GT.0) CALL PCGN2AR(IUNIT(59),IFREFM,MXITER,IGRID)
IF(IUNIT(50).GT.0) CALL GWF2MNW27AR(IUNIT(50),IGRID)
IF(IUNIT(51).GT.0) CALL GWF2MNW2I7AR(IUNIT(51),IUNIT(50),IGRID)
IF(IUNIT(52).GT.0) CALL GWF2MNW17AR(IUNIT(52),IUNIT(9),&
&IUNIT(10),IUNIT(63),0,IUNIT(13),&
&0,IUNIT(42),FNAME,IGRID)
IF(IUNIT(57).GT.0) CALL GWF2SWT7AR(IUNIT(57),IGRID)
IF(IUNIT(64).GT.0) CALL GWF2SWR7AR(IUNIT(64),&
&IUNIT(1),IUNIT(23),IUNIT(37),&
&IUNIT(62),IUNIT(44),IUNIT(63),IGRID)  !SWR  - JDH
IF(IUNIT(65).GT.0) CALL GWF2SWI2AR(IUNIT(65),&
&IUNIT(1),IUNIT(23),IUNIT(37),IUNIT(62),&
&IGRID)   !SWI2  - JDH
   !     IF(IUNIT(66).GT.0) CALL GWF2GFB7AR(IUNIT(66),IGRID)
IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7AR(IUNIT(43),IGRID)
IF(IUNIT(43).GT.0 .AND. IUNIT(19).GT.0)&
&CALL GWF2HYD7IBS7AR(IUNIT(43),IGRID)
IF(IUNIT(43).GT.0 .AND. IUNIT(54).GT.0)&
&CALL GWF2HYD7SUB7AR(IUNIT(43),IGRID)
IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)&
&CALL GWF2HYD7STR7AR(IUNIT(43),IGRID)
IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)&
&CALL GWF2HYD7SFR7AR(IUNIT(43),IGRID)
IF(IUNIT(49).GT.0) CALL LMT8BAS7AR(INUNIT,CUNIT,IGRID)
IF(IUNIT(66).GT.0) CALL GWF2AG7AR(IUNIT(66),IUNIT(44),IUNIT(63))
   !      IF(IUNIT(61).GT.0) THEN
   !        CALL FMP2AR(
   !     1  IUNIT(61),IUNIT(44),IUNIT(52),IUNIT(55),IGRID)                  !FMP2AR CALL ADDED BY SCHMID
   !        CALL FMP2RQ(IUNIT(61),IUNIT(44),IUNIT(52),IGRID)                !FMP2RQ CALL ADDED BY SCHMID
   !      ENDIF
   !
   !  Observation allocate and read
CALL OBS2BAS7AR(IUNIT(28),IGRID)
IF(IUNIT(33).GT.0) CALL OBS2DRN7AR(IUNIT(33),IUNIT(3),IGRID)
IF(IUNIT(34).GT.0) CALL OBS2RIV7AR(IUNIT(34),IUNIT(4),IGRID)
IF(IUNIT(35).GT.0) CALL OBS2GHB7AR(IUNIT(35),IUNIT(7),IGRID)
IF(IUNIT(36).GT.0) CALL OBS2STR7AR(IUNIT(36),IUNIT(18),IGRID)
IF(IUNIT(38).GT.0) CALL OBS2CHD7AR(IUNIT(38),IGRID)
   ! Modify conductance for HFB when using UPW.
 !IF ( IUNIT(62).GT.0 ) THEN
 !  IF(IUNIT(21).GT.0) CALL GWF2HFB7UPW(IGRID)
 !END IF
   !
   !7------SIMULATE EACH STRESS PERIOD.
DO 100 KPER = 1, NPER
   KKPER = KPER
   IF(IUNIT(62).GT.0 ) CALL GWF2UPWUPDATE(1,Igrid)
   CALL GWF2BAS7ST(KKPER,IGRID)
   IF(IUNIT(19).GT.0) CALL GWF2IBS7ST(KKPER,IGRID)
   IF(IUNIT(54).GT.0) CALL GWF2SUB7ST(KKPER,IGRID)
   IF(IUNIT(57).GT.0) CALL GWF2SWT7ST(KKPER,IGRID)
   !
   !7B-----READ AND PREPARE INFORMATION FOR STRESS PERIOD.
   !----------READ USING PACKAGE READ AND PREPARE MODULES.
   IF(IUNIT(2).GT.0) CALL GWF2WEL7RP(IUNIT(2),KPER,IGRID)
   IF(IUNIT(3).GT.0) CALL GWF2DRN7RP(IUNIT(3),IGRID)
   IF(IUNIT(4).GT.0) CALL GWF2RIV7RP(IUNIT(4),IGRID)
   IF(IUNIT(5).GT.0) CALL GWF2EVT7RP(IUNIT(5),IGRID)
   IF(IUNIT(7).GT.0) CALL GWF2GHB7RP(IUNIT(7),IGRID)
   !        IF(IUNIT(8).GT.0) CALL GWF2RCH7RP(IUNIT(8),IUNIT(44),IGRID)
   IF(IUNIT(8).GT.0) CALL GWF2RCH7RP(IUNIT(8),IGRID)
   IF(IUNIT(17).GT.0) CALL GWF2RES7RP(IUNIT(17),IGRID)
   IF(IUNIT(18).GT.0) CALL GWF2STR7RP(IUNIT(18),IGRID)
   IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)&
   &CALL GWF2HYD7STR7RP(IUNIT(43),KKPER,IGRID)
   IF(IUNIT(20).GT.0) CALL GWF2CHD7RP(IUNIT(20),IGRID)
   IF(IUNIT(44).GT.0) CALL GWF2SFR7RP(IUNIT(44),IUNIT(15),&
   &IUNIT(22),KKPER,KKSTP,&
   &NSOL,IOUTS,IUNIT(55),IGRID)
   IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)&
   &CALL GWF2HYD7SFR7RP(IUNIT(43),KKPER,IGRID)
   IF(IUNIT(55).GT.0) CALL GWF2UZF1RP(IUNIT(55),KKPER,KKSTP,&
   &IUNIT(44),IGRID)
   IF(IUNIT(22).GT.0) CALL GWF2LAK7RP(IUNIT(22),IUNIT(1),&
   &IUNIT(15),IUNIT(23),IUNIT(37),IUNIT(44),IUNIT(55),&
   &IUNIT(62),KKPER,NSOL,IOUTS,IGRID)
   IF(IUNIT(46).GT.0.AND.KKPER.EQ.1) CALL GWF2GAG7RP(IUNIT(15),&
   &IUNIT(22),IUNIT(55),NSOL,IGRID)
   IF(IUNIT(39).GT.0) CALL GWF2ETS7RP(IUNIT(39),IGRID)
   IF(IUNIT(40).GT.0) CALL GWF2DRT7RP(IUNIT(40),IGRID)
   IF(IUNIT(50).GT.0) CALL GWF2MNW27RP(IUNIT(50),kper,IUNIT(9),&
   &IUNIT(10),0,IUNIT(13),&
   &IUNIT(15),IUNIT(63),IGRID)
   IF(IUNIT(51).GT.0.AND.KKPER.EQ.1) CALL GWF2MNW2I7RP(IUNIT(51),&
   &0,IGRID)
   IF(IUNIT(52).GT.0) CALL GWF2MNW17RP(IUNIT(52),IUNIT(1),&
   &IUNIT(23),IUNIT(37),IUNIT(62),KKPER,&
   &IGRID)
   IF ( IUNIT(66).GT.0 ) CALL GWF2AG7RP(IUNIT(66),IUNIT(44),KKPER)
   !        IF(IUNIT(61).GT.0) CALL FMP2RP(IUNIT(61),ISTARTFL,KKPER,        !FMP2AR CALL ADDED BY SCHMID
   !     1                          IUNIT(44),IUNIT(52),IGRID)
   IF(IUNIT(64).GT.0) CALL GWF2SWR7RP(IUNIT(64),KKPER,IGRID)  !SWR - JDH
   !IF(IUNIT(66).GT.0 ) CALL GWF2AG7AD(IUNIT(66),KKPER)      not sure why this is here.
   !
   !7C-----SIMULATE EACH TIME STEP.
   DO 90 KSTP = 1, NSTP(KPER)
      KKSTP = KSTP
   !
   !7C1----CALCULATE TIME STEP LENGTH. SET HOLD=HNEW.
      IF(IUNIT(62).GT.0 ) CALL GWF2UPWUPDATE(1,Igrid)
      CALL GWF2BAS7AD(KKPER,KKSTP,IGRID)
      IF(IUNIT(62).GT.0) CALL GWF2UPW1AD(IGRID)
      IF(IUNIT(20).GT.0) CALL GWF2CHD7AD(KKPER,IGRID)
      IF(IUNIT(1).GT.0) CALL GWF2BCF7AD(KKPER,IGRID)
      IF(IUNIT(17).GT.0) CALL GWF2RES7AD(KKSTP,KKPER,IGRID)
      IF(IUNIT(23).GT.0) CALL GWF2LPF7AD(KKPER,IGRID)
      IF(IUNIT(37).GT.0) CALL GWF2HUF7AD(KKPER,IGRID)
      IF(IUNIT(16).GT.0) CALL GWF2FHB7AD(IGRID)
      IF(IUNIT(22).GT.0) CALL GWF2LAK7AD(KKPER,KKSTP,IUNIT(15),&
      &IGRID)
      IF(IUNIT(55).GT.0) CALL GWF2UZF1AD(IUNIT(55), KKPER, KKSTP,&
      &Igrid)
      IF(IUNIT(65).GT.0) CALL GWF2SWI2AD(KKSTP,KKPER,IGRID)  !SWI2
      IF( IUNIT(44).GT.0 ) CALL GWF2SFR7AD(IUNIT(44),IUNIT(22),&
      &KKSTP,KKPER,IGRID)  !rgn 6/12/12
      IF(IUNIT(50).GT.0) THEN
         IF (IUNIT(1).GT.0) THEN
            CALL GWF2MNW27BCF(KPER,IGRID)
         ELSE IF (IUNIT(23).GT.0) THEN
            CALL GWF2MNW27LPF(KPER,IGRID)
         ELSE IF(IUNIT(37).GT.0) THEN
            CALL GWF2MNW27HUF(KPER,IGRID)
         ELSE IF(IUNIT(62).GT.0) THEN
            CALL GWF2MNW27UPW(KPER,IGRID)
         ELSE
            WRITE(IOUT,1000)
1000        FORMAT(/1X,&
            &'***ERROR: MNW2 PACKAGE DOES NOT SUPPORT',/,&
            &' SELECTED FLOW PACKAGE',/,&
            &' (MNW2 DOES FULLY SUPPORT BCF, LPF, HUF, AND UPW PACKAGES)',/,&
            &' -- STOP EXECUTION')
            CALL USTOP('MNW2 error-flow package')
         END IF
         CALL GWF2MNW27AD(KKSTP,KKPER,IUNIT(62),IGRID)
      END IF
      IF(IUNIT(52).GT.0) CALL GWF2MNW17AD(IUNIT(1),IUNIT(23),&
      &IUNIT(37),IUNIT(62),IGRID)
   !          IF(IUNIT(61).GT.0) THEN                                       !FMP2AD CALL ADDED BY SCHMID
   !             IF(IRTFL.EQ.3.OR.ICUFL.EQ.3.OR.IPFL.EQ.3
   !     1                           .OR.IEBFL.EQ.1.OR.IEBFL.EQ.3)
   !     2       CALL FMP2AD(ISTARTFL,KKPER,IGRID)
   !          ENDIF
      IF(IUNIT(64).GT.0) CALL GWF2SWR7AD(KKPER,KKSTP,&
      &IGRID,IUNIT(54))  !SWR - JDH
      IF(IUNIT(66).GT.0) CALL GWF2AG7AD(IUNIT(66),KKPER)
   !
   !---------INDICATE IN PRINTOUT THAT SOLUTION IS FOR HEADS
      CALL UMESPR('SOLVING FOR HEAD',' ',IOUT)
      WRITE(*,25)KPER,KSTP
25    FORMAT(' Solving:  Stress period: ',i5,4x,&
      &'Time step: ',i5,4x,'Groundwater-Flow Eqn.')
26    FORMAT('Skipping:  Stress period: ',i5,4x,&
      &'Time step: ',i5)
   !
   !7C2----ITERATIVELY FORMULATE AND SOLVE THE FLOW EQUATIONS.
   !          DO 30 KITER = 1, MXITER
      KITER = 0
      ITREAL2 = 0
      IF ( IUNIT(63).GT.0 ) ITREAL = 0
      DO WHILE (ITREAL2.LT.MXITER)
         KITER = KITER + 1
         KKITER = KITER
         IF ( IUNIT(63).EQ.0 ) ITREAL2 = KITER
         IF(IUNIT(62).GT.0) CALL GWF2UPWUPDATE(2,Igrid)

   !
   !7C2A---FORMULATE THE FINITE DIFFERENCE EQUATIONS.
         CALL GWF2BAS7FM(IGRID)
         IF(IUNIT(1).GT.0) CALL GWF2BCF7FM(KKITER,KKSTP,&
         &KKPER,IGRID)
         IF(IUNIT(62).GT.0) CALL GWF2UPWFMS(KKITER,KKSTP,KKPER,IGRID)
         IF(IUNIT(23).GT.0) CALL GWF2LPF7FM(KKITER,KKSTP,KKPER,IGRID)
         IF(IUNIT(37).GT.0) CALL GWF2HUF7FM(KKITER,&
         &KKSTP,KKPER,IUNIT(47),IGRID)
         IF ( IUNIT(62).EQ.0 ) THEN
            IF(IUNIT(21).GT.0) CALL GWF2HFB7FM(IGRID)
         END IF
         IF(IUNIT(2).GT.0) CALL GWF2WEL7FM(IUNIT(63),IGRID)
         IF(IUNIT(3).GT.0) CALL GWF2DRN7FM(IGRID)
         IF(IUNIT(4).GT.0) CALL GWF2RIV7FM(IGRID)
         IF(IUNIT(5).GT.0) THEN
            IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(&
            &0,IGRID)
            CALL GWF2EVT7FM(IGRID)
            IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(&
            &1,IGRID)
         END IF
         IF(IUNIT(7).GT.0) CALL GWF2GHB7FM(IGRID)
         IF(IUNIT(8).GT.0) THEN
            IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(&
            &0,IGRID)
            CALL GWF2RCH7FM(IGRID)
            IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(&
            &1,IGRID)
         END IF
         IF(IUNIT(16).GT.0) CALL GWF2FHB7FM(IGRID)
         IF(IUNIT(17).GT.0) CALL GWF2RES7FM(IGRID)
         IF(IUNIT(18).GT.0) CALL GWF2STR7FM(IGRID)
         IF(IUNIT(19).GT.0) CALL GWF2IBS7FM(KKPER,IGRID)
         IF(IUNIT(39).GT.0) CALL GWF2ETS7FM(IGRID)
         IF(IUNIT(40).GT.0) CALL GWF2DRT7FM(IGRID)
   !            IF(IUNIT(61).GT.0) CALL FMP2FM(KKITER,KKPER,KKSTP,ISTARTFL, !FMP2FM CALL ADDED BY SCHMID
   !     1                              IUNIT(44),IUNIT(52),IUNIT(55),IGRID)
         IF(IUNIT(55).GT.0) CALL GWF2UZF1FM(KKPER,KKSTP,KKITER,&
         &IUNIT(44),IUNIT(22),IUNIT(63),&
         &IUNIT(64),IGRID)  !SWR - JDH ADDED IUNIT(64)
         IF(IUNIT(44).GT.0) CALL GWF2SFR7FM(KKITER,KKPER,KKSTP,&
         &IUNIT(22),IUNIT(63),IUNIT(8),&
         &IUNIT(55),IGRID)   !cjm (added IUNIT(8))
         IF(IUNIT(22).GT.0) CALL GWF2LAK7FM(KKITER,KKPER,KKSTP,&
         &IUNIT(44),IUNIT(55),IGRID)
         IF(IUNIT(50).GT.0) THEN
            IF (IUNIT(1).GT.0) THEN
               CALL GWF2MNW27BCF(KPER,IGRID)
            ELSE IF (IUNIT(23).GT.0) THEN
               CALL GWF2MNW27LPF(KPER,IGRID)
            ELSE IF(IUNIT(37).GT.0) THEN
               CALL GWF2MNW27HUF(KPER,IGRID)
            ELSE IF(IUNIT(62).GT.0) THEN
               CALL GWF2MNW27UPW(KPER,IGRID)
            END IF
            CALL GWF2MNW27FM(KKITER,IUNIT(62),kkstp,kkper,IGRID)
         END IF
         IF(IUNIT(52).GT.0) CALL GWF2MNW17FM(KKITER,IUNIT(1),&
         &IUNIT(23),IUNIT(37),IUNIT(62),&
         &IGRID)
         IF(IUNIT(54).GT.0) CALL GWF2SUB7FM(KKPER,KKITER,&
         &IUNIT(9),IGRID)
         IF(IUNIT(57).GT.0) CALL GWF2SWT7FM(KKPER,IGRID)
         IF(IUNIT(64).GT.0) CALL GWF2SWR7FM(KKITER,KKPER,KKSTP,IGRID)  !SWR - JDH
         IF(IUNIT(66).GT.0 )&
         &CALL GWF2AG7FM(Kkper,Kkstp,Kkiter,IUNIT(63),AGCONVERGE)
   !            IF(IUNIT(66).GT.0) CALL GWF2GFB7FM(IGRID)
   !-------------SWI2 FORMULATE (GWF2SWI2FM) NEEDS TO BE THE LAST PACKAGE
   !             ENTRY SINCE SWI2 SAVES THE RHS (RHSFRESH) PRIOR TO ADDING SWI TERMS
   !             RHSFRESH IS USED TO CALCULATE BOUNDARY CONDITION FLUXES
         IF(IUNIT(65).GT.0) CALL GWF2SWI2FM(KKSTP,KKPER,KKITER,IGRID)  !SWI2 - JDH
   !
   !7C2B---MAKE ONE CUT AT AN APPROXIMATE SOLUTION.
         IF (IUNIT(9).GT.0) THEN
            CALL SIP7PNT(IGRID)
            CALL SIP7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,EL,FL,GL,&
            &V,W,HDCG,LRCH,NPARM,KKITER,HCLOSE,ACCL,ICNVG,&
            &KKSTP,KKPER,IPCALC,IPRSIP,MXITER,NSTP(KKPER),&
            &NCOL,NROW,NLAY,NODES,IOUT,0,IERR)
         END IF
         IF (IUNIT(10).GT.0) THEN
            CALL DE47PNT(IGRID)
            CALL DE47AP(HNEW,IBOUND,AU,AL,IUPPNT,IEQPNT,D4B,MXUP,&
            &MXLOW,MXEQ,MXBW,CR,CC,CV,HCOF,RHS,ACCLDE4,KITER,&
            &ITMX,MXITER,NITERDE4,HCLOSEDE4,IPRD4,ICNVG,NCOL,&
            &NROW,NLAY,IOUT,LRCHDE4,HDCGDE4,IFREQ,KKSTP,KKPER,&
            &DELT,NSTP(KKPER),ID4DIR,ID4DIM,MUTD4,&
            &DELTL,NBWL,NUPL,NLOWL,NLOW,NEQ,NUP,NBW,IERR)
         END IF
         IF (IUNIT(13).GT.0) THEN
            CALL PCG7PNT(IGRID)
            CALL PCG7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,VPCG,SS,&
            &P,CD,HCHG,LHCH,RCHG,LRCHPCG,KKITER,NITER,&
            &HCLOSEPCG,RCLOSEPCG,ICNVG,KKSTP,KKPER,IPRPCG,&
            &MXITER,ITER1,NPCOND,NBPOL,NSTP(KKPER),NCOL,NROW,&
            &NLAY,NODES,RELAXPCG,IOUT,MUTPCG,IT1,DAMPPCG,BUFF,&
            &HCSV,IERR,HPCG,DAMPPCGT,ISSFLG(KKPER),HDRY,&
            &IHCOFADD)
         END IF
         IF (IUNIT(42).GT.0) THEN
            CALL GMG7PNT(IGRID)
            CALL GMG7AP(HNEW,RHS,CR,CC,CV,HCOF,HNOFLO,IBOUND,&
            &IITER,MXITER,RCLOSEGMG,HCLOSEGMG,&
            &KKITER,KKSTP,KKPER,NCOL,NROW,NLAY,ICNVG,&
            &SITER,TSITER,DAMPGMG,IADAMPGMG,IOUTGMG,&
            &IOUT,GMGID,&
            &IUNITMHC,DUP,DLOW,CHGLIMIT,&
            &BIGHEADCHG,HNEWLAST)
         ENDIF
   !            IF (IUNIT(59).GT.0) THEN
   !              CALL PCGN2AP(HNEW,RHS,CR,CC,CV,HCOF,IBOUND,
   !     1              KKITER,KKSTP,KKPER,ICNVG,HNOFLO,IGRID)
   !            ENDIF
   ! Calculate new heads using Newton solver
         IF(IUNIT(63).GT.0 )&
         &CALL GWF2NWT1FM(KKITER,ICNVG,KSTP,KPER,Mxiter,&
         &IUNIT(22),AGCONVERGE,IGRID)
         IF ( IUNIT(63).GT.0 )ITREAL2 = ITREAL
   !
   !-------ENSURE CONVERGENCE OF SWR - BASEFLOW CHANGES LESS THAN TOLF - JDH
         IF(IUNIT(64).GT.0) THEN
            CALL GWF2SWR7CV(KKITER,IGRID,ICNVG,MXITER)
         END IF
   !
   !7C2C---IF CONVERGENCE CRITERION HAS BEEN MET STOP ITERATING.
         IF (ICNVG.EQ.1) GOTO 33
   !  30      CONTINUE
      END DO
      KITER = MXITER
   !
33    CONTINUE
      IF(IUNIT(62).GT.0 ) CALL GWF2UPWUPDATE(2,Igrid)
   !
   !7C3----DETERMINE WHICH OUTPUT IS NEEDED.
      CALL GWF2BAS7OC(KKSTP,KKPER,ICNVG,IUNIT(12),IGRID)
   !
   !7C4----CALCULATE BUDGET TERMS. SAVE CELL-BY-CELL FLOW TERMS.
      MSUM = 1
      IF (IUNIT(1).GT.0) THEN
         CALL GWF2BCF7BDS(KKSTP,KKPER,IGRID)
         CALL GWF2BCF7BDCH(KKSTP,KKPER,IGRID)
         IBDRET=0
         IC1=1
         IC2=NCOL
         IR1=1
         IR2=NROW
         IL1=1
         IL2=NLAY
         DO 37 IDIR = 1, 3
            CALL GWF2BCF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,&
            &IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
37       CONTINUE
      ENDIF
      IF(IUNIT(23).GT.0) THEN
         CALL GWF2LPF7BDS(KKSTP,KKPER,IGRID)
         CALL GWF2LPF7BDCH(KKSTP,KKPER,IGRID)
         IBDRET=0
         IC1=1
         IC2=NCOL
         IR1=1
         IR2=NROW
         IL1=1
         IL2=NLAY
         DO 157 IDIR=1,3
            CALL GWF2LPF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,&
            &IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
157      CONTINUE
      ENDIF
   ! ADDED CALLS FOR FLOW BUDGETS TO UPW PACKAGE.
      IF(IUNIT(62).GT.0) THEN
         CALL GWF2UPWBDS(KKSTP,KKPER,IGRID)
         CALL GWF2UPWBDCH(KKSTP,KKPER,IGRID)
         IBDRET=0
         IC1=1
         IC2=NCOL
         IR1=1
         IR2=NROW
         IL1=1
         IL2=NLAY
         DO 158 IDIR=1,3
            CALL GWF2UPWBDADJ(KKSTP,KKPER,IDIR,IBDRET,&
            &IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
158      CONTINUE
      ENDIF
      IF(IUNIT(37).GT.0) THEN
         CALL GWF2HUF7BDS(KKSTP,KKPER,IGRID)
         CALL GWF2HUF7BDCH(KKSTP,KKPER,IUNIT(47),IGRID)
         IBDRET=0
         IC1=1
         IC2=NCOL
         IR1=1
         IR2=NROW
         IL1=1
         IL2=NLAY
         DO 159 IDIR=1,3
            CALL GWF2HUF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,&
            &IC1,IC2,IR1,IR2,IL1,IL2,IUNIT(47),IGRID)
159      CONTINUE
      ENDIF
      IF(IUNIT(2).GT.0) CALL GWF2WEL7BD(KKSTP,KKPER,IUNIT(63),IGRID)
      IF(IUNIT(3).GT.0) CALL GWF2DRN7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(4).GT.0) CALL GWF2RIV7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(5).GT.0) THEN
         IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(&
         &0,IGRID)
         CALL GWF2EVT7BD(KKSTP,KKPER,IGRID)
         IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(&
         &1,IGRID)
      END IF
      IF(IUNIT(7).GT.0) CALL GWF2GHB7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(8).GT.0) THEN
         IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(&
         &0,IGRID)
   !           CALL GWF2RCH7BD(KKSTP,KKPER,IUNIT(44),IGRID)
         CALL GWF2RCH7BD(KKSTP,KKPER,IGRID)
         IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(&
         &1,IGRID)
      END IF
      IF(IUNIT(16).GT.0) CALL GWF2FHB7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(17).GT.0) CALL GWF2RES7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(18).GT.0) CALL GWF2STR7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(19).GT.0) CALL GWF2IBS7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(39).GT.0) CALL GWF2ETS7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(40).GT.0) CALL GWF2DRT7BD(KKSTP,KKPER,IGRID)
   ! (CJM) Added RCH unit number for RCH->SFR.
      IF(IUNIT(44).GT.0) CALL GWF2SFR7BD(KKSTP,KKPER,IUNIT(15),&
      &IUNIT(22),IUNIT(46),IUNIT(55),NSOL,&
      &IUNIT(8),IGRID)
      IF(IUNIT(22).GT.0) CALL GWF2LAK7BD(KKSTP,KKPER,IUNIT(15),&
      &IUNIT(46),IUNIT(44),IUNIT(55),NSOL,IGRID)
   ! Moved call to UZF1BD to follow SFR7BD for printing net recharge in UZF.
      IF(IUNIT(55).GT.0) CALL GWF2UZF1BD(KKSTP,KKPER,IUNIT(22),&
      &IUNIT(44),IGRID)
      IF(IUNIT(50).GT.0) CALL GWF2MNW27BD(KKSTP,KKPER,IUNIT(62),&
      &IGRID)
      IF(IUNIT(52).GT.0) CALL GWF2MNW17BD(NSTP(KPER),KKSTP,KKPER,&
      &IGRID)
      IF(IUNIT(54).GT.0) CALL GWF2SUB7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(57).GT.0) CALL GWF2SWT7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(64).GT.0) CALL GWF2SWR7BD(KKSTP,KKPER,IGRID)  !SWR - JDH
   !         IF(IUNIT(66).GT.0) CALL GWF2GFB7BD(KKSTP,KKPER,IGRID)
      IF(IUNIT(65).GT.0) CALL GWF2SWI2BD(KKSTP,KKPER,IGRID)  !SWI2 - JDH
      IF(IUNIT(66).GT.0) CALL GWF2AG7BD(KKSTP,KKPER,IUNIT(63))
   !LMT
   !LMT----CALL LINK-MT3DMS SUBROUTINES TO SAVE FLOW-TRANSPORT LINK FILE
   !LMT----FOR USE BY MT3DMS FOR TRANSPORT SIMULATION
   !LMT
      IF(IUNIT(49).GT.0) CALL LMT8BD(KKSTP,KKPER,IGRID)
   !LMT
   !
   !
   !  Set NWT heads to Hdry when head is below bottom.
      IF(IUNIT(63).GT.0)CALL GWF2NWT1BD()
   !  Observation simulated equivalents
      CALL OBS2BAS7SE(IUNIT(28),IGRID)
      IF(IUNIT(33).GT.0) CALL OBS2DRN7SE(IGRID)
      IF(IUNIT(34).GT.0) CALL OBS2RIV7SE(IGRID)
      IF(IUNIT(35).GT.0) CALL OBS2GHB7SE(IGRID)
      IF(IUNIT(36).GT.0) CALL OBS2STR7SE(IGRID)
      IF(IUNIT(38).GT.0) CALL OBS2CHD7SE(KKPER,IUNIT(62),IGRID)
      IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7SE(1,IGRID)
      IF(IUNIT(43).GT.0 .AND. IUNIT(19).GT.0)&
      &CALL GWF2HYD7IBS7SE(1,IGRID)
      IF(IUNIT(43).GT.0 .AND. IUNIT(54).GT.0)&
      &CALL GWF2HYD7SUB7SE(1,IGRID)
      IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)&
      &CALL GWF2HYD7STR7SE(1,IGRID)
      IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)&
      &CALL GWF2HYD7SFR7SE(1,IGRID)
   !
   !7C5---PRINT AND/OR SAVE DATA.
      CALL GWF2BAS7OT(KKSTP,KKPER,ICNVG,1,IGRID,BUDPERC)
      IF(IUNIT(19).GT.0) CALL GWF2IBS7OT(KKSTP,KKPER,IUNIT(19),&
      &IGRID)
      IF(IUNIT(37).GT.0)THEN
         IF(IOHUFHDS .NE.0 .OR.IOHUFFLWS .NE.0)&
         &CALL GWF2HUF7OT(KKSTP,KKPER,ICNVG,1,IGRID)
      ENDIF
      IF(IUNIT(51).NE.0) CALL GWF2MNW2I7OT(NSTP(KKPER),KKSTP,&
      &KKPER,IGRID)
      IF(IUNIT(54).GT.0) CALL GWF2SUB7OT(KKSTP,KKPER,IUNIT(54),&
      &IGRID)
      IF(IUNIT(57).GT.0) CALL GWF2SWT7OT(KKSTP,KKPER,IGRID)
      IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7OT(KKSTP,KKPER,IGRID)
   !
   !7C6---JUMP TO END OF PROGRAM IF CONVERGENCE WAS NOT ACHIEVED.
      IF ( IUNIT(63).GT.0 ) THEN
         IF(ICNVG.EQ.0) THEN
            NCVGERR=NCVGERR+1
            IF ( ICNVGFLG.EQ.0 ) THEN
               WRITE(IOUT,87) BUDPERC
               WRITE(IOUT,*) 'STOPPING SIMULATION'
               GO TO 110
            ELSE
               WRITE(IOUT,*) 'CONTINUING EXECUTION'
            END IF
         END IF
      ELSE
         IF(ICNVG.EQ.0) THEN
            NCVGERR=NCVGERR+1
            WRITE(IOUT,87) BUDPERC
87          FORMAT(1X,'FAILURE TO MEET SOLVER CONVERGENCE CRITERIA',/&
            &1X,'BUDGET PERCENT DISCREPANCY IS',F10.4)
            IF(ABS(BUDPERC).GT.STOPER) THEN
               WRITE(IOUT,*) 'STOPPING SIMULATION'
               GO TO 110
            ELSE
               WRITE(IOUT,*) 'CONTINUING EXECUTION'
            END IF
         END IF
      END IF
   !
   !-----END OF TIME STEP (KSTP) AND STRESS PERIOD (KPER) LOOPS
90 CONTINUE
100 CONTINUE
   !
   !
IF(IUNIT(52).NE.0) CALL GWF2MNW17OT(IGRID)
   !
   !8------END OF SIMULATION
   !-------SAVE RESTART RECORDS FOR SUB PACKAGE
110 IF(IUNIT(54).GT.0) CALL GWF2SUB7SV(IGRID)
   !
   !  Observation output
IF(IUNIT(28).GT.0) CALL OBS2BAS7OT(IUNIT(28),IGRID)
IF(IUNIT(33).GT.0) CALL OBS2DRN7OT(IGRID)
IF(IUNIT(34).GT.0) CALL OBS2RIV7OT(IGRID)
IF(IUNIT(35).GT.0) CALL OBS2GHB7OT(IGRID)
IF(IUNIT(36).GT.0) CALL OBS2STR7OT(IGRID)
IF(IUNIT(38).GT.0) CALL OBS2CHD7OT(IGRID)
CALL GLO1BAS6ET(IOUT,IBDT,1)
   !
   !-------OUTPUT RESULTS OF SWR TIMER
IF(IUNIT(64).GT.0) CALL GWF2SWR7OT(IGRID)
   !
   !9------CLOSE FILES AND DEALLOCATE MEMORY.  GWF2BAS7DA MUST BE CALLED
   !9------LAST BECAUSE IT DEALLOCATES IUNIT.
CALL SGWF2BAS7PNT(IGRID)
IF(IUNIT(1).GT.0) CALL GWF2BCF7DA(IGRID)
IF(IUNIT(2).GT.0) CALL GWF2WEL7DA(IGRID)
IF(IUNIT(3).GT.0) CALL GWF2DRN7DA(IGRID)
IF(IUNIT(4).GT.0) CALL GWF2RIV7DA(IGRID)
IF(IUNIT(5).GT.0) CALL GWF2EVT7DA(IGRID)
IF(IUNIT(7).GT.0) CALL GWF2GHB7DA(IGRID)
IF(IUNIT(8).GT.0) CALL GWF2RCH7DA(IGRID)
IF(IUNIT(9).GT.0) CALL SIP7DA(IGRID)
IF(IUNIT(10).GT.0) CALL DE47DA(IGRID)
IF(IUNIT(13).GT.0) CALL PCG7DA(IGRID)
   !      IF(IUNIT(59).GT.0) CALL PCGN2DA(IGRID)
IF(IUNIT(63).GT.0) THEN
   IF(LINMETH.EQ.1) THEN
      CALL GMRES7DA(IGRID)
   ELSEIF(LINMETH.EQ.2) THEN
      CALL XMD7DA(IGRID)
   !        ELSEIF(LINMETH.EQ.3) THEN
   !          CALL SAMG7DA(IGRID)
   END IF
   CALL GWF2NWT1DA(IGRID)
END IF
IF(IUNIT(62).GT.0) CALL GWF2UPW1DA(IGRID)
IF(IUNIT(16).GT.0) CALL GWF2FHB7DA(IGRID)
IF(IUNIT(17).GT.0) CALL GWF2RES7DA(IGRID)
IF(IUNIT(18).GT.0) CALL GWF2STR7DA(IGRID)
IF(IUNIT(19).GT.0) CALL GWF2IBS7DA(IGRID)
IF(IUNIT(20).GT.0) CALL GWF2CHD7DA(IGRID)
   !      IF(IUNIT(43).GT.0) CALL GWF2HYD7DA(IGRID)
IF(IUNIT(21).GT.0) CALL GWF2HFB7DA(IGRID)
IF(IUNIT(22).GT.0 .OR. IUNIT(44).GT.0)CALL GWF2LAK7DA(IUNIT(22),&
&IGRID)
IF(IUNIT(23).GT.0) CALL GWF2LPF7DA(IGRID)
IF(IUNIT(37).GT.0) CALL GWF2HUF7DA(IGRID)
IF(IUNIT(39).GT.0) CALL GWF2ETS7DA(IGRID)
IF(IUNIT(40).GT.0) CALL GWF2DRT7DA(IGRID)
IF(IUNIT(42).GT.0) CALL GMG7DA(IGRID)
IF(IUNIT(44).GT.0) CALL GWF2SFR7DA(IGRID)
IF(IUNIT(46).GT.0) CALL GWF2GAG7DA(IGRID)
IF(IUNIT(50).GT.0) CALL GWF2MNW27DA(IGRID)
IF(IUNIT(51).GT.0) CALL GWF2MNW2I7DA(IGRID)
IF(IUNIT(52).GT.0) CALL GWF2MNW17DA(IGRID)
IF(IUNIT(54).GT.0) CALL GWF2SUB7DA(IGRID)
IF(IUNIT(55).GT.0) CALL GWF2UZF1DA(IGRID)
IF(IUNIT(57).GT.0) CALL GWF2SWT7DA(IGRID)
IF(IUNIT(64).GT.0) CALL GWF2SWR7DA(IGRID)  !SWR - JDH
IF(IUNIT(65).GT.0) CALL GWF2SWI2DA(IGRID)  !SW12 - JDH
IF(IUNIT(66).GT.0) CALL GWF2AG7DA()
   !      IF(IUNIT(66).GT.0) CALL GWF2GFB7DA(IGRID)
CALL OBS2BAS7DA(IUNIT(28),IGRID)
IF(IUNIT(33).GT.0) CALL OBS2DRN7DA(IGRID)
IF(IUNIT(34).GT.0) CALL OBS2RIV7DA(IGRID)
IF(IUNIT(35).GT.0) CALL OBS2GHB7DA(IGRID)
IF(IUNIT(36).GT.0) CALL OBS2STR7DA(IGRID)
IF(IUNIT(38).GT.0) CALL OBS2CHD7DA(IGRID)
IF(IUNIT(43).GT.0) CALL GWF2HYD7DA(IGRID)
IF(IUNIT(49).GT.0) CALL LMT8DA(IGRID)
   !      IF(IUNIT(61).GT.0) CALL FMP2DA(IGRID)
CALL GWF2BAS7DA(IGRID)
   !
   !10-----END OF PROGRAM.
IF(NCVGERR.GT.0) THEN
   WRITE(*,*) 'FAILED TO MEET SOLVER CONVERGENCE CRITERIA ',&
   &NCVGERR,' TIME(S)'
ELSE
   WRITE(*,*) ' Normal termination of simulation'
END IF
CALL USTOP(' ')
   !
END
SUBROUTINE GETNAMFIL(FNAME)
   !     ******************************************************************
   !     GET THE NAME OF THE NAME FILE
   !     ******************************************************************
   !        SPECIFICATIONS:
   !
   !     ------------------------------------------------------------------
   CHARACTER*(*) FNAME
   CHARACTER*200 COMLIN
   LOGICAL EXISTS
   !     ------------------------------------------------------------------
   !
   ! Get name file from command line or user interaction.
   FNAME=' '
   COMLIN=' '
   ! *** Subroutines GETARG and GETCL are extensions to Fortran 90/95 that
   ! *** allow a program to retrieve command-line arguments.  To enable
   ! *** Modflow-2000 to read the name of a Name file from the command
   ! *** line, either GETARG or GETCL must be called, but not both.  As
   ! *** distributed, the call to GETARG is uncommented.  For compilers
   ! *** that support GETCL but not GETARG, comment out the call to GETARG
   ! *** and uncomment the call to GETCL.  The calls to both GETARG and
   ! *** GETCL may be commented out for compilers that do not support
   ! *** either extension.
   CALL GETARG(1,COMLIN)
   !        CALL GETCL(COMLIN)
   ICOL = 1
   IF(COMLIN.NE.' ') THEN
      FNAME=COMLIN
   ELSE
15    WRITE (*,*) ' Enter the name of the NAME FILE: '
      READ (*,'(A)') FNAME
      CALL URWORD(FNAME,ICOL,ISTART,ISTOP,0,N,R,0,0)
      FNAME=FNAME(ISTART:ISTOP)
      IF (FNAME.EQ.' ') GOTO 15
   ENDIF
   INQUIRE (FILE=FNAME,EXIST=EXISTS)
   IF(.NOT.EXISTS) THEN
      NC=INDEX(FNAME,' ')
      FNAME(NC:NC+3)='.nam'
      INQUIRE (FILE=FNAME,EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
         WRITE (*,480) FNAME(1:NC-1),FNAME(1:NC+3)
480      FORMAT(1X,'Can''t find name file ',A,' or ',A)
         CALL USTOP(' ')
      ENDIF
   ENDIF
   !
   RETURN
END
SUBROUTINE GLO1BAS6ET(IOUT,IBDT,IPRTIM)
   !     ******************************************************************
   !     Get end time and calculate elapsed time
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   INTEGER IBDT(8), IEDT(8), IDPM(12)
   DATA IDPM/31,28,31,30,31,30,31,31,30,31,30,31/ ! Days per month
   DATA NSPD/86400/  ! Seconds per day
   !     ------------------------------------------------------------------
   !
   !     Get current date and time, assign to IEDT, and write.
   CALL DATE_AND_TIME(VALUES=IEDT)
   WRITE(*,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)
1000 FORMAT(1X,'Run end date and time (yyyy/mm/dd hh:mm:ss): ',&
   &I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)
   IF(IPRTIM.GT.0) THEN
      WRITE(IOUT,'(1X)')
      WRITE(IOUT,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)
   END IF
   !
   !     Calculate elapsed time in days and seconds
   NDAYS=0
   LEAP=0
   IF (MOD(IEDT(1),4).EQ.0) LEAP = 1
   IBD = IBDT(3)            ! BEGIN DAY
   IED = IEDT(3)            ! END DAY
   !     FIND DAYS
   IF (IBDT(2).NE.IEDT(2)) THEN
   !       MONTHS DIFFER
      MB = IBDT(2)             ! BEGIN MONTH
      ME = IEDT(2)             ! END MONTH
      NM = ME-MB+1             ! NUMBER OF MONTHS TO LOOK AT
      IF (MB.GT.ME) NM = NM+12
      MC=MB-1
      DO 10 M=1,NM
         MC=MC+1                ! MC IS CURRENT MONTH
         IF (MC.EQ.13) MC = 1
         IF (MC.EQ.MB) THEN
            NDAYS = NDAYS+IDPM(MC)-IBD
            IF (MC.EQ.2) NDAYS = NDAYS + LEAP
         ELSEIF (MC.EQ.ME) THEN
            NDAYS = NDAYS+IED
         ELSE
            NDAYS = NDAYS+IDPM(MC)
            IF (MC.EQ.2) NDAYS = NDAYS + LEAP
         ENDIF
10    CONTINUE
   ELSEIF (IBD.LT.IED) THEN
   !       START AND END IN SAME MONTH, ONLY ACCOUNT FOR DAYS
      NDAYS = IED-IBD
   ENDIF
   ELSEC=NDAYS*NSPD
   !
   !     ADD OR SUBTRACT SECONDS
   ELSEC = ELSEC+(IEDT(5)-IBDT(5))*3600.0
   ELSEC = ELSEC+(IEDT(6)-IBDT(6))*60.0
   ELSEC = ELSEC+(IEDT(7)-IBDT(7))
   ELSEC = ELSEC+(IEDT(8)-IBDT(8))*0.001
   !
   !     CONVERT SECONDS TO DAYS, HOURS, MINUTES, AND SECONDS
   NDAYS = ELSEC/NSPD
   RSECS = MOD(ELSEC,86400.0)
   NHOURS = RSECS/3600.0
   RSECS = MOD(RSECS,3600.0)
   NMINS = RSECS/60.0
   RSECS = MOD(RSECS,60.0)
   NSECS = RSECS
   RSECS = MOD(RSECS,1.0)
   MSECS = NINT(RSECS*1000.0)
   NRSECS = NSECS
   IF (RSECS.GE.0.5) NRSECS=NRSECS+1
   !
   !     Write elapsed time to screen
   IF (NDAYS.GT.0) THEN
      WRITE(*,1010) NDAYS,NHOURS,NMINS,NRSECS
1010  FORMAT(1X,'Elapsed run time: ',I3,' Days, ',I2,' Hours, ',I2,&
      &' Minutes, ',I2,' Seconds',/)
   ELSEIF (NHOURS.GT.0) THEN
      WRITE(*,1020) NHOURS,NMINS,NRSECS
1020  FORMAT(1X,'Elapsed run time: ',I2,' Hours, ',I2,&
      &' Minutes, ',I2,' Seconds',/)
   ELSEIF (NMINS.GT.0) THEN
      WRITE(*,1030) NMINS,NSECS,MSECS
1030  FORMAT(1X,'Elapsed run time: ',I2,' Minutes, ',&
      &I2,'.',I3.3,' Seconds',/)
   ELSE
      WRITE(*,1040) NSECS,MSECS
1040  FORMAT(1X,'Elapsed run time: ',I2,'.',I3.3,' Seconds',/)
   ENDIF
   !
   !     Write times to file if requested
   IF(IPRTIM.GT.0) THEN
      IF (NDAYS.GT.0) THEN
         WRITE(IOUT,1010) NDAYS,NHOURS,NMINS,NRSECS
      ELSEIF (NHOURS.GT.0) THEN
         WRITE(IOUT,1020) NHOURS,NMINS,NRSECS
      ELSEIF (NMINS.GT.0) THEN
         WRITE(IOUT,1030) NMINS,NSECS,MSECS
      ELSE
         WRITE(IOUT,1040) NSECS,MSECS
      ENDIF
   ENDIF
   !
   RETURN
END
