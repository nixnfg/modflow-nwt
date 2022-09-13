MODULE OBSBASMODULE
   INTEGER, SAVE,POINTER ::ITS,NH,MAXM,MOBS,IUHOBSV,IDRY,JDRY
   INTEGER, SAVE,POINTER ::IPRT, OBSTART
   REAL,    SAVE,POINTER ::HOBDRY
   INTEGER, SAVE, DIMENSION(:,:), POINTER ::NDER
   INTEGER, SAVE, DIMENSION(:,:), POINTER ::MLAY
   INTEGER, SAVE, DIMENSION(:),   POINTER ::IOFF
   INTEGER, SAVE, DIMENSION(:),   POINTER ::JOFF
   INTEGER, SAVE, DIMENSION(:),   POINTER ::IHOBWET
   REAL,    SAVE, DIMENSION(:),   POINTER ::H
   REAL,    SAVE, DIMENSION(:),   POINTER ::HOBS
   REAL,    SAVE, DIMENSION(:),   POINTER ::TOFF
   REAL,    SAVE, DIMENSION(:),   POINTER ::ROFF
   REAL,    SAVE, DIMENSION(:),   POINTER ::COFF
   REAL,    SAVE, DIMENSION(:),   POINTER ::OTIME
   REAL,    SAVE, DIMENSION(:,:), POINTER ::PR
   REAL,    SAVE, DIMENSION(:,:), POINTER ::RINT
   CHARACTER*12,SAVE,DIMENSION(:),POINTER ::OBSNAM
   TYPE OBSBASTYPE
      INTEGER,  POINTER    ::ITS,NH,MAXM,MOBS,IUHOBSV,IDRY,JDRY
      INTEGER,  POINTER    ::IPRT, OBSTART
      REAL,     POINTER    ::HOBDRY
      INTEGER,  DIMENSION(:,:), POINTER ::NDER
      INTEGER,  DIMENSION(:,:), POINTER ::MLAY
      INTEGER,  DIMENSION(:),   POINTER ::IOFF
      INTEGER,  DIMENSION(:),   POINTER ::JOFF
      INTEGER,  DIMENSION(:),   POINTER ::IHOBWET
      REAL,     DIMENSION(:),   POINTER ::H
      REAL,     DIMENSION(:),   POINTER ::HOBS
      REAL,     DIMENSION(:),   POINTER ::TOFF
      REAL,     DIMENSION(:),   POINTER ::ROFF
      REAL,     DIMENSION(:),   POINTER ::COFF
      REAL,     DIMENSION(:),   POINTER ::OTIME
      REAL,     DIMENSION(:,:), POINTER ::PR
      REAL,     DIMENSION(:,:), POINTER ::RINT
      CHARACTER*12,DIMENSION(:),POINTER ::OBSNAM
   END TYPE
   TYPE(OBSBASTYPE), SAVE :: OBSBASDAT(10)
END MODULE OBSBASMODULE
   !  NDER(1,n) -- Observation layer
   !  NDER(2,n) -- Observation row
   !  NDER(3,n) -- Observation column
   !  NDER(4,n) -- Observation time step
   !  NDER(5,n) -- Observation number for computing observation as a head change
   !  MLAY(MAXM,MOBS) -- Layer numbers for multilayer observations
   !  IOFF(NH) -- Row offset for neighboring cell for interpolation
   !  JOFF(NH) -- Column offset for neighboring cell for interpolation
   !  IHOBWET(NH) -- Flag for observation -- 1 for wet, and -1 for dry
   !  H(NH) -- Simulated value
   !  HOBS(NH) -- Observed value
   !  TOFF(NH) -- Fractional offset between time steps
   !  ROFF(NH) -- Fractional offset from center of cell in Y direction (between rows)
   !  COFF(NH) -- Fractional offset from center of cell in X direction (between columns)
   !  OTIME(NH) -- Observation time in model time units
   !  PR(MAXM,MOBS) -- Fractional value for each layer of multilayer observations
   !  RINT(4,NH) -- Interpolation coefficients for the 4 nodes surrounding observation
   !  OBSNAM(NH) -- Observation name


SUBROUTINE OBS2BAS7AR(IUHDOB,IGRID)
   !     ******************************************************************
   !     INITIALIZE AND READ VARIABLES FOR HEAD OBSERVATIONS
   !     ******************************************************************
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY: NCOL,NROW,NLAY,DELR,DELC,&
   &NPER,NSTP,PERLEN,TSMULT,ISSFLG,IOUT,ITRSS
   !      USE GSFMODFLOW, ONLY : Modflow_skip_time_step
   USE OBSBASMODULE
   !
   CHARACTER*200 LINE
   !     ------------------------------------------------------------------
   !
   !1------ALLOCATE AND INITIALIZE TIME STEP COUNTER FOR USE BY ANY
   !1------OBSERVATION PACKAGE.
   ALLOCATE(ITS,OBSTART)
   ITS=0.0
   OBSTART = 0
   IF(IUHDOB.LE.0) GO TO 700
   !
   !2------ALLOCATE OTHER SCALARS IF HEAD OBSERVATIONS ARE BEING SPECIFIED.
   ALLOCATE(NH,MAXM,MOBS,IUHOBSV,IDRY,JDRY,IPRT)
   ALLOCATE(HOBDRY)
   !
   !3------DEFINE CONSTANTS
   IDRY = 0
   JDRY = 0
   ZERO=0.0
   ONEN=-1.0
   ML=0
   IERR=0
   !
   !4------WRITE VERSION.
   WRITE (IOUT,14) IUHDOB
14 FORMAT (/,' OBS2BAS7 -- HEAD OBSERVATIONS, ',&
   &'VERSION 2.0, 2/28/2006',/,' INPUT READ FROM UNIT ',I3)
   !
   !5------READ & PRINT ITEM 1 OF THE HOB INPUT FILE
   CALL URDCOM(IUHDOB,IOUT,LINE)
   LLOC = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NH,DUM,IOUT,IUHDOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MOBS,DUM,IOUT,IUHDOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MAXM,DUM,IOUT,IUHDOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUHOBSV,DUM,IOUT,IUHDOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,HOBDRY,IOUT,IUHDOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,IDUM,DUM,IOUT,IUHDOB)
   IPRT=1
   IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
      IPRT=0
      WRITE(IOUT,*) 'NOPRINT option for HEAD OBSERVATIONS'
   END IF
   IF (MAXM.EQ.1) THEN
      WRITE (IOUT,17)
17    FORMAT (/,' MAXM CAN NOT EQUAL 1 -- STOP EXECUTION')
      CALL USTOP(' ')
   ENDIF
   WRITE (IOUT,19) NH, MOBS, MAXM
19 FORMAT (/,&
   &' NUMBER OF HEADS....................................:',I5,/,&
   &'   NUMBER OF MULTILAYER HEADS.......................:',I5,/,&
   &'   MAXIMUM NUMBER OF LAYERS FOR MULTILAYER HEADS....:',I5)
   IF(NH.LE.0) THEN
      WRITE(IOUT,*) ' NH LESS THAN OR EQUAL TO 0'
      CALL USTOP(' ')
   END IF
   IF(IUHOBSV.GT.0) THEN
      WRITE(IOUT,21) IUHOBSV
21    FORMAT(1X,&
      &'HEAD OBSERVATIONS WILL BE SAVED ON UNIT............:',I5)
   ELSE
      WRITE(IOUT,22)
22    FORMAT(1X,'HEAD OBSERVATIONS WILL NOT BE SAVED IN A FILE')
   END IF
   WRITE(IOUT,23) HOBDRY
23 FORMAT(1X,'SIMULATED EQUIVALENT HEAD AT DRY CELLS WILL BE:',&
   &1P,1G15.6)
   !
   !6------ALLOCATE ARRAY DATA.
   ALLOCATE(NDER(5,NH))
   ALLOCATE(IOFF(NH))
   ALLOCATE(JOFF(NH))
   ALLOCATE(IHOBWET(NH))
   ALLOCATE(OBSNAM(NH))
   ALLOCATE(H(NH))
   ALLOCATE(HOBS(NH))
   ALLOCATE(TOFF(NH))
   ALLOCATE(OTIME(NH))
   ALLOCATE(ROFF(NH))
   ALLOCATE(COFF(NH))
   ALLOCATE(RINT(4,NH))
   IF(MOBS.GT.0 .AND. MAXM.GT.1) THEN
      ALLOCATE(MLAY(MAXM,MOBS))
      ALLOCATE(PR(MAXM,MOBS))
   ELSE
      ALLOCATE(MLAY(1,1))
      ALLOCATE(PR(1,1))
   END IF
   !
   !7------INITIALIZE OTIME, SIMULATED EQUIVALENT HEAD, NDER(5,n), and IHOBWET.
   DO 50 N = 1, NH
      OTIME(N) = ONEN
      H(N) = ZERO
      NDER(5,N)=0
      IHOBWET(N)=1
50 CONTINUE
   !
   !8------READ ITEM 2
   READ (IUHDOB,*) TOMULTH
   !
   !9------WRITE ITEM 2 AND TITLE FOR OBSERVATION TIMES.
   IF(IPRT.NE.0) WRITE (IOUT,530) TOMULTH
530 FORMAT (/,' OBSERVED HEAD DATA -- TIME OFFSETS ARE',&
   &' MULTIPLIED BY: ',G12.5,//,&
   &20X,'REFER.',/,&
   &7X,'OBSERVATION',2X,'STRESS',4X,'TIME',/,&
   &2X,'OBS#    NAME',6X,'PERIOD',3X,'OFFSET    OBSERVATION')
   !
   !10-----INITIALIZE N, WHICH IS THE COUNT OF THE NUMBER OF OBSERVATIONS
   !10-----AS THEY ARE READ.
   N=0
   !
   !11-----READ NAME, LOCATION, TIME, AND OBSERVED VALUE (ITEM 3)
60 N=N+1
   READ (IUHDOB,*) OBSNAM(N), (NDER(I,N),I=1,3), IREFSP, TOFFSET,&
   &ROFF(N), COFF(N), HOBS(N)
   IF(IPRT.NE.0) WRITE (IOUT,535) N,OBSNAM(N),IREFSP,TOFFSET,HOBS(N)
535 FORMAT (1X,I5,1X,A12,2X,I4,2X,G11.4,1X,G11.4)
   !
   !11A----FOR SINGLE-TIME OBSERVATION (IREFSP>0), CALL UOBSTI TO DETERMINE
   !11A----WHEN OBSERVATION OCCURS.
   IF (IREFSP.GE.0) THEN
      CALL UOBSTI(OBSNAM(N),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,&
      &NDER(4,N),PERLEN,TOFF(N),TOFFSET,TOMULTH,TSMULT,&
      &0,OTIME(N))
   END IF
   !
   !12-----CHECK ROW AND COLUMN LOCATION.
   I = NDER(2,N)
   J = NDER(3,N)
   IF (J.LE.0 .OR. J.GT.NCOL .OR. I.LE.0 .OR. I.GT.NROW) THEN
      WRITE (IOUT,550) N
550   FORMAT (' FOR OBS',I5,' ROW OR COLUMN NUMBER INVALID -- ',&
      &'STOP EXECUTION (OBS2BAS7HRP)',/)
      IERR = 1
   ENDIF
   !
   !13-----Check if multi-layer
   IF(NDER(1,N).GE.0) THEN
   !
   !13A----SINGLE LAYER -- CHECK FOR VALID LAYER.
      IF (NDER(1,N).LE.0 .OR. NDER(1,N).GT.NLAY) THEN
         WRITE (IOUT,265) NDER(1,N)
265      FORMAT (' FOR OBS ',I5,' LAYER INVALID -- STOP EXECUTION',&
         &' (OBS2BAS7AR)',/)
         IERR = 1
      END IF
   ELSE
   !
   !13B----MULTI-LAYER -- CHECK LIMITS AND READ THE LAYERS AND PROPORTIONS.
      NL=-NDER(1,N)
      ML = ML + 1
      IF(ML.GT.MOBS) THEN
         WRITE (IOUT,565)
565      FORMAT (/,' NUMBER OF MULTILAYER OBSERVATIONS EXCEEDS MOBS -- ',&
         &'STOP EXECUTION (OBS2BAS7AR)',/)
         CALL USTOP(' ')
      END IF
      IF(NL.GT.MAXM) THEN
         WRITE(IOUT,620) NL
620      FORMAT(/,1X,'ERROR: VALUE ENTERED FOR MAXM IN HOB FILE IS',&
         &' SMALLER THAN THE MAXIMUM NUMBER',/,' OF LAYERS',&
         &' IN A MULTILAYER HEAD OBSERVATION, WHICH IS ',&
         &I3,' -- INCREASE MAXM.')
         CALL USTOP(' ')
      END IF
      DO 268 M=1,MAXM
         MLAY(M,ML) = 0
268   CONTINUE
      READ (IUHDOB,*) (MLAY(M,ML),PR(M,ML),M=1,NL)
      IF(IPRT.NE.0) WRITE(IOUT,540) (MLAY(M,ML),PR(M,ML),M=1,NL)
540   FORMAT (5X,'MULTIPLE LAYERS AND PROPORTIONS :',5(I5,',',F5.2,3X))
   !
   !
   !13C----CHECK LAYER NUMBERS AND ADD PROPORTIONS FOR MULTILAYER
   !13C----OBSERVATION WELLS.
      TPR=ZERO
      DO 270 K = 1,NL
         KK = MLAY(K,ML)
         TPR = TPR + PR(K,ML)
         IF (KK.LE.0 .OR. KK.GT.NLAY) THEN
            WRITE (IOUT,265) N
            IERR = 1
         ENDIF
270   CONTINUE
   !
   !13D----CHECK SUM OF PROPORTIONS FOR MULTILAYER OBS WELLS
      IF (ABS(1.-TPR).GT..02) THEN
         WRITE (IOUT,560) N
560      FORMAT (/,' FOR OBS',I5,' MULTILAYER PROPORTIONS DO NOT SUM ',&
         &'TO 1.0 -- STOP EXECUTION (OBS2BAS7AR)',/)
         IERR = 1
      ENDIF
   END IF
   !
   !14-----CALCULATE INTERPOLATION COEFFICIENTS FOR THE LOCATION.
   CALL SOBS2BAS7HIA(N,ML)
   !
   !15-----CHECK FOR MULTI-TIME
   NT=-IREFSP
   IF(NT.GT.0) THEN
   !
   !15A----READ FLAG FOR USING TEMPORAL CHANGES IN HEAD (ITEM 5)
      READ (IUHDOB,*) ITT
      IF(IPRT.NE.0) WRITE (IOUT,515) ITT
515   FORMAT (2X,'TRANSIENT DATA AT THIS LOCATION, ITT =',I4)
      IF (ITT.NE.1 .AND. ITT.NE.2) THEN
         WRITE (IOUT,575) N
575      FORMAT (' FOR OBS',I5,&
         &' ITT MUST = 1 OR 2 -- STOP EXECUTION (OBS2BAS7HRP)',/)
         CALL USTOP(' ')
      ENDIF
   !
   !15B----LOOP THROUGH THE TIMES
      NBASE=N
      DO 200 J=1,NT
         IF(J.NE.1) THEN
   !
   !15B1---DUPLICATE THE LOCATION INFORMATION FOR THE OBSERVATIONS AT THE
   !15B1---SAME LOCATION.
            N=N+1
            IF(N.GT.NH) THEN
               WRITE(IOUT,127)
127            FORMAT(1X,/,1X,'ABORTING BECAUSE THERE ARE MORE HEAD',&
               &' OBSERVATIONS THAN SPECIFIED BY NH')
               CALL USTOP(' ')
            END IF
            DO 140 I = 1, 3
               NDER(I,N) = NDER(I,N-1)
140         CONTINUE
            ROFF(N) = ROFF(N-1)
            COFF(N) = COFF(N-1)
            IOFF(N) = IOFF(N-1)
            JOFF(N) = JOFF(N-1)
            DO 150 I = 1, 4
               RINT(I,N) = RINT(I,N-1)
150         CONTINUE
            IF (NDER(1,N-1).LT.0) THEN
               ML1 = ML
               ML = ML + 1
               DO 160 M = 1, MAXM
                  PR(M,ML) = PR(M,ML1)
                  MLAY(M,ML) = MLAY(M,ML1)
160            CONTINUE
            ENDIF
         END IF
   !
   !15B2---READ ONE MULTI-TIME OBSERVATION.
         READ (IUHDOB,*) OBSNAM(N), IREFSP, TOFFSET, HOBS(N)
         IF(ITT.EQ.2 .AND. J.NE.1) THEN
            HOBS(N)=HOBS(N)-HOBS(NBASE)
            NDER(5,N)=NBASE
         END IF
   !
   !15B3---WRITE ONE MULTI-TIME OBSERVATION.
         IF(IPRT.NE.0) WRITE (IOUT,535) N,OBSNAM(N),IREFSP,TOFFSET,&
         &HOBS(N)
         CALL UOBSTI(OBSNAM(N),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,&
         &NDER(4,N),PERLEN,TOFF(N),TOFFSET,TOMULTH,&
         &TSMULT,0,OTIME(N))
         IF(J.EQ.NT .AND. IPRT.NE.0) WRITE (IOUT,570)
570      FORMAT (' ')
200   CONTINUE
   END IF
   !
   !16-----READ ANOTHER OBSERVATION (ITEM 3) IF THERE ARE STILL MORE OBSERVATIONS.
   IF(N.LT.NH) GO TO 60
   !
   !17-----DONE READING HEAD OBSERVATIONS.
   !17-----PRINT TABLE SHOWING LOCATION OF OBSERVATIONS.
   IF(IPRT.NE.0) THEN
      WRITE (IOUT,590)
590   FORMAT (/,53X,'HEAD CHANGE',/,54X,'REFERENCE',/,&
      &8X,'OBSERVATION',19X,'ROW',5X,'COL    OBSERVATION',/,&
      &2X,'OBS#',5X,'NAME',7X,'LAY  ROW  COL  OFFSET  OFFSET',3X,&
      &'(IF > 0)')
      DO 450 N = 1, NH
         WRITE (IOUT,600) N, OBSNAM(N), (NDER(I,N),I=1,3), ROFF(N),&
         &COFF(N), NDER(5,N)
600      FORMAT (1X,I5,2X,A12,2X,I3,2(1X,I4),2(2X,F6.3),3X,I6)
450   CONTINUE
   END IF
   !
   !18-----IF ERROR OCCURRED ABOVE, PRINT MESSAGE AND STOP.
   IF (IERR.GT.0) THEN
      WRITE(IOUT,610)
610   FORMAT(/,1X,'ERROR: SEARCH ABOVE FOR ERROR MESSAGE(S)',/,1X,&
      &'STOP EXECUTION -- (OBS2BAS7AR)')
      CALL USTOP(' ')
   ENDIF
   !
   !19-----RETURN.
700 CALL SOBS2BAS7PSV(IUHDOB,IGRID)
   RETURN
END
SUBROUTINE OBS2BAS7SE(IUHDOB,IGRID)
   !     ******************************************************************
   !     INTERPOLATE HEADS.  ACCOUNT FOR DRY CELLS, IF NEEDED.
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY: NCOL,NROW,NLAY,DELR,DELC,IBOUND,HNEW,&
   &NPER,NSTP,PERLEN,TSMULT,ISSFLG,IOUT
   USE OBSBASMODULE
   DOUBLE PRECISION V
   !     ------------------------------------------------------------------
   CALL SOBS2BAS7PNT(IUHDOB,IGRID)
   !
   !1------INCREMENT TIME STEP COUNTER, WHICH IS USED BY OTHER OBSERVATION
   !1------PACKAGES.  RETURN IF THERE ARE NO HEAD OBSERVATIONS.
   ITS=ITS+1
   IF(IUHDOB.LE.0) RETURN
   ZERO = 0.0
   !
   !2------CHECK FOR NODES USED TO INTERPOLATE HEADS THAT HAVE GONE DRY OR
   !2------ARE OTHERWISE INACTIVE.
   !2------ELIMINATE OBSERVATIONS OR RECALC. INTERPOLATION COEFFICIENTS.
   !2 -----CHECK FOR OBSERVATIONS THAT NEED TO BE OMITTED OR NEED TO HAVE
   !2------THE INTERPOLATION RECALCULATED
   !2------IDRY = # OBS OMITTED; JDRY = # INTERPOLATIONS CHANGED
   ML = 0
   DO 30 N = 1, NH
      K = NDER(1,N)
      MM = 1
      IF (K.LT.0) THEN
         ML = ML + 1
         MM = -K
      ENDIF
      IF((IHOBWET(N).LT. 0) .OR.&
      &(ITS.NE.NDER(4,N) .AND. ITS.NE.NDER(4,N)+1)) GO TO 30
      II = NDER(2,N)
      JJ = NDER(3,N)
   !  Moved the following to inside the DO 20 loop so updates of IOFF and
   !  JOFF by SOBS2BAS7HIB will be applied.
   !        IO = IOFF(N)
   !        JO = JOFF(N)
   !
   !4------CHECK FOR DRY OBSERVATIONS OR INTERPOLATIONS AFFECTED BY DRY
   !4------CELLS
      DO 20 M = 1, MM
         IO = IOFF(N)
         JO = JOFF(N)
         KK = K
         IF (K.LT.0) KK = MLAY(M,ML)
         IF (KK.EQ.0) GOTO 30
         IF (IBOUND(JJ,II,KK).EQ.0) THEN
            IDRY = IDRY + 1
            IHOBWET(N)=-1
            WRITE (IOUT,495) N, OBSNAM(N)
495         FORMAT (/,' HEAD OBS#',I5,', ID ',A,' IS DRY -- OMIT',&
            &' (OBS2BAS7SE)')
            GOTO 30
   !
   !5------CHECK TO SEE IF A CELL USED IN INTERPOLATION IS INACTIVE
         ELSEIF ((RINT(2,N).NE.ZERO.AND.IBOUND(JJ+JO,II,KK)&
         &.EQ.0) .OR.&
         &(RINT(3,N).NE.ZERO.AND.IBOUND(JJ,II+IO,KK)&
         &.EQ.0) .OR.&
         &(RINT(4,N).NE.ZERO.AND.IBOUND(JJ+JO,II+IO,KK)&
         &.EQ.0)) THEN
            IF(M.GT.1) THEN
               IDRY = IDRY + 1
               IHOBWET(N)=-1
   !              WRITE (IOUT,500) N, OBSNAM(N)
   !  500 FORMAT (/,' HEAD OBS#',I5,', ID ',A,
   !     &' OMITTED BECAUSE IBOUND=0 FOR CELL(S)',/,' REQUIRED FOR',
   !     &' MULTILAYER INTERPOLATION (OBS2BAS7SE)')
               WRITE (IOUT,500) N, OBSNAM(N),KK
500            FORMAT (/,' HEAD OBS#',I5,', ID ',A,&
               &' OMITTED BECAUSE IBOUND=0 FOR CELL(S) IN MODEL LAYER',I10,/,&
               &' REQUIRED FOR MULTILAYER INTERPOLATION (OBS2BAS7SE)')
   !
               GOTO 30
            ENDIF
            WRITE (IOUT,505) N, OBSNAM(N)
505         FORMAT (/,' INTERPOLATION FOR HEAD OBS#',I5,', ID ',A,' CHANGED',&
            &' BECAUSE AT LEAST ONE',/,&
            &' NEIGHBORING CELL REQUIRED FOR INTERPOLATION IS DRY',&
            &' OR INACTIVE (OBS2BAS7SE)')
            MLL = 0
            IF (NDER(1,N).LT.0) MLL = MLAY(1,ML)
            CALL SOBS2BAS7HIB(NDER(:,N),COFF(N),ROFF(N),DELR,DELC,&
            &IBOUND,NCOL,NROW,NLAY,RINT(:,N),&
            &JOFF(N),IOFF(N),MLL)
            JDRY = JDRY + 1
   !
   !6------COULD INSERT ELSEIF TO SEE IF A NEIGHBORING CELL HAS REWET, IF SO,
   !6------RECALCULATE RINT
         ENDIF
20    CONTINUE
30 CONTINUE
   !
   !7------INTERPOLATION
   ML = 0
   DO 60 N = 1, NH
   !
   !8------UPDATE COUNTER FOR MULTILAYER WELLS
      K = NDER(1,N)
      MM = 1
      IF (K.LT.0) THEN
         ML = ML + 1
         MM = -K
      ENDIF
   !
   !9------OBSERVATION AT THIS TIME STEP?
      IF((IHOBWET(N).LT.0) .OR.&
      &(NDER(4,N).NE.ITS .AND. NDER(4,N)+1.NE.ITS)) GO TO 60
   !
      II = NDER(2,N)
      JJ = NDER(3,N)
      IO = IOFF(N)
      JO = JOFF(N)
      V = 0.0
      DO 40 M = 1, MM
         KK = K
         PROP = 1.
         IF (K.LT.0) THEN
            KK = MLAY(M,ML)
            PROP = PR(M,ML)
         ENDIF
         IF (KK.EQ.0) GOTO 50
   !
   !10-----CALCULATE CONTRIBUTION FROM THIS LAYER TO HEADS
         V = V + PROP*(RINT(1,N)*HNEW(JJ,II,KK)+&
         &RINT(2,N)*HNEW(JJ+JO,II,KK)+&
         &RINT(3,N)*HNEW(JJ,II+IO,KK)+&
         &RINT(4,N)*HNEW(JJ+JO,II+IO,KK))
40    CONTINUE
   !
   !11-----INDEX WHICH, IF NOT ZERO, IDENTIFIES THE HEAD USED TO
   !11-----CALCULATE DRAWDOWN
50    N1 = NDER(5,N)
   !
   !12-----INTERPOLATE OVER TIME AND COMPUTE DRAWDOWN IF INDICATED.
      IF(ITS.EQ.NDER(4,N)) THEN
         H(N)=V
   !
      ELSE IF(NDER(4,N)+1.EQ.ITS) THEN
         IF(ITS.EQ.1) THEN
   !  For observations in first time, H(N) will not have been initialized
   !  because this routine is not called at the beginning of the simulation.
   !  Set H(N)=V because observations in the first time step must be at
   !  the end -- i.e. no interpolation between time steps 0 and 1.
            H(N)=V
         ELSE
            H(N) = H(N) + TOFF(N)*(V-H(N))
         END IF
         IF(N1.GT.0) THEN
            IF(IHOBWET(N1).LT.0) THEN
               IHOBWET(N)=-1
            ELSE
               H(N) = H(N) - H(N1)
            END IF
         END IF
      END IF
60 CONTINUE
   !
   !13------RETURN.
   RETURN
END
SUBROUTINE OBS2BAS7OT(IUHDOB,IGRID)
   !     ******************************************************************
   !     WRITE HEAD OBSERVATIONS
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY: IOUT
   USE OBSBASMODULE
   DOUBLE PRECISION  SQ,SUMSQ
   !     ------------------------------------------------------------------
   CALL SOBS2BAS7PNT(IUHDOB,IGRID)
   !
   !1------WRITE OBSERVATIONS TO LISTING FILE.
   IF(IPRT.NE.0) WRITE(IOUT,17)
17 FORMAT(1X,/,1X,'HEAD AND DRAWDOWN OBSERVATIONS',/,&
   &1X,'OBSERVATION       OBSERVED           SIMULATED',/&
   &1X,'  NAME              VALUE              VALUE',&
   &'             DIFFERENCE',/&
   &1X,'-----------------------------------------------',&
   &'---------------------')
   SUMSQ=0.
   DO 100 N=1,NH
      IF(IHOBWET(N).LT.0) THEN
         H(N)=HOBDRY
         IF(IPRT.NE.0) WRITE(IOUT,27) OBSNAM(N),HOBS(N),H(N)
      ELSE
         DIFF=HOBS(N)-H(N)
         SQ=DIFF*DIFF
         SUMSQ=SUMSQ+SQ
         IF(IPRT.NE.0) WRITE(IOUT,27) OBSNAM(N),HOBS(N),H(N),DIFF
      END IF
27    FORMAT(1X,A,1P,3G20.11)
100 CONTINUE
   WRITE(IOUT,28) SUMSQ
28 FORMAT(1X,/,1X,'HEAD/DRAWDOWN SUM OF SQUARED DIFFERENCE:',&
   &1P,E15.5)
   !
   !2------WRITE OBSERVATIONS TO SEPARATE FILE.
   IF(IUHOBSV.GT.0) CALL UOBSSV(IUHOBSV,NH,H,HOBS,OBSNAM,1)
   !
   !3------RETURN.
   RETURN
END
SUBROUTINE SOBS2BAS7HIA(N,ML)
   !     ******************************************************************
   !     CALCULATE INTERPOLATION COEFFICIENTS FOR LOCATING OBSERVED HEADS
   !     ASSUMING ALL CELLS ARE ACTIVE.
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY: NCOL,NROW,NLAY,DELR,DELC,IBOUND,HNEW,STRT,&
   &NPER,NSTP,PERLEN,TSMULT,ISSFLG,IOUT
   USE OBSBASMODULE
   !     ------------------------------------------------------------------
   !
   K = NDER(1,N)
   IF (K.LT.0) K = MLAY(1,ML)
   I = NDER(2,N)
   J = NDER(3,N)
   I1 = I + 1
   J1 = J + 1
   IOFF(N) = 1
   JOFF(N) = 1
   IF (ROFF(N).LT.0.) THEN
      I1 = I - 1
      IOFF(N) = -1
   ENDIF
   IF (COFF(N).LT.0.) THEN
      J1 = J - 1
      JOFF(N) = -1
   ENDIF
   IF (I1.GE.1 .AND. I1.LE.NROW) IBI = 1
   IF (J1.GE.1 .AND. J1.LE.NCOL) IBJ = 1
   IF (I1.GE.1 .AND. I1.LE.NROW .AND. J1.GE.1 .AND. J1.LE.NCOL)&
   &IBIJ = 1
   IF (I1.LT.1 .OR. I1.GT.NROW) THEN
      ROFF(N) = 0.
      IBI = 0
   ENDIF
   IF (J1.LT.1 .OR. J1.GT.NCOL) THEN
      COFF(N) = 0.
      IBJ = 0
   ENDIF
   IF (I1.LT.1 .OR. I1.GT.NROW .OR. J1.LT.1 .OR. J1.GT.NCOL) IBIJ = 0
   !
   CALL SOBS2BAS7HBF(COFF(N),DELC,DELR,I,I1,IBI,IBIJ,IBJ,IOFF(N),&
   &J,J1,JOFF(N),NCOL,NROW,RINT(:,N),ROFF(N))
   !
   RETURN
END
SUBROUTINE SOBS2BAS7HIB(NDER,COFF,ROFF,DELR,DELC,IBOUND,NCOL,NROW,&
&NLAY,RINT,JOFF,IOFF,MLAY)
   !     ******************************************************************
   !     CALCULATE INTERPOLATION COEFFICIENTS FOR LOCATING OBSERVED HEADS
   !     USING CURRENT IBOUND VALUES.
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   DIMENSION NDER(5), DELR(NCOL), DELC(NROW), IBOUND(NCOL,NROW,NLAY),&
   &RINT(4)
   !     ------------------------------------------------------------------
   !
   K = NDER(1)
   IF (K.LT.0) K = MLAY
   I = NDER(2)
   J = NDER(3)
   I1 = I + 1
   J1 = J + 1
   IOFF = 1
   JOFF = 1
   IF (ROFF.LT.0.) THEN
      I1 = I - 1
      IOFF = -1
   ENDIF
   IF (COFF.LT.0.) THEN
      J1 = J - 1
      JOFF = -1
   ENDIF
   IF (I1.GE.1 .AND. I1.LE.NROW) IBI = IBOUND(J,I1,K)
   IF (J1.GE.1 .AND. J1.LE.NCOL) IBJ = IBOUND(J1,I,K)
   IF (I1.GE.1 .AND. I1.LE.NROW .AND. J1.GE.1 .AND. J1.LE.NCOL)&
   &IBIJ = IBOUND(J1,I1,K)
   IF (I1.LT.1 .OR. I1.GT.NROW) THEN
      ROFF = 0.
      IBI = 0
   ENDIF
   IF (J1.LT.1 .OR. J1.GT.NCOL) THEN
      COFF = 0.
      IBJ = 0
   ENDIF
   IF (I1.LT.1 .OR. I1.GT.NROW .OR. J1.LT.1 .OR. J1.GT.NCOL) IBIJ = 0
   !
   CALL SOBS2BAS7HBF(COFF,DELC,DELR,I,I1,IBI,IBIJ,IBJ,IOFF,J,J1,JOFF,&
   &NCOL,NROW,RINT,ROFF)
   !
   RETURN
END
SUBROUTINE SOBS2BAS7HBF(COFF,DELC,DELR,I,I1,IBI,IBIJ,IBJ,IOFF,J,&
&J1,JOFF,NCOL,NROW,RINT,ROFF)
   !     ******************************************************************
   !     CALCULATE BASIS FUNCTIONS FOR INTERPOLATING OBSERVED HEADS.
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   DIMENSION DELR(NCOL), DELC(NROW), RINT(4)
   !     ------------------------------------------------------------------
   A=0.
   !
   !1------MOVE OBSERVATION TO NODE IF CLOSE TO NODE OR IF NEIGHBORS ARE
   !1------NO FLOW
   IF ((ABS(ROFF).LT..001.AND.ABS(COFF).LT..001) .OR.&
   &(ABS(ROFF).LT..001.AND.IBJ.EQ.0) .OR.&
   &(ABS(COFF).LT..001.AND.IBI.EQ.0) .OR. (IBI.EQ.0.AND.IBJ.EQ.0))&
   &THEN
      IOFF = 0
      JOFF = 0
      DO 10 IR = 1, 4
         RINT(IR) = .25
10    CONTINUE
      RETURN
   ENDIF
   !
   !2------CALCULATE CONSTANTS
   IF (ABS(ROFF).GE..001) THEN
      DC = (DELC(I)+DELC(I1))/2.
      DCF = ABS(ROFF)*DELC(I)
   ENDIF
   IF (ABS(COFF).GE..001) THEN
      DR = (DELR(J)+DELR(J1))/2.
      DRF = ABS(COFF)*DELR(J)
   ENDIF
   IF (ABS(ROFF).GE..001 .AND. ABS(COFF).GE..001) A = 1/(DC*DR)
   !
   !3------LINEAR INTERPOLATION
   IF (ABS(ROFF).LT..001 .OR. (IBI.EQ.0.AND.IBIJ.EQ.0)) THEN
      IOFF = 0
      RINT(1) = 0.5*(1.-DRF/DR)
      RINT(2) = 0.5*DRF/DR
      RINT(3) = RINT(1)
      RINT(4) = RINT(2)
   !
   ELSEIF (ABS(COFF).LT..001 .OR. (IBJ.EQ.0.AND.IBIJ.EQ.0)) THEN
      JOFF = 0
      RINT(1) = 0.5*(1.-DCF/DC)
      RINT(2) = RINT(1)
      RINT(3) = 0.5*DCF/DC
      RINT(4) = RINT(3)
   !
   !4------CALCULATE BASIS FUNCTIONS FOR INTERPOLATION ON A RECTANGLE
   ELSEIF (IBJ.NE.0 .AND. IBI.NE.0 .AND. IBIJ.NE.0) THEN
      RINT(3) = A*(DR-DRF)*DCF
      RINT(4) = A*DRF*DCF
      RINT(2) = A*DRF*(DC-DCF)
      RINT(1) = A*(DR-DRF)*(DC-DCF)
   !
   !5------CALCULATE BASIS FUNCTIONS FOR INTERPOLATION ON A TRIANGLE
   ELSEIF (IBJ.EQ.0) THEN
      RINT(1) = A*(DR*DC-DR*DCF)
      RINT(2) = 0.0
      RINT(3) = A*(DR*DCF-DC*DRF)
      RINT(4) = A*(DC*DRF)
   !
   ELSEIF (IBI.EQ.0) THEN
      RINT(1) = A*(DR*DC-DC*DRF)
      RINT(4) = A*(DR*DCF)
      RINT(2) = A*(DC*DRF-DR*DCF)
      RINT(3) = 0.0
   !
   ELSEIF (IBIJ.EQ.0) THEN
      RINT(1) = A*(DR*DC-DC*DRF-DR*DCF)
      RINT(3) = A*(DR*DCF)
      RINT(2) = A*(DC*DRF)
      RINT(4) = 0.0
   ENDIF
   !
   !6------
   RETURN
END
SUBROUTINE UOBSTI(ID,IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,NUMTS,&
&PERLEN,TOFF1,TOFFSET,TOMULT,TSMULT,ITR1ST,&
&OBSTIME)
   !     ******************************************************************
   !     ASSIGN OBSERVATION TIME STEP (NUMTS) AND TOFF GIVEN REFERENCE
   !     STRESS PERIOD (IREFSP), OBSERVATION-TIME OFFSET (TOFFSET), AND
   !     TIME-OFFSET MULTIPLIER (TOMULT)
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   CHARACTER*(*) ID
   INTEGER IREFSP, ISSFLG, ITR1ST, NPER, NSTP, NUMTS
   REAL DELT, ENDTIME, PERLEN, TIME, TOFF1, TOFFMULT, TOFFSET,&
   &TOMULT, TSMULT
   DIMENSION NSTP(NPER), PERLEN(NPER), TSMULT(NPER), ISSFLG(NPER)
   !     ------------------------------------------------------------------
   ZERO=0.
   !
   !1------ENSURE THAT SPECIFIED REFERENCE STRESS PERIOD IS VALID
   IF (IREFSP.LT.1 .OR. IREFSP.GT.NPER) THEN
      WRITE(IOUT,505) IREFSP
505   FORMAT(/,' REFERENCE STRESS PERIOD (IREFSP) WAS SPECIFIED AS ',&
      &I5,', BUT IT MUST BE',/,&
      &' BETWEEN 1 AND NPER (OF THE DISCRETIZATION INPUT FILE)',/,&
      &' -- STOP EXECUTION (UOBSTI)')
      CALL USTOP(' ')
   ENDIF
   !
   !2------ENSURE THAT TOFFSET IS NOT NEGATIVE
   IF (TOFFSET.LT.ZERO) THEN
      WRITE(IOUT,510) TRIM(ID)
510   FORMAT(/,' TOFFSET IS NEGATIVE FOR OBSERVATION "',A,&
      &'" -- STOP EXECUTION (UOBSTI)')
      CALL USTOP(' ')
   ENDIF
   !
   !3------FIND NUMBER OF TIME STEPS PRECEDING REFERENCE STRESS PERIOD
   NUMTS = 0
   OBSTIME = ZERO
   IF (IREFSP.GT.1) THEN
      DO 10 I = 1, IREFSP-1
         NUMTS = NUMTS + NSTP(I)
         OBSTIME = OBSTIME + PERLEN(I)
10    CONTINUE
   ENDIF
   !
   ! Note that variables TIME and ENDTIME are relative to the reference stress
   !  period. Variable OBSTIME is relative to the start of the simulation.
   TIME = ZERO
   !
   !4------USE TOMULT TO CONVERT TOFFSET TO MODEL-TIME UNITS (ASSUMES THAT
   !4------USER HAS DEFINED TOMULT CORRECTLY)
   TOFFMULT = TOFFSET*TOMULT
   OBSTIME = OBSTIME + TOFFMULT
   !
   !5------FIND STRESS PERIOD IN WHICH OBSERVATION TIME FALLS.
   !5------LOOP THROUGH STRESS PERIODS STARTING AT REFERENCE STRESS PERIOD.
   !5------TOFF1 IS OBSERVATION TIME IN TIME STEP, AS A FRACTION OF THE TIME
   !5------STEP. NUMTS IS THE NUMBER OF THE TIME STEP PRECEDING THE TIME STEP
   !5------IN WHICH THE OBSERVATION TIME OCCURS.
   DO 60 I = IREFSP, NPER
      ENDTIME = TIME+PERLEN(I)
      IF (ENDTIME.GE.TOFFMULT) THEN
   !
   !6------FIND TIME STEP PRECEDING OBSERVATION TIME
   !         CALCULATE LENGTH OF FIRST TIME STEP IN CURRENT STRESS PERIOD
         DELT = PERLEN(I)/FLOAT(NSTP(I))
         IF (TSMULT(I).NE.1.) DELT = PERLEN(I)*(1.-TSMULT(I))/&
         &(1.-TSMULT(I)**NSTP(I))
   !
   !7------LOOP THROUGH TIME STEPS
         DO 40 J = 1, NSTP(I)
            ENDTIME = TIME+DELT
            IF (ENDTIME.GE.TOFFMULT) THEN
               IF (ISSFLG(I).NE.0 .OR. ITRSS.EQ.0) THEN
   !
   !8------STEADY-STATE TIME STEP.
   !8------SET NUMTS AS THE START OF THE NEXT TIME STEP.  TELL USER UNLESS
   !8------THE OBSERVATION TIME IS THE END OF THE TIME STEP
                  IF(TOFFMULT.LT.ENDTIME) THEN
                     WRITE(IOUT,33)
33                   FORMAT(1X,'Observation within a steady-state time step has',&
                     &' been moved to the end of the time step.')
                  END IF
                  TOFF1 = 1.0
               ELSE
   !
   !9------Transient time step.
   !9A-----CALCULATE TOFF1 AS FRACTION OF TIME-STEP DURATION
                  TOFF1 = (TOFFMULT-TIME)/DELT
   !
   !9B-----CHECK FOR INITIAL TRANSIENT TIME STEP
                  IF (NUMTS.EQ.0) THEN
                     IF(ITR1ST.EQ.1) THEN
   !
   !9C-----ITR1ST IS A FLAG THAT INDICATES IF, FOR A CERTAIN
   !9C-----OBSERVATION TYPE, THE OBSERVATION TIME CAN BE BEFORE THE
   !9C-----END OF AN INITIAL TRANSIENT TIME STEP.  ITR1ST=0
   !9C-----INDICATES THAT THE OBS. TIME CAN BE IN INITIAL TRANSIENT
   !9C-----TIME STEP.  ITR1ST=1 INDICATES IT CAN'T.
                        WRITE(IOUT,37)
37                      FORMAT(1X,'The observation is in the first time step of the',&
                        &' simulation, but the observation type does not allow this.',/&
                        &1X,'The observation is being moved to the end of the time step.')
                        TOFF1 = 1.0
                     ELSE
   !
   !9D-----STOP IF THE OBSERVATION IS AT THE BEGINNING OF AN INITIAL TRANSIENT
   !9D-----TIME STEP.
                        IF(TOFFMULT.EQ.ZERO) THEN
                           WRITE(IOUT,38)
38                         FORMAT(1X,'An observation cannot be placed at the very',&
                           &' beginning of the simulation if the first period is transient.')
                           CALL USTOP(' ')
                        END IF
                     END IF
                  ENDIF
               ENDIF
               GOTO 80
            ENDIF
            TIME = TIME+DELT
            DELT = DELT*TSMULT(I)
            NUMTS = NUMTS+1
40       CONTINUE
      ELSE
         NUMTS = NUMTS+NSTP(I)
         TIME = TIME+PERLEN(I)
      ENDIF
60 CONTINUE
   !
   !10-----ALLOW FOR ROUND-OFF ERROR, SO THAT OBSERVATION TIMES SPECIFIED
   !10-----AT THE EXACT END OF THE SIMULATION ARE NOT FLAGGED AS ERRORS
   TOLERANCE = 1.0E-6*PERLEN(NPER)
   TDIFF = TOFFMULT-TIME
   IF (TDIFF.LT.TOLERANCE) THEN
      TOFF1 = 1.0
      NUMTS=NUMTS-1
   ELSE
      WRITE(IOUT,500) ID
500   FORMAT(/,' TIME SPECIFIED FOR OBSERVATION "',A,&
      &'" IS AFTER END OF SIMULATION',/,' -- STOP EXECUTION (UOBSTI)')
      CALL USTOP(' ')
   ENDIF
   !
   !11-----The Time step and interpolation coefficient have been determined.
80 CONTINUE
   RETURN
END
SUBROUTINE UOBSSV(IUOBSSV,NOBS,H,HOBS,OBSNAM,LABEL)
   !     ******************************************************************
   !     SAVE OBSERVATIONS TO A DISK FILE
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   DIMENSION H(NOBS),HOBS(NOBS)
   CHARACTER*(*) OBSNAM(NOBS)
   !     ------------------------------------------------------------------
   !
   IF(IUOBSSV.GT.0) THEN
   !
   !1------WRITE LABEL IF "LABEL" IS NOT 0
      IF(LABEL.NE.0) WRITE(IUOBSSV,18)
18    FORMAT('"SIMULATED EQUIVALENT"',3X,'"OBSERVED VALUE"',&
      &4X,'"OBSERVATION NAME"')
   !
   !2------WRITE OBSERVATIONS
      DO 100 N=1,NOBS
         WRITE(IUOBSSV,28) H(N),HOBS(N),OBSNAM(N)
28       FORMAT(1X,1P,E19.11,E20.11,2X,A)
100   CONTINUE
   END IF
   !
   !3------RETURN
   RETURN
END
SUBROUTINE OBS2BAS7DA(IUHDOB,IGRID)
   !  Deallocate OBSBAS memory
   USE OBSBASMODULE
   !
   CALL SOBS2BAS7PNT(IUHDOB,IGRID)
   DEALLOCATE(ITS)
   IF(IUHDOB.LE.0) RETURN
   !
   DEALLOCATE(NH)
   DEALLOCATE(MAXM)
   DEALLOCATE(MOBS)
   DEALLOCATE(IUHOBSV)
   DEALLOCATE(IDRY)
   DEALLOCATE(JDRY)
   DEALLOCATE(IPRT)
   DEALLOCATE(HOBDRY)
   DEALLOCATE(NDER)
   DEALLOCATE(MLAY)
   DEALLOCATE(IOFF)
   DEALLOCATE(JOFF)
   DEALLOCATE(IHOBWET)
   DEALLOCATE(H)
   DEALLOCATE(HOBS)
   DEALLOCATE(TOFF)
   DEALLOCATE(ROFF)
   DEALLOCATE(COFF)
   DEALLOCATE(OTIME)
   DEALLOCATE(PR)
   DEALLOCATE(RINT)
   DEALLOCATE(OBSNAM)
   !
   RETURN
END
SUBROUTINE SOBS2BAS7PNT(IUHDOB,IGRID)
   !  Change OBSBAS data to a different grid.
   USE OBSBASMODULE
   !
   ITS=>OBSBASDAT(IGRID)%ITS
   IF(IUHDOB.LE.0) RETURN
   !
   NH=>OBSBASDAT(IGRID)%NH
   MAXM=>OBSBASDAT(IGRID)%MAXM
   MOBS=>OBSBASDAT(IGRID)%MOBS
   IUHOBSV=>OBSBASDAT(IGRID)%IUHOBSV
   IDRY=>OBSBASDAT(IGRID)%IDRY
   JDRY=>OBSBASDAT(IGRID)%JDRY
   IPRT=>OBSBASDAT(IGRID)%IPRT
   HOBDRY=>OBSBASDAT(IGRID)%HOBDRY
   NDER=>OBSBASDAT(IGRID)%NDER
   MLAY=>OBSBASDAT(IGRID)%MLAY
   IOFF=>OBSBASDAT(IGRID)%IOFF
   JOFF=>OBSBASDAT(IGRID)%JOFF
   IHOBWET=>OBSBASDAT(IGRID)%IHOBWET
   H=>OBSBASDAT(IGRID)%H
   HOBS=>OBSBASDAT(IGRID)%HOBS
   TOFF=>OBSBASDAT(IGRID)%TOFF
   ROFF=>OBSBASDAT(IGRID)%ROFF
   COFF=>OBSBASDAT(IGRID)%COFF
   OTIME=>OBSBASDAT(IGRID)%OTIME
   PR=>OBSBASDAT(IGRID)%PR
   RINT=>OBSBASDAT(IGRID)%RINT
   OBSNAM=>OBSBASDAT(IGRID)%OBSNAM
   !
   RETURN
END
SUBROUTINE SOBS2BAS7PSV(IUHDOB,IGRID)
   !  Save OBSBAS data for a grid.
   USE OBSBASMODULE
   !
   !
   OBSBASDAT(IGRID)%ITS=>ITS
   IF(IUHDOB.LE.0) RETURN
   !
   OBSBASDAT(IGRID)%NH=>NH
   OBSBASDAT(IGRID)%MAXM=>MAXM
   OBSBASDAT(IGRID)%MOBS=>MOBS
   OBSBASDAT(IGRID)%IUHOBSV=>IUHOBSV
   OBSBASDAT(IGRID)%IDRY=>IDRY
   OBSBASDAT(IGRID)%JDRY=>JDRY
   OBSBASDAT(IGRID)%IPRT=>IPRT
   OBSBASDAT(IGRID)%HOBDRY=>HOBDRY
   OBSBASDAT(IGRID)%NDER=>NDER
   OBSBASDAT(IGRID)%MLAY=>MLAY
   OBSBASDAT(IGRID)%IOFF=>IOFF
   OBSBASDAT(IGRID)%JOFF=>JOFF
   OBSBASDAT(IGRID)%IHOBWET=>IHOBWET
   OBSBASDAT(IGRID)%H=>H
   OBSBASDAT(IGRID)%HOBS=>HOBS
   OBSBASDAT(IGRID)%TOFF=>TOFF
   OBSBASDAT(IGRID)%ROFF=>ROFF
   OBSBASDAT(IGRID)%COFF=>COFF
   OBSBASDAT(IGRID)%OTIME=>OTIME
   OBSBASDAT(IGRID)%PR=>PR
   OBSBASDAT(IGRID)%RINT=>RINT
   OBSBASDAT(IGRID)%OBSNAM=>OBSNAM
   !
   RETURN
END
