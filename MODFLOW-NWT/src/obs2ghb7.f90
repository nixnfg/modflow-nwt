MODULE OBSGHBMODULE
   INTEGER, SAVE, POINTER  ::NQGB,NQCGB,NQTGB,IUGBOBSV,IPRT
   INTEGER, SAVE, DIMENSION(:),   POINTER ::NQOBGB
   INTEGER, SAVE, DIMENSION(:),   POINTER ::NQCLGB
   INTEGER, SAVE, DIMENSION(:),   POINTER ::IOBTS
   REAL,    SAVE, DIMENSION(:),   POINTER ::FLWSIM
   REAL,    SAVE, DIMENSION(:),   POINTER ::FLWOBS
   REAL,    SAVE, DIMENSION(:),   POINTER ::TOFF
   REAL,    SAVE, DIMENSION(:),   POINTER ::OTIME
   REAL,    SAVE, DIMENSION(:,:), POINTER ::QCELL
   CHARACTER*12,SAVE,DIMENSION(:),POINTER ::OBSNAM
   TYPE OBSGHBTYPE
      INTEGER, POINTER  ::NQGB,NQCGB,NQTGB,IUGBOBSV,IPRT
      INTEGER,     DIMENSION(:),   POINTER ::NQOBGB
      INTEGER,     DIMENSION(:),   POINTER ::NQCLGB
      INTEGER,     DIMENSION(:),   POINTER ::IOBTS
      REAL,        DIMENSION(:),   POINTER ::FLWSIM
      REAL,        DIMENSION(:),   POINTER ::FLWOBS
      REAL,        DIMENSION(:),   POINTER ::TOFF
      REAL,        DIMENSION(:),   POINTER ::OTIME
      REAL,        DIMENSION(:,:), POINTER ::QCELL
      CHARACTER*12,DIMENSION(:),POINTER ::OBSNAM
   END TYPE
   TYPE(OBSGHBTYPE),  SAVE   ::OBSGHBDAT(10)
END MODULE
   !  NQGB -- number of cell groups
   !  NQCGB -- total number of cells in all groups
   !  NQTGB -- total number of observations -- sum of the number of times for each group
   !  NQOBGB(NQGB) -- The number of observations in each observation group
   !  NQCLGB(NQGB) -- The number of cells in each observation group
   !  IOBTS(NQTGB) -- Observation time step
   !  FLWSIM(NQTGB) -- Simulated value
   !  FLWOBS(NQTGB) -- Observed value
   !  TOFF(NQTGB) -- Fractional offset between time steps
   !  OTIME(NQTGB) -- Observation time in model time units
   !  QCELL(4,NQCGB) -- Location and proportion factor for each observation cell


SUBROUTINE OBS2GHB7AR(IUGBOB,IUGB,IGRID)
   !     ******************************************************************
   !     ALLOCATE MEMORY AND READ FLOW OBSERVATIONS AT GHB CELLS
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,  ONLY:IOUT,NPER,NSTP,PERLEN,TSMULT,ISSFLG,&
   &NCOL,NROW,NLAY,ITRSS
   USE OBSGHBMODULE
   CHARACTER*200 LINE
   !     ------------------------------------------------------------------
   ALLOCATE(NQGB,NQCGB,NQTGB,IUGBOBSV,IPRT)
   !
   ZERO=0.0
   IERR=0
   !  NT is the observation counter.
   NT=0
   !  NC is the cell counter.
   NC=0
   !
   !     IDENTIFY PROCESS AND PACKAGE
   WRITE(IOUT,7) IUGBOB
7  FORMAT(/,' OBS2GHB7 -- OBSERVATION PROCESS (GHB FLOW ',&
   &'OBSERVATIONS)',/,' VERSION 2, 02/28/2006',/,&
   &' INPUT READ FROM UNIT ',I4)
   !
   !x------Turn off observation package if GWFGHB is not active
   IF (IUGB.EQ.0) THEN
      WRITE (IOUT,29 )
29    FORMAT (/,' GHB PACKAGE OF GWF IS NOT OPEN')
      CALL USTOP(' ')
   ENDIF
   !
   !x------Read items 0 and 1.
   CALL URDCOM(IUGBOB,IOUT,LINE)
   LLOC = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQGB,DUM,IOUT,IUGBOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQCGB,DUM,IOUT,IUGBOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQTGB,DUM,IOUT,IUGBOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUGBOBSV,DUM,IOUT,IUGBOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,IDUM,DUM,IOUT,IUGBOB)
   IPRT=1
   IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
      IPRT=0
      WRITE(IOUT,*) 'NOPRINT option for GHB OBSERVATIONS'
   END IF
   WRITE (IOUT,9) NQGB, NQCGB, NQTGB
9  FORMAT (/,&
   &' NUMBER OF FLOW-OBSERVATION GHB-CELL GROUPS.....: ',I6,/,&
   &'   NUMBER OF CELLS IN GHB-CELL GROUPS...........: ',I6,/,&
   &'   NUMBER OF GHB-CELL FLOWS.....................: ',I6)
   IF(NQTGB.LE.0) THEN
      WRITE(IOUT,*) ' NUMBER OF OBSERVATIONS LESS THAN OR EQUAL TO 0'
      CALL USTOP(' ')
   END IF
   IF(IUGBOBSV.GT.0) THEN
      WRITE(IOUT,21) IUGBOBSV
21    FORMAT(1X,&
      &'GHB OBSERVATIONS WILL BE SAVED ON UNIT.........:',I7)
   ELSE
      WRITE(IOUT,22)
22    FORMAT(1X,'GHB OBSERVATIONS WILL NOT BE SAVED IN A FILE')
   END IF
   !
   !x------Allocate memory
   ALLOCATE(NQOBGB(NQGB))
   ALLOCATE(NQCLGB(NQGB))
   ALLOCATE(IOBTS(NQTGB))
   ALLOCATE(FLWSIM(NQTGB))
   ALLOCATE(FLWOBS(NQTGB))
   ALLOCATE(TOFF(NQTGB))
   ALLOCATE(OTIME(NQTGB))
   ALLOCATE(QCELL(4,NQCGB))
   ALLOCATE(OBSNAM(NQTGB))
   !
   !x------Initialize simulated equivalents
   DO 15 N=1,NQTGB
      FLWSIM(N)=ZERO
15 CONTINUE
   !
   !x------READ AND WRITE TIME-OFFSET MULTIPLIER FOR FLOW-OBSERVATION TIMES.
   READ(IUGBOB,*) TOMULTRV
   IF(IPRT.NE.0) WRITE (IOUT,20) TOMULTRV
20 FORMAT (/,' OBSERVED GHB-CELL FLOW DATA',/,' -- TIME OFFSETS',&
   &' ARE MULTIPLIED BY: ',G12.5)
   !
   !x------LOOP THROUGH CELL GROUPS.
   DO 200 IQ = 1,NQGB
   !
   !x------READ NUMBER OF OBSERVATINS AND NUMBER OF CELLS FOR ONE GROUP
   !x------(ITEM 3).
      READ (IUGBOB,*) NQOBGB(IQ), NQCLGB(IQ)
      IF(IPRT.NE.0) WRITE (IOUT,25) IQ, 'GHB', NQCLGB(IQ), NQOBGB(IQ)
25    FORMAT (/,'   GROUP NUMBER: ',I6,'   BOUNDARY TYPE: ',A,&
      &'   NUMBER OF CELLS IN GROUP: ',I6,/,&
      &'   NUMBER OF FLOW OBSERVATIONS: ',I6,//,&
      &40X,'OBSERVED',/,&
      &20X,'REFER.',13X,'GHB FLOW',/,&
      &7X,'OBSERVATION',2X,'STRESS',4X,'TIME',5X,'GAIN (-) OR',14X,/,&
      &2X,'OBS#    NAME',6X,'PERIOD   OFFSET',5X,'LOSS (+)')
   !
   !x------SET FLAG FOR SETTING ALL PORTION FACTORS TO 1
      IFCTFLG = 0
      IF (NQCLGB(IQ).LT.0) THEN
         IFCTFLG = 1
         NQCLGB(IQ) = -NQCLGB(IQ)
      ENDIF
   !
   !
   !x------READ THE OBSERVATION NAMES, TIMES, AND MEASURED VALUES FOR
   !x------ONE CELL GROUP (ITEM 4)
      NT1 = NT + 1
      NT2 = NT + NQOBGB(IQ)
      DO 30 J = NT1, NT2
         READ (IUGBOB,*) OBSNAM(J), IREFSP, TOFFSET, FLWOBS(J)
         IF(IPRT.NE.0) WRITE(IOUT,27 ) J,OBSNAM(J),IREFSP,TOFFSET,&
         &FLWOBS(J)
27       FORMAT (I6,1X,A12,2X,I4,2X,G11.4,1X,G11.4)
         CALL UOBSTI(OBSNAM(J),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,&
         &IOBTS(J),PERLEN,TOFF(J),TOFFSET,TOMULTRV,TSMULT,1,&
         &OTIME(J))
30    CONTINUE
   !
   !x------READ LAYER, ROW, COLUMN, AND FACTOR (ITEM 5) FOR EACH CELL IN
   !x------THE CELL GROUP.
      NC1 = NC + 1
      NC2 = NC + NQCLGB(IQ)
      IF(IPRT.NE.0) WRITE (IOUT,54)
54    FORMAT (/,'       LAYER  ROW  COLUMN    FACTOR')
      DO 100 L = NC1, NC2
         READ (IUGBOB,*) (QCELL(I,L),I=1,4)
         IF(IFCTFLG.EQ.1) QCELL(4,L) = 1.
         IF(IPRT.NE.0) WRITE (IOUT,55) (QCELL(I,L),I=1,4)
55       FORMAT (4X,F8.0,F6.0,F7.0,F9.2)
         I = QCELL(2,L)
         J = QCELL(3,L)
         IF (J.LE.0 .OR. J.GT.NCOL .OR. I.LE.0 .OR. I.GT.NROW) THEN
            WRITE (IOUT,59)
59          FORMAT (/,' ROW OR COLUMN NUMBER INVALID',&
            &' -- STOP EXECUTION (OBS2GHB7AR)',/)
            IERR = 1
         ENDIF
100   CONTINUE
   !
   !x------END OF INPUT FOR ONE CELL GROUP -- UPDATE COUNTERS.
      NC = NC2
      NT = NT2
200 CONTINUE
   !
   !
   IF (IERR.GT.0) THEN
      WRITE(IOUT,620)
620   FORMAT (/,1X,'ERROR: SEARCH ABOVE FOR ERROR MESSAGE(S)',/,&
      &' -- STOP EXECUTION (OBS2GHB7AR)')
      CALL USTOP(' ')
   ENDIF
   !
   !x------RETURN.
   CALL SOBS2GHB7PSV(IGRID)
   RETURN
END
SUBROUTINE OBS2GHB7SE(IGRID)
   !     ******************************************************************
   !     CALCULATE SIMULATED EQUIVALENTS TO OBSERVED FLOWS FOR THE GHB
   !     PACKAGE
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:IOUT,HNEW,IBOUND
   USE GWFGHBMODULE, ONLY:NBOUND,BNDS
   USE OBSBASMODULE,ONLY:ITS
   USE OBSGHBMODULE
   DOUBLE PRECISION HHNEW, C, HB , RBOT
   !     ------------------------------------------------------------------
   CALL SGWF2GHB7PNT(IGRID)
   CALL SOBS2GHB7PNT(IGRID)
   !
   !-------INITIALIZE VARIABLES
   ZERO = 0.0
   NC = 0
   NT1 = 1
   !
   !x------JRBOT IS FLAG FOR PRINTING THE HEADING FOR CELLS THAT ARE NOT
   !x------HEAD DEPENDENT.  THE FLAG IS USED TO PRINT THIS HEADING ONLY ONCE.
   JRBOT = 0
   !
   !-------LOOP THROUGH CELL GROUPS
   DO 800 IQ = 1, NQGB
      NT2 = NT1 + NQOBGB(IQ) - 1
   !
   !x------LOOK THROUGH ALL OBSERVATIONS FOR THIS GROUP TO SEE FIND OBSERVATIONS
   !x------FOR THE CURRENT TIME STEP.
      DO 600 NT = NT1, NT2
         IF (IOBTS(NT).EQ.ITS .OR.&
         &(IOBTS(NT).EQ.ITS-1.AND.TOFF(NT).GT.ZERO)) THEN
   !
   !x------FOUND AN OBSERVATION FOR CURRENT TIME STEP.
   !x------INITIALIZE NUMBER OF DRY CELLS IN OBSERVATION (KRBOT)
            KRBOT = 0
   !
   !x------LOOP THROUGH CELLS IN THE CELL GROUP
            NC1 = NC + 1
            NC2 = NC + NQCLGB(IQ)
            NB = 0
            DO 400 N = NC1, NC2
               K = QCELL(1,N)
               I = QCELL(2,N)
               J = QCELL(3,N)
   !
   !x------LOOP THROUGH ACTIVE GHB CELLS TO FIND A MATCH.
               DO 100 MNB = 1, NBOUND
                  NB = NB + 1
                  IF (NB.GT.NBOUND) NB = 1
                  KK = BNDS(1,NB)
                  II = BNDS(2,NB)
                  JJ = BNDS(3,NB)
   !
   !x------DO SIMULATED EQUIVALENT CALCULATIONS IF THIS IS A MATCH.
                  IF (I.EQ.II.AND.J.EQ.JJ.AND.K.EQ.KK) THEN
   !
   !x------CHECK IF THE MATCHED REACH IS IN A DRY CELL.
                     IF (IBOUND(J,I,K).EQ.0) THEN
                        KRBOT = KRBOT + 1
                        GOTO 400
                     ENDIF
   !
   !x------COMPUTE FLOW FOR THE BOUNDARY.
                     HHNEW = HNEW(J,I,K)
                     HB = BNDS(4,NB)
                     C = BNDS(5,NB)
                     HH = C*(HB-HHNEW)
   !
   !x------CALCULATE THE FACTOR FOR TEMPORAL INTERPOLATION.
                     TFACT = 1.0
                     IF (TOFF(NT).GT.ZERO) THEN
                        IF (IOBTS(NT).EQ.ITS) TFACT = 1. - TOFF(NT)
                        IF (IOBTS(NT).EQ.ITS-1) TFACT = TOFF(NT)
                     ENDIF
   !
   !x------ADD FLOW FOR THE REACH TO THE SIMULATED EQUIVALENT.
   !x------QCELL(4,N) IS THE PORTION FACTOR.
                     FLWSIM(NT) = FLWSIM(NT) + HH*TFACT*QCELL(4,N)
                     GO TO 400
                  ENDIF
100            CONTINUE
   !
   !x------LOOKED THROUGH ENTIRE LIST OF ACTIVE GHB CELLS WITHOUT
   !x------FINDING OBSERVATION CELL.  STOP.
               WRITE (IOUT,140) N, IQ, OBSNAM(NT),K,I,J
140            FORMAT  (' CELL ',I6,&
               &' OF GHB OBSERVATION CELL GROUP',I5,/,&
               &' NOT FOUND IN CELLS LISTED FOR GHB PACKAGE',/,&
               &' OBSERVATION NAME:',A,/,&
               &' CELL LAYER, ROW, AND COLUMN:',3I8,/,&
               &'  -- STOP EXECUTION (OBS2GHB7SE)')
               CALL USTOP(' ')
   !
   !x------END OF LOOP FOR THE CELLS IN ONE CELL GROUP FOR ONE OBSERVATION TIME..
400         CONTINUE
   !
   !-------CHECK FOR ALL CELLS IN OBSERVATION BEING DRY.
            IF(KRBOT.EQ.NQCLGB(IQ)) THEN
               WRITE (IOUT,535)
535            FORMAT(' ALL CELLS INCLUDED IN THIS OBSERVATION ARE DRY')
            ENDIF
         ENDIF
   !
   !x------END OF LOOP FOR OBSERVATION TIMES IN ONE CELL GROUP
600   CONTINUE
   !
   !-------UPDATE COUNTERS
700   NC = NC + NQCLGB(IQ)
      NT1 = NT2 + 1
   !
   !x------END OF LOOP FOR ALL CELL GROUPS.
800 CONTINUE
   !
   !x------RETURN
   RETURN
END
SUBROUTINE OBS2GHB7OT(IGRID)
   !     ******************************************************************
   !     WRITE ALL OBSERVATIONS TO LISTING FILE.
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY: IOUT
   USE OBSGHBMODULE
   DOUBLE PRECISION SQ,SUMSQ
   !     ------------------------------------------------------------------
   CALL SOBS2GHB7PNT(IGRID)
   !
   !1------WRITE OBSERVATIONS TO LISTING FILE.
   IF(IPRT.NE.0) WRITE(IOUT,17)
17 FORMAT(1X,/,1X,'GHB FLOW OBSERVATIONS',/,&
   &1X,'OBSERVATION       OBSERVED           SIMULATED',/&
   &1X,'  NAME              VALUE              VALUE',&
   &'             DIFFERENCE',/&
   &1X,'----------------------------------------------',&
   &'----------------------')
   SUMSQ=0.
   DO 100 N=1,NQTGB
      DIFF=FLWOBS(N)-FLWSIM(N)
      SQ=DIFF*DIFF
      SUMSQ=SUMSQ+SQ
      IF(IPRT.NE.0) WRITE(IOUT,27) OBSNAM(N),FLWOBS(N),FLWSIM(N),DIFF
27    FORMAT(1X,A,1P,3G20.11)
100 CONTINUE
   WRITE(IOUT,28) SUMSQ
28 FORMAT(1X,/,1X,'GHB FLOW SUM OF SQUARED DIFFERENCE:',1P,E15.5)
   !
   !2------WRITE OBSERVATIONS TO SEPARATE FILE.
   IF(IUGBOBSV.GT.0) CALL UOBSSV(IUGBOBSV,NQTGB,FLWSIM,&
   &FLWOBS,OBSNAM,0)
   !
   !3------RETURN.
   RETURN
END
SUBROUTINE OBS2GHB7DA(IGRID)
   !  Deallocate OBSGHB memory
   USE OBSGHBMODULE
   !
   CALL SOBS2GHB7PNT(IGRID)
   DEALLOCATE(NQGB)
   DEALLOCATE(NQCGB)
   DEALLOCATE(NQTGB)
   DEALLOCATE(IUGBOBSV)
   DEALLOCATE(IPRT)
   DEALLOCATE(NQOBGB)
   DEALLOCATE(NQCLGB)
   DEALLOCATE(IOBTS)
   DEALLOCATE(FLWSIM)
   DEALLOCATE(FLWOBS)
   DEALLOCATE(TOFF)
   DEALLOCATE(OTIME)
   DEALLOCATE(QCELL)
   DEALLOCATE(OBSNAM)
   !
   RETURN
END
SUBROUTINE SOBS2GHB7PNT(IGRID)
   !  Change OBSGHB data to a different grid.
   USE OBSGHBMODULE
   !
   NQGB=>OBSGHBDAT(IGRID)%NQGB
   NQCGB=>OBSGHBDAT(IGRID)%NQCGB
   NQTGB=>OBSGHBDAT(IGRID)%NQTGB
   IUGBOBSV=>OBSGHBDAT(IGRID)%IUGBOBSV
   IPRT=>OBSGHBDAT(IGRID)%IPRT
   NQOBGB=>OBSGHBDAT(IGRID)%NQOBGB
   NQCLGB=>OBSGHBDAT(IGRID)%NQCLGB
   IOBTS=>OBSGHBDAT(IGRID)%IOBTS
   FLWSIM=>OBSGHBDAT(IGRID)%FLWSIM
   FLWOBS=>OBSGHBDAT(IGRID)%FLWOBS
   TOFF=>OBSGHBDAT(IGRID)%TOFF
   OTIME=>OBSGHBDAT(IGRID)%OTIME
   QCELL=>OBSGHBDAT(IGRID)%QCELL
   OBSNAM=>OBSGHBDAT(IGRID)%OBSNAM
   !
   RETURN
END
SUBROUTINE SOBS2GHB7PSV(IGRID)
   !  Save OBSGHB data for a grid.
   USE OBSGHBMODULE
   !
   OBSGHBDAT(IGRID)%NQGB=>NQGB
   OBSGHBDAT(IGRID)%NQCGB=>NQCGB
   OBSGHBDAT(IGRID)%NQTGB=>NQTGB
   OBSGHBDAT(IGRID)%IUGBOBSV=>IUGBOBSV
   OBSGHBDAT(IGRID)%IPRT=>IPRT
   OBSGHBDAT(IGRID)%NQOBGB=>NQOBGB
   OBSGHBDAT(IGRID)%NQCLGB=>NQCLGB
   OBSGHBDAT(IGRID)%IOBTS=>IOBTS
   OBSGHBDAT(IGRID)%FLWSIM=>FLWSIM
   OBSGHBDAT(IGRID)%FLWOBS=>FLWOBS
   OBSGHBDAT(IGRID)%TOFF=>TOFF
   OBSGHBDAT(IGRID)%OTIME=>OTIME
   OBSGHBDAT(IGRID)%QCELL=>QCELL
   OBSGHBDAT(IGRID)%OBSNAM=>OBSNAM
   !
   RETURN
END
