MODULE OBSDRNMODULE
   INTEGER, SAVE, POINTER  ::NQDR,NQCDR,NQTDR,IUDROBSV,IPRT
   INTEGER, SAVE, DIMENSION(:),   POINTER ::NQOBDR
   INTEGER, SAVE, DIMENSION(:),   POINTER ::NQCLDR
   INTEGER, SAVE, DIMENSION(:),   POINTER ::IOBTS
   REAL,    SAVE, DIMENSION(:),   POINTER ::FLWSIM
   REAL,    SAVE, DIMENSION(:),   POINTER ::FLWOBS
   REAL,    SAVE, DIMENSION(:),   POINTER ::TOFF
   REAL,    SAVE, DIMENSION(:),   POINTER ::OTIME
   REAL,    SAVE, DIMENSION(:,:), POINTER ::QCELL
   CHARACTER*12,SAVE,DIMENSION(:),POINTER ::OBSNAM
   TYPE OBSDRNTYPE
      INTEGER, POINTER  ::NQDR,NQCDR,NQTDR,IUDROBSV,IPRT
      INTEGER,     DIMENSION(:),   POINTER ::NQOBDR
      INTEGER,     DIMENSION(:),   POINTER ::NQCLDR
      INTEGER,     DIMENSION(:),   POINTER ::IOBTS
      REAL,        DIMENSION(:),   POINTER ::FLWSIM
      REAL,        DIMENSION(:),   POINTER ::FLWOBS
      REAL,        DIMENSION(:),   POINTER ::TOFF
      REAL,        DIMENSION(:),   POINTER ::OTIME
      REAL,        DIMENSION(:,:), POINTER ::QCELL
      CHARACTER*12,DIMENSION(:),POINTER ::OBSNAM
   END TYPE
   TYPE(OBSDRNTYPE),  SAVE   ::OBSDRNDAT(10)
END MODULE
   !  NQDR -- number of cell groups
   !  NQCDR -- total number of cells in all groups
   !  NQTDR -- total number of observations -- sum of the number of times for each group
   !  NQOBDR(NQDR) -- The number of observations in each observation group
   !  NQCLDR(NQDR) -- The number of cells in each observation group
   !  IOBTS(NQTDR) -- Observation time step
   !  FLWSIM(NQTDR) -- Simulated value
   !  FLWOBS(NQTDR) -- Observed value
   !  TOFF(NQTDR) -- Fractional offset between time steps
   !  OTIME(NQTDR) --
   !  QCELL(4,NQCDR) -- Location and proportion factor for each observation cell


SUBROUTINE OBS2DRN7AR(IUDROB,IUDRN,IGRID)
   !     ******************************************************************
   !     ALLOCATE MEMORY AND READ FLOW OBSERVATIONS AT DRAIN CELLS
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,  ONLY:IOUT,NPER,NSTP,PERLEN,TSMULT,ISSFLG,&
   &NCOL,NROW,NLAY,ITRSS
   USE OBSDRNMODULE
   CHARACTER*200 LINE
   !     ------------------------------------------------------------------
   !
   ALLOCATE(NQDR,NQCDR,NQTDR,IUDROBSV,IPRT)
   !
   ZERO=0.0
   IERR=0
   !  NT is the observation counter.
   NT=0
   !  NC is the cell counter.
   NC=0
   !
   !     IDENTIFY PROCESS AND PACKAGE
   WRITE(IOUT,7) IUDROB
7  FORMAT(/,' OBS2DRN7 -- OBSERVATION PROCESS (DRAIN FLOW ',&
   &'OBSERVATIONS)',/,' VERSION 2, 02/28/2006',/,&
   &' INPUT READ FROM UNIT ',I4)
   !
   !x------Stop if GWFDRN is not active
   IF (IUDRN.EQ.0) THEN
      WRITE (IOUT,29 )
29    FORMAT (/,' DRAIN PACKAGE OF GWF IS NOT OPEN')
      CALL USTOP(' ')
   ENDIF
   !
   !x------Read items 0 and 1.
   CALL URDCOM(IUDROB,IOUT,LINE)
   LLOC = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQDR,DUM,IOUT,IUDROB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQCDR,DUM,IOUT,IUDROB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQTDR,DUM,IOUT,IUDROB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUDROBSV,DUM,IOUT,IUDROB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,IDUM,DUM,IOUT,IUDROB)
   IPRT=1
   IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
      IPRT=0
      WRITE(IOUT,*) 'NOPRINT option for DRAIN OBSERVATIONS'
   END IF
   WRITE (IOUT,9) NQDR, NQCDR, NQTDR
9  FORMAT (/,&
   &' NUMBER OF FLOW-OBSERVATION DRAIN-CELL GROUPS.....: ',I6,/,&
   &'   NUMBER OF CELLS IN DRAIN-CELL GROUPS...........: ',I6,/,&
   &'   NUMBER OF DRAIN-CELL FLOWS.....................: ',I6)
   IF(NQTDR.LE.0) THEN
      WRITE(IOUT,*) ' NUMBER OF OBSERVATIONS LESS THAN OR EQUAL TO 0'
      CALL USTOP(' ')
   END IF
   IF(IUDROBSV.GT.0) THEN
      WRITE(IOUT,21) IUDROBSV
21    FORMAT(1X,&
      &'DRAIN OBSERVATIONS WILL BE SAVED ON UNIT.........:',I7)
   ELSE
      WRITE(IOUT,22)
22    FORMAT(1X,'DRAIN OBSERVATIONS WILL NOT BE SAVED IN A FILE')
   END IF
   !
   !x------Allocate memory
   ALLOCATE(NQOBDR(NQDR))
   ALLOCATE(NQCLDR(NQDR))
   ALLOCATE(IOBTS(NQTDR))
   ALLOCATE(FLWSIM(NQTDR))
   ALLOCATE(FLWOBS(NQTDR))
   ALLOCATE(TOFF(NQTDR))
   ALLOCATE(OTIME(NQTDR))
   ALLOCATE(QCELL(4,NQCDR))
   ALLOCATE(OBSNAM(NQTDR))
   !
   !x------Initialize simulated equivalents
   DO 15 N=1,NQTDR
      FLWSIM(N)=ZERO
15 CONTINUE
   !
   !x------READ AND WRITE TIME-OFFSET MULTIPLIER FOR FLOW-OBSERVATION TIMES.
   READ(IUDROB,*) TOMULTDR
   IF(IPRT.NE.0) WRITE (IOUT,20) TOMULTDR
20 FORMAT (/,' OBSERVED DRAIN-CELL FLOW DATA',/,' -- TIME OFFSETS',&
   &' ARE MULTIPLIED BY: ',G12.5)
   !
   !x------LOOP THROUGH CELL GROUPS.
   DO 200 IQ = 1,NQDR
   !
   !x------READ NUMBER OF OBSERVATIONS AND NUMBER OF CELLS FOR ONE GROUP
   !x------(ITEM 3).
      READ (IUDROB,*) NQOBDR(IQ), NQCLDR(IQ)
      IF(IPRT.NE.0) WRITE (IOUT,25) IQ, 'DRN', NQCLDR(IQ), NQOBDR(IQ)
25    FORMAT (/,'   GROUP NUMBER: ',I6,'   BOUNDARY TYPE: ',A,&
      &'   NUMBER OF CELLS IN GROUP: ',I6,/,&
      &'   NUMBER OF FLOW OBSERVATIONS: ',I6,//,&
      &40X,'OBSERVED',/,&
      &20X,'REFER.',13X,'DRAIN FLOW',/,&
      &7X,'OBSERVATION',2X,'STRESS',4X,'TIME',5X,'GAIN (-) OR',14X,/,&
      &2X,'OBS#    NAME',6X,'PERIOD   OFFSET',5X,'LOSS (+)')
   !
   !x------SET FLAG FOR SETTING ALL PORTION FACTORS TO 1
      IFCTFLG = 0
      IF (NQCLDR(IQ).LT.0) THEN
         IFCTFLG = 1
         NQCLDR(IQ) = -NQCLDR(IQ)
      ENDIF
   !
   !
   !x------READ THE OBSERVATION NAMES, TIMES, AND MEASURED VALUES FOR
   !x------ONE CELL GROUP (ITEM 4)
      NT1 = NT + 1
      NT2 = NT + NQOBDR(IQ)
      DO 30 J = NT1, NT2
         READ (IUDROB,*) OBSNAM(J), IREFSP, TOFFSET, FLWOBS(J)
         IF(IPRT.NE.0) WRITE(IOUT,27 ) J,OBSNAM(J),IREFSP,TOFFSET,&
         &FLWOBS(J)
27       FORMAT (I6,1X,A12,2X,I4,2X,G11.4,1X,G11.4)
         CALL UOBSTI(OBSNAM(J),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,&
         &IOBTS(J),PERLEN,TOFF(J),TOFFSET,TOMULTDR,TSMULT,1,&
         &OTIME(J))
30    CONTINUE
   !
   !x------READ LAYER, ROW, COLUMN, AND FACTOR (ITEM 5) FOR EACH CELL IN
   !x------THE CELL GROUP.
      NC1 = NC + 1
      NC2 = NC + NQCLDR(IQ)
      IF(IPRT.NE.0) WRITE (IOUT,54)
54    FORMAT (/,'       LAYER  ROW  COLUMN    FACTOR')
      DO 100 L = NC1, NC2
         READ (IUDROB,*) (QCELL(I,L),I=1,4)
         IF(IFCTFLG.EQ.1) QCELL(4,L) = 1.
         IF(IPRT.NE.0) WRITE (IOUT,55) (QCELL(I,L),I=1,4)
55       FORMAT (4X,F8.0,F6.0,F7.0,F9.2)
         I = QCELL(2,L)
         J = QCELL(3,L)
         IF (J.LE.0 .OR. J.GT.NCOL .OR. I.LE.0 .OR. I.GT.NROW) THEN
            WRITE (IOUT,59)
59          FORMAT (/,' ROW OR COLUMN NUMBER INVALID',&
            &' -- STOP EXECUTION (OBS2DRN7RP)',/)
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
      &' -- STOP EXECUTION (OBS2DRN7RP)')
      CALL USTOP(' ')
   ENDIF
   !
   !x------RETURN.
   CALL SOBS2DRN7PSV(IGRID)
   RETURN
END
SUBROUTINE OBS2DRN7SE(IGRID)
   !     ******************************************************************
   !     CALCULATE SIMULATED EQUIVALENTS TO OBSERVED FLOWS FOR THE DRAIN
   !     PACKAGE
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:IOUT,HNEW,IBOUND
   USE GWFDRNMODULE, ONLY:NDRAIN,DRAI
   USE OBSBASMODULE,ONLY:ITS
   USE OBSDRNMODULE
   DOUBLE PRECISION HHNEW, HB, C
   !     ------------------------------------------------------------------
   CALL SGWF2DRN7PNT(IGRID)
   CALL SOBS2DRN7PNT(IGRID)
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
   DO 800 IQ = 1, NQDR
      NT2 = NT1 + NQOBDR(IQ) - 1
   !
   !x------LOOK THROUGH ALL OBSERVATIONS FOR THIS GROUP TO SEE FIND OBSERVATIONS
   !x------FOR THE CURRENT TIME STEP.
      DO 600 NT = NT1, NT2
         IF (IOBTS(NT).EQ.ITS .OR.&
         &(IOBTS(NT).EQ.ITS-1.AND.TOFF(NT).GT.ZERO)) THEN
   !
   !x------FOUND AN OBSERVATION FOR CURRENT TIME STEP.
   !x------INITIALIZE NUMBER OF DRY CELLS IN OBSERVATION (KRBOT) AND
   !x------NUMBER OF NON-HEAD-DEPENDENT CELLS (IRBOT).
            IRBOT = 0
            KRBOT = 0
   !
   !x------LOOP THROUGH CELLS IN THE CELL GROUP
            NC1 = NC + 1
            NC2 = NC + NQCLDR(IQ)
            NB = 0
            DO 400 N = NC1, NC2
               K = QCELL(1,N)
               I = QCELL(2,N)
               J = QCELL(3,N)
   !
   !x------LOOP THROUGH ACTIVE DRAIN REACHES TO FIND A MATCH.
               DO 100 MNB = 1, NDRAIN
                  NB = NB + 1
                  IF (NB.GT.NDRAIN) NB = 1
                  KK = DRAI(1,NB)
                  II = DRAI(2,NB)
                  JJ = DRAI(3,NB)
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
   !x------COMPUTE FLOW FOR THE REACH.
                     HHNEW = HNEW(J,I,K)
                     HB = DRAI(4,NB)
                     C = DRAI(5,NB)
                     HH = C*(HB-HHNEW)
                     IF(HHNEW.LE.HB) THEN
                        HH = ZERO
                        IF (JRBOT.EQ.0) WRITE (IOUT,83 )
83                      FORMAT (/,&
                        &' HEADS AT DRAIN CELLS ARE BELOW THE',&
                        &' BOTTOM OF THE DRAIN BED AT THE CELLS LISTED',/,&
                        &' BELOW.  THESE CONDITIONS DIMINISH THE IMPACT',&
                        &' OF THE OBSERVATION ON ESTIMATES OF',/,&
                        &' ALL PARAMETERS EXCEPT THOSE THAT CONTROL THE HYDRAULIC',&
                        &' CONDUCTIVITY OF THE',/,&
                        &' DRAIN BED.  (SEE TEXT FOR MORE INFORMATION).')
                        JRBOT = 1
                        IF (IRBOT.EQ.0) THEN
                           WRITE (IOUT,92 ) NT, OBSNAM(NT), ITS
92                         FORMAT (/,' OBS# ',I6,', ID ',A,', TIME STEP ',I5)
                           WRITE (IOUT,93 )
93                         FORMAT ('    LAYER   ROW  COLUMN')
                        ENDIF
                        IRBOT = IRBOT + 1
                        WRITE (IOUT,97 ) K, I, J
97                      FORMAT(3I7)
                     ENDIF
   !
   !x------CALCULATE THE FACTOR FOR TEMPORAL INTERPOLATION.
                     FACT = 1.0
                     IF (TOFF(NT).GT.ZERO) THEN
                        IF (IOBTS(NT).EQ.ITS) FACT = 1. - TOFF(NT)
                        IF (IOBTS(NT).EQ.ITS-1) FACT = TOFF(NT)
                     ENDIF
   !
   !x------ADD FLOW FOR THE REACH TO THE SIMULATED EQUIVALENT.
   !x------QCELL(4,N) IS THE PORTION FACTOR.
                     FLWSIM(NT) = FLWSIM(NT) + HH*FACT*QCELL(4,N)
                     GO TO 400
                  ENDIF
100            CONTINUE
   !
   !x------LOOKED THROUGH ENTIRE LIST OF ACTIVE DRAIN REACHES WITHOUT
   !x------FINDING OBSERVATION CELL.  STOP.
               WRITE (IOUT,140) N, IQ, OBSNAM(NT),K,I,J
140            FORMAT  (' CELL ',I6,&
               &' OF DRAIN OBSERVATION CELL GROUP',I5,/,&
               &' NOT FOUND IN CELLS LISTED FOR DRAIN PACKAGE',/,&
               &' OBSERVATION NAME:',A,/,&
               &' CELL LAYER, ROW, AND COLUMN:',3I8,/,&
               &'  -- STOP EXECUTION (OBS2DRN7SE)')
               CALL USTOP(' ')
   !
   !x------END OF LOOP FOR THE CELLS IN ONE CELL GROUP FOR ONE OBSERVATION TIME..
400         CONTINUE
   !
   !-------PRINT NUMBER OF CELLS AT WHICH HEAD IS BELOW THE BOTTOM OF THE
   !-------DRAIN BED; CHECK FOR ALL CELLS IN OBSERVATION BEING DRY.
            IF(IRBOT.GT.0) WRITE (IOUT,530) IRBOT, NQCLDR(IQ)
530         FORMAT (I7,' OF THE',I7,' CELLS USED TO SIMULATE THE',&
            &' GAIN OR LOSS ARE',/,22X,'AFFECTED.')
            IF(KRBOT.EQ.NQCLDR(IQ)) THEN
               WRITE (IOUT,535)
535            FORMAT(' ALL CELLS INCLUDED IN THIS OBSERVATION ARE DRY')
            ENDIF
         ENDIF
   !
   !x------END OF LOOP FOR OBSERVATION TIMES IN ONE CELL GROUP
600   CONTINUE
   !
   !-------UPDATE COUNTERS
700   NC = NC + NQCLDR(IQ)
      NT1 = NT2 + 1
   !
   !x------END OF LOOP FOR ALL CELL GROUPS.
800 CONTINUE
   !
   !x------RETURN
   RETURN
END
SUBROUTINE OBS2DRN7OT(IGRID)
   !     ******************************************************************
   !     WRITE ALL OBSERVATIONS TO LISTING FILE.
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY: IOUT
   USE OBSDRNMODULE
   DOUBLE PRECISION SQ,SUMSQ
   !     ------------------------------------------------------------------
   CALL SOBS2DRN7PNT(IGRID)
   !
   !1------WRITE OBSERVATIONS TO LISTING FILE.
   IF(IPRT.NE.0) WRITE(IOUT,17)
17 FORMAT(1X,/,1X,'DRAIN FLOW OBSERVATIONS',/,&
   &1X,'OBSERVATION       OBSERVED           SIMULATED',/&
   &1X,'  NAME              VALUE              VALUE',&
   &'             DIFFERENCE',/&
   &1X,'----------------------------------------------',&
   &'----------------------')
   SUMSQ=0.
   DO 100 N=1,NQTDR
      DIFF=FLWOBS(N)-FLWSIM(N)
      SQ=DIFF*DIFF
      SUMSQ=SUMSQ+SQ
      IF(IPRT.NE.0) WRITE(IOUT,27) OBSNAM(N),FLWOBS(N),FLWSIM(N),DIFF
27    FORMAT(1X,A,1P,3G20.11)
100 CONTINUE
   WRITE(IOUT,28) SUMSQ
28 FORMAT(1X,/,1X,'DRN FLOW SUM OF SQUARED DIFFERENCE:',1P,E15.5)
   !
   !2------WRITE OBSERVATIONS TO SEPARATE FILE.
   IF(IUDROBSV.GT.0) CALL UOBSSV(IUDROBSV,NQTDR,FLWSIM,&
   &FLWOBS,OBSNAM,0)
   !
   !3------RETURN.
   RETURN
END
SUBROUTINE OBS2DRN7DA(IGRID)
   !  Deallocate OBSDRN memory
   USE OBSDRNMODULE
   !
   CALL SOBS2DRN7PNT(IGRID)
   DEALLOCATE(NQDR)
   DEALLOCATE(NQCDR)
   DEALLOCATE(NQTDR)
   DEALLOCATE(IUDROBSV)
   DEALLOCATE(IPRT)
   DEALLOCATE(NQOBDR)
   DEALLOCATE(NQCLDR)
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
SUBROUTINE SOBS2DRN7PNT(IGRID)
   !  Change OBSDRN data to a different grid.
   USE OBSDRNMODULE
   !
   NQDR=>OBSDRNDAT(IGRID)%NQDR
   NQCDR=>OBSDRNDAT(IGRID)%NQCDR
   NQTDR=>OBSDRNDAT(IGRID)%NQTDR
   IUDROBSV=>OBSDRNDAT(IGRID)%IUDROBSV
   IPRT=>OBSDRNDAT(IGRID)%IPRT
   NQOBDR=>OBSDRNDAT(IGRID)%NQOBDR
   NQCLDR=>OBSDRNDAT(IGRID)%NQCLDR
   IOBTS=>OBSDRNDAT(IGRID)%IOBTS
   FLWSIM=>OBSDRNDAT(IGRID)%FLWSIM
   FLWOBS=>OBSDRNDAT(IGRID)%FLWOBS
   TOFF=>OBSDRNDAT(IGRID)%TOFF
   OTIME=>OBSDRNDAT(IGRID)%OTIME
   QCELL=>OBSDRNDAT(IGRID)%QCELL
   OBSNAM=>OBSDRNDAT(IGRID)%OBSNAM
   !
   RETURN
END
SUBROUTINE SOBS2DRN7PSV(IGRID)
   !  Save OBSDRN data for a grid.
   USE OBSDRNMODULE
   !
   OBSDRNDAT(IGRID)%NQDR=>NQDR
   OBSDRNDAT(IGRID)%NQCDR=>NQCDR
   OBSDRNDAT(IGRID)%NQTDR=>NQTDR
   OBSDRNDAT(IGRID)%IUDROBSV=>IUDROBSV
   OBSDRNDAT(IGRID)%IPRT=>IPRT
   OBSDRNDAT(IGRID)%NQOBDR=>NQOBDR
   OBSDRNDAT(IGRID)%NQCLDR=>NQCLDR
   OBSDRNDAT(IGRID)%IOBTS=>IOBTS
   OBSDRNDAT(IGRID)%FLWSIM=>FLWSIM
   OBSDRNDAT(IGRID)%FLWOBS=>FLWOBS
   OBSDRNDAT(IGRID)%TOFF=>TOFF
   OBSDRNDAT(IGRID)%OTIME=>OTIME
   OBSDRNDAT(IGRID)%QCELL=>QCELL
   OBSDRNDAT(IGRID)%OBSNAM=>OBSNAM
   !
   RETURN
END
