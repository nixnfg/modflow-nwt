MODULE OBSSTRMODULE
   INTEGER, SAVE, POINTER  ::NQST,NQCST,NQTST,IUSTOBSV,IPRT
   INTEGER, SAVE, DIMENSION(:),   POINTER ::NQOBST
   INTEGER, SAVE, DIMENSION(:),   POINTER ::NQCLST
   INTEGER, SAVE, DIMENSION(:),   POINTER ::IOBTS
   REAL,    SAVE, DIMENSION(:),   POINTER ::FLWSIM
   REAL,    SAVE, DIMENSION(:),   POINTER ::FLWOBS
   REAL,    SAVE, DIMENSION(:),   POINTER ::TOFF
   REAL,    SAVE, DIMENSION(:),   POINTER ::OTIME
   REAL,    SAVE, DIMENSION(:,:), POINTER ::QCELL
   CHARACTER*12,SAVE,DIMENSION(:),POINTER ::OBSNAM
   TYPE OBSSTRTYPE
      INTEGER, POINTER  ::NQST,NQCST,NQTST,IUSTOBSV,IPRT
      INTEGER,     DIMENSION(:),   POINTER ::NQOBST
      INTEGER,     DIMENSION(:),   POINTER ::NQCLST
      INTEGER,     DIMENSION(:),   POINTER ::IOBTS
      REAL,        DIMENSION(:),   POINTER ::FLWSIM
      REAL,        DIMENSION(:),   POINTER ::FLWOBS
      REAL,        DIMENSION(:),   POINTER ::TOFF
      REAL,        DIMENSION(:),   POINTER ::OTIME
      REAL,        DIMENSION(:,:), POINTER ::QCELL
      CHARACTER*12,DIMENSION(:),POINTER ::OBSNAM
   END TYPE
   TYPE(OBSSTRTYPE),  SAVE   ::OBSSTRDAT(10)
END MODULE
   !  NQST -- number of cell groups
   !  NQCST -- total number of cells in all groups
   !  NQTST -- total number of observations -- sum of the number of times for each group
   !  NQOBST(NQST) -- The number of observations in each observation group
   !  NQCLST(NQST) -- The number of cells in each observation group
   !  IOBTS(NQTST) -- observation time step
   !  FLWSIM(NQTST) -- Simulated value
   !  FLWOBS(NQTST) -- observed value
   !  TOFF(NQTST) -- Fractional offset between time steps
   !  OTIME(NQTST) -- observation time in model time units
   !  QCELL(4,NQCST) -- Data for each observation cell -- Segment, Reach, Unused,
   !                    Proportion factor


SUBROUTINE OBS2STR7AR(IUSTOB,IUST,IGRID)
   !     ******************************************************************
   !     ALLOCATE MEMORY AND READ FLOW OBSERVATIONS AT STREAM CELLS
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,  ONLY:IOUT,NPER,NSTP,PERLEN,TSMULT,ISSFLG,ITRSS
   USE OBSSTRMODULE
   CHARACTER*200 LINE
   !     ------------------------------------------------------------------
   ALLOCATE(NQST,NQCST,NQTST,IUSTOBSV,IPRT)
   !
   ZERO=0.0
   IERR=0
   !  NT is the observation counter.
   NT=0
   !  NC is the cell counter.
   NC=0
   !
   !     IDENTIFY PROCESS AND PACKAGE
   WRITE(IOUT,7) IUSTOB
7  FORMAT(/,' OBS2STR7 -- OBSERVATION PROCESS (STREAMFLOW ',&
   &'OBSERVATIONS)',/,' VERSION 2, 08/05/2009',/,&
   &' INPUT READ FROM UNIT ',I4)
   !
   !x------Stop if GWFSTR is not active
   IF (IUST.EQ.0) THEN
      WRITE (IOUT,29 )
29    FORMAT (/,' STREAMFLOW PACKAGE OF GWF IS NOT OPEN')
      CALL USTOP(' ')
      RETURN
   ENDIF
   !
   !x------Read items 0 and 1.
   CALL URDCOM(IUSTOB,IOUT,LINE)
   LLOC = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQST,DUM,IOUT,IUSTOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQCST,DUM,IOUT,IUSTOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQTST,DUM,IOUT,IUSTOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUSTOBSV,DUM,IOUT,IUSTOB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,IDUM,DUM,IOUT,IUSTOB)
   IPRT=1
   IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
      IPRT=0
      WRITE(IOUT,*) 'NOPRINT option for STREAM OBSERVATIONS'
   END IF
   WRITE (IOUT,9) NQST, NQCST, NQTST
9  FORMAT (/,&
   &' NUMBER OF FLOW-OBSERVATION STREAM-CELL GROUPS.....: ',I6,/,&
   &'   NUMBER OF CELLS IN STREAM-CELL GROUPS...........: ',I6,/,&
   &'   NUMBER OF STREAM-CELL FLOWS.....................: ',I6)
   IF(NQTST.LE.0) THEN
      WRITE(IOUT,*) ' NUMBER OF OBSERVATIONS LESS THAN OR EQUAL TO 0'
      CALL USTOP(' ')
   END IF
   IF(IUSTOBSV.GT.0) THEN
      WRITE(IOUT,21) IUSTOBSV
21    FORMAT(1X,&
      &'STREAM OBSERVATIONS WILL BE SAVED ON UNIT.........:',I7)
   ELSE
      WRITE(IOUT,22)
22    FORMAT(1X,'STREAM OBSERVATIONS WILL NOT BE SAVED IN A FILE')
   END IF
   !
   !x------Allocate memory
   ALLOCATE(NQOBST(NQST))
   ALLOCATE(NQCLST(NQST))
   ALLOCATE(IOBTS(NQTST))
   ALLOCATE(FLWSIM(NQTST))
   ALLOCATE(FLWOBS(NQTST))
   ALLOCATE(TOFF(NQTST))
   ALLOCATE(OTIME(NQTST))
   ALLOCATE(QCELL(4,NQCST))
   ALLOCATE(OBSNAM(NQTST))
   !
   !x------Initialize simulated equivalents
   DO 15 N=1,NQTST
      FLWSIM(N)=ZERO
15 CONTINUE
   !
   !x------READ AND WRITE TIME-OFFSET MULTIPLIER FOR FLOW-OBSERVATION TIMES.
   READ(IUSTOB,*) TOMULTST
   IF(IPRT.NE.0) WRITE (IOUT,20) TOMULTST
20 FORMAT (/,' OBSERVED STREAM-CELL FLOW DATA',/,' -- TIME OFFSETS',&
   &' ARE MULTIPLIED BY: ',G12.5)
   !
   !x------LOOP THROUGH CELL GROUPS.
   DO 200 IQ = 1,NQST
   !
   !x------READ NUMBER OF OBSERVATIONS AND NUMBER OF CELLS FOR ONE GROUP
   !x------(ITEM 3).
      READ (IUSTOB,*) NQOBST(IQ), NQCLST(IQ)
      IF(IPRT.NE.0) WRITE (IOUT,25) IQ, 'STR', NQCLST(IQ), NQOBST(IQ)
25    FORMAT (/,'   GROUP NUMBER: ',I6,'   BOUNDARY TYPE: ',A,&
      &'   NUMBER OF CELLS IN GROUP: ',I6,/,&
      &'   NUMBER OF FLOW OBSERVATIONS: ',I6,//,&
      &40X,'OBSERVED',/,&
      &20X,'REFER.',13X,'STREAMFLOW',/,&
      &7X,'OBSERVATION',2X,'STRESS',4X,'TIME',5X,'GAIN (-) OR',14X,/,&
      &2X,'OBS#    NAME',6X,'PERIOD   OFFSET',5X,'LOSS (+)')
   !
   !x------SET FLAG FOR SETTING ALL PORTION FACTORS TO 1
      IFCTFLG = 0
      IF (NQCLST(IQ).LT.0) THEN
         IFCTFLG = 1
         NQCLST(IQ) = -NQCLST(IQ)
      ENDIF
   !
   !
   !x------READ THE OBSERVATION NAMES, TIMES, AND MEASURED VALUES FOR
   !x------ONE CELL GROUP (ITEM 4)
      NT1 = NT + 1
      NT2 = NT + NQOBST(IQ)
      DO 30 J = NT1, NT2
         READ (IUSTOB,*) OBSNAM(J), IREFSP, TOFFSET, FLWOBS(J)
         IF(IPRT.NE.0) WRITE(IOUT,27 ) J,OBSNAM(J),IREFSP,TOFFSET,&
         &FLWOBS(J)
27       FORMAT (I6,1X,A12,2X,I4,2X,G11.4,1X,G11.4)
         CALL UOBSTI(OBSNAM(J),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,&
         &IOBTS(J),PERLEN,TOFF(J),TOFFSET,TOMULTST,TSMULT,1,&
         &OTIME(J))
30    CONTINUE
   !
   !x------READ SEGMENT, REACH, AND FACTOR (ITEM 5) FOR EACH CELL IN
   !x------THE CELL GROUP.
      NC1 = NC + 1
      NC2 = NC + NQCLST(IQ)
      IF(IPRT.NE.0) WRITE (IOUT,54)
54    FORMAT (/,'     SEGMENT  REACH    FACTOR')
      DO 100 L = NC1, NC2
         READ (IUSTOB,*) (QCELL(I,L),I=1,2),QCELL(4,L)
         IF(IFCTFLG.EQ.1) QCELL(4,L) = 1.
         IF(IPRT.NE.0) WRITE (IOUT,55) (QCELL(I,L),I=1,2),QCELL(4,L)
55       FORMAT (4X,F8.0,F6.0,F9.2)
         QCELL(3,L)=0.0
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
      &' -- STOP EXECUTION (OBS2STR7)')
      CALL USTOP(' ')
   ENDIF
   !
   !x------RETURN.
   CALL SOBS2STR7PSV(IGRID)
   RETURN
END
SUBROUTINE OBS2STR7SE(IGRID)
   !     ******************************************************************
   !     CALCULATE SIMULATED EQUIVALENTS TO OBSERVED FLOWS FOR THE STREAM
   !     PACKAGE
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:IOUT,HNEW,IBOUND
   USE GWFSTRMODULE, ONLY:NSTREM,STRM,ISTRM
   USE OBSBASMODULE,ONLY:ITS
   USE OBSSTRMODULE
   DOUBLE PRECISION HHNEW, C, HB , RBOT
   !     ------------------------------------------------------------------
   CALL SGWF2STR7PNT(IGRID)
   CALL SOBS2STR7PNT(IGRID)
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
   DO 800 IQ = 1, NQST
      NT2 = NT1 + NQOBST(IQ) - 1
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
            NC2 = NC + NQCLST(IQ)
            NB = 0
            DO 400 N = NC1, NC2
               K = QCELL(1,N)
               I = QCELL(2,N)
   !
   !x------LOOP THROUGH ACTIVE STREAM REACHES TO FIND A MATCH.
               DO 100 MNB = 1, NSTREM
                  NB = NB + 1
                  IF (NB.GT.NSTREM) NB = 1
                  IS=ISTRM(4,NB)
                  IR=ISTRM(5,NB)
   !
   !x------DO SIMULATED EQUIVALENT CALCULATIONS IF THIS IS A MATCH.
                  IF (I.EQ.IR .AND. K.EQ.IS) THEN
                     KK = ISTRM(1,NB)
                     II = ISTRM(2,NB)
                     JJ = ISTRM(3,NB)
   !
   !x------CHECK IF THE MATCHED REACH IS IN A DRY CELL.
                     IF (IBOUND(JJ,II,KK).EQ.0) THEN
                        KRBOT = KRBOT + 1
                        GOTO 400
                     ENDIF
   !
   !x------COMPUTE FLOW FOR THE REACH.
                     HHNEW = HNEW(JJ,II,KK)
                     HB = STRM(2,NB)
                     IF(STRM(10,NB).LE.ZERO) HB=STRM(5,NB)
                     C = STRM(3,NB)
                     RBOT = STRM(4,NB)
                     FLWCEL = C*(HB-HHNEW)
                     ISTRF=0
                     IF(STRM(10,NB).LE.STRM(11,NB)) ISTRF=1
                     IF(HHNEW.LE.RBOT .OR. ISTRF.EQ.1) THEN
                        FLWCEL = C*(HB-RBOT)
                        IF(ISTRF.EQ.1) FLWCEL=STRM(10,NB)
                        IF (JRBOT.EQ.0) WRITE (IOUT,41 )
41                      FORMAT(/,'For the cells listed below, one of two',&
                        &' conditions exist, as indicated by the',/,&
                        &' absence or presence of an asterisk (*). If the',&
                        &' observation is being used for',/,&
                        &' parameter estimation, these conditions can',&
                        &' diminish the impact of the',/,&
                        &' observation for some parameters.',/,&
                        &'  No *: Aquifer head at the stream cell is',&
                        &' below the bottom of the streambed.',/,&
                        &'  *   : Seepage to the aquifer is limited by',&
                        &' available streamflow.')
                        JRBOT = 1
                        IF (IRBOT.EQ.0) THEN
                           WRITE (IOUT,42 ) NT, OBSNAM(NT), ITS
42                         FORMAT (/,' OBS# ',I6,', ID ',A,', TIME STEP ',I5)
                           WRITE (IOUT,43 )
43                         FORMAT ('   SEGMENT REACH')
                        ENDIF
                        IRBOT = IRBOT + 1
                        IF(ISTRF.EQ.1) THEN
                           WRITE (IOUT,44 ) IS,IR
44                         FORMAT(' *',I5,I7)
                        ELSE
                           WRITE (IOUT,45 ) IS,IR
45                         FORMAT(2I7)
                        END IF
                     ENDIF
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
                     FLWSIM(NT) = FLWSIM(NT) + FLWCEL*TFACT*QCELL(4,N)
                     GO TO 400
                  ENDIF
100            CONTINUE
   !
   !x------LOOKED THROUGH ENTIRE LIST OF ACTIVE STREAM REACHES WITHOUT
   !x------FINDING OBSERVATION CELL.  STOP.
               WRITE (IOUT,140) N, IQ, OBSNAM(NT),K,I
140            FORMAT  (' CELL ',I6,&
               &' OF STREAM OBSERVATION CELL GROUP',I5,/,&
               &' NOT FOUND IN CELLS LISTED FOR STREAM PACKAGE',/,&
               &' OBSERVATION NAME:',A,/,&
               &' STREAM SEGMENT AND REACH:',2I8,/,&
               &'  -- STOP EXECUTION (OBS2STR7SE)')
               CALL USTOP(' ')
   !
   !x------END OF LOOP FOR THE CELLS IN ONE CELL GROUP FOR ONE OBSERVATION TIME..
400         CONTINUE
   !
   !-------PRINT NUMBER OF CELLS AT WHICH HEAD IS BELOW THE BOTTOM OF THE
   !-------STREAMBED; CHECK FOR ALL CELLS IN OBSERVATION BEING DRY.
            IF(IRBOT.GT.0) WRITE (IOUT,530) IRBOT, NQCLST(IQ)
530         FORMAT (I7,' OF THE',I7,' CELLS USED TO SIMULATE THE',&
            &' GAIN OR LOSS ARE',/,22X,'AFFECTED.')
            IF(KRBOT.EQ.NQCLST(IQ)) THEN
               WRITE (IOUT,535)
535            FORMAT(' ALL CELLS INCLUDED IN THIS OBSERVATION ARE DRY')
            ENDIF
         ENDIF
   !
   !x------END OF LOOP FOR OBSERVATION TIMES IN ONE CELL GROUP
600   CONTINUE
   !
   !-------UPDATE COUNTERS
700   NC = NC + NQCLST(IQ)
      NT1 = NT2 + 1
   !
   !x------END OF LOOP FOR ALL CELL GROUPS.
800 CONTINUE
   !
   !x------RETURN
   RETURN
END
SUBROUTINE OBS2STR7OT(IGRID)
   !     ******************************************************************
   !     WRITE ALL OBSERVATIONS TO LISTING FILE.
   !     ******************************************************************
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY: IOUT
   USE OBSSTRMODULE
   DOUBLE PRECISION SQ,SUMSQ
   !     ------------------------------------------------------------------
   CALL SOBS2STR7PNT(IGRID)
   !
   !1------WRITE OBSERVATIONS TO LISTING FILE.
   IF(IPRT.NE.0) WRITE(IOUT,17)
17 FORMAT(1X,/,1X,'STREAMFLOW OBSERVATIONS',/,&
   &1X,'OBSERVATION       OBSERVED           SIMULATED',/&
   &1X,'  NAME              VALUE              VALUE',&
   &'             DIFFERENCE',/&
   &1X,'----------------------------------------------',&
   &'----------------------')
   SUMSQ=0.
   DO 100 N=1,NQTST
      DIFF=FLWOBS(N)-FLWSIM(N)
      SQ=DIFF*DIFF
      SUMSQ=SUMSQ+SQ
      IF(IPRT.NE.0) WRITE(IOUT,27) OBSNAM(N),FLWOBS(N),FLWSIM(N),DIFF
27    FORMAT(1X,A,1P,3G20.11)
100 CONTINUE
   WRITE(IOUT,28) SUMSQ
28 FORMAT(1X,/,1X,'STR FLOW SUM OF SQUARED DIFFERENCE:',1P,E15.5)
   !
   !2------WRITE OBSERVATIONS TO SEPARATE FILE.
   IF(IUSTOBSV.GT.0) CALL UOBSSV(IUSTOBSV,NQTST,FLWSIM,&
   &FLWOBS,OBSNAM,0)
   !
   !3------RETURN.
   RETURN
END
SUBROUTINE OBS2STR7DA(IGRID)
   !  Deallocate OBSSTR memory
   USE OBSSTRMODULE
   !
   CALL SOBS2STR7PNT(IGRID)
   DEALLOCATE(NQST)
   DEALLOCATE(NQCST)
   DEALLOCATE(NQTST)
   DEALLOCATE(IUSTOBSV)
   DEALLOCATE(IPRT)
   DEALLOCATE(NQOBST)
   DEALLOCATE(NQCLST)
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
SUBROUTINE SOBS2STR7PNT(IGRID)
   !  Change OBSSTR data to a different grid.
   USE OBSSTRMODULE
   !
   NQST=>OBSSTRDAT(IGRID)%NQST
   NQCST=>OBSSTRDAT(IGRID)%NQCST
   NQTST=>OBSSTRDAT(IGRID)%NQTST
   IUSTOBSV=>OBSSTRDAT(IGRID)%IUSTOBSV
   IPRT=>OBSSTRDAT(IGRID)%IPRT
   NQOBST=>OBSSTRDAT(IGRID)%NQOBST
   NQCLST=>OBSSTRDAT(IGRID)%NQCLST
   IOBTS=>OBSSTRDAT(IGRID)%IOBTS
   FLWSIM=>OBSSTRDAT(IGRID)%FLWSIM
   FLWOBS=>OBSSTRDAT(IGRID)%FLWOBS
   TOFF=>OBSSTRDAT(IGRID)%TOFF
   OTIME=>OBSSTRDAT(IGRID)%OTIME
   QCELL=>OBSSTRDAT(IGRID)%QCELL
   OBSNAM=>OBSSTRDAT(IGRID)%OBSNAM
   !
   RETURN
END
SUBROUTINE SOBS2STR7PSV(IGRID)
   !  Save OBSSTR data for a grid.
   USE OBSSTRMODULE
   !
   OBSSTRDAT(IGRID)%NQST=>NQST
   OBSSTRDAT(IGRID)%NQCST=>NQCST
   OBSSTRDAT(IGRID)%NQTST=>NQTST
   OBSSTRDAT(IGRID)%IUSTOBSV=>IUSTOBSV
   OBSSTRDAT(IGRID)%IPRT=>IPRT
   OBSSTRDAT(IGRID)%NQOBST=>NQOBST
   OBSSTRDAT(IGRID)%NQCLST=>NQCLST
   OBSSTRDAT(IGRID)%IOBTS=>IOBTS
   OBSSTRDAT(IGRID)%FLWSIM=>FLWSIM
   OBSSTRDAT(IGRID)%FLWOBS=>FLWOBS
   OBSSTRDAT(IGRID)%TOFF=>TOFF
   OBSSTRDAT(IGRID)%OTIME=>OTIME
   OBSSTRDAT(IGRID)%QCELL=>QCELL
   OBSSTRDAT(IGRID)%OBSNAM=>OBSNAM
   !
   RETURN
END
