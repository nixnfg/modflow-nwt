      MODULE GWFGHBMODULE
        INTEGER,SAVE,POINTER  ::NBOUND,MXBND,NGHBVL,IGHBCB,IPRGHB
        INTEGER,SAVE,POINTER  ::NPGHB,IGHBPB,NNPGHB
        CHARACTER(LEN=16),SAVE, DIMENSION(:),   POINTER     ::GHBAUX
        REAL,             SAVE, DIMENSION(:,:), POINTER     ::BNDS
      TYPE GWFGHBTYPE
        INTEGER,POINTER  ::NBOUND,MXBND,NGHBVL,IGHBCB,IPRGHB
        INTEGER,POINTER  ::NPGHB,IGHBPB,NNPGHB
        CHARACTER(LEN=16), DIMENSION(:),   POINTER     ::GHBAUX
        REAL,              DIMENSION(:,:), POINTER     ::BNDS
      END TYPE
      TYPE(GWFGHBTYPE), SAVE:: GWFGHBDAT(10)
      END MODULE GWFGHBMODULE


      SUBROUTINE GWF2GHB7AR(IN,IGRID)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE AND READ PARAMETER DEFINITIONS FOR GHB
C     PACKAGE
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM
      USE GWFGHBMODULE, ONLY:NBOUND,MXBND,NGHBVL,IGHBCB,IPRGHB,NPGHB,
     1                       IGHBPB,NNPGHB,GHBAUX,BNDS
C
      CHARACTER*200 LINE
C     ------------------------------------------------------------------
      ALLOCATE(NBOUND,MXBND,NGHBVL,IGHBCB,IPRGHB)
      ALLOCATE(NPGHB,IGHBPB,NNPGHB)
C
C1------IDENTIFY PACKAGE AND INITIALIZE NBOUND.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'GHB -- GENERAL-HEAD BOUNDARY PACKAGE, VERSION 7',
     1   ', 5/2/2005',/,9X,'INPUT READ FROM UNIT ',I4)
      NBOUND=0
      NNPGHB=0
C
C2------READ MAXIMUM NUMBER OF GHB'S AND UNIT OR FLAG FOR
C2------CELL-BY-CELL FLOW TERMS.
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARLSTAL(IN,IOUT,LINE,NPGHB,MXPB)
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(2I10)') MXACTB,IGHBCB
         LLOC=21
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXACTB,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IGHBCB,R,IOUT,IN)
      END IF
      WRITE(IOUT,3) MXACTB
    3 FORMAT(1X,'MAXIMUM OF ',I6,' ACTIVE GHB CELLS AT ONE TIME')
      IF(IGHBCB.LT.0) WRITE(IOUT,7)
    7 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
      IF(IGHBCB.GT.0) WRITE(IOUT,8) IGHBCB
    8 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
C
C3------READ AUXILIARY VARIABLES AND PRINT OPTION.
      ALLOCATE (GHBAUX(20))
      NAUX=0
      IPRGHB=1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(NAUX.LT.5) THEN
            NAUX=NAUX+1
            GHBAUX(NAUX)=LINE(ISTART:ISTOP)
            WRITE(IOUT,12) GHBAUX(NAUX)
   12       FORMAT(1X,'AUXILIARY GHB VARIABLE: ',A)
         END IF
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
         WRITE(IOUT,13)
   13    FORMAT(1X,'LISTS OF GENERAL-HEAD BOUNDARY CELLS WILL NOT BE',
     &          ' PRINTED')
         IPRGHB = 0
         GO TO 10
      END IF
C3A-----THERE ARE FIVE INPUT DATA VALUES PLUS ONE LOCATION FOR
C3A-----CELL-BY-CELL FLOW.
      NGHBVL=6+NAUX
C
C4------ALLOCATE SPACE FOR THE BNDS ARRAY.
      IGHBPB=MXACTB+1
      MXBND=MXACTB+MXPB
      ALLOCATE (BNDS(NGHBVL,MXBND))
C
C-------READ NAMED PARAMETERS.
      WRITE(IOUT,1000) NPGHB
 1000 FORMAT(1X,//1X,I5,' GHB parameters')
      IF(NPGHB.GT.0) THEN
        NAUX=NGHBVL-6
        LSTSUM=IGHBPB
        DO 120 K=1,NPGHB
          LSTBEG=LSTSUM
          CALL UPARLSTRP(LSTSUM,MXBND,IN,IOUT,IP,'GHB','GHB',1,
     &                  NUMINST)
          NLST=LSTSUM-LSTBEG
          IF (NUMINST.EQ.0) THEN
C5A-----READ LIST OF CELLS WITHOUT INSTANCES.
            CALL ULSTRD(NLST,BNDS,LSTBEG,NGHBVL,MXBND,1,IN,IOUT,
     &      'BOUND. NO. LAYER   ROW   COL     STAGE    STRESS FACTOR',
     &      GHBAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRGHB)
          ELSE
C5B-----READ INSTANCES
            NINLST=NLST/NUMINST
            DO 110 I=1,NUMINST
            CALL UINSRP(I,IN,IOUT,IP,IPRGHB)
            CALL ULSTRD(NINLST,BNDS,LSTBEG,NGHBVL,MXBND,1,IN,IOUT,
     &      'BOUND. NO. LAYER   ROW   COL     STAGE    STRESS FACTOR',
     &      GHBAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRGHB)
            LSTBEG=LSTBEG+NINLST
  110       CONTINUE
          END IF
  120   CONTINUE
      END IF
C
C6------RETURN
      CALL SGWF2GHB7PSV(IGRID)
      RETURN
      END
      SUBROUTINE GWF2GHB7RP(IN,IGRID)
C     ******************************************************************
C     READ GHB HEAD, CONDUCTANCE AND BOTTOM ELEVATION
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,BOTM,LBOTM,
     1                       IBOUND
      USE GWFGHBMODULE, ONLY:NBOUND,MXBND,NGHBVL,IPRGHB,NPGHB,
     1                       IGHBPB,NNPGHB,GHBAUX,BNDS
C     ------------------------------------------------------------------
      CALL SGWF2GHB7PNT(IGRID)
C
C1------READ ITMP (NUMBER OF GHB'S OR FLAG TO REUSE DATA) AND
C1------NUMBER OF PARAMETERS.
      IF(NPGHB.GT.0) THEN
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(2I10)') ITMP,NP
         ELSE
            READ(IN,*) ITMP,NP
         END IF
      ELSE
         NP=0
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10)') ITMP
         ELSE
            READ(IN,*) ITMP
         END IF
      END IF
C
C------CALCULATE SOME CONSTANTS
      NAUX=NGHBVL-6
      IOUTU = IOUT
      IF (IPRGHB.EQ.0) IOUTU=-IOUT
C
C2------DETERMINE THE NUMBER OF NON-PARAMETER GHB'S.
      IF(ITMP.LT.0) THEN
         WRITE(IOUT,7)
    7    FORMAT(1X,/1X,
     1   'REUSING NON-PARAMETER GHB CELLS FROM LAST STRESS PERIOD')
      ELSE
         NNPGHB=ITMP
      END IF
C
C3------IF THERE ARE NEW NON-PARAMETER GHB'S, READ THEM.
      MXACTB=IGHBPB-1
      IF(ITMP.GT.0) THEN
         IF(NNPGHB.GT.MXACTB) THEN
            WRITE(IOUT,99) NNPGHB,MXACTB
   99       FORMAT(1X,/1X,'THE NUMBER OF ACTIVE GHB CELLS (',I6,
     1                     ') IS GREATER THAN MXACTB(',I6,')')
            CALL USTOP(' ')
         END IF
         CALL ULSTRD(NNPGHB,BNDS,1,NGHBVL,MXBND,1,IN,IOUT,
     1      'BOUND. NO. LAYER   ROW   COL     STAGE      CONDUCTANCE',
     2      GHBAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRGHB)
      END IF
      NBOUND=NNPGHB
C
C1C-----IF THERE ARE ACTIVE GHB PARAMETERS, READ THEM AND SUBSTITUTE
      CALL PRESET('GHB')
      IF(NP.GT.0) THEN
         NREAD=NGHBVL-1
         DO 30 N=1,NP
         CALL UPARLSTSUB(IN,'GHB',IOUTU,'GHB',BNDS,NGHBVL,MXBND,NREAD,
     1                MXACTB,NBOUND,5,5,
     2      'BOUND. NO. LAYER   ROW   COL     STAGE      CONDUCTANCE',
     3            GHBAUX,20,NAUX)
   30    CONTINUE
      END IF
C4
      ierr = 0
      DO 100 L=1,NBOUND
C
C5------GET COLUMN, ROW AND LAYER OF CELL CONTAINING BOUNDARY.
      IL=BNDS(1,L)
      IR=BNDS(2,L)
      IC=BNDS(3,L)
C
C6------IF THE CELL IS EXTERNAL THEN SKIP IT.
      IF(IBOUND(IC,IR,IL).LE.0) GO TO 100
C
C7------SINCE THE CELL IS INTERNAL GET THE BOUNDARY DATA.
      HB=BNDS(4,L)
      BOT=BOTM(IC,IR,LBOTM(IL))
      IF ( HB.LT.BOT ) THEN
        WRITE(IOUT,103)IC,IR,IL,HB,BOT
        ierr = 1
      END IF
  103 FORMAT('GHB HEAD SET TO BELOW CELL BOTTOM. MODEL STOPPING. ',
     +                'CELL WITH ERROR (IC,IR,IL): ',3I5,/,
     +                'HEAD:',F10.3,'BOTTOM:',F10.3)

  100 CONTINUE
      IF ( ierr==1 ) CALL USTOP(' ')
C
C3------PRINT NUMBER OF GHB'S IN CURRENT STRESS PERIOD.
      WRITE (IOUT,101) NBOUND
  101 FORMAT(1X,/1X,I6,' GHB CELLS')
C
C8------RETURN.
  260 RETURN
      END
      SUBROUTINE GWF2GHB7FM(IGRID)
C     ******************************************************************
C     ADD GHB TERMS TO RHS AND HCOF
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,HCOF
      USE GWFGHBMODULE, ONLY:NBOUND,BNDS
C     ------------------------------------------------------------------
      CALL SGWF2GHB7PNT(IGRID)
C
C1------IF NBOUND<=0 THEN THERE ARE NO GENERAL HEAD BOUNDS. RETURN.
      IF(NBOUND.LE.0) RETURN
C
C2------PROCESS EACH ENTRY IN THE GENERAL HEAD BOUND LIST (BNDS).
      DO 100 L=1,NBOUND
C
C3------GET COLUMN, ROW AND LAYER OF CELL CONTAINING BOUNDARY.
      IL=BNDS(1,L)
      IR=BNDS(2,L)
      IC=BNDS(3,L)
C
C4------IF THE CELL IS EXTERNAL THEN SKIP IT.
      IF(IBOUND(IC,IR,IL).LE.0) GO TO 100
C
C5------SINCE THE CELL IS INTERNAL GET THE BOUNDARY DATA.
      HB=BNDS(4,L)
      C=BNDS(5,L)
C
C6------ADD TERMS TO RHS AND HCOF.
      HCOF(IC,IR,IL)=HCOF(IC,IR,IL)-C
      RHS(IC,IR,IL)=RHS(IC,IR,IL)-C*HB
  100 CONTINUE
C
C7------RETURN.
      RETURN
      END
      SUBROUTINE GWF2GHB7BD(KSTP,KPER,IGRID)
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR GHB
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,HNEW,BUFF
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      USE GWFGHBMODULE,ONLY:NBOUND,IGHBCB,BNDS,NGHBVL,GHBAUX
C
      DOUBLE PRECISION CCGHB,CHB,RATIN,RATOUT,RRATE
      CHARACTER*16 TEXT
      DATA TEXT /' HEAD DEP BOUNDS'/
C     ------------------------------------------------------------------
      CALL SGWF2GHB7PNT(IGRID)
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD) AND
C1------ACCUMULATORS (RATIN AND RATOUT).
      ZERO=0.
      RATOUT=ZERO
      RATIN=ZERO
      IBD=0
      IF(IGHBCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IGHBCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IF(IBD.EQ.2) THEN
         NAUX=NGHBVL-6
         IF(IAUXSV.EQ.0) NAUX=0
         CALL UBDSV4(KSTP,KPER,TEXT,NAUX,GHBAUX,IGHBCB,NCOL,NROW,NLAY,
     1          NBOUND,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 IL=1,NLAY
      DO 50 IR=1,NROW
      DO 50 IC=1,NCOL
      BUFF(IC,IR,IL)=ZERO
50    CONTINUE
C
C4------IF NO BOUNDARIES, SKIP FLOW CALCULATIONS.
      IF(NBOUND.EQ.0) GO TO 200
C
C5------LOOP THROUGH EACH BOUNDARY CALCULATING FLOW.
      DO 100 L=1,NBOUND
C
C5A-----GET LAYER, ROW AND COLUMN OF EACH GENERAL HEAD BOUNDARY.
      IL=BNDS(1,L)
      IR=BNDS(2,L)
      IC=BNDS(3,L)
      RATE=ZERO
C
C5B-----IF CELL IS NO-FLOW OR CONSTANT-HEAD, THEN IGNORE IT.
      IF(IBOUND(IC,IR,IL).LE.0) GO TO 99
C
C5C-----GET PARAMETERS FROM BOUNDARY LIST.
      HB=BNDS(4,L)
      C=BNDS(5,L)
      CCGHB=C
C
C5D-----CALCULATE THE FOW RATE INTO THE CELL.
      CHB=C*HB
      RRATE=CHB - CCGHB*HNEW(IC,IR,IL)
      RATE=RRATE
C
C5E-----PRINT THE INDIVIDUAL RATES IF REQUESTED(IGHBCB<0).
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
         WRITE(IOUT,62) L,IL,IR,IC,RATE
   62    FORMAT(1X,'BOUNDARY ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5,'   RATE ',1PG15.6)
         IBDLBL=1
      END IF
C
C5F-----ADD RATE TO BUFFER.
      BUFF(IC,IR,IL)=BUFF(IC,IR,IL)+RATE
C
C5G-----SEE IF FLOW IS INTO AQUIFER OR OUT OF AQUIFER.
      IF(RATE.LT.ZERO) THEN
C
C5H------FLOW IS OUT OF AQUIFER SUBTRACT RATE FROM RATOUT.
        RATOUT=RATOUT-RRATE
      ELSE
C
C5I-----FLOW IS INTO AQIFER; ADD RATE TO RATIN.
        RATIN=RATIN+RRATE
      END IF
C
C5J-----IF SAVING CELL-BY-CELL FLOWS IN LIST, WRITE FLOW.  ALSO
C5J-----FLOW TO BNDS.
   99 IF(IBD.EQ.2) CALL UBDSVB(IGHBCB,NCOL,NROW,IC,IR,IL,RATE,
     1                  BNDS(:,L),NGHBVL,NAUX,6,IBOUND,NLAY)
      BNDS(NGHBVL,L)=RATE
  100 CONTINUE
C
C6------IF CELL-BY-CELL TERMS WILL BE SAVED AS A 3-D ARRAY, THEN CALL
C6------UTILITY MODULE UBUDSV TO SAVE THEM.
      IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IGHBCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
C
C7------MOVE RATES, VOLUMES AND LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVL(3,MSUM)=RIN
      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
      VBVL(4,MSUM)=ROUT
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBNM(MSUM)=TEXT
C
C8------INCREMENT THE BUDGET TERM COUNTER.
      MSUM=MSUM+1
C
C9------RETURN.
      RETURN
      END
      SUBROUTINE GWF2GHB7DA(IGRID)
C  Deallocate GHB MEMORY
      USE GWFGHBMODULE
C
        CALL SGWF2GHB7PNT(IGRID)
        DEALLOCATE(NBOUND)
        DEALLOCATE(MXBND)
        DEALLOCATE(NGHBVL)
        DEALLOCATE(IGHBCB)
        DEALLOCATE(IPRGHB)
        DEALLOCATE(NPGHB)
        DEALLOCATE(IGHBPB)
        DEALLOCATE(NNPGHB)
        DEALLOCATE(GHBAUX)
        DEALLOCATE(BNDS)
C
      RETURN
      END
      SUBROUTINE SGWF2GHB7PNT(IGRID)
C  Change GHB data to a different grid.
      USE GWFGHBMODULE
C
        NBOUND=>GWFGHBDAT(IGRID)%NBOUND
        MXBND=>GWFGHBDAT(IGRID)%MXBND
        NGHBVL=>GWFGHBDAT(IGRID)%NGHBVL
        IGHBCB=>GWFGHBDAT(IGRID)%IGHBCB
        IPRGHB=>GWFGHBDAT(IGRID)%IPRGHB
        NPGHB=>GWFGHBDAT(IGRID)%NPGHB
        IGHBPB=>GWFGHBDAT(IGRID)%IGHBPB
        NNPGHB=>GWFGHBDAT(IGRID)%NNPGHB
        GHBAUX=>GWFGHBDAT(IGRID)%GHBAUX
        BNDS=>GWFGHBDAT(IGRID)%BNDS
C
      RETURN
      END
      SUBROUTINE SGWF2GHB7PSV(IGRID)
C  Save GHB data for a grid.
      USE GWFGHBMODULE
C
        GWFGHBDAT(IGRID)%NBOUND=>NBOUND
        GWFGHBDAT(IGRID)%MXBND=>MXBND
        GWFGHBDAT(IGRID)%NGHBVL=>NGHBVL
        GWFGHBDAT(IGRID)%IGHBCB=>IGHBCB
        GWFGHBDAT(IGRID)%IPRGHB=>IPRGHB
        GWFGHBDAT(IGRID)%NPGHB=>NPGHB
        GWFGHBDAT(IGRID)%IGHBPB=>IGHBPB
        GWFGHBDAT(IGRID)%NNPGHB=>NNPGHB
        GWFGHBDAT(IGRID)%GHBAUX=>GHBAUX
        GWFGHBDAT(IGRID)%BNDS=>BNDS
C
      RETURN
      END
