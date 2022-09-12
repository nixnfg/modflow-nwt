MODULE GWFETSMODULE
   INTEGER,SAVE,POINTER   ::NETSOP,IETSCB,NPETS,IETSPF,NETSEG
   INTEGER,      SAVE, DIMENSION(:,:),   POINTER ::IETS
   REAL,         SAVE, DIMENSION(:,:),   POINTER ::ETSR
   REAL,         SAVE, DIMENSION(:,:),   POINTER ::ETSX
   REAL,         SAVE, DIMENSION(:,:),   POINTER ::ETSS
   REAL,         SAVE, DIMENSION(:,:,:), POINTER ::PXDP
   REAL,         SAVE, DIMENSION(:,:,:), POINTER ::PETM
   TYPE GWFETSTYPE
      INTEGER, POINTER   ::NETSOP,IETSCB,NPETS,IETSPF,NETSEG
      INTEGER,       DIMENSION(:,:),   POINTER ::IETS
      REAL,          DIMENSION(:,:),   POINTER ::ETSR
      REAL,          DIMENSION(:,:),   POINTER ::ETSX
      REAL,          DIMENSION(:,:),   POINTER ::ETSS
      REAL,          DIMENSION(:,:,:), POINTER ::PXDP
      REAL,          DIMENSION(:,:,:), POINTER ::PETM
   END TYPE
   TYPE(GWFETSTYPE), SAVE :: GWFETSDAT(10)
END MODULE GWFETSMODULE

SUBROUTINE GWF2ETS7AR(IN,IGRID)
   !     ******************************************************************
   !     ALLOCATE ARRAY STORAGE FOR EVAPOTRANSPIRATION SEGMENTS AND READ
   !     PARAMETER DEFINITIONS
   !     Modified 11/21/2001 to support parameter instances - ERB
   !     Modified 8/17/2009 to support NETSOP=3 - ERB
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:IOUT,NCOL,NROW,IFREFM
   USE GWFETSMODULE, ONLY:NETSOP,IETSCB,NPETS,IETSPF,NETSEG,&
   &IETS,ETSR,ETSX,ETSS,PXDP,PETM
   !
   CHARACTER*4 PTYP
   CHARACTER*200 LINE
   !     ------------------------------------------------------------------
500 FORMAT(1X,/&
   &1X,'ETS7 -- EVAPOTRANSPIRATION SEGMENTS PACKAGE, VERSION 7,',&
   &' 2/28/2006',/,9X,'INPUT READ FROM UNIT ',I4)
510 FORMAT(&
   &1X,I5,' SEGMENTS DEFINE EVAPOTRANSPIRATION RATE FUNCTION')
520 FORMAT(' EVAPOTRANSPIRATION RATE FUNCTION IS LINEAR')
530 FORMAT(&
   &' ERROR: EVAPOTRANSPIRATION RATE FUNCTION MUST CONTAIN AT',/,&
   &' LEAST ONE SEGMENT -- STOP EXECUTION (GWF2ETS7ALP)')
540 FORMAT(1X,'ILLEGAL ET OPTION CODE. SIMULATION ABORTING')
550 FORMAT(1X,'OPTION 1 -- EVAPOTRANSPIRATION FROM TOP LAYER')
560 FORMAT(1X,'OPTION 2 -- EVAPOTRANSPIRATION FROM ONE SPECIFIED',&
   &' NODE IN EACH VERTICAL COLUMN')
564 FORMAT(1X,'OPTION 3 -- EVAPOTRANSPIRATION FROM UPPERMOST ACTIVE ',&
   &'CELL')
570 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
580 FORMAT(1X,I10,' ELEMENTS IN RX ARRAY ARE USED BY ETS')
590 FORMAT(1X,I10,' ELEMENTS IN IR ARRAY ARE USED BY ETS')
   !
   ALLOCATE (NETSOP,IETSCB,NPETS,IETSPF,NETSEG)
   !
   !1------IDENTIFY PACKAGE.
   IETSPF=20
   WRITE(IOUT,500)IN
   !
   !     READ COMMENT LINE(S) (ITEM 0)
   CALL URDCOM(IN,IOUT,LINE)
   !
   !2------READ ET OPTION (NETSOP), UNIT OR FLAG FOR CELL-BY-CELL FLOW
   !       TERMS (IETSCB), NUMBER OF PARAMETERS (NPETS), AND NUMBER OF
   !       SEGMENTS (NETSEG) (ITEM 1)
   IF (IFREFM.EQ.0) THEN
      READ(LINE,'(4I10)') NETSOP,IETSCB,NPETS,NETSEG
   ELSE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NETSOP,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IETSCB,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPETS,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NETSEG,R,IOUT,IN)
   ENDIF
   !
   !3------CHECK TO SEE THAT ET OPTION IS LEGAL.
   IF (NETSOP.GE.1 .AND. NETSOP.LE.3) GO TO 10
   !
   !3A-----OPTION IS ILLEGAL -- PRINT A MESSAGE & ABORT SIMULATION.
   WRITE(IOUT,540)
   CALL USTOP(' ')
   !
   !4------OPTION IS LEGAL -- PRINT THE OPTION CODE.
10 CONTINUE
   IF (NETSOP.EQ.1) WRITE(IOUT,550)
   IF (NETSOP.EQ.2) WRITE(IOUT,560)
   IF (NETSOP.EQ.3) WRITE(IOUT,564) ! Add option 3 ERB 5/8/2009
   !
   !5------IF CELL-BY-CELL FLOWS ARE TO BE SAVED, THEN PRINT UNIT NUMBER.
   IF (IETSCB.GT.0) WRITE(IOUT,570) IETSCB
   !
   !-----PRINT NUMBER OF PARAMETERS TO BE USED
   CALL UPARARRAL(-1,IOUT,LINE,NPETS)
   !
   !     PRINT MESSAGE IDENTIFYING NUMBER OF SEGMENTS IN ET VS. HEAD CURVE
   IF(NETSEG.GT.1) THEN
      WRITE(IOUT,510) NETSEG
   ELSEIF (NETSEG.EQ.1) THEN
      WRITE(IOUT,520)
   ELSE
      WRITE(IOUT,530)
      CALL USTOP(' ')
   ENDIF
   !
   !6------ALLOCATE SPACE FOR THE ARRAYS ETSR, ETSX, ETSS, PXDP, AND PETM.
   ALLOCATE (ETSR(NCOL,NROW))
   ALLOCATE (ETSX(NCOL,NROW))
   ALLOCATE (ETSS(NCOL,NROW))
   IF( NETSEG.GT.1) THEN
      ALLOCATE (PXDP(NCOL,NROW,NETSEG))
      ALLOCATE (PETM(NCOL,NROW,NETSEG))
   ELSE
      ALLOCATE (PXDP(1,1,1))
      ALLOCATE (PETM(1,1,1))
   END IF
   !
   !7------ALLOCATE SPACE FOR LAYER INDICATOR ARRAY (IETS) EVEN IF ET
   !7------OPTION IS NOT 2.
   ALLOCATE (IETS(NCOL,NROW))
   !
   !-------READ NAMED PARAMETERS
   WRITE(IOUT,50) NPETS
50 FORMAT(1X,//1X,I5,' Evapotranspiration segments parameters')
   IF (NPETS.GT.0) THEN
      DO 100 K=1,NPETS
   !         UPARARRRP READS PARAMETER NAME AND DEFINITION (ITEMS 2 AND 3)
         CALL UPARARRRP(IN,IOUT,N,0,PTYP,1,1,0)
         IF(PTYP.NE.'ETS') THEN
            WRITE(IOUT,57)
57          FORMAT(1X,'Parameter type must be ETS')
            CALL USTOP(' ')
         ENDIF
100   CONTINUE
   ENDIF
   !
   !8------RETURN
   CALL SGWF2ETS7PSV(IGRID)
   RETURN
END
SUBROUTINE GWF2ETS7RP(IN,IGRID)
   !     ******************************************************************
   !     READ EVAPOTRANSPIRATION DATA, AND PERFORM SUBSTITUTION USING
   !     PARAMETER VALUES IF ETS PARAMETERS ARE DEFINED
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:NCOL,NROW,IOUT,DELR,DELC,IFREFM
   USE GWFETSMODULE, ONLY:NETSOP,NETSEG,NPETS,IETSPF,&
   &IETS,ETSR,ETSX,ETSS,PXDP,PETM
   !
   CHARACTER*24 ANAME(6)
   DATA ANAME(1) /'   ET LAYER INDEX (IETS)'/
   DATA ANAME(2) /'       ET SURFACE (ETSS)'/
   DATA ANAME(3) /' EVAPOTRANS. RATE (ETSR)'/
   DATA ANAME(4) /' EXTINCTION DEPTH (ETSX)'/
   DATA ANAME(5) /'EXTINCT. DEP. PROPORTION'/
   DATA ANAME(6) /'      ET RATE PROPORTION'/
   !     ------------------------------------------------------------------
   CALL SGWF2ETS7PNT(IGRID)
   !
   !1------READ FLAGS SHOWING WHETHER DATA FROM PREVIOUS STRESS PERIOS ARE
   !       TO BE REUSED.
   IF (NETSEG.GT.1) THEN
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(5I10)') INETSS,INETSR,INETSX,INIETS,INSGDF
      ELSE
         READ(IN,*) INETSS,INETSR,INETSX,INIETS,INSGDF
      ENDIF
   ELSE
      IF(NETSOP.EQ.2) THEN
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(4I10)') INETSS,INETSR,INETSX,INIETS
         ELSE
            READ(IN,*) INETSS,INETSR,INETSX,INIETS
         ENDIF
      ELSE
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(3I10)') INETSS,INETSR,INETSX
         ELSE
            READ(IN,*) INETSS,INETSR,INETSX
         ENDIF
      ENDIF
   ENDIF
   !
   !2------TEST INETSS TO SEE WHERE SURFACE ELEVATION COMES FROM.
   IF (INETSS.LT.0) THEN
   !2A------IF INETSS<0 THEN REUSE SURFACE ARRAY FROM LAST STRESS PERIOD
      WRITE(IOUT,10)
10    FORMAT(1X,/1X,'REUSING ETSS FROM LAST STRESS PERIOD')
   ELSE
   !3-------IF INETSS=>0 THEN CALL MODULE U2DREL TO READ SURFACE.
      CALL U2DREL(ETSS,ANAME(2),NROW,NCOL,0,IN,IOUT)
   ENDIF
   !
   !4------TEST INETSR TO SEE WHERE MAX ET RATE COMES FROM.
   IF (INETSR.LT.0) THEN
   !4A-----IF INETSR<0 THEN REUSE MAX ET RATE.
      WRITE(IOUT,20)
20    FORMAT(1X,/1X,'REUSING ETSR FROM LAST STRESS PERIOD')
   ELSE
   !5------IF INETSR=>0 CALL MODULE U2DREL TO READ MAX ET RATE.
      IF(NPETS.EQ.0) THEN
         CALL U2DREL(ETSR,ANAME(3),NROW,NCOL,0,IN,IOUT)
      ELSE
   !    INETSR is the number of parameters to use this stress period
         CALL PRESET('ETS')
         WRITE(IOUT,30)
30       FORMAT(1X,///1X,&
         &'ETSR array defined by the following parameters:')
         IF (INETSR.EQ.0) THEN
            WRITE(IOUT,35)
35          FORMAT(' ERROR: When parameters are defined for the ETS',&
            &' Package, at least one parameter',/,' must be specified',&
            &' each stress period -- STOP EXECUTION (GWF2ETS7RPSS)')
            CALL USTOP(' ')
         ENDIF
         CALL UPARARRSUB2(ETSR,NCOL,NROW,0,INETSR,IN,IOUT,'ETS',&
         &ANAME(3),'ETS',IETSPF)
      ENDIF
   !
   !6------MULTIPLY MAX ET RATE BY CELL AREA TO GET VOLUMETRIC RATE
      DO 50 IR=1,NROW
         DO 40 IC=1,NCOL
            ETSR(IC,IR)=ETSR(IC,IR)*DELR(IC)*DELC(IR)
40       CONTINUE
50    CONTINUE
   ENDIF
   !
   !7------TEST INETSX TO SEE WHERE EXTINCTION DEPTH COMES FROM
   IF (INETSX.LT.0) THEN
   !7A------IF INETSX<0 REUSE EXTINCTION DEPTH FROM LAST STRESS PERIOD
      WRITE(IOUT,60)
60    FORMAT(1X,/1X,'REUSING ETSX FROM LAST STRESS PERIOD')
   ELSE
   !8-------IF INETSX=>0 CALL MODULE U2DREL TO READ EXTINCTION DEPTH
      CALL U2DREL(ETSX,ANAME(4),NROW,NCOL,0,IN,IOUT)
   ENDIF
   !
   !9------IF OPTION(NETSOP) IS 2 THEN WE NEED AN INDICATOR ARRAY.
   IF (NETSOP.EQ.2) THEN
   !10------IF INIETS<0 THEN REUSE LAYER INDICATOR ARRAY.
      IF (INIETS.LT.0) THEN
         WRITE(IOUT,70)
70       FORMAT(1X,/1X,'REUSING IETS FROM LAST STRESS PERIOD')
      ELSE
   !11------IF INIETS=>0 THEN CALL MODULE U2DINT TO READ INDICATOR ARRAY.
         CALL U2DINT(IETS,ANAME(1),NROW,NCOL,0,IN,IOUT)
      ENDIF
   ENDIF
   !
   !12------IF ET FUNCTION IS SEGMENTED PXDP AND PETM ARRAYS ARE NEEDED.
   IF (NETSEG.GT.1) THEN
   !13------IF INSGDF<0 THEN REUSE PXDP AND PETM ARRAYS.
      IF (INSGDF.LT.0) THEN
         WRITE(IOUT,80)
80       FORMAT(1X,/1X,&
         &'REUSING PXDP AND PETM FROM LAST STRESS PERIOD')
   !14------IF INSGDF=>0 THEN CALL MODULE U2DREL TO READ PXDP AND PETM
   !        ARRAYS.
      ELSE
         DO 90 I = 1,NETSEG-1
            WRITE(IOUT,100) I
            CALL U2DREL(PXDP(:,:,I),ANAME(5),NROW,NCOL,0,IN,IOUT)
            CALL U2DREL(PETM(:,:,I),ANAME(6),NROW,NCOL,0,IN,IOUT)
90       CONTINUE
      ENDIF
   ENDIF
100 FORMAT(/,' PXDP AND PETM ARRAYS FOR INTERSECTION ',I4,&
   &' OF HEAD/ET RELATION:')
   !
   !15-----RETURN
   RETURN
END
SUBROUTINE GWF2ETS7FM(IGRID)
   !     ******************************************************************
   !        ADD EVAPOTRANSPIRATION TO RHS AND HCOF
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,       ONLY:NCOL,NROW,NLAY,RHS,HCOF,IBOUND,HNEW
   USE GWFETSMODULE, ONLY:NETSOP,NETSEG,IETS,ETSR,ETSX,ETSS,PXDP,PETM
   !
   DOUBLE PRECISION HH, SS, XX, DD, PXDP1, PXDP2
   !     ------------------------------------------------------------------
   CALL SGWF2ETS7PNT(IGRID)
   !
   !1------PROCESS EACH HORIZONTAL CELL LOCATION
   DO 30 IR=1,NROW
      DO 20 IC=1,NCOL
   !
         IF (NETSOP.EQ.1) THEN
   !2----------SET THE LAYER INDEX EQUAL TO 1      .
            IL=1
         ELSEIF (NETSOP.EQ.2) THEN
   !3----------GET LAYER INDEX FROM IETS ARRAY
            IL=IETS(IC,IR)
            IF (IL.EQ.0) GO TO 20
         ELSEIF (NETSOP.EQ.3) THEN
   !3A---------FIND UPPERMOST ACTIVE CELL
            IL = 1 ! If stack is inactive, this is the default
            FINDFIRST: DO ILQ=1,NLAY
               IF (IBOUND(IC,IR,ILQ).NE.0) THEN
                  IL = ILQ
                  EXIT FINDFIRST
               ENDIF
            ENDDO FINDFIRST
         ENDIF
   !
   !4------IF THE CELL IS EXTERNAL IGNORE IT.
         IF (IBOUND(IC,IR,IL).GT.0) THEN
            C=ETSR(IC,IR)
            S=ETSS(IC,IR)
            SS=S
            HH=HNEW(IC,IR,IL)
   !
   !5------IF HEAD IN CELL IS GREATER THAN OR EQUAL TO ETSS, ET IS CONSTANT
            IF(HH.GE.SS) THEN
   !
   !5A-----SUBTRACT -ETSR FROM RHS
               RHS(IC,IR,IL)=RHS(IC,IR,IL) + C
            ELSE
   !
   !6------IF DEPTH TO WATER>=EXTINCTION DEPTH THEN ET IS 0
               DD=SS-HH
               X=ETSX(IC,IR)
               XX=X
               IF (DD.LT.XX) THEN
   !7------VARIABLE RANGE. ADD ET TERMS TO BOTH RHS AND HCOF.
   !
                  IF (NETSEG.GT.1) THEN
   !                 DETERMINE WHICH SEGMENT APPLIES BASED ON HEAD, AND
   !                 CALCULATE TERMS TO ADD TO RHS AND HCOF
   !
   !                 SET PROPORTIONS CORRESPONDING TO ETSS ELEVATION
                     PXDP1 = 0.0
                     PETM1 = 1.0
                     DO 10 ISEG = 1,NETSEG
   !                   SET PROPORTIONS CORRESPONDING TO LOWER END OF
   !                   SEGMENT
                        IF (ISEG.LT.NETSEG) THEN
                           PXDP2 = PXDP(IC,IR,ISEG)
                           PETM2 = PETM(IC,IR,ISEG)
                        ELSE
                           PXDP2 = 1.0
                           PETM2 = 0.0
                        ENDIF
                        IF (DD.LE.PXDP2*XX) THEN
   !                     HEAD IS IN DOMAIN OF THIS SEGMENT
                           GOTO 15
                        ENDIF
   !                   PROPORTIONS AT LOWER END OF SEGMENT WILL BE FOR
   !                   UPPER END OF SEGMENT NEXT TIME THROUGH LOOP
                        PXDP1 = PXDP2
                        PETM1 = PETM2
10                   CONTINUE
15                   CONTINUE
   !                 CALCULATE TERMS TO ADD TO RHS AND HCOF BASED ON
   !                 SEGMENT THAT APPLIES AT HEAD ELEVATION
                     THCOF = -(PETM1-PETM2)*C/((PXDP2-PXDP1)*X)
                     TRHS = THCOF*(S-PXDP1*X) + PETM1*C
                  ELSE
   !                 CALCULATE TERMS TO ADD TO RHS AND HCOF BASED ON SIMPLE
   !                 LINEAR RELATION OF ET VS. HEAD
                     TRHS = C-C*S/X
                     THCOF = -C/X
                  ENDIF
                  RHS(IC,IR,IL)=RHS(IC,IR,IL)+TRHS
                  HCOF(IC,IR,IL)=HCOF(IC,IR,IL)+THCOF
               ENDIF
            ENDIF
         ENDIF
20    CONTINUE
30 CONTINUE
   !
   !8------RETURN
   RETURN
END
SUBROUTINE GWF2ETS7BD(KSTP,KPER,IGRID)
   !     ******************************************************************
   !     CALCULATE VOLUMETRIC BUDGET FOR EVAPOTRANSPIRATION SEGMENTS
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL ,      ONLY: IOUT,HNEW,IBOUND,BUFF,NCOL,NROW,NLAY
   USE GWFBASMODULE, ONLY: MSUM,VBNM,VBVL,PERTIM,TOTIM,DELT,ICBCFL
   USE GWFETSMODULE, ONLY: NETSOP,IETSCB,NETSEG,IETS,ETSR,ETSX,ETSS,&
   &PXDP,PETM
   !
   DOUBLE PRECISION RATOUT, QQ, HH, SS, DD, XX, HHCOF, RRHS,&
   &PXDP1, PXDP2
   CHARACTER*16 TEXT
   DATA TEXT /'     ET SEGMENTS'/
   !     ------------------------------------------------------------------
   CALL SGWF2ETS7PNT(IGRID)
   !
   !1------CLEAR THE RATE ACCUMULATOR.
   ZERO=0.
   RATOUT=ZERO
   !
   !2------SET CELL-BY-CELL BUDGET SAVE FLAG (IBD) AND CLEAR THE BUFFER.
   IBD=0
   IF(IETSCB.GT.0) IBD=ICBCFL
   DO 30 IL=1,NLAY
      DO 20 IR=1,NROW
         DO 10 IC=1,NCOL
            BUFF(IC,IR,IL)=ZERO
10       CONTINUE
20    CONTINUE
30 CONTINUE
   !
   !3------PROCESS EACH HORIZONTAL CELL LOCATION.
   DO 70 IR=1,NROW
      DO 60 IC=1,NCOL
   !
         IF (NETSOP.EQ.1) THEN
   !2----------SET THE LAYER INDEX EQUAL TO 1      .
            IL=1
         ELSEIF (NETSOP.EQ.2) THEN
   !3----------GET LAYER INDEX FROM IETS ARRAY
            IL=IETS(IC,IR)
            IF (IL.EQ.0) GO TO 60
         ELSEIF (NETSOP.EQ.3) THEN
   !3A---------FIND UPPERMOST ACTIVE CELL
            IL = 1 ! If stack is inactive, this is the default
            FINDFIRST: DO ILQ=1,NLAY
               IF (IBOUND(IC,IR,ILQ).NE.0) THEN
                  IL = ILQ
                  EXIT FINDFIRST
               ENDIF
            ENDDO FINDFIRST
            IETS(IC,IR) = IL
         ENDIF
   !
   !6------IF CELL IS EXTERNAL THEN IGNORE IT.
         IF (IBOUND(IC,IR,IL).GT.0) THEN
            C=ETSR(IC,IR)
            S=ETSS(IC,IR)
            SS=S
            HH=HNEW(IC,IR,IL)
   !
   !7------IF HEAD IN CELL => ETSS,SET Q=MAX ET RATE.
            IF (HH.GE.SS) THEN
               QQ=-C
            ELSE
   !
   !8------IF DEPTH=>EXTINCTION DEPTH, ET IS 0.
               X=ETSX(IC,IR)
               XX=X
               DD=SS-HH
               IF (DD.LT.XX) THEN
   !9------VARIABLE RANGE.  CALCULATE Q DEPENDING ON NUMBER OF SEGMENTS
   !
                  IF (NETSEG.GT.1) THEN
   !                 DETERMINE WHICH SEGMENT APPLIES BASED ON HEAD, AND
   !                 CALCULATE TERMS TO ADD TO RHS AND HCOF
   !
   !                 SET PROPORTIONS CORRESPONDING TO ETSS ELEVATION
                     PXDP1 = 0.0
                     PETM1 = 1.0
                     DO 40 ISEG = 1,NETSEG
   !                   SET PROPORTIONS CORRESPONDING TO LOWER END OF
   !                   SEGMENT
                        IF (ISEG.LT.NETSEG) THEN
                           PXDP2 = PXDP(IC,IR,ISEG)
                           PETM2 = PETM(IC,IR,ISEG)
                        ELSE
                           PXDP2 = 1.0
                           PETM2 = 0.0
                        ENDIF
                        IF (DD.LE.PXDP2*XX) THEN
   !                     HEAD IS IN DOMAIN OF THIS SEGMENT
                           GOTO 50
                        ENDIF
   !                   PROPORTIONS AT LOWER END OF SEGMENT WILL BE FOR
   !                   UPPER END OF SEGMENT NEXT TIME THROUGH LOOP
                        PXDP1 = PXDP2
                        PETM1 = PETM2
40                   CONTINUE
50                   CONTINUE
   !9------CALCULATE ET RATE BASED ON SEGMENT THAT APPLIES AT HEAD
   !9------ELEVATION
                     HHCOF = -(PETM1-PETM2)*C/((PXDP2-PXDP1)*X)
                     RRHS = -HHCOF*(S-PXDP1*X) - PETM1*C
                  ELSE
   !10-----SIMPLE LINEAR RELATION.  Q=-ETSR*(HNEW-(ETSS-ETSX))/ETSX, WHICH
   !10-----IS FORMULATED AS Q= -HNEW*ETSR/ETSX + (ETSR*ETSS/ETSX -ETSR).
                     HHCOF = -C/X
                     RRHS = (C*S/X) - C
                  ENDIF
                  QQ = HH*HHCOF + RRHS
               ELSE
                  QQ = 0.0
               ENDIF
            ENDIF
   !
   !10-----ACCUMULATE TOTAL FLOW RATE.
            Q=QQ
            RATOUT=RATOUT-QQ
   !
   !11-----ADD Q TO BUFFER.
            BUFF(IC,IR,IL)=Q
         ENDIF
60    CONTINUE
70 CONTINUE
   !
   !12-----IF CELL-BY-CELL FLOW TO BE SAVED, CALL APPROPRIATE UTILITY
   !12-----MODULE SAVE THEM.
   IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IETSCB,BUFF,NCOL,NROW,&
   &NLAY,IOUT)
   IF(IBD.EQ.2) CALL UBDSV3(KSTP,KPER,TEXT,IETSCB,BUFF,IETS,NETSOP,&
   &NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
   !
   !13-----MOVE TOTAL ET RATE INTO VBVL FOR PRINTING BY BAS1OT.
   ROUT=RATOUT
   VBVL(3,MSUM)=ZERO
   VBVL(4,MSUM)=ROUT
   !
   !14-----ADD ET(ET_RATE TIMES STEP LENGTH) TO VBVL.
   VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
   !
   !15-----MOVE BUDGET TERM LABELS TO VBNM FOR PRINT BY MODULE BAS1OT.
   VBNM(MSUM)=TEXT
   !
   !16-----INCREMENT BUDGET TERM COUNTER.
   MSUM=MSUM+1
   !
   !17-----RETURN.
   RETURN
END
SUBROUTINE GWF2ETS7DA(IGRID)
   !  Deallocate ETS MEMORY
   USE GWFETSMODULE
   !
   CALL SGWF2ETS7PNT(IGRID)
   DEALLOCATE(NETSOP)
   DEALLOCATE(IETSCB)
   DEALLOCATE(NPETS)
   DEALLOCATE(IETSPF)
   DEALLOCATE(NETSEG)
   DEALLOCATE(IETS)
   DEALLOCATE(ETSR)
   DEALLOCATE(ETSX)
   DEALLOCATE(ETSS)
   DEALLOCATE(PXDP)
   DEALLOCATE(PETM)
   !
   RETURN
END
SUBROUTINE SGWF2ETS7PNT(IGRID)
   !  Change ETS data to a different grid.
   USE GWFETSMODULE
   !
   NETSOP=>GWFETSDAT(IGRID)%NETSOP
   IETSCB=>GWFETSDAT(IGRID)%IETSCB
   NPETS=>GWFETSDAT(IGRID)%NPETS
   IETSPF=>GWFETSDAT(IGRID)%IETSPF
   NETSEG=>GWFETSDAT(IGRID)%NETSEG
   IETS=>GWFETSDAT(IGRID)%IETS
   ETSR=>GWFETSDAT(IGRID)%ETSR
   ETSX=>GWFETSDAT(IGRID)%ETSX
   ETSS=>GWFETSDAT(IGRID)%ETSS
   PXDP=>GWFETSDAT(IGRID)%PXDP
   PETM=>GWFETSDAT(IGRID)%PETM
   !
   RETURN
END
SUBROUTINE SGWF2ETS7PSV(IGRID)
   !  Save ETS data for a grid.
   USE GWFETSMODULE
   !
   GWFETSDAT(IGRID)%NETSOP=>NETSOP
   GWFETSDAT(IGRID)%IETSCB=>IETSCB
   GWFETSDAT(IGRID)%NPETS=>NPETS
   GWFETSDAT(IGRID)%IETSPF=>IETSPF
   GWFETSDAT(IGRID)%NETSEG=>NETSEG
   GWFETSDAT(IGRID)%IETS=>IETS
   GWFETSDAT(IGRID)%ETSR=>ETSR
   GWFETSDAT(IGRID)%ETSX=>ETSX
   GWFETSDAT(IGRID)%ETSS=>ETSS
   GWFETSDAT(IGRID)%PXDP=>PXDP
   GWFETSDAT(IGRID)%PETM=>PETM
   !
   RETURN
END
