MODULE GWFRCHMODULE
   INTEGER, SAVE, POINTER                 ::NRCHOP,IRCHCB
   INTEGER, SAVE, POINTER                 ::NPRCH,IRCHPF
   REAL,    SAVE,   DIMENSION(:,:),  POINTER      ::RECH
   INTEGER, SAVE,   DIMENSION(:,:),  POINTER      ::IRCH
   TYPE GWFRCHTYPE
      INTEGER,  POINTER                 ::NRCHOP,IRCHCB
      INTEGER,  POINTER                 ::NPRCH,IRCHPF
      REAL,       DIMENSION(:,:),  POINTER      ::RECH
      INTEGER,    DIMENSION(:,:),  POINTER      ::IRCH
   END TYPE
   TYPE(GWFRCHTYPE), SAVE ::GWFRCHDAT(10)
END MODULE GWFRCHMODULE


SUBROUTINE GWF2RCH7AR(IN,IGRID)
   !     ******************************************************************
   !     ALLOCATE ARRAY STORAGE FOR RECHARGE
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:IOUT,NCOL,NROW,IFREFM
   USE GWFRCHMODULE,ONLY:NRCHOP,IRCHCB,NPRCH,IRCHPF,RECH,IRCH
   !
   CHARACTER*200 LINE
   CHARACTER*4 PTYP
   !     ------------------------------------------------------------------
   !
   !1-------ALLOCATE SCALAR VARIABLES.
   ALLOCATE(NRCHOP,IRCHCB)
   ALLOCATE(NPRCH,IRCHPF)
   !
   !2------IDENTIFY PACKAGE.
   IRCHPF=0
   WRITE(IOUT,1)IN
1  FORMAT(1X,/1X,'RCH -- RECHARGE PACKAGE, VERSION 7, 5/2/2005',&
   &' INPUT READ FROM UNIT ',I4)
   !
   !3------READ NRCHOP AND IRCHCB.
   CALL URDCOM(IN,IOUT,LINE)
   CALL UPARARRAL(IN,IOUT,LINE,NPRCH)
   IF(IFREFM.EQ.0) THEN
      READ(LINE,'(2I10)') NRCHOP,IRCHCB
   ELSE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NRCHOP,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IRCHCB,R,IOUT,IN)
   END IF
   !
   !4------CHECK TO SEE THAT OPTION IS LEGAL.
   IF(NRCHOP.LT.1.OR.NRCHOP.GT.3) THEN
      WRITE(IOUT,8) NRCHOP
8     FORMAT(1X,'ILLEGAL RECHARGE OPTION CODE (NRCHOP = ',I5,&
      &') -- SIMULATION ABORTING')
      CALL USTOP(' ')
   END IF
   !
   !5------OPTION IS LEGAL -- PRINT OPTION CODE.
   IF(NRCHOP.EQ.1) WRITE(IOUT,201)
201 FORMAT(1X,'OPTION 1 -- RECHARGE TO TOP LAYER')
   IF(NRCHOP.EQ.2) WRITE(IOUT,202)
202 FORMAT(1X,'OPTION 2 -- RECHARGE TO ONE SPECIFIED NODE IN EACH',&
   &' VERTICAL COLUMN')
   IF(NRCHOP.EQ.3) WRITE(IOUT,203)
203 FORMAT(1X,'OPTION 3 -- RECHARGE TO HIGHEST ACTIVE NODE IN',&
   &' EACH VERTICAL COLUMN')
   !
   !6------IF CELL-BY-CELL FLOWS ARE TO BE SAVED, THEN PRINT UNIT NUMBER.
   IF(IRCHCB.GT.0) WRITE(IOUT,204) IRCHCB
204 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
   !
   !7------ALLOCATE SPACE FOR THE RECHARGE (RECH) AND INDICATOR (IRCH)
   !7------ARRAYS.
   ALLOCATE (RECH(NCOL,NROW))
   ALLOCATE (IRCH(NCOL,NROW))
   !
   !8------READ NAMED PARAMETERS
   WRITE(IOUT,5) NPRCH
5  FORMAT(1X,//1X,I5,' Recharge parameters')
   IF(NPRCH.GT.0) THEN
      DO 20 K=1,NPRCH
         CALL UPARARRRP(IN,IOUT,N,0,PTYP,1,1,0)
         IF(PTYP.NE.'RCH') THEN
            WRITE(IOUT,7)
7           FORMAT(1X,'Parameter type must be RCH')
            CALL USTOP(' ')
         END IF
20    CONTINUE
   END IF
   !
   !9------RETURN
   CALL SGWF2RCH7PSV(IGRID)
   RETURN
END
SUBROUTINE GWF2RCH7RP(IN,IGRID)
   !     ******************************************************************
   !     READ RECHARGE DATA FOR STRESS PERIOD
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,DELR,DELC
   USE GWFRCHMODULE,ONLY:NRCHOP,NPRCH,IRCHPF,RECH,IRCH
   !
   CHARACTER*24 ANAME(2)
   !
   DATA ANAME(1) /'    RECHARGE LAYER INDEX'/
   DATA ANAME(2) /'                RECHARGE'/
   !     ------------------------------------------------------------------
   !
   !1------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2RCH7PNT(IGRID)
   !
   !2------READ FLAGS SHOWING WHETHER DATA IS TO BE REUSED.
   IF(NRCHOP.EQ.2) THEN
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(2I10)') INRECH,INIRCH
      ELSE
         READ(IN,*) INRECH,INIRCH
      END IF
   ELSE
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(I10)') INRECH
      ELSE
         READ(IN,*) INRECH
      END IF
   END IF
   !
   !3------TEST INRECH TO SEE HOW TO DEFINE RECH.
   IF(INRECH.LT.0) THEN
   !
   !3A-----INRECH<0, SO REUSE RECHARGE ARRAY FROM LAST STRESS PERIOD.
      WRITE(IOUT,3)
3     FORMAT(1X,/1X,'REUSING RECH FROM LAST STRESS PERIOD')
   ELSE
   !
   !3B-----INRECH=>0, SO READ RECHARGE RATE.
      IF(NPRCH.EQ.0) THEN
   !
   !3B1----THERE ARE NO PARAMETERS, SO READ RECH USING U2DREL.
         CALL U2DREL(RECH,ANAME(2),NROW,NCOL,0,IN,IOUT)
      ELSE
   !
   !3B2----DEFINE RECH USING PARAMETERS.  INRECH IS THE NUMBER OF
   !3B2----PARAMETERS TO USE THIS STRESS PERIOD.
         CALL PRESET('RCH')
         WRITE(IOUT,33)
33       FORMAT(1X,///1X,&
         &'RECH array defined by the following parameters:')
         IF(INRECH.EQ.0) THEN
            WRITE(IOUT,34)
34          FORMAT(' ERROR: When parameters are defined for the RCH',&
            &' Package, at least one parameter',/,' must be specified',&
            &' each stress period -- STOP EXECUTION (GWF2RCH7RPLL)')
            CALL USTOP(' ')
         END IF
         CALL UPARARRSUB2(RECH,NCOL,NROW,0,INRECH,IN,IOUT,'RCH',&
         &ANAME(2),'RCH',IRCHPF)
      END IF
   !
   !4------MULTIPLY RECHARGE RATE BY CELL AREA TO GET VOLUMETRIC RATE.
      DO 50 IR=1,NROW
         DO 50 IC=1,NCOL
            RECH(IC,IR)=RECH(IC,IR)*DELR(IC)*DELC(IR)
50    CONTINUE
   END IF
   !
   !5------IF NRCHOP=2 THEN A LAYER INDICATOR ARRAY IS NEEDED.  TEST INIRCH
   !5------TO SEE HOW TO DEFINE IRCH.
   IF(NRCHOP.EQ.2) THEN
      IF(INIRCH.LT.0) THEN
   !
   !5A-----INIRCH<0, SO REUSE LAYER INDICATOR ARRAY FROM LAST STRESS PERIOD.
         WRITE(IOUT,2)
2        FORMAT(1X,/1X,'REUSING IRCH FROM LAST STRESS PERIOD')
      ELSE
   !
   !5B-----INIRCH=>0, SO CALL U2DINT TO READ LAYER INDICATOR ARRAY(IRCH)
         CALL U2DINT(IRCH,ANAME(1),NROW,NCOL,0,IN,IOUT)
         DO 57 IR=1,NROW
            DO 57 IC=1,NCOL
               IF(IRCH(IC,IR).LT.1 .OR. IRCH(IC,IR).GT.NLAY) THEN
                  WRITE(IOUT,56) IC,IR,IRCH(IC,IR)
56                FORMAT(1X,/1X,'INVALID LAYER NUMBER IN IRCH FOR COLUMN',I4,&
                  &'  ROW',I4,'  :',I4)
                  CALL USTOP(' ')
               END IF
57       CONTINUE
      END IF
   END IF
   !
   !6------RETURN
   RETURN
END
SUBROUTINE GWF2RCH7FM(IGRID)
   !     ******************************************************************
   !     SUBTRACT RECHARGE FROM RHS
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,RHS
   USE GWFRCHMODULE,ONLY:NRCHOP,RECH,IRCH
   !     ------------------------------------------------------------------
   !
   !1------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2RCH7PNT(IGRID)
   !
   !2------DETERMINE WHICH RECHARGE OPTION.
   IF(NRCHOP.EQ.1) THEN
   !
   !3------NRCHOP IS 1, SO RECHARGE IS IN TOP LAYER. LAYER INDEX IS 1.
      DO 10 IR=1,NROW
         DO 10 IC=1,NCOL
   !
   !3A-----IF CELL IS VARIABLE HEAD, SUBTRACT RECHARGE RATE FROM
   !3A-----RIGHT-HAND-SIDE.
            IF(IBOUND(IC,IR,1).GT.0) RHS(IC,IR,1)=RHS(IC,IR,1)-RECH(IC,IR)
10    CONTINUE
   ELSE IF(NRCHOP.EQ.2) THEN
   !
   !4------NRCHOP IS 2, SO RECHARGE IS INTO LAYER IN INDICATOR ARRAY
      DO 20 IR=1,NROW
         DO 20 IC=1,NCOL
            IL=IRCH(IC,IR)
   !
   !4A-----IF THE CELL IS VARIABLE HEAD, SUBTRACT RECHARGE FROM
   !4A-----RIGHT-HAND-SIDE.
            IF(IL.EQ.0) GO TO 20
            IF(IBOUND(IC,IR,IL).GT.0)&
            &RHS(IC,IR,IL)=RHS(IC,IR,IL)-RECH(IC,IR)
20    CONTINUE
   ELSE
   !
   !5------NRCHOP IS 3, RECHARGE IS INTO HIGHEST VARIABLE-HEAD CELL, EXCEPT
   !5------CANNOT PASS THROUGH CONSTANT HEAD NODE
      DO 30 IR=1,NROW
         DO 30 IC=1,NCOL
            DO 28 IL=1,NLAY
   !
   !5A-----IF CELL IS CONSTANT HEAD MOVE ON TO NEXT HORIZONTAL LOCATION.
               IF(IBOUND(IC,IR,IL).LT.0) GO TO 30
   !
   !5B-----IF THE CELL IS VARIABLE HEAD, SUBTRACT RECHARGE FROM
   !5B-----RIGHT-HAND-SIDE AND MOVE TO NEXT HORIZONTAL LOCATION.
               IF(IBOUND(IC,IR,IL).GT.0) THEN
                  RHS(IC,IR,IL)=RHS(IC,IR,IL)-RECH(IC,IR)
                  GO TO 30
               END IF
28          CONTINUE
30    CONTINUE
   END IF
   !
   !6------RETURN
   RETURN
END
SUBROUTINE GWF2RCH7BD(KSTP,KPER,IGRID)
   !     ******************************************************************
   !     CALCULATE VOLUMETRIC BUDGET FOR RECHARGE
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF
   USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,ICBCFL,DELT,PERTIM,TOTIM
   USE GWFRCHMODULE,ONLY:NRCHOP,IRCHCB,RECH,IRCH
   !
   DOUBLE PRECISION RATIN,RATOUT,QQ
   CHARACTER*16 TEXT
   DATA TEXT /'        RECHARGE'/
   !     ------------------------------------------------------------------
   !
   !1------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2RCH7PNT(IGRID)
   !
   !2------CLEAR THE RATE ACCUMULATORS.
   ZERO=0.
   RATIN=ZERO
   RATOUT=ZERO
   !
   !3------CLEAR THE BUFFER & SET FLAG FOR SAVING CELL-BY-CELL FLOW TERMS.
   DO 2 IL=1,NLAY
      DO 2 IR=1,NROW
         DO 2 IC=1,NCOL
            BUFF(IC,IR,IL)=ZERO
2  CONTINUE
   IBD=0
   IF(IRCHCB.GT.0) IBD=ICBCFL
   !
   !4------DETERMINE THE RECHARGE OPTION.
   IF(NRCHOP.EQ.1) THEN
   !
   !5------NRCHOP=1, SO RECH GOES INTO LAYER 1. PROCESS EACH HORIZONTAL
   !5------CELL LOCATION.
      DO 10 IR=1,NROW
         DO 10 IC=1,NCOL
   !
   !5A-----IF CELL IS VARIABLE HEAD, THEN DO BUDGET FOR IT.
            IF(IBOUND(IC,IR,1).GT.0) THEN
               Q=RECH(IC,IR)
               QQ=Q
   !
   !5B-----ADD RECH TO BUFF.
               BUFF(IC,IR,1)=Q
   !
   !5C-----IF RECH POSITIVE ADD IT TO RATIN, ELSE ADD IT TO RATOUT.
               IF(Q.GE.ZERO) THEN
                  RATIN=RATIN+QQ
               ELSE
                  RATOUT=RATOUT-QQ
               END IF
            END IF
10    CONTINUE
   ELSE IF(NRCHOP.EQ.2) THEN
   !
   !6------NRCHOP=2, RECH IS IN LAYER SPECIFIED IN INDICATOR ARRAY(IRCH).
   !6------PROCESS EACH HORIZONTAL CELL LOCATION.
      DO 20 IR=1,NROW
         DO 20 IC=1,NCOL
   !
   !6A-----GET LAYER INDEX FROM INDICATOR ARRAY(IRCH).
            IL=IRCH(IC,IR)
   !
   !6B-----IF CELL IS VARIABLE HEAD, THEN DO BUDGET FOR IT.
            IF(IL.EQ.0) GO TO 20
            IF(IBOUND(IC,IR,IL).GT.0) THEN
               Q=RECH(IC,IR)
               QQ=Q
   !
   !6C-----ADD RECHARGE TO BUFF.
               BUFF(IC,IR,IL)=Q
   !
   !6D-----IF RECHARGE IS POSITIVE ADD TO RATIN, ELSE ADD IT TO RATOUT.
               IF(Q.GE.ZERO) THEN
                  RATIN=RATIN+QQ
               ELSE
                  RATOUT=RATOUT-QQ
               END IF
            END IF
20    CONTINUE
   ELSE
   !
   !7------NRCHOP=3; RECHARGE IS INTO HIGHEST CELL IN A VERTICAL COLUMN
   !7------THAT IS NOT NO FLOW.  PROCESS EACH HORIZONTAL CELL LOCATION.
      DO 30 IR=1,NROW
         DO 29 IC=1,NCOL
   !
   !7A-----INITIALIZE IRCH TO 1, AND LOOP THROUGH CELLS IN A VERTICAL
   !7A-----COLUMN TO FIND WHERE TO PLACE RECHARGE.
            IRCH(IC,IR)=1
            DO 28 IL=1,NLAY
   !
   !7B-----IF CELL IS CONSTANT HEAD, MOVE ON TO NEXT HORIZONTAL LOCATION.
               IF(IBOUND(IC,IR,IL).LT.0) GO TO 29
   !
   !7C-----IF CELL IS VARIABLE HEAD, THEN DO BUDGET FOR IT.
               IF (IBOUND(IC,IR,IL).GT.0) THEN
                  Q=RECH(IC,IR)
                  QQ=Q
   !
   !7D-----ADD RECHARGE TO BUFFER, AND STORE LAYER NUMBER IN IRCH.
                  BUFF(IC,IR,IL)=Q
                  IRCH(IC,IR)=IL
   !
   !7E-----IF RECH IS POSITIVE ADD IT TO RATIN, ELSE ADD IT TO RATOUT.
                  IF(Q.GE.ZERO) THEN
                     RATIN=RATIN+QQ
                  ELSE
                     RATOUT=RATOUT-QQ
                  END IF
                  GO TO 29
               END IF
28          CONTINUE
29       CONTINUE
30    CONTINUE
   !
   END IF
   !
   !8------IF CELL-BY-CELL FLOW TERMS SHOULD BE SAVED, CALL APPROPRIATE
   !8------UTILITY MODULE TO WRITE THEM.
100 IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IRCHCB,BUFF,NCOL,NROW,&
   &NLAY,IOUT)
   IF(IBD.EQ.2) CALL UBDSV3(KSTP,KPER,TEXT,IRCHCB,BUFF,IRCH,NRCHOP,&
   &NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
   !
   !9------MOVE TOTAL RECHARGE RATE INTO VBVL FOR PRINTING BY BAS1OT.
   ROUT=RATOUT
   RIN=RATIN
   VBVL(4,MSUM)=ROUT
   VBVL(3,MSUM)=RIN
   !
   !10-----ADD RECHARGE FOR TIME STEP TO RECHARGE ACCUMULATOR IN VBVL.
   VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
   VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
   !
   !11-----MOVE BUDGET TERM LABELS TO VBNM FOR PRINT BY MODULE BAS_OT.
   VBNM(MSUM)=TEXT
   !
   !12-----INCREMENT BUDGET TERM COUNTER.
   MSUM=MSUM+1
   !
   !13-----RETURN
   RETURN
END
SUBROUTINE GWF2RCH7DA(IGRID)
   !  Deallocate RCH DATA
   USE GWFRCHMODULE
   !
   DEALLOCATE(GWFRCHDAT(IGRID)%NRCHOP)
   DEALLOCATE(GWFRCHDAT(IGRID)%IRCHCB)
   DEALLOCATE(GWFRCHDAT(IGRID)%NPRCH)
   DEALLOCATE(GWFRCHDAT(IGRID)%IRCHPF)
   DEALLOCATE(GWFRCHDAT(IGRID)%RECH)
   DEALLOCATE(GWFRCHDAT(IGRID)%IRCH)
   !
   RETURN
END
SUBROUTINE SGWF2RCH7PNT(IGRID)
   !  Set RCH pointers for grid.
   USE GWFRCHMODULE
   !
   NRCHOP=>GWFRCHDAT(IGRID)%NRCHOP
   IRCHCB=>GWFRCHDAT(IGRID)%IRCHCB
   NPRCH=>GWFRCHDAT(IGRID)%NPRCH
   IRCHPF=>GWFRCHDAT(IGRID)%IRCHPF
   RECH=>GWFRCHDAT(IGRID)%RECH
   IRCH=>GWFRCHDAT(IGRID)%IRCH
   !
   RETURN
END
SUBROUTINE SGWF2RCH7PSV(IGRID)
   !  Save RCH pointers for grid.
   USE GWFRCHMODULE
   !
   GWFRCHDAT(IGRID)%NRCHOP=>NRCHOP
   GWFRCHDAT(IGRID)%IRCHCB=>IRCHCB
   GWFRCHDAT(IGRID)%NPRCH=>NPRCH
   GWFRCHDAT(IGRID)%IRCHPF=>IRCHPF
   GWFRCHDAT(IGRID)%RECH=>RECH
   GWFRCHDAT(IGRID)%IRCH=>IRCH
   !
   RETURN
END
