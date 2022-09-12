MODULE GWFHFBMODULE
   INTEGER, POINTER  ::MXHFB,NHFB,IPRHFB,NHFBNP,NPHFB,IHFBPB
   REAL,    DIMENSION(:,:), POINTER   ::HFB
   TYPE GWFHFBTYPE
      INTEGER, POINTER  ::MXHFB,NHFB,IPRHFB,NHFBNP,NPHFB,IHFBPB
      REAL,    DIMENSION(:,:), POINTER   ::HFB
   END TYPE
   TYPE(GWFHFBTYPE), SAVE    ::GWFHFBDAT(10)
END MODULE GWFHFBMODULE


SUBROUTINE GWF2HFB7AR(INHFB,IGRID)
   !     ******************************************************************
   !     ALLOCATE ARRAY STORAGE FOR HORIZONTAL FLOW BARRIER PACKAGE
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,LAYHDT,CR,CC,BOTM,LBOTM,&
   &DELR,DELC,IOUT
   USE GWFHFBMODULE,ONLY:MXHFB,NHFB,IPRHFB,NHFBNP,NPHFB,IHFBPB,HFB
   !
   INTEGER INHFB, MXACTFB
   CHARACTER*16 AUX(1)
   CHARACTER*200 LINE
   !     ------------------------------------------------------------------
   !
   !1------Allocate scalar data.
   ALLOCATE(MXHFB,NHFB,IPRHFB,NHFBNP,NPHFB,IHFBPB)
   !
   !2------IDENTIFY PACKAGE.
   WRITE(IOUT,1) INHFB
1  FORMAT(1X,/1X,'HFB -- HORIZONTAL-FLOW BARRIER',&
   &' PACKAGE, NWT VERSION 1.3.0, 7/01/2022.',/,&
   &'   INPUT READ FROM UNIT ',I4)
   !
   !3------READ AND PRINT NPHFB, MXFB, NHFBNP
   CALL URDCOM(INHFB,IOUT,LINE)
   LLOC = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPHFB,DUM,IOUT,INHFB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXFBP,DUM,IOUT,INHFB)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NHFBNP,DUM,IOUT,INHFB)
   WRITE(IOUT,500) NPHFB,MXFBP
500 FORMAT(1X,I5,' PARAMETERS DEFINE A MAXIMUM OF ',I6,&
   &' HORIZONTAL FLOW BARRIERS')
   WRITE(IOUT,530) NHFBNP
530 FORMAT(1X,I6,' HORIZONTAL FLOW BARRIERS NOT DEFINED BY',&
   &' PARAMETERS')
   !
   !4------LOOK FOR NOPRINT OPTION.
   IPRHFB = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INHFB)
   IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
      WRITE(IOUT,3)
3     FORMAT(1X,&
      &'LISTS OF HORIZONTAL FLOW BARRIER CELLS WILL NOT BE PRINTED')
      IPRHFB = 0
   END IF
   !
   !5------CALCULATE AMOUNT OF SPACE USED BY HFB PACKAGE AND ALLOCATE HFB.
   MXACTFB = NHFBNP+MXFBP
   IHFBPB = MXACTFB + 1
   MXHFB = MXACTFB + MXFBP
   ALLOCATE (HFB(7,MXHFB))
   !
   !6------CHECK THAT THE FLOW PACKAGE IS A KIND THAT HFB CAN SUPPORT.
   !6------LAYHDT IS -1 UNLESS THE FLOW PACKAGE CHANGES IT.  IF LAYHDT
   !6------IS STILL NEGATIVE, IT IS ASSUMED THAT HFB WILL NOT WORK.
   IF (LAYHDT(1).LT.0) THEN
      WRITE(IOUT,550)
550   FORMAT(/,&
      &' ERROR: SELECTED FLOW PACKAGE DOES NOT SUPPORT HFB PACKAGE',/,&
      &' -- STOP EXECUTION (GWF2HFB7AR)')
      CALL USTOP(' ')
   ENDIF
   !
   !7------READ PARAMETER DEFINITIONS (ITEMS 2 AND 3)
   WRITE(IOUT,600) NPHFB
600 FORMAT(//,1X,I5,' HFB parameters')
   IF (NPHFB.GT.0) THEN
      LSTSUM = IHFBPB
      DO 20 K = 1,NPHFB
         LSTBEG = LSTSUM
         CALL UPARLSTRP(LSTSUM,MXHFB,INHFB,IOUT,IP,'HFB ','HFB ',&
         &1,NUMINST)
         IF(NUMINST.GT.0) THEN
            WRITE(IOUT,*) ' INSTANCES ARE NOT SUPPORTED FOR HFB'
            CALL USTOP(' ')
         END IF
         NLST=LSTSUM-LSTBEG
         CALL SGWF2HFB7RL(NLST,HFB,LSTBEG,MXHFB,INHFB,IOUT,&
         &'BARRIER  LAYER  IROW1  ICOL1  IROW2  ICOL2     FACTOR',&
         &NCOL,NROW,NLAY,IPRHFB)
         CALL SGWF2HFB7CK(LSTBEG,LSTSUM-1)
20    CONTINUE
   ENDIF
   !
   !8------READ BARRIERS NOT DEFINED BY PARAMETERS (ITEM 4)
   NHFB = 0
   WRITE(IOUT,610) NHFBNP
610 FORMAT(/,1X,I6,' BARRIERS NOT DEFINED BY PARAMETERS')
   IF(NHFBNP.GT.0) THEN
      CALL SGWF2HFB7RL(NHFBNP,HFB,1,MXHFB,INHFB,IOUT,&
      &'BARRIER  LAYER  IROW1  ICOL1  IROW2  ICOL2    HYDCHR',&
      &NCOL,NROW,NLAY,IPRHFB)
      NHFB = NHFB + NHFBNP
      CALL SGWF2HFB7CK(1,NHFBNP)
   ENDIF
   !
   !9------SUBSTITUTE DATA FOR PARAMETERIZED BARRIERS INTO ACTIVE SECTION
   !9------OF HFB ARRAY
   IOUTU = IOUT
   IF (IPRHFB.EQ.0) IOUTU = -IOUT
   MXACTFB= IHFBPB-1
   CALL PRESET('HFB ')
   IF(NPHFB.GT.0) THEN
   !
   !10-----READ NUMBER OF ACTIVE HFB PARAMETERS (ITEM 5)
      READ(INHFB,*) NACTHFB
      IF (NACTHFB.GT.0) THEN
         DO 650 I = 1,NACTHFB
   !
   !11-----READ AND ACTIVATE AN HFB PARAMETER (ITEM 6)
            CALL SGWF2HFB7SUB(INHFB,'HFB ',IOUTU,'HFB ',HFB,7,MXHFB,&
            &MXACTFB,NHFB,&
            &'BARRIER  LAYER  IROW1  ICOL1  IROW2  ICOL2     HYDCHR')
650      CONTINUE
      ENDIF
   ENDIF
   !
   !12-----MODIFY HORIZONTAL BRANCH CONDUCTANCES FOR CONSTANT T LAYERS.
   CALL SGWF2HFB7MC()
   WRITE (IOUT,660) NHFB
660 FORMAT(/,1X,1I6,' HFB BARRIERS')
   !
   !13-----SAVE POINTERS TO GRID AND RETURN.
   CALL SGWF2HFB7PSV(IGRID)
   RETURN
END
SUBROUTINE GWF2HFB7FM(IGRID)
   !     ******************************************************************
   !     MODIFY HORIZONTAL BRANCH CONDUCTANCES IN VARIABLE-TRANSMISSIVITY
   !     LAYERS TO ACCOUNT FOR HORIZONTAL FLOW BARRIERS. STORE UNMODIFIED
   !     HORIZONTAL CONDUCTANCE IN HFB(7,#) TO ALLOW CALCULATION OF
   !     SENSITIVITIES.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,HNEW,LAYHDT,CR,CC,BOTM,LBOTM,&
   &DELR,DELC
   USE GWFHFBMODULE,ONLY:NHFB,HFB
   !     ------------------------------------------------------------------
   !
   !1------Set pointers to the specified grid.
   CALL SGWF2HFB7PNT(IGRID)
   !
   !2------FOR EACH BARRIER, MODIFY HORIZONTAL BRANCH CONDUCTANCES IF LAYER
   !2------IS CONVERTIBLE.
   DO 10 II=1,NHFB
      K = HFB(1,II)
   !
   !3-----IF LAYHDT=0, THICKNESS AND CONDUCTANCE DO NOT VARY, AND
   !3-----MODIFICATION OF CONDUCTANCE DUE TO BARRIER WAS DONE IN
   !3-----SGWF1HFBMC
      IF (LAYHDT(K).GT.0) THEN
   !
   !4------CELL (J1,I1,K) IS THE ONE WHOSE HORIZONTAL BRANCH
   !4------CONDUCTANCES ARE TO BE MODIFIED.
         I1 = HFB(2,II)
         J1 = HFB(3,II)
   !
   !5------CELL (J2,I2,K) IS THE CELL NEXT TO CELL (J1,I1,K) AND
   !5------SEPARATED FROM IT BY THE BARRIER.
         I2 = HFB(4,II)
         J2 = HFB(5,II)
         HCDW = HFB(6,II)
   !
   !6------IF I1=I2, MODIFY HORIZONTAL BRANCH CONDUCTANCES ALONG ROW
   !6------DIRECTION.
         IF (I1.EQ.I2) THEN
   !
   !7------IF CR(J1,I1,K) NOT = 0, CELLS ON EITHER SIDE OF BARRIER ARE
   !7------ACTIVE
            IF (CR(J1,I1,K).NE.0.) THEN
   !
   !8------CALCULATE AVERAGE SATURATED THICKNESS BETWEEN CELLS
   !8------(J1,I1,K) AND (J2,I2,K).  NOTE: NEGATIVE SATURATED
   !8------THICKNESS DOES NOT OCCUR; OTHERWISE, CR(J1,I1,K) WOULD BE
   !8------ZERO AND THE FOLLOWING CALCULATION FOR SATURATED THICKNESS
   !8------WOULD BE SKIPPED.
               HD1 = HNEW(J1,I1,K)
               HD2 = HNEW(J2,I2,K)
               IF (HD1.GT.BOTM(J1,I1,LBOTM(K)-1)) HD1 =&
               &BOTM(J1,I1,LBOTM(K)-1)
               IF (HD2.GT.BOTM(J2,I2,LBOTM(K)-1)) HD2 =&
               &BOTM(J2,I2,LBOTM(K)-1)
               THKAVG = ((HD1-BOTM(J1,I1,LBOTM(K))) +&
               &(HD2-BOTM(J2,I2,LBOTM(K))))/2.
   !
   !9------STORE UNMODIFIED CR FOR CALCULATING SENSITIVITIES
               HFB(7,II) = CR(J1,I1,K)
   !
   !10-----MODIFY CR(J1,I1,K) TO ACCOUNT FOR BARRIER.
               TDW = THKAVG*HCDW
               CR(J1,I1,K) = TDW*CR(J1,I1,K)*DELC(I1)/&
               &(TDW*DELC(I1)+CR(J1,I1,K))
            ENDIF
   !
   !11-----CASE OF J1=J2. MODIFY HORIZONTAL BRANCH CONDUCTANCES ALONG
   !11-----COLUMN DIRECTION.
         ELSE
   !
   !12-----IF CC(J1,I1,K) NOT = 0, CELLS ON EITHER SIDE OF BARRIER ARE
   !12-----ACTIVE
            IF (CC(J1,I1,K).NE.0.) THEN
   !
   !13-----CALCULATE AVERAGE SATURATED THICKNESS BETWEEN CELLS
   !13-----(J1,I1,K) AND (J2,I2,K).  NEGATIVE SATURATED THICKNESS
   !13-----DOES NOT OCCUR FOR THE SAME REASON AS DESCRIBED ABOVE.
               HD1 = HNEW(J1,I1,K)
               HD2 = HNEW(J2,I2,K)
               IF (HD1.GT.BOTM(J1,I1,LBOTM(K)-1)) HD1 =&
               &BOTM(J1,I1,LBOTM(K)-1)
               IF (HD2.GT.BOTM(J2,I2,LBOTM(K)-1)) HD2 =&
               &BOTM(J2,I2,LBOTM(K)-1)
               THKAVG = ((HD1-BOTM(J1,I1,LBOTM(K))) +&
               &(HD2-BOTM(J2,I2,LBOTM(K))))/2.
   !
   !14-----STORE UNMODIFIED CC FOR CALCULATING SENSITIVITIES.
               HFB(7,II) = CC(J1,I1,K)
   !
   !15-----MODIFY CC(J1,I1,K) TO ACCOUNT FOR BARRIER.
               TDW = THKAVG*HCDW
               CC(J1,I1,K) = TDW*CC(J1,I1,K)*DELR(J1)/&
               &(TDW*DELR(J1)+CC(J1,I1,K))
            ENDIF
         ENDIF
      ENDIF
10 CONTINUE
   !
   !16-----RETURN
   RETURN
END
SUBROUTINE SGWF2HFB7MC()
   !     ******************************************************************
   !     MODIFY HORIZONTAL CONDUCTANCES (CR AND CC) FOR CONFINED LAYERS TO
   !     ACCOUNT FOR HORIZONTAL FLOW BARRIERS.  STORE UNMODIFIED HORIZONTAL
   !     CONDUCTANCES IN HFB(7,#) TO ALLOW CALCULATION OF SENSITIVITIES.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:IOUT,BOTM,LBOTM,DELR,DELC,CR,CC,LAYHDT
   USE GWFHFBMODULE,ONLY:NHFB,HFB
   !     ------------------------------------------------------------------
   !
   !1------INITIALIZE ERROR FLAG TO ZERO.
   IERFLG=0
   !
   !2----DO FOR EACH BARRIER IN RANGE.
   DO 10 II = 1,NHFB
      K = HFB(1,II)
   !
   !3------FIND ROW AND COLUMN NUMBERS OF THE TWO CELLS ON BOTH SIDES
   !3------OF THE BARRIER.
      I1 = HFB(2,II)
      J1 = HFB(3,II)
      I2 = HFB(4,II)
      J2 = HFB(5,II)
      TH0 = BOTM(J1,I1,LBOTM(K)-1) - BOTM(J1,I1,LBOTM(K))
      TH1 = BOTM(J2,I2,LBOTM(K)-1) - BOTM(J2,I2,LBOTM(K))
      THKAVG = (TH0+TH1)/2.0
      TDW = THKAVG*HFB(6,II)
   !
   !4------IF I1=I2, BARRIER IS BETWEEN TWO CELLS ON THE SAME ROW.
      IF (I1.EQ.I2) THEN
   !
   !5------IF J2-J1=1, THE TWO CELLS ARE NEXT TO ONE ANOTHER (DATA OK).
         IF ((J2-J1).EQ.1) THEN
   !
   !6------BARRIER CELLS ARE ADJACENT.
   !6------IF LAYER IS CONFINED AND BOTH CELLS ARE ACTIVE, SAVE
   !-------ORIGINAL CR FOR COMPUTING SENSITIVITIES AND MODIFY CR
            IF (LAYHDT(K).EQ.0) THEN
   !
   !7------IF CR(J1,I1,K) NOT 0, BOTH CELLS ARE ACTIVE.
               IF (CR(J1,I1,K).NE.0.) THEN
                  HFB(7,II) = CR(J1,I1,K)
   !
   !8------MODIFY CR(J1,I1,K) TO ACCOUNT FOR BARRIER.
                  CR(J1,I1,K) = TDW*CR(J1,I1,K)*DELC(I1)/&
                  &(TDW*DELC(I1)+CR(J1,I1,K))
               ENDIF
            ENDIF
         ENDIF
   !
   !9------IF J1=J2, BARRIER IS BETWEEN TWO CELLS ON THE SAME COLUMN.
      ELSEIF (J1.EQ.J2) THEN
   !
   !10-----IF I2-I1=1, THE TWO CELLS ARE NEXT TO ONE ANOTHER (DATA OK).
         IF ((I2-I1).EQ.1) THEN
   !
   !11-----BARRIER CELLS ARE ADJACENT.
   !11-----IF LAYER IS CONFINED AND BOTH CELLS ARE ACTIVE, SAVE
   !11-----ORIGINAL CC FOR COMPUTING SENSITIVITIES AND MODIFY CC
            IF (LAYHDT(K).EQ.0) THEN
   !
   !12-----IF CC(J1,I1,K) NOT 0, BOTH CELLS ARE ACTIVE.
               IF (CC(J1,I1,K).NE.0.) THEN
                  HFB(7,II) = CC(J1,I1,K)
   !
   !13-----MODIFY CC(J1,I1,K) TO ACCOUNT FOR BARRIER.
                  CC(J1,I1,K) = TDW*CC(J1,I1,K)*DELR(J1)/&
                  &(TDW*DELR(J1)+CC(J1,I1,K))
               ENDIF
            ENDIF
         ENDIF
      ENDIF
10 CONTINUE
   !
   !14-----RETURN
   RETURN
END
SUBROUTINE SGWF2HFB7CK(IB1,IB2)
   !     ******************************************************************
   !     CHECK HFB CELL LOCATIONS
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:IOUT
   USE GWFHFBMODULE,ONLY:HFB
   !     ------------------------------------------------------------------
   !
   !1----INITIALIZE ERROR FLAG TO ZERO.
   IERFLG=0
   !
   !2----CHECK EACH BARRIER IN RANGE.
   DO 10 II = IB1,IB2
   !
   !3------FIND ROW AND COLUMN NUMBERS OF THE TWO CELLS ON BOTH SIDES
   !3------OF THE BARRIER AND REARRANGE HFB ARRAY.
      I1 = MIN(HFB(2,II),HFB(4,II))
      J1 = MIN(HFB(3,II),HFB(5,II))
      I2 = MAX(HFB(2,II),HFB(4,II))
      J2 = MAX(HFB(3,II),HFB(5,II))
      HFB(2,II) = I1
      HFB(3,II) = J1
      HFB(4,II) = I2
      HFB(5,II) = J2
      ID = I2 - I1
      JD = J2 - J1
      IF (ID.LT.0 .OR. ID.GT.1 .OR. JD.LT.0 .OR. JD.GT.1 .OR.&
      &ID.EQ.JD) THEN
   !
   !4------CELLS ARE NOT ADJACENT. PRINT ERROR MESSAGE AND SET ERROR FLAG.
80       WRITE (IOUT,1) II-IB1+1
1        FORMAT (1X,'ERROR DETECTED IN LOCATION DATA OF BARRIER NO. ',&
         &I6)
         IERFLG=1
      ENDIF
10 CONTINUE
   !
   !5------HALT EXECUTION IF ERRORS ARE DETECTED.
   IF (IERFLG.EQ.1) CALL USTOP(' ')
   !
   !6------RETURN
   RETURN
END
SUBROUTINE SGWF2HFB7RL(NLIST,HFB,LSTBEG,MXHFB,INPACK,&
&IOUT,LABEL,NCOL,NROW,NLAY,IPRFLG)
   !     ******************************************************************
   !     Read and print a list of HFB barriers.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   use openspec
   CHARACTER*(*) LABEL
   DIMENSION HFB(7,MXHFB)
   CHARACTER*200 LINE,FNAME
   CHARACTER*1 DASH(120)
   DATA DASH/120*'-'/
   DATA NUNOPN/99/
   !     ------------------------------------------------------------------
   !
   !1------Check for and decode EXTERNAL and SFAC records.
   IN = INPACK
   ICLOSE = 0
   READ(IN,'(A)') LINE
   SFAC = 1.
   LLOC = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
   IF(LINE(ISTART:ISTOP).EQ.'EXTERNAL') THEN
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
      IN = I
      IF (IPRFLG.EQ.1) WRITE(IOUT,111) IN
111   FORMAT(1X,'Reading list on unit ',I4)
      READ(IN,'(A)') LINE
   ELSE IF(LINE(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
      FNAME = LINE(ISTART:ISTOP)
      IN = NUNOPN
      IF (IPRFLG.EQ.1) WRITE(IOUT,115) IN,FNAME
115   FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
      OPEN(UNIT=IN,FILE=FNAME,ACTION=ACTION(1))
      ICLOSE = 1
      READ(IN,'(A)') LINE
   END IF
   LLOC = 1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
   IF(LINE(ISTART:ISTOP).EQ.'SFAC') THEN
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SFAC,IOUT,IN)
      IF (IPRFLG.EQ.1) THEN
         WRITE(IOUT,116) SFAC
116      FORMAT(1X,'LIST SCALING FACTOR= ',1PG12.5)
      ENDIF
      READ(IN,'(A)') LINE
   END IF
   !
   !2------Define label for printout.
   NBUF = LEN(LABEL)+3
   IF (IPRFLG.EQ.1) THEN
      WRITE(IOUT,103) LABEL
      WRITE(IOUT,104) (DASH(J),J=1,NBUF)
103   FORMAT(1X,/1X,A)
104   FORMAT(1X,400A)
   ENDIF
   !
   !3------Loop through the number of cells to read.
   N = NLIST+LSTBEG-1
   DO 250 II=LSTBEG,N
   !
   !4------Read a line into the buffer.  (The first line has already been read
   !4------in order to scan for EXTERNAL and SFAC records.)
      IF(II.NE.LSTBEG) READ(IN,'(A)') LINE
   !
   !5------Read the non-optional values from a line.
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,K,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I1,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J1,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I2,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J2,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,FACTOR,IOUT,IN)
      HFB(1,II) = K
      HFB(2,II) = I1
      HFB(3,II) = J1
      HFB(4,II) = I2
      HFB(5,II) = J2
      HFB(6,II) = FACTOR*SFAC
      HFB(7,II) = 0.0
   !
   !6------Write the values that were read.
      NN = II-LSTBEG+1
      IF (IPRFLG.EQ.1) WRITE(IOUT,205) NN,K,I1,J1,I2,J2,HFB(6,II)
205   FORMAT(1X,I6,2X,I5,1X,4(2X,I5),2X,1PG11.4)
   !
   !7------Check for illegal grid location.
      IF(K.LT.1 .OR. K.GT.NLAY) THEN
         WRITE(IOUT,*) ' Layer number in list is outside of the grid'
         CALL USTOP(' ')
      END IF
      IF(I1.LT.1 .OR. I1.GT.NROW .OR. I2.LT.1 .OR. I2.GT.NROW) THEN
         WRITE(IOUT,*) ' Row number in list is outside of the grid'
         CALL USTOP(' ')
      END IF
      IF(J1.LT.1 .OR. J1.GT.NCOL .OR. J2.LT.1 .OR. J2.GT.NCOL) THEN
         WRITE(IOUT,*) ' Column number in list is outside of the grid'
         CALL USTOP(' ')
      END IF
250 CONTINUE
   IF(ICLOSE.NE.0) CLOSE(UNIT=IN)
   !
   !8------Return.
   RETURN
END
SUBROUTINE SGWF2HFB7SUB(IN,PACK,IOUTU,PTYP,HFB,LSTVL,MXHFB,&
&MXACTFB,NHFB,LABEL)
   !     ******************************************************************
   !     Read a parameter name, look it up in the list of parameters,
   !     and substitute values into active part of HFB array.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE PARAMMODULE
   CHARACTER*(*) PACK,PTYP
   DIMENSION HFB(LSTVL,MXHFB)
   CHARACTER*(*) LABEL
   CHARACTER*200 LINE
   CHARACTER*10 CTMP1,CTMP2,CTMP3,CTMP4
   !     ------------------------------------------------------------------
   !
   !1------The Listing File file unit is the absolute value of IOUTU.
   !1------Read the parameter name.
   IOUT = ABS(IOUTU)
   READ(IN,'(A)') LINE
   LLOC=1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
   WRITE(IOUT,1) LINE(ISTART:ISTOP)
1  FORMAT(/,' Parameter:  ',A)
   IF(LINE(ISTART:ISTOP).EQ.' ') THEN
      WRITE(IOUT,*) ' Blank parameter name in the ',PACK,' file.'
      CALL USTOP(' ')
   END IF
   !
   !2------Find the parameter in the list of parameters.
   CTMP1=LINE(ISTART:ISTOP)
   CALL UPCASE(CTMP1)
   DO 100 IP=1,IPSUM
      CTMP2=PARNAM(IP)
      CALL UPCASE(CTMP2)
      IF(CTMP1.EQ.CTMP2) THEN
         IF(PARTYP(IP).NE.PTYP) THEN
            WRITE(IOUT,11) PARNAM(IP),PARTYP(IP),PACK,PTYP
11          FORMAT(1X,'Parameter type conflict:',/&
            &1X,'Named parameter:',A,' was defined as type:',A,/&
            &1X,'However, this parameter is used in the ',A,&
            &' file, so it should be type:',A)
            CALL USTOP(' ')
         END IF
   !
   !3------Set indices to point to the barriers that correspond to the
   !3------specified parameter.
         NLST=IPLOC(2,IP)-IPLOC(1,IP)+1
         NI=1
   !
   !4------Check that the parameter is not already active.
         IF (IACTIVE(IP).GT.0) THEN
            WRITE(IOUT,73) PARNAM(IP)
73          FORMAT(/,1X,'*** ERROR: PARAMETER "',A,&
            &'" HAS ALREADY BEEN ACTIVATED THIS STRESS PERIOD',/,&
            &' -- STOP EXECUTION (UPARLSTSUB)')
            CALL USTOP(' ')
         ENDIF
   !
   !5------Set the active flag.
         IACTIVE(IP)=NI
   !
   !6------Accumulate the total number of active barriers in the list.
         NHFB=NHFB+NLST
         IF(NHFB.GT.MXACTFB) THEN
            WRITE(IOUT,83) NHFB,MXACTFB
83          FORMAT(1X,/1X,'THE NUMBER OF ACTIVE LIST ENTRIES (',I6,&
            &')',/1X,'IS GREATER THAN THE MAXIMUM ALLOWED (',I6,')')
            CALL USTOP(' ')
         END IF
   !
   !7------Write label for barrier values if IOUTU is positive.
         IF(IOUTU.GT.0) THEN
            WRITE(IOUT,'(1X,A)') LABEL
            WRITE(IOUT,84)
84          FORMAT(1X,56('-'))
         END IF
   !
   !8------Copy the values from the paramter location into the front part
   !8------of the list where the currently active list is kept.
         DO 90 I=1,NLST
            II=NHFB-NLST+I
            III=I-1+IPLOC(1,IP)+(NI-1)*NLST
            DO 85 J=1,7
               HFB(J,II)=HFB(J,III)
85          CONTINUE
   !
   !8A-----Scale HYDCHR by the parameter value.
            HFB(6,II)=HFB(6,II)*B(IP)
            IL=HFB(1,II)
            IR1=HFB(2,II)
            IC1=HFB(3,II)
            IR2=HFB(4,II)
            IC2=HFB(5,II)
            IF(IOUTU.GT.0) WRITE(IOUT,89) II,IL,IR1,IC1,IR2,IC2,&
            &HFB(6,II)
89          FORMAT(1X,I6,2X,I5,1X,4(2X,I5),2X,1PG11.4)
90       CONTINUE
   !
   !8B------After moving the data, return.
         RETURN
      END IF
100 CONTINUE
   !
   !9------All parameter names have been checked without finding the
   !9------parameter. Write an error message and stop.
   WRITE(IOUT,*) ' The ',PACK,&
   &' file specifies an undefined parameter:',LINE(ISTART:ISTOP)
   CALL USTOP(' ')
   !
END
   !
SUBROUTINE GWF2HFB7UPW(IGRID)
   !     ******************************************************************
   !     MODIFY HORIZONTAL BRANCH CONDUCTANCES IN VARIABLE-TRANSMISSIVITY
   !     LAYERS TO ACCOUNT FOR HORIZONTAL FLOW BARRIERS IN MODFLOW-NWT.
   !     STORE UNMODIFIED HORIZONTAL CONDUCTANCE IN HFB(7,#) TO ALLOW
   !     CALCULATION OF SENSITIVITIES.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,HNEW,LAYHDT,CR,CC,BOTM,LBOTM,&
   &DELR,DELC
   USE GWFHFBMODULE,ONLY:NHFB,HFB
   !     ------------------------------------------------------------------
   !
   !1------Set pointers to the specified grid.
   CALL SGWF2HFB7PNT(IGRID)
   !
   !2------FOR EACH BARRIER, MODIFY HORIZONTAL BRANCH CONDUCTANCES IF LAYER
   !2------IS CONVERTIBLE.
   DO 10 II=1,NHFB
      K = HFB(1,II)
   !
   !3-----IF LAYHDT=0, THICKNESS AND CONDUCTANCE DO NOT VARY, AND
   !3-----MODIFICATION OF CONDUCTANCE DUE TO BARRIER WAS DONE IN
   !3-----SGWF1HFBMC
      IF (LAYHDT(K).GT.0) THEN
   !
   !4------CELL (J1,I1,K) IS THE ONE WHOSE HORIZONTAL BRANCH
   !4------CONDUCTANCES ARE TO BE MODIFIED.
         I1 = HFB(2,II)
         J1 = HFB(3,II)
   !
   !5------CELL (J2,I2,K) IS THE CELL NEXT TO CELL (J1,I1,K) AND
   !5------SEPARATED FROM IT BY THE BARRIER.
         I2 = HFB(4,II)
         J2 = HFB(5,II)
         HCDW = HFB(6,II)
   !
   !6------IF I1=I2, MODIFY HORIZONTAL BRANCH CONDUCTANCES ALONG ROW
   !6------DIRECTION.
         IF (I1.EQ.I2) THEN
   !
   !7------IF CR(J1,I1,K) NOT = 0, CELLS ON EITHER SIDE OF BARRIER ARE
   !7------ACTIVE
            IF (CR(J1,I1,K).NE.0.) THEN
   !
   !
   !9------STORE UNMODIFIED CR FOR CALCULATING SENSITIVITIES
               HFB(7,II) = CR(J1,I1,K)
   !
   !10-----MODIFY CR(J1,I1,K) TO ACCOUNT FOR BARRIER.
               TDW = HCDW
               CR(J1,I1,K) = TDW*CR(J1,I1,K)*DELC(I1)/&
               &(TDW*DELC(I1)+CR(J1,I1,K))
            ENDIF
   !
   !11-----CASE OF J1=J2. MODIFY HORIZONTAL BRANCH CONDUCTANCES ALONG
   !11-----COLUMN DIRECTION.
         ELSE
   !
   !12-----IF CC(J1,I1,K) NOT = 0, CELLS ON EITHER SIDE OF BARRIER ARE
   !12-----ACTIVE
            IF (CC(J1,I1,K).NE.0.) THEN
   !
   !14-----STORE UNMODIFIED CC FOR CALCULATING SENSITIVITIES.
               HFB(7,II) = CC(J1,I1,K)
   !
   !15-----MODIFY CC(J1,I1,K) TO ACCOUNT FOR BARRIER.
               TDW = HCDW
               CC(J1,I1,K) = TDW*CC(J1,I1,K)*DELR(J1)/&
               &(TDW*DELR(J1)+CC(J1,I1,K))
            ENDIF
         ENDIF
      ENDIF
10 CONTINUE
   !
   !16-----RETURN
   RETURN
END
SUBROUTINE GWF2HFB7DA(IGRID)
   !  Deallocate HFB data for a grid.
   USE GWFHFBMODULE
   !
   DEALLOCATE(GWFHFBDAT(IGRID)%MXHFB)
   DEALLOCATE(GWFHFBDAT(IGRID)%NHFB)
   DEALLOCATE(GWFHFBDAT(IGRID)%IPRHFB)
   DEALLOCATE(GWFHFBDAT(IGRID)%NHFBNP)
   DEALLOCATE(GWFHFBDAT(IGRID)%NPHFB)
   DEALLOCATE(GWFHFBDAT(IGRID)%IHFBPB)
   DEALLOCATE(GWFHFBDAT(IGRID)%HFB)
   !
   RETURN
END
SUBROUTINE SGWF2HFB7PNT(IGRID)
   !  Set pointers to HFB data for a grid.
   USE GWFHFBMODULE
   !
   MXHFB=>GWFHFBDAT(IGRID)%MXHFB
   NHFB=>GWFHFBDAT(IGRID)%NHFB
   IPRHFB=>GWFHFBDAT(IGRID)%IPRHFB
   NHFBNP=>GWFHFBDAT(IGRID)%NHFBNP
   NPHFB=>GWFHFBDAT(IGRID)%NPHFB
   IHFBPB=>GWFHFBDAT(IGRID)%IHFBPB
   HFB=>GWFHFBDAT(IGRID)%HFB
   !
   RETURN
END
SUBROUTINE SGWF2HFB7PSV(IGRID)
   !  Save pointers to HFB data for a grid.
   USE GWFHFBMODULE
   !
   GWFHFBDAT(IGRID)%MXHFB=>MXHFB
   GWFHFBDAT(IGRID)%NHFB=>NHFB
   GWFHFBDAT(IGRID)%IPRHFB=>IPRHFB
   GWFHFBDAT(IGRID)%NHFBNP=>NHFBNP
   GWFHFBDAT(IGRID)%NPHFB=>NPHFB
   GWFHFBDAT(IGRID)%IHFBPB=>IHFBPB
   GWFHFBDAT(IGRID)%HFB=>HFB
   !
   RETURN
END
