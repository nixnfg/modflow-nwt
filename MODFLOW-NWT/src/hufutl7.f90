   !     VERSION 2.3.1, 05/13/2003
SUBROUTINE UHUF7RMLT(RMLT0,J,I,NZ,NM,ICL)
   !
   !     ******************************************************************
   !     Calculate RMLT for specified cell.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE PARAMMODULE, ONLY:IPCLST,IZON,RMLT
   !
   RMLT0=1.0
   IF (NZ.GT.0) THEN
      RMLT0=0.
      DO 30 JJ = 5,IPCLST(4,ICL)
         IF(IZON(J,I,NZ).EQ.IPCLST(JJ,ICL)) THEN
            IF(NM.GT.0) THEN
               RMLT0=RMLT(J,I,NM)
            ELSE
               RMLT0=1.0
            ENDIF
         END IF
30    CONTINUE
   ELSEIF(NM.GT.0) THEN
      RMLT0=RMLT(J,I,NM)
   ENDIF
   !
   !
   !4------RETURN
   RETURN
END
   !======================================================================
SUBROUTINE UHUF7PARRP(IN,IOUT,NP,PTYP,ITERP,NHUF)
   !
   !     ******************************************************************
   !     Read and store array parameter definition information for HUF package
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE PARAMMODULE
   USE GWFHUFMODULE, ONLY:HGUNAM
   CHARACTER*(*) PTYP
   CHARACTER*200 LINE
   CHARACTER*10 PN,CTMP1,CTMP2
   !     ------------------------------------------------------------------
   !
   ILFLG=1
   !  Read a parameter definition line and decode the parameter name, type,
   !  and value
   READ(IN,'(A)') LINE
   LLOC=1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
   PN=LINE(ISTART:ISTOP)
   CTMP1=PN
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
   PTYP=LINE(ISTART:ISTOP)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,N,PV,IOUT,IN)
   !
   !  Look for the parameter name in the parameter list
   DO 10 NP=1,MXPAR
      CTMP2=PARNAM(NP)
      CALL UPCASE(CTMP2)
      IF(CTMP1.EQ.CTMP2) THEN
   !
   !  If found, determine if it is an illegal duplicate or if it was
   !  predefined.
   !swm              IF(PARTYP(NP).NE.' ' .AND. IDEFPAR.EQ.0) THEN
         IF(PARTYP(NP).NE.' ' .AND. ITERP.EQ.1) THEN
   !  Illegal duplicate
            WRITE(IOUT,*) ' Duplicate parameter name'
            CALL USTOP(' ')
         END IF
   !  Parameter was predefined -- leave its value alone (i.e. ignore PV).
         GO TO 100
      ELSE IF(PARNAM(NP).EQ.' ') THEN
   !  Parameter was not found in the list, so it is a new definition.
   !  Put values in the list.
         PARNAM(NP)=PN
         B(NP)=PV
         IPSUM=IPSUM+1
         GO TO 100
      END IF
10 CONTINUE
   !  Too many parameters
   WRITE(IOUT,11)
11 FORMAT(1X,'The number of parameters has exceeded the maximum')
   CALL USTOP(' ')
   !
   !  Parse the rest of the parameter definition.
100 PARTYP(NP)=PTYP
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NDHUF,R,IOUT,IN)
   IF(IPLOC(1,NP).EQ.0) THEN
      ICLSUM=ICLSUM+1
      IPLOC(1,NP)=ICLSUM
      ICLSUM=ICLSUM+NDHUF-1
      IPLOC(2,NP)=ICLSUM
   END IF
   IACTIVE(NP)=0
   !
   IF(IPLOC(2,NP).GT.MXCLST) THEN
      WRITE(IOUT,117) IPLOC(2,NP),MXCLST
117   FORMAT(1X,I5,&
      &' CLUSTERS WERE SPECIFIED, BUT THERE IS SPACE FOR ONLY',I5)
      WRITE(IOUT,*) NP,NDHUF
      WRITE(IOUT,'(A)') PARNAM(NP)
      WRITE(IOUT,'(2I10)') IPLOC
      CALL USTOP(' ')
   END IF
   WRITE(IOUT,121) PARNAM(NP),PARTYP(NP),NDHUF
121 FORMAT(1X/,1X,'PARAMETER NAME:',A,'   TYPE:',A,' UNITS:',I4)
   WRITE(IOUT,122) PV
122 FORMAT(1X,'The parameter value from the package file is:',1PG13.5)
   IF(B(NP).NE.PV) WRITE(IOUT,123) B(NP)
123 FORMAT(1X,'This parameter value has been replaced by the',&
   &' value from the',/1X,'Parameter Value file:',1PG13.5)
   !
   !  Read clusters
   DO 200 I=IPLOC(1,NP),IPLOC(2,NP)
      READ(IN,'(A)') LINE
      LLOC=1
   !
   !  Store layer number for LVDA
      IF(PARTYP(NP).EQ.'LVDA') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPCLST(1,I),R,-1,IN)
      ELSEIF(PARTYP(NP).EQ.'SYTP')THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPCLST(1,I),R,-1,IN)
         IPCLST(1,I)=1
      ELSE
   !  Find hydrogeologic-unit number
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         PN=LINE(ISTART:ISTOP)
         CTMP1=PN
   !
   !  Look for the unit name in the list of unit names
         DO 220 NU=1,NHUF
            CTMP2=HGUNAM(NU)
            CALL UPCASE(CTMP2)
            IF(CTMP1.EQ.CTMP2) THEN
               IPCLST(1,I)=NU
               WRITE(IOUT,38) CTMP1,NU
38             FORMAT('UNIT ',A10,'CORRESPONDS TO UNIT NO. ',I5)
               GO TO 221
            END IF
220      CONTINUE
221      CONTINUE
      ENDIF
   !
   ! Parse multiplication and zone array information
      CALL URWORD(LINE,LLOC,IM1,IM2,1,N,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,IZ1,IZ2,1,N,R,IOUT,IN)
      DO 30 J=5,14
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPCLST(J,I),R,-1,IN)
         IF(IPCLST(J,I).EQ.0) THEN
            IPCLST(4,I)=J-1
            GO TO 32
         END IF
30    CONTINUE
      IPCLST(4,I)=14
32    IF(ILFLG.NE.0) THEN
         WRITE(IOUT,36) IPCLST(1,I),LINE(IM1:IM2),LINE(IZ1:IZ2)
36       FORMAT(1X,'               UNIT:',I3,'   MULTIPLIER:',A,&
         &'   ZONE:',A)
      ELSE
         WRITE(IOUT,37) LINE(IM1:IM2),LINE(IZ1:IZ2)
37       FORMAT(1X,'               MULTIPLIER:',A,'   ZONE:',A)
      END IF
   !
   !  Find the multiplier array number
      CTMP1=LINE(IM1:IM2)
      IF(CTMP1.EQ.'NONE') THEN
         IPCLST(2,I)=0
      ELSE
         DO 40 J=1,NMLTAR
            CTMP2=MLTNAM(J)
            CALL UPCASE(CTMP2)
            IF(CTMP1.EQ.CTMP2) GO TO 45
40       CONTINUE
         WRITE(IOUT,'(A)') ' Multiplier array has not been defined'
         CALL USTOP(' ')
45       IPCLST(2,I)=J
      END IF
   !
   !  Find the zone array number
      CTMP1=LINE(IZ1:IZ2)
      IF(CTMP1.EQ.'ALL') THEN
         IPCLST(3,I)=0
      ELSE
         IF(IPCLST(4,I).EQ.4) THEN
            WRITE(IOUT,47)
47          FORMAT(&
            &1X,'There were no zone values specified in the cluster',/&
            &1X,'At least one zone must be specified')
            CALL USTOP(' ')
         END IF
         WRITE(IOUT,48) (IPCLST(J,I),J=5,IPCLST(4,I))
48       FORMAT(1X,'               ZONE VALUES:',10I5)
         DO 50 J=1,NZONAR
            CTMP2=ZONNAM(J)
            CALL UPCASE(CTMP2)
            IF(CTMP1.EQ.CTMP2) GO TO 55
50       CONTINUE
         WRITE(IOUT,'(A)') ' Zone array has not been defined'
         CALL USTOP(' ')
55       IPCLST(3,I)=J
      END IF
200 CONTINUE
   !
   RETURN
END
   !======================================================================
SUBROUTINE UHUF7POPL(VDHD,NCOL,NROW,NLAY,I,J)
   !
   !     ******************************************************************
   !     Populate VDHD array.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE PARAMMODULE
   DIMENSION VDHD(NCOL,NROW,NLAY)
   !     ------------------------------------------------------------------

   !
   ! Loop through parameters
   DO 250 NP=1,IPSUM
      IF(PARTYP(NP).EQ.'LVDA') THEN
   ! Loop through layers that apply to this parameter
         DO 300 ND=IPLOC(1,NP),IPLOC(2,NP)
            NL=IPCLST(1,ND)
            NM=IPCLST(2,ND)
            NZ=IPCLST(3,ND)
   !
            CALL UHUF7RMLT(RMLT0,J,I,NZ,NM,ND)
   !
   !---Populate LVDA array
            VDHD(J,I,NL)=VDHD(J,I,NL)+RMLT0*B(NP)
300      CONTINUE
      ENDIF
250 CONTINUE
   !
   !
   !4------RETURN
   RETURN
END
   !======================================================================
SUBROUTINE UHUF7POP(HUFARRAY,PTYPE,I,J,NNU,IOUT)
   !
   !     ******************************************************************
   !     Populate HUF arrays.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE PARAMMODULE
   USE GWFHUFMODULE, ONLY:NHUF,HUFTHK
   CHARACTER*4 PTYPE,PTEMP
   DIMENSION HUFARRAY(NHUF)
   !
   ! Loop through parameters
   DO 250 NP=1,IPSUM
      PTEMP=PARTYP(NP)
      IF(PTEMP.EQ.PTYPE) THEN
   ! Loop through units that apply to this parameter
         IP1=IPLOC(1,NP)
         IP2=IPLOC(2,NP)
         DO 300 ND=IP1,IP2
            NU=IPCLST(1,ND)
            IF(NNU.GT.0.AND.NNU.NE.NU) GOTO 300
            NM=IPCLST(2,ND)
            NZ=IPCLST(3,ND)
   !
   ! First, skip this unit if thickness if zero
            THCKU=HUFTHK(J,I,NU,2)
            CALL UHUF7RMLT(RMLT0,J,I,NZ,NM,ND)
            THCKU=RMLT0*THCKU
            IF(THCKU.LE.0) GOTO 300
   !
   !---Populate HUF array
            IF(PTYPE.EQ.'VANI') THEN
               IF(RMLT0.NE.0..AND.HUFARRAY(NU).NE.0) THEN
                  WRITE(IOUT,100)
100               FORMAT(//,'Additive VANI parameters not allowed! ',&
                  &'STOP EXECUTION(UHUF7POP)')
                  CALL USTOP(' ')
               ENDIF
            ENDIF
            HUFARRAY(NU)=HUFARRAY(NU)+RMLT0*B(NP)
300      CONTINUE
      ENDIF
250 CONTINUE
   !
   !
   !4------RETURN
   RETURN
END

   !======================================================================
SUBROUTINE UHUF7THK(TOP,BOT,TOPU,THKU,THCK,ATOP,ABOT)
   !
   !     ******************************************************************
   !     Determine contributing thicknesses of hydrogeologic units.
   !     Return adjusted top and bottom of unit in ATOP & ABOT
   !     ******************************************************************
   !
   ABOT=0.0
   ATOP=0.0
   TOPL=TOP
   BOTU=TOPU-THKU
   IF(TOPU.LE.BOT.OR.BOTU.GE.TOPL) THEN
      THCK=0
   ELSE
      ATOP=TOPU
      ABOT=BOTU
      IF(TOPU.GT.TOPL) ATOP=TOPL
      IF(BOTU.LT.BOT) ABOT=BOT
      THCK=ATOP-ABOT
   ENDIF
   IF(ABOT.NE.0.0) THEN
      IF(ABS(THCK/ABOT).LT.1E-4) THCK=0.
   ELSEIF(ATOP.NE.0.0) THEN
      IF(ABS(THCK/ATOP).LT.1E-4) THCK=0.
   ELSE
      IF(ABS(THCK).LT.1E-4) THCK=0.
   ENDIF
   !
   !4------RETURN
   RETURN
END
   !======================================================================
SUBROUTINE UHUF7MMTH(&
&ARRAY1,M1,ARRAY2,M2,ARRAY3,CHAR,NCOL,NROW)
   !
   !     ******************************************************************
   !     Perform matrix math on ARRAY1 & ARRAY2, put results in ARRAY3
   !       ARRAY3 = ARRAY1 'CHAR' ARRAY2, where CHAR = +,-,*, or /
   !       M1 & M2 are constants that can replace the arrays
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   CHARACTER*1 CHAR
   REAL M1, M2
   DIMENSION ARRAY1(NCOL,NROW),ARRAY2(NCOL,NROW),ARRAY3(NCOL,NROW)
   !
   AM1 = M1
   AM2 = M2
   IF(CHAR.EQ.'+') THEN
      DO 10 I=1,NROW
         DO 20 J=1,NCOL
            IF(M1.EQ.0.0) AM1 = ARRAY1(J,I)
            IF(M2.EQ.0.0) AM2 = ARRAY2(J,I)
            ARRAY3(J,I) = AM1 + AM2
20       CONTINUE
10    CONTINUE
   ENDIF
   IF(CHAR.EQ.'-') THEN
      DO 12 I=1,NROW
         DO 22 J=1,NCOL
            IF(M1.EQ.0.0) AM1 = ARRAY1(J,I)
            IF(M2.EQ.0.0) AM2 = ARRAY2(J,I)
            ARRAY3(J,I) = AM1 - AM2
22       CONTINUE
12    CONTINUE
   ENDIF
   IF(CHAR.EQ.'*') THEN
      DO 14 I=1,NROW
         DO 24 J=1,NCOL
            IF(M1.EQ.0.0) AM1 = ARRAY1(J,I)
            IF(M2.EQ.0.0) AM2 = ARRAY2(J,I)
            ARRAY3(J,I) = AM1 * AM2
24       CONTINUE
14    CONTINUE
   ENDIF
   IF(CHAR.EQ.'/') THEN
      DO 16 I=1,NROW
         DO 26 J=1,NCOL
            IF(M1.EQ.0.0) AM1 = ARRAY1(J,I)
            IF(M2.EQ.0.0) AM2 = ARRAY2(J,I)
            IF(AM2.EQ.0.0) THEN
               ARRAY3(J,I) = 0.0
            ELSE
               ARRAY3(J,I) = AM1 / AM2
            ENDIF
26       CONTINUE
16    CONTINUE
   ENDIF
   !
   !4------RETURN
   RETURN
END
   !======================================================================
SUBROUTINE UHUFPRWC(A,NCOL,NROW,HGUNAM,IOUT,IPRN,ANAME)
   !
   !     ******************************************************************
   !     CHECK TO SEE IF AN ARRAY IS CONSTANT, AND PRINT IT APPROPRIATELY
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   DIMENSION A(NCOL,NROW)
   CHARACTER*(*) ANAME
   CHARACTER*(*) HGUNAM
   !     ------------------------------------------------------------------
   !
   !  Check to see if entire array is a constant.
   ILAY=1
   TMP=0.0
   DO 300 I=1,NROW
      DO 300 J=1,NCOL
         IF(A(J,I).NE.TMP) THEN
            IF(TMP.EQ.0.0) THEN
               TMP = A(J,I)
            ELSE
               GO TO 400
            ENDIF
         ENDIF
300 CONTINUE
   IF(ILAY.GT.0) THEN
      WRITE(IOUT,302) ANAME,TMP,HGUNAM
302   FORMAT(1X,/1X,A,' =',1P,G14.6,' FOR HYDROGEOLOGIC UNIT ',A)
   ELSE IF(ILAY.EQ.0) THEN
      WRITE(IOUT,303) ANAME,TMP
303   FORMAT(1X,/1X,A,' =',1P,G14.6)
   ENDIF
   RETURN
   !
   !  Print the array.
400 IF(ILAY.GT.0) THEN
      WRITE(IOUT,494) ANAME,HGUNAM
494   FORMAT(1X,//11X,A,' FOR HYDROGEOLOGIC UNIT ',A)
   ELSE IF(ILAY.EQ.0) THEN
      WRITE(IOUT,495) ANAME
495   FORMAT(1X,//11X,A)
   END IF
   IF(IPRN.GE.0) CALL ULAPRW(A,ANAME,0,0,NCOL,NROW,0,IPRN,IOUT)
   !
   RETURN
END

