MODULE GWFUPWMODULE
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER :: HEPS = 1.0E-7
   DOUBLE PRECISION, PARAMETER :: CLOSEZERO = 1.0E-15
   DOUBLE PRECISION,PARAMETER :: BIG = 1.0D20
   DOUBLE PRECISION,PARAMETER :: SMALL = 1.0D-5
   DOUBLE PRECISION, SAVE, DIMENSION(:), POINTER :: Sn, So
   INTEGER, SAVE,   POINTER :: Iuupw
   ! Cell property data
   INTEGER, SAVE,   POINTER ::IUPWCB,IWDFLG,IWETIT,IHDWET,IPHDRY
   INTEGER, SAVE,   POINTER ::ISFAC,ICONCV,ITHFLG,NOCVCO,NOVFC
   REAL,    SAVE,   POINTER ::WETFCT
   INTEGER, SAVE,   POINTER, DIMENSION(:)     ::LAYTYPUPW
   INTEGER, SAVE,   POINTER, DIMENSION(:)     ::LAYAVG
   REAL,    SAVE,   POINTER, DIMENSION(:)     ::CHANI
   INTEGER, SAVE,   POINTER, DIMENSION(:)     ::LAYVKAUPW
   INTEGER, SAVE,   POINTER, DIMENSION(:)     ::LAYWET
   INTEGER, SAVE,   POINTER, DIMENSION(:)     ::LAYSTRT
   INTEGER, SAVE,   POINTER, DIMENSION(:,:)   ::LAYFLG
   INTEGER, SAVE,    DIMENSION(:,:,:), POINTER ::IBOUND2
   REAL,    SAVE,   POINTER, DIMENSION(:,:,:) ::VKAUPW
   REAL,    SAVE,   POINTER, DIMENSION(:,:,:) ::VKCB
   REAL,    SAVE,   POINTER, DIMENSION(:,:,:) ::SC1
   REAL,    SAVE,   POINTER, DIMENSION(:,:,:) ::SC2UPW
   REAL,    SAVE,   POINTER, DIMENSION(:,:,:) ::HANI
   REAL,    SAVE,   POINTER, DIMENSION(:,:,:) ::WETDRY
   REAL,    SAVE,   POINTER, DIMENSION(:,:,:) ::HKUPW
   TYPE GWFUPWTYPE
      INTEGER, POINTER :: Iuupw
   ! Cell property data
      INTEGER, POINTER ::IUPWCB,IWDFLG,IWETIT,IHDWET,IPHDRY
      INTEGER, POINTER ::ISFAC,ICONCV,ITHFLG,NOCVCO,NOVFC
      REAL, POINTER    ::WETFCT
      DOUBLE PRECISION, DIMENSION(:), POINTER :: Sn, So
      INTEGER,   POINTER, DIMENSION(:)     ::LAYTYPUPW
      INTEGER,   POINTER, DIMENSION(:)     ::LAYAVG
      REAL,      POINTER, DIMENSION(:)     ::CHANI
      INTEGER,   POINTER, DIMENSION(:)     ::LAYVKAUPW
      INTEGER,   POINTER, DIMENSION(:)     ::LAYWET
      INTEGER,   POINTER, DIMENSION(:)     ::LAYSTRT
      INTEGER,   POINTER, DIMENSION(:,:)   ::LAYFLG
      INTEGER,   POINTER, DIMENSION(:,:,:) ::IBOUND2
      REAL,      POINTER, DIMENSION(:,:,:) ::VKAUPW
      REAL,      POINTER, DIMENSION(:,:,:) ::VKCB
      REAL,      POINTER, DIMENSION(:,:,:) ::SC1
      REAL,      POINTER, DIMENSION(:,:,:) ::SC2UPW
      REAL,      POINTER, DIMENSION(:,:,:) ::HANI
      REAL,      POINTER, DIMENSION(:,:,:) ::WETDRY
      REAL,      POINTER, DIMENSION(:,:,:) ::HKUPW
   END TYPE GWFUPWTYPE
   TYPE (GWFUPWTYPE) , SAVE::Gwfupwdat(10)
END MODULE GWFUPWMODULE
   !

   !
   !-------SUBROUTINE GWF2UPW1AR
   !
SUBROUTINE GWF2UPW1AR(In, Igrid)

   USE GLOBAL,     ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD,&
   &NCNFBD,IBOUND,BUFF,BOTM,NBOTM,DELR,DELC,IOUT,&
   &LBOTM,HNEW,IUNIT
   USE GWFBASMODULE,ONLY:HDRY
   USE GWFNWTMODULE, ONLY: Numcell
   USE GWFUPWMODULE
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   INTRINSIC INT
   EXTERNAL URDCOM, URWORD
   EXTERNAL SGWF2UPW1PSV
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     ------------------------------------------------------------------
   INTEGER In, Igrid
   !     ------------------------------------------------------------------
   !     LOCAL VARIABLES
   !     ------------------------------------------------------------------
   INTEGER lloc, istart, istop, i, ic, ir, il, jj
   CHARACTER(LEN=200) line
   INTEGER NPHK,NPVKCB,NPVK,NPVANI,NPSS,NPSY,NPHANI
   INTEGER IANAME,KHANI,N,KK,j,k,NCNVRT,NHANI,NWETD
   !     LOCAL VARIABLES FOR DEFINING CELL PROPERTES (FROM LPF)
   INTEGER NPUPW, NOPCHK
   REAL ZERO, R
   !
   CHARACTER*14 LAYPRN(5),AVGNAM(3),TYPNAM(2),VKANAM(2),WETNAM(2),&
   &HANNAM
   CHARACTER*100 WARN1,WARN2,WARN3,WARN4,WARN5,WARN6
   DATA AVGNAM/'      HARMONIC','   LOGARITHMIC','     LOG-ARITH'/
   DATA TYPNAM/'      CONFINED','   CONVERTIBLE'/
   DATA VKANAM/'    VERTICAL K','    ANISOTROPY'/
   DATA WETNAM/'  NON-WETTABLE','      WETTABLE'/
   DATA HANNAM/'      VARIABLE'/
   CHARACTER*24 ANAME(9),STOTXT
   CHARACTER*4 PTYP
   LOGICAL :: OBACTIVE
   !
   DATA ANAME(1) /'   HYD. COND. ALONG ROWS'/
   DATA ANAME(2) /'  HORIZ. ANI. (COL./ROW)'/
   DATA ANAME(3) /'     VERTICAL HYD. COND.'/
   DATA ANAME(4) /' HORIZ. TO VERTICAL ANI.'/
   DATA ANAME(5) /'QUASI3D VERT. HYD. COND.'/
   DATA ANAME(6) /'        SPECIFIC STORAGE'/
   DATA ANAME(7) /'          SPECIFIC YIELD'/
   DATA ANAME(9) /'     STORAGE COEFFICIENT'/
   DATA WARN1 /' ***WARNING*** IPHDRY IS SET GREATER THAN ZERO. '/
   DATA WARN2 /'GROUNDWATER HEAD IN DRY CELLS WILL BE SET TO HDRY. '/
   DATA WARN3 /'OBSERVATION PACKAGES ARE ACTIVE AND HDRY COULD BE '/
   DATA WARN4 /'USED TO CALCULATE OBSERVATION VALUES. IPHDRY '/
   DATA WARN5 /'SHOULD BE SET TO ZERO IF OBSERVATION PACKAGES '/
   DATA WARN6 /'ARE ACTIVE.'/
   !     ------------------------------------------------------------------
   !1------IDENTIFY PACKAGE AND INITIALIZE.
   WRITE (Iout, 9001) In
9001 FORMAT (1X, /' UPW1 -- UPSTREAM WEIGHTING FLOW PACKAGE, ',&
   &'VERSION 1.3.0, 7/01/2022', /, 9X, 'INPUT READ FROM UNIT',&
   &I3,/)
   !  ALLOCATE, READ AND SET DATA FOR CELL PROPERTIES (FROM LPF)
   !1------Allocate scalar data.
   ALLOCATE(IUPWCB,NOVFC,Iuupw)
   !  STORE UPW UNIT NUMBER IN MODULE VARIABLE
   Iuupw = In
   ALLOCATE(ISFAC,ICONCV,ITHFLG,NOCVCO,IPHDRY)
   ZERO=0.
   !
   !3------READ COMMENTS AND ITEM 1.
   CALL URDCOM(IN,IOUT,LINE)
   LLOC=1
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUPWCB,R,IOUT,IN)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,IOUT,IN)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPUPW,R,IOUT,IN)
   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPHDRY,R,IOUT,IN)

   !
   !
   !------CHECK IF ANY OBS PACKAGES ARE ACTIVE AND IPHDRY>0
   OBACTIVE = .FALSE.
   IF ( IPHDRY > 0 ) THEN
      IF ( IUNIT(28)>0 ) OBACTIVE = .TRUE.
      IF ( IUNIT(33)>0 ) OBACTIVE = .TRUE.
      IF ( IUNIT(34)>0 ) OBACTIVE = .TRUE.
      IF ( IUNIT(35)>0 ) OBACTIVE = .TRUE.
      IF ( IUNIT(36)>0 ) OBACTIVE = .TRUE.
      IF ( IUNIT(38)>0 ) OBACTIVE = .TRUE.
   END IF
   IF ( OBACTIVE ) THEN
      WRITE(IOUT,111) WARN1
      WRITE(IOUT,111) WARN2
      WRITE(IOUT,111) WARN3
      WRITE(IOUT,111) WARN4
      WRITE(IOUT,111) WARN5
      WRITE(IOUT,111) WARN6
      WRITE(IOUT,*)
      WRITE(*,111) WARN1
      WRITE(*,111) WARN2
      WRITE(*,111) WARN3
      WRITE(*,111) WARN4
      WRITE(*,111) WARN5
      WRITE(*,111) WARN6
   END IF
111 FORMAT(A100)
   !
   !3A-----WRITE ITEM 1
   IF(IUPWCB.LT.0) WRITE(IOUT,8)
8  FORMAT(1X,'CONSTANT-HEAD CELL-BY-CELL FLOWS WILL BE PRINTED',&
   &' WHEN ICBCFL IS NOT 0')
   IF(IUPWCB.GT.0) WRITE(IOUT,9) IUPWCB
9  FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
   IF(NPUPW.GT.0) THEN
      WRITE(IOUT,15) NPUPW
15    FORMAT(1X,I5,' Named Parameters     ')
   ELSE
      NPUPW=0
      WRITE(IOUT,'(A)') ' No named parameters'
   END IF
   !
   !3B-----GET OPTIONS.
   ISFAC=0
   ICONCV=0
   ITHFLG=0
   NOCVCO=0
   NOVFC=0
   NOPCHK=0
   STOTXT=ANAME(6)
20 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
   IF(LINE(ISTART:ISTOP).EQ.'STORAGECOEFFICIENT') THEN
      ISFAC=1
      STOTXT=ANAME(9)
      WRITE(IOUT,21)
21    FORMAT(1X,'STORAGECOEFFICIENT OPTION:',/,&
      &1X,'Read storage coefficient rather than specific storage'&
      &1X,'Option not supported in UPW Package for MODFLOW-NWT')
   ELSE IF(LINE(ISTART:ISTOP).EQ.'CONSTANTCV') THEN
      ICONCV=1
      WRITE(IOUT,23)
23    FORMAT(1X,'CONSTANTCV OPTION:',/,1X,'Constant vertical',&
      &' conductance for convertible layers'&
      &1X,'Option not supported in UPW Package for MODFLOW-NWT')
   ELSE IF(LINE(ISTART:ISTOP).EQ.'THICKSTRT') THEN
      ITHFLG=1
      WRITE(IOUT,25)
25    FORMAT(1X,'THICKSTRT OPTION:',/,1X,'Negative LAYTYP indicates',&
      &' confined layer with thickness computed from STRT-BOT'&
      &1X,'Option not supported in UPW Package for MODFLOW-NWT')
   ELSE IF(LINE(ISTART:ISTOP).EQ.'NOCVCORRECTION') THEN
      NOCVCO=1
      WRITE(IOUT,27)
27    FORMAT(1X,'NOCVCORRECTION OPTION:',/,1X,&
      &'Do not adjust vertical conductance when applying',&
      &' the vertical flow correction'&
      &1X,'Option not supported, no vertical flow correction '&
      &' in UPW Package')
   ELSE IF(LINE(ISTART:ISTOP).EQ.'NOVFC') THEN
      NOVFC=1
      NOCVCO=1
      WRITE(IOUT,29)
29    FORMAT(1X,'NOVFC OPTION:',/,1X,&
      &'vertical flow correction does not apply in UPW Package')
   ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPARCHECK') THEN
      NOPCHK=1
      WRITE(IOUT,30)
30    FORMAT(1X,'NOPARCHECK  OPTION:',/,1X,&
      &'For data defined by parameters, do not check to see if ',&
      &'parameters define data at all cells')
   END IF
   IF(LLOC.LT.200) GO TO 20
   !
   !4------ALLOCATE AND READ LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET, LAYSTRT.
   ALLOCATE (Sn(Numcell), So(Numcell))
   Sn = 0.0D0
   So = 0.0D0
   ALLOCATE(LAYTYPUPW(NLAY))
   ALLOCATE(LAYAVG(NLAY))
   ALLOCATE(CHANI(NLAY))
   ALLOCATE(LAYVKAUPW(NLAY))
   ALLOCATE(LAYWET(NLAY))
   ALLOCATE(LAYSTRT(NLAY))
   ALLOCATE(IBOUND2(NCOL,NROW,NLAY))
   IBOUND2 = 0
   READ(IN,*) (LAYTYPUPW(K),K=1,NLAY)
   READ(IN,*) (LAYAVG(K),K=1,NLAY)
   READ(IN,*) (CHANI(K),K=1,NLAY)
   READ(IN,*) (LAYVKAUPW(K),K=1,NLAY)
   READ(IN,*) (LAYWET(K),K=1,NLAY)
   !
   !4A-----PRINT A TABLE OF VALUES FOR LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET.
   WRITE(IOUT,47)
47 FORMAT(1X,/3X,'LAYER FLAGS:',/1X,&
   &'LAYER       LAYTYP          LAYAVG    CHANI    ',&
   &'       LAYVKA           LAYWET',/1X,75('-'))
   DO 50 K=1,NLAY
      WRITE(IOUT,48) K,LAYTYPUPW(K),LAYAVG(K),CHANI(K),LAYVKAUPW(K),&
      &LAYWET(K)
48    FORMAT(1X,I4,2I14,1PE14.3,2I14)
   !
   !4A1----SET GLOBAL HEAD-DEPENDENT TRANSMISSIVITY AND STORAGE FLAGS.
      IF (LAYTYPUPW(K).GT.0) THEN
         LAYHDT(K)=1
         LAYHDS(K)=1
      ELSE
         LAYHDT(K)=0
         LAYHDS(K)=0
      ENDIF
50 CONTINUE
   !
   !4A2----SET LAYSTRT AND RESET LAYTYP IF THICKSTRT OPTION IS ACTIVE.
   DO 60 K=1,NLAY
      LAYSTRT(K)=0
      IF(LAYTYPUPW(K).LT.0 .AND. ITHFLG.NE.0) THEN
         LAYSTRT(K)=1
         LAYTYPUPW(K)=0
         LAYHDT(K)=0
         LAYHDS(K)=0
         WRITE(IOUT,57) K
57       FORMAT(1X,'Layer',I5,&
         &' is confined because LAYTYP<0 and THICKSTRT option is active')
      END IF
60 CONTINUE
   !
   !4B-----BASED ON LAYTYP, LAYAVG, CHANI, LAYWET, COUNT THE NUMBER OF EACH
   !4B-----TYPE OF 2-D ARRAY; CHECK VALUES FOR CONSISTENCY; AND SETUP
   !4B-----POINTERS IN LAYTYP, CHANI, AND LAYWET FOR CONVENIENT ACCESS
   !4B-----TO SC2, HANI, and WETDRY.  PRINT INTERPRETED VALUES OF FLAGS.
   NCNVRT=0
   NHANI=0
   NWETD=0
   WRITE(IOUT,67)
67 FORMAT(1X,/3X,'INTERPRETATION OF LAYER FLAGS:',/1X,&
   &'                       INTERBLOCK     HORIZONTAL',&
   &'    DATA IN',/1X,&
   &'        LAYER TYPE   TRANSMISSIVITY   ANISOTROPY',&
   &'   ARRAY VKA   WETTABILITY',/1X,&
   &'LAYER   (LAYTYP)        (LAYAVG)      (CHANI)   ',&
   &'  (LAYVKA)       (LAYWET)',/1X,75('-'))
   DO 100 K=1,NLAY
      IF(LAYTYPUPW(K).GT.0) THEN
         NCNVRT=NCNVRT+1
         LAYTYPUPW(K)=NCNVRT
      END IF
      IF(CHANI(K).LE.ZERO) THEN
         NHANI=NHANI+1
         CHANI(K)=-NHANI
      END IF
      IF(LAYWET(K).NE.0) THEN
   !         IF(LAYTYPUPW(K).EQ.0) THEN
         WRITE(IOUT,*)
         WRITE(IOUT,*)&
         &' LAYWET is not 0 and wetting does not apply in UPW '
         WRITE(IOUT,*) ' LAYWET must be 0 when using the UPW Package'
         CALL USTOP(' ')
   !         ELSE
   !            NWETD=NWETD+1
   !            LAYWET(K)=NWETD
   !         END IF
      END IF
      IF(LAYAVG(K).LT.0 .OR. LAYAVG(K).GT.2) THEN
         WRITE(IOUT,74) LAYAVG(K)
74       FORMAT(1X,I8,&
         &' IS AN INVALID LAYAVG VALUE -- MUST BE 0, 1, or 2')
         CALL USTOP(' ')
      END IF
      LAYPRN(1)=TYPNAM(1)
      IF(LAYTYPUPW(K).GT.0) LAYPRN(1)=TYPNAM(2)
      LAYPRN(2)=AVGNAM(LAYAVG(K)+1)
      IF(CHANI(K).LE.0) THEN
         LAYPRN(3)=HANNAM
      ELSE
         WRITE(LAYPRN(3),'(1PE14.3)') CHANI(K)
      END IF
      LAYPRN(4)=VKANAM(1)
      IF(LAYVKAUPW(K).NE.0) LAYPRN(4)=VKANAM(2)
      LAYPRN(5)=WETNAM(1)
      IF(LAYWET(K).NE.0) LAYPRN(5)=WETNAM(2)
      WRITE(IOUT,78) K,(LAYPRN(I),I=1,5)
78    FORMAT(1X,I4,5A)
100 CONTINUE
   !
   !4C-----PRINT WETTING INFORMATION.  RGN commented out because this does not apply
   !      IF(NWETD.EQ.0) THEN
   !         WRITE(IOUT,13)
   !   13    FORMAT(1X,/,1X,'WETTING CAPABILITY IS NOT ACTIVE IN ANY LAYER')
   !         IWDFLG=0
   !      ELSE
   !         WRITE(IOUT,12) NWETD
   !   12    FORMAT(1X,/,1X,'WETTING CAPABILITY IS ACTIVE IN',I4,' LAYERS')
   !         IWDFLG=1
   !         READ(IN,*) WETFCT,IWETIT,IHDWET
   !         IF(IWETIT.LE.0) IWETIT=1
   !         WRITE(IOUT,*) ' WETTING FACTOR=',WETFCT
   !         WRITE(IOUT,*) ' WETTING ITERATION INTERVAL=',IWETIT
   !         WRITE(IOUT,*) ' IHDWET=',IHDWET
   !      END IF
   ! ALLOCATE THESE BECAUSE THEY ARE LEFT IN THE CODE (DEALLOCATED)
   ALLOCATE(WETFCT,IWETIT,IHDWET,IWDFLG)
   !
   !5------ALLOCATE MEMORY FOR ARRAYS.
   ALLOCATE(LAYFLG(6,NLAY))
   ALLOCATE(HKUPW(NCOL,NROW,NLAY))
   ALLOCATE(VKAUPW(NCOL,NROW,NLAY))
   IF(NCNFBD.GT.0) THEN
      ALLOCATE(VKCB(NCOL,NROW,NCNFBD))
   ELSE
      ALLOCATE(VKCB(1,1,1))
   END IF
   IF(ITRSS.NE.0) THEN
      ALLOCATE(SC1(NCOL,NROW,NLAY))
   ELSE
      ALLOCATE(SC1(1,1,1))
   END IF
   SC1 = 0.0
   ! RGN 6/25/09      IF(ITRSS.NE.0 .AND. NCNVRT.GT.0) THEN
   IF(ITRSS.NE.0) THEN
   ! RGN 6/25/09        ALLOCATE(SC2UPW(NCOL,NROW,NCNVRT))
      ALLOCATE(SC2UPW(NCOL,NROW,NLAY))
   ELSE
      ALLOCATE(SC2UPW(1,1,1))
   END IF
   SC2UPW = 0.0
   IF(NHANI.GT.0) THEN
      ALLOCATE(HANI(NCOL,NROW,NHANI))
   ELSE
      ALLOCATE(HANI(1,1,1))
   END IF
   IF(NWETD.GT.0) THEN
      ALLOCATE(WETDRY(NCOL,NROW,NWETD))
   ELSE
      ALLOCATE(WETDRY(1,1,1))
   END IF
   !
   !6------READ PARAMETER DEFINITIONS
   NPHK=0
   NPVKCB=0
   NPVK=0
   NPVANI=0
   NPSS=0
   NPSY=0
   NPHANI=0
   IF(NPUPW.GT.0) THEN
      WRITE(IOUT,115)
115   FORMAT(/,' PARAMETERS DEFINED IN THE UPW PACKAGE')
      DO 120 K=1,NPUPW
         CALL UPARARRRP(IN,IOUT,N,1,PTYP,1,0,-1)
   !   Note that NPHK and the other NP variables in
   !   this group are used only as flags, not counts
         IF(PTYP.EQ.'HK') THEN
            NPHK=1
         ELSE IF(PTYP.EQ.'HANI') THEN
   !6A-----WHEN A HANI PARAMETER IS USED, THEN ALL HORIZONTAL ANISOTROPY
   !6A-----MUST BE DEFINED USING PARAMETERS.  ENSURE THAT ALL CHANI <= 0
            DO 118 I = 1, NLAY
               IF (CHANI(I).GT.0.0) THEN
                  WRITE(IOUT,117)
117               FORMAT(/,&
                  &'ERROR: WHEN A HANI PARAMETER IS USED, CHANI FOR ALL LAYERS',/,&
                  &'MUST BE LESS THAN OR EQUAL TO 0.0 -- STOP EXECUTION',&
                  &' (GWF2UPW1AR)')
                  CALL USTOP(' ')
               ENDIF
118         CONTINUE
            NPHANI=1
         ELSE IF(PTYP.EQ.'VKCB') THEN
            NPVKCB=1
         ELSE IF(PTYP.EQ.'VK') THEN
            NPVK=1
            CALL SGWF2UPWCK(IOUT,N,'VK  ')
         ELSE IF(PTYP.EQ.'VANI') THEN
            NPVANI=1
            CALL SGWF2UPWCK(IOUT,N,'VANI')
         ELSE IF(PTYP.EQ.'SS') THEN
            NPSS=1
         ELSE IF(PTYP.EQ.'SY') THEN
            NPSY=1
         ELSE
            WRITE(IOUT,*) ' Invalid parameter type for UPW Package'
            CALL USTOP(' ')
         END IF
120   CONTINUE
   END IF
   !
   !7------DEFINE DATA FOR EACH LAYER -- VIA READING OR NAMED PARAMETERS.
   DO 200 K=1,NLAY
      KK=K
   !
   !7A-----DEFINE HORIZONTAL HYDRAULIC CONDUCTIVITY (HK)
      IF(NPHK.EQ.0) THEN
         CALL U2DREL(HKUPW(:,:,KK),ANAME(1),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(1,K)
         WRITE(IOUT,121) ANAME(1),K,LAYFLG(1,K)
121      FORMAT(1X,/1X,A,' FOR LAYER',I4,&
         &' WILL BE DEFINED BY PARAMETERS',/1X,'(PRINT FLAG=',I4,')')
         CALL UPARARRSUB1(HKUPW(:,:,KK),NCOL,NROW,KK,'HK',&
         &IOUT,ANAME(1),LAYFLG(1,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,NLAY,&
         &NROW,'HK  ')
      END IF
   !
   !7B-----READ HORIZONTAL ANISOTROPY IF CHANI IS NON-ZERO
      IF(CHANI(K).LE.ZERO) THEN
         KHANI=-CHANI(K)
         IF(NPHANI.EQ.0) THEN
            CALL U2DREL(HANI(:,:,KHANI),ANAME(2),NROW,NCOL,KK,IN,IOUT)
         ELSE
            READ(IN,*) LAYFLG(6,K)
            WRITE(IOUT,121) ANAME(2),K,LAYFLG(6,K)
            CALL UPARARRSUB1(HANI(:,:,KHANI),NCOL,NROW,KK,'HANI',&
            &IOUT,ANAME(2),LAYFLG(6,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,&
            &NLAY,NROW,'HANI')
         END IF
      END IF
   !
   !7C-----DEFINE VERTICAL HYDRAULIC CONDUCTIVITY OR HORIZONTAL TO VERTICAL
   !7C-----ANISOTROPY (VKA).
      IANAME=3
      PTYP='VK'
      IF(LAYVKAUPW(K).NE.0) THEN
         IANAME=4
         PTYP='VANI'
      END IF
      IF(NPVK.EQ.0 .AND. NPVANI.EQ.0) THEN
         CALL U2DREL(VKAUPW(:,:,KK),ANAME(IANAME),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(2,K)
         WRITE(IOUT,121) ANAME(IANAME),K,LAYFLG(2,K)
         CALL UPARARRSUB1(VKAUPW(:,:,KK),NCOL,NROW,KK,PTYP,IOUT,&
         &ANAME(IANAME),LAYFLG(2,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,NLAY,&
         &NROW,PTYP)
      END IF
   !
   !7D-----DEFINE SPECIFIC STORAGE OR STORAGE COEFFICIENT IN ARRAY SC1 IF TRANSIENT.
      IF(ITRSS.NE.0) THEN
         IF(NPSS.EQ.0) THEN
            CALL U2DREL(SC1(:,:,KK),STOTXT,NROW,NCOL,KK,IN,IOUT)
         ELSE
            READ(IN,*) LAYFLG(3,K)
            WRITE(IOUT,121) STOTXT,K,LAYFLG(3,K)
            CALL UPARARRSUB1(SC1(:,:,KK),NCOL,NROW,KK,'SS',&
            &IOUT,STOTXT,LAYFLG(3,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,&
            &NLAY,NROW,'SS  ')
         END IF
         IF(ISFAC.EQ.0) THEN
            CALL SGWF2UPWSC(SC1(:,:,KK),KK,1)
         ELSE
            CALL SGWF2UPWSC(SC1(:,:,KK),KK,0)
         END IF
      END IF
   !
   !7E-----DEFINE SPECIFIC YIELD IN ARRAY SC2 IF TRANSIENT AND LAYER IS
   !7E-----IS CONVERTIBLE.
      IF(LAYTYPUPW(K).GT.0) THEN
         IF(ITRSS.NE.0) THEN
            IF(NPSY.EQ.0) THEN
               CALL U2DREL(SC2UPW(:,:,LAYTYPUPW(K)),ANAME(7),NROW,NCOL,KK,&
               &IN,IOUT)
            ELSE
               READ(IN,*) LAYFLG(4,K)
               WRITE(IOUT,121) ANAME(7),K,LAYFLG(4,K)
               CALL UPARARRSUB1(SC2UPW(:,:,LAYTYPUPW(K)),NCOL,&
               &NROW,KK,'SY',IOUT,ANAME(7),LAYFLG(4,KK))
               IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,&
               &NCOL,NLAY,NROW,'SY  ')
            END IF
            CALL SGWF2UPWSC(SC2UPW(:,:,LAYTYPUPW(K)),KK,0)
         END IF
      ELSE
         IF(ITRSS.NE.0) THEN
            DO J=1,NROW
               DO I=1,NCOL
                  SC2UPW(I,J,KK) = 0.0D0    !SC1(I,J,KK)
               END DO
            END DO
         END IF
      END IF
   !
   !7F-----READ CONFINING BED VERTICAL HYDRAULIC CONDUCTIVITY (VKCB)
      IF(LAYCBD(K).NE.0) THEN
         IF(NPVKCB.EQ.0) THEN
            CALL U2DREL(VKCB(:,:,LAYCBD(K)),ANAME(5),NROW,NCOL,KK,IN,&
            &IOUT)
         ELSE
            READ(IN,*) LAYFLG(5,K)
            WRITE(IOUT,121) ANAME(5),K,LAYFLG(5,K)
            CALL UPARARRSUB1(VKCB(:,:,LAYCBD(K)),NCOL,NROW,KK,&
            &'VKCB',IOUT,ANAME(5),LAYFLG(5,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,&
            &NLAY,NROW,'VKCB')
         END IF
      END IF
   !
   !7G-----READ WETDRY CODES IF WETTING CAPABILITY HAS BEEN INVOKED
   !7G-----(LAYWET NOT 0).
      IF(LAYWET(K).NE.0) THEN
         CALL U2DREL(WETDRY(:,:,LAYWET(K)),ANAME(8),NROW,NCOL,KK,IN,&
         &IOUT)
      END IF
200 CONTINUE
   !
   !8------PREPARE AND CHECK LPF DATA.
   CALL SGWF2UPWN()
   !9------Calculate constant part of conductance. Conductance includes
   !       cell thickness for confined conditions.
   DO K=1,NLAY
      IF(LAYAVG(K).EQ.0) THEN
         IF ( LAYTYPUPW(K).GT.0 ) THEN
            CALL SGWF2UPW1HHARM(K)
         ELSE
            CALL SGWF2UPW1HHARMCON(K)
         END IF
      ELSE IF(LAYAVG(K).EQ.1) THEN
         IF ( LAYTYPUPW(K).GT.0 ) THEN
            CALL SGWF2UPW1HLOG(K)
         ELSE
            CALL SGWF2UPW1HLOGCON(K)
         END IF
      ELSE IF(LAYAVG(K).EQ.2) THEN
         IF ( LAYTYPUPW(K).GT.0 ) THEN
            CALL SGWF2UPW1HUNCNF(K)
         ELSE
            CALL SGWF2UPW1HUNCNFCON(K)
         END IF
      END IF
      CALL SGWF2UPW1VCOND(K)
   END DO
   !10-----SAVE POINTERS FOR GRID AND RETURN
   CALL SGWF2UPW1PSV(Igrid)
   !
   !11-----RETURN
END SUBROUTINE GWF2UPW1AR
   !     -----------------------------------------------------------------
SUBROUTINE SGWF2UPWN()
   !     ******************************************************************
   !     INITIALIZE AND CHECK UPW DATA
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,LAYCBD,CV,&
   &BOTM,NBOTM,DELR,DELC,IOUT
   USE GWFBASMODULE,ONLY:HNOFLO
   USE GWFUPWMODULE,ONLY:HKUPW,VKAUPW
   !     ------------------------------------------------------------------
   !
   !1------DEFINE CONSTANTS.
   ZERO=0.
   HCNV=HNOFLO
   !
   !2------INSURE THAT EACH ACTIVE CELL HAS AT LEAST ONE NON-ZERO
   !2------TRANSMISSIVE PARAMETER.
   DO 60 K=1,NLAY
      DO 40 I=1,NROW
         DO 40 J=1,NCOL
            IF(IBOUND(J,I,K).EQ.0 ) GO TO 40
   !
   !3A-----CHECK HORIZONTAL HYDRAULIC CONDUCTIVITY (HK).
            IF(HKUPW(J,I,K).NE.ZERO) GO TO 40
   !
   !3B-----CHECK VERTICAL HYDRAULIC CONDUCTIVITY AND CONFINING BED
   !3B-----VERTICAL HYDRAULIC CONDUCTIVITY.
            IF(NLAY.GT.1) THEN

               IF(VKAUPW(J,I,K).NE.ZERO) THEN
                  IF(K.NE.NLAY) THEN
                     IF (VKAUPW(J,I,K+1).NE.ZERO) GO TO 40
                  END IF
                  IF(K.NE.1) THEN
                     IF (VKAUPW(J,I,K-1).NE.ZERO) GO TO 40
                  END IF
               END IF
            END IF
   !
   !3C-----ALL TRANSMISSIVE TERMS ARE ALL 0, SO CONVERT CELL TO NO FLOW.
            IBOUND(J,I,K)=0
            HNEW(J,I,K)=HCNV
            WRITE(IOUT,43) K,I,J
40    CONTINUE
43    FORMAT(1X,'NODE (LAYER,ROW,COL) ',I3,2(1X,I5),&
      &' ELIMINATED BECAUSE ALL HYDRAULIC',/,&
      &' CONDUCTIVITIES TO NODE ARE 0')
   !
60 CONTINUE
   !
   !7------RETURN.
   RETURN
END
   !
SUBROUTINE GWF2UPWFMS(KITER,KSTP,KPER,IGRID)
   !     ******************************************************************
   !     ADD LEAKAGE CORRECTION AND STORAGE TO HCOF AND RHS, AND CALCULATE
   !     CONDUCTANCE AS REQUIRED.
   !     ******************************************************************
   !
   !        SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,BOTM,NBOTM,DELR,DELC,&
   &LBOTM,CV,HNEW,RHS,HCOF,HOLD,ISSFLG,IOUT
   USE GWFBASMODULE,ONLY:DELT
   USE GWFNWTMODULE,ONLY:A, ICELL, IA
   USE GWFUPWMODULE
   DOUBLE PRECISION, EXTERNAL :: DHORIZUPW
   DOUBLE PRECISION ZERO
   DOUBLE PRECISION HTMP, TP, BT, TLED, ONE, SOLD, SNEW, STRG
   DOUBLE PRECISION RHO1, RHO2, HLD, THICK, dS
   !     ------------------------------------------------------------------
   !
   !1------SET POINTERS TO DATA, GET STEADY-STATE FLAG FOR STRESS PERIOD,
   !1------DEFINE CONSTANT.
   CALL SGWF2UPW1PNT(IGRID)
   ISS=ISSFLG(KPER)
   ONE=1.0D0
   ZERO = 0.0D0
   !
   !2------IF THE STRESS PERIOD IS TRANSIENT, ADD STORAGE TO HCOF AND RHS
   IF(ISS.EQ.0) THEN
      TLED=ONE/DBLE(DELT)
      DO 200 K=1,NLAY
   !3------CHECK OLD AND NEW HEADS TO DETERMINE
   !3------WHEN TO USE PRIMARY AND SECONDARY STORAGE
         DO 180 I=1,NROW
            DO 180 J=1,NCOL
   !
   !4A-----IF THE CELL IS EXTERNAL THEN SKIP IT.
               IF(IBOUND(J,I,K).LE.0) GO TO 180
               TP=dble(BOTM(J,I,LBOTM(K)-1))
               BT=dble(BOTM(J,I,LBOTM(K)))
               THICK = (TP-BT)
               HTMP = HNEW(J,I,K)
               HLD = DBLE(HOLD(J,I,K))
               RHO2 = SC2UPW(J,I,K)*TLED
   !            IF ( HTMP.GT.TP ) RHO2 = SC1(J,I,K)*TLED
               RHO1 = SC1(J,I,K)*TLED
   !            IF ( HLD.GT.TP ) RHO1 = SC1(J,I,K)*TLED
               dS = DHORIZUPW(HTMP, TP, BT, K)

   ! Add derivatives to jacobian
               ij = Icell(J,I,K)
   ! Derivative for HCOF
               A(IA(ij)) = A(IA(ij)) - dS*THICK*RHO1*HTMP
   ! Derivative for RHS
               A(IA(ij)) = A(IA(ij)) - THICK*RHO2*dS&
               &+ THICK*RHO1*dS*HLD

   !
   !6------ADD STORAGE TERMS TO RHS AND HCOF.
               HCOF(J,I,K) = HCOF(J,I,K) - Sn(ij)*THICK*RHO1
               RHS(J,I,K) = RHS(J,I,K) + THICK*RHO2*(Sn(ij)-So(ij)) -&
               &Sn(ij)*THICK*RHO1*HLD
   !
180      CONTINUE
   !
200   CONTINUE
   END IF
   !
   !10-----RETURN
   RETURN
END
   !
SUBROUTINE GWF2UPWBDADJ(KSTP,KPER,IDIR,IBDRET,&
&IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
   !     ******************************************************************
   !     COMPUTE FLOW BETWEEN ADJACENT CELLS IN A SUBREGION OF THE GRID
   !     ******************************************************************
   !
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,CR,CC,CV,&
   &BOTM,LBOTM,IOUT
   USE GWFBASMODULE,ONLY:ICBCFL,DELT,PERTIM,TOTIM,ICHFLG
   USE GWFUPWMODULE,ONLY:IUPWCB,Sn,LAYTYPUPW
   USE GWFNWTMODULE,ONLY:ICELL
   !
   CHARACTER*16 TEXT(3)
   DOUBLE PRECISION HD, ttop, bbot, HDIFF, TMP, TOP
   INTEGER iltyp
   !
   DATA TEXT(1),TEXT(2),TEXT(3)&
   &/'FLOW RIGHT FACE ','FLOW FRONT FACE ','FLOW LOWER FACE '/
   !     ------------------------------------------------------------------
   !
   CALL SGWF2UPW1PNT(IGRID)
   !
   !1------IF CELL-BY-CELL FLOWS WILL BE SAVED IN A FILE, SET FLAG IBD.
   !1------RETURN IF FLOWS ARE NOT BEING SAVED OR RETURNED.
   ZERO=0.
   IBD=0
   IF(IUPWCB.GT.0) IBD=ICBCFL
   IF(IBD.EQ.0 .AND. IBDRET.EQ.0) RETURN
   !
   !2------SET THE SUBREGION EQUAL TO THE ENTIRE GRID IF VALUES ARE BEING
   !2------SAVED IN A FILE.
   IF(IBD.NE.0) THEN
      K1=1
      K2=NLAY
      I1=1
      I2=NROW
      J1=1
      J2=NCOL
   END IF
   !
   !3------TEST FOR DIRECTION OF CALCULATION;  IF NOT ACROSS COLUMNS, GO TO
   !3------STEP 4.  IF ONLY 1 COLUMN, RETURN.
   IF(IDIR.NE.1) GO TO 405
   IF(NCOL.EQ.1) RETURN
   !
   !3A-----CALCULATE FLOW ACROSS COLUMNS (THROUGH RIGHT FACE).  IF NOT
   !3A-----SAVING IN A FILE, SET THE SUBREGION.  CLEAR THE BUFFER.
   IF(IBD.EQ.0) THEN
      K1=IL1
      K2=IL2
      I1=IR1
      I2=IR2
      J1=IC1-1
      IF(J1.LT.1) J1=1
      J2=IC2
   END IF
   DO 310 K=K1,K2
      DO 310 I=I1,I2
         DO 310 J=J1,J2
            BUFF(J,I,K)=ZERO
310 CONTINUE
   !
   !3B-----FOR EACH CELL CALCULATE FLOW THRU RIGHT FACE & STORE IN BUFFER.
   IF(J2.EQ.NCOL) J2=J2-1
   DO 400 K=K1,K2
      iltyp = LAYTYPUPW(K)
      DO 400 I=I1,I2
         DO 400 J=J1,J2
            IF(ICHFLG.EQ.0) THEN
               IF((IBOUND(J,I,K).LE.0) .AND. (IBOUND(J+1,I,K).LE.0)) GO TO 400
               IF((IBOUND(J,I,K).EQ.0) .OR. (IBOUND(J+1,I,K).EQ.0)) GO TO 400
            ELSE
               IF((IBOUND(J,I,K).EQ.0) .OR. (IBOUND(J+1,I,K).EQ.0)) GO TO 400
            END IF
            HDIFF=HNEW(J,I,K)-HNEW(J+1,I,K)
            IF ( HDIFF.GT.0.0 ) THEN
               IF ( iltyp.GT.0 ) THEN
                  TTOP = dble(BOTM(J,I,LBOTM(K)-1))
                  BBOT = dble(BOTM(J,I,LBOTM(K)))
                  ij = Icell(J,I,K)
                  BUFF(J,I,K)=HDIFF*CR(J,I,K)*&
                  &(TTOP-BBOT)*Sn(ij)
               ELSE
                  BUFF(J,I,K)=HDIFF*CR(J,I,K)
               END IF
            ELSE
               IF ( iltyp.GT.0 ) THEN
                  TTOP = dble(BOTM(J+1,I,LBOTM(K)-1))
                  BBOT = dble(BOTM(J+1,I,LBOTM(K)))
                  ij = Icell(J+1,I,K)
                  BUFF(J,I,K)=HDIFF*CR(J,I,K)*(TTOP-BBOT)*Sn(ij)
               ELSE
                  BUFF(J,I,K)=HDIFF*CR(J,I,K)
               END IF
            END IF
400 CONTINUE
   !
   !3C-----RECORD CONTENTS OF BUFFER AND RETURN.
   IF(IBD.EQ.1)&
   &CALL UBUDSV(KSTP,KPER,TEXT(1),IUPWCB,BUFF,NCOL,NROW,NLAY,IOUT)
   IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT(1),IUPWCB,BUFF,NCOL,NROW,&
   &NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
   RETURN
   !
   !4------TEST FOR DIRECTION OF CALCULATION;  IF NOT ACROSS ROWS, GO TO
   !4------STEP 5.  IF ONLY 1 ROW, RETURN.
405 IF(IDIR.NE.2) GO TO 505
   IF(NROW.EQ.1) RETURN
   !
   !4A-----CALCULATE FLOW ACROSS ROWS (THROUGH FRONT FACE).  IF NOT SAVING
   !4A-----IN A FILE, SET THE SUBREGION.  CLEAR THE BUFFER.
   IF(IBD.EQ.0) THEN
      K1=IL1
      K2=IL2
      I1=IR1-1
      IF(I1.LT.1) I1=1
      I2=IR2
      J1=IC1
      J2=IC2
   END IF
   DO 410 K=K1,K2
      DO 410 I=I1,I2
         DO 410 J=J1,J2
            BUFF(J,I,K)=ZERO
410 CONTINUE
   !
   !4B-----FOR EACH CELL CALCULATE FLOW THRU FRONT FACE & STORE IN BUFFER.
   IF(I2.EQ.NROW) I2=I2-1
   DO 500 K=K1,K2
      iltyp = LAYTYPUPW(K)
      DO 500 I=I1,I2
         DO 500 J=J1,J2
            IF(ICHFLG.EQ.0) THEN
               IF((IBOUND(J,I,K).LE.0) .AND. (IBOUND(J,I+1,K).LE.0)) GO TO 500
               IF((IBOUND(J,I,K).EQ.0) .OR. (IBOUND(J,I+1,K).EQ.0)) GO TO 500
            ELSE
               IF((IBOUND(J,I,K).EQ.0) .OR. (IBOUND(J,I+1,K).EQ.0)) GO TO 500
            END IF
            HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
            IF ( HDIFF.GT.0.0 ) THEN
               IF ( iltyp.GT.0 ) THEN
                  TTOP = dble(BOTM(J,I,LBOTM(K)-1))
                  BBOT = dble(BOTM(J,I,LBOTM(K)))
                  ij = Icell(J,I,K)
                  BUFF(J,I,K)=HDIFF*CC(J,I,K)*(TTOP-BBOT)*Sn(ij)
               ELSE
                  BUFF(J,I,K)=HDIFF*CC(J,I,K)
               END IF
            ELSE
               IF ( iltyp.GT.0 ) THEN
                  TTOP = dble(BOTM(J,I+1,LBOTM(K)-1))
                  BBOT = dble(BOTM(J,I+1,LBOTM(K)))
                  ij = Icell(J,I+1,K)
                  BUFF(J,I,K)=HDIFF*CC(J,I,K)*(TTOP-BBOT)*Sn(ij)
               ELSE
                  BUFF(J,I,K)=HDIFF*CC(J,I,K)
               END IF
            END IF
500 CONTINUE
   !
   !4C-----RECORD CONTENTS OF BUFFER AND RETURN.
   IF(IBD.EQ.1)&
   &CALL UBUDSV(KSTP,KPER,TEXT(2),IUPWCB,BUFF,NCOL,NROW,NLAY,IOUT)
   IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT(2),IUPWCB,BUFF,NCOL,NROW,&
   &NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
   RETURN
   !
   !5------DIRECTION OF CALCULATION IS ACROSS LAYERS BY ELIMINATION.  IF
   !5------ONLY 1 LAYER, RETURN.
505 IF(NLAY.EQ.1) RETURN
   !
   !5A-----CALCULATE FLOW ACROSS LAYERS (THROUGH LOWER FACE).  IF NOT
   !5A-----SAVING IN A FILE, SET THE SUBREGION.  CLEAR THE BUFFER.
   IF(IBD.EQ.0) THEN
      K1=IL1-1
      IF(K1.LT.1) K1=1
      K2=IL2
      I1=IR1
      I2=IR2
      J1=IC1
      J2=IC2
   END IF
   DO 510 K=K1,K2
      DO 510 I=I1,I2
         DO 510 J=J1,J2
            BUFF(J,I,K)=ZERO
510 CONTINUE
   !
   !5B-----FOR EACH CELL CALCULATE FLOW THRU LOWER FACE & STORE IN BUFFER.
   IF(K2.EQ.NLAY) K2=K2-1
   DO 600 K=1,K2
      IF(K.LT.K1) GO TO 600
      DO 590 I=I1,I2
         DO 590 J=J1,J2
            IF(ICHFLG.EQ.0) THEN
               IF((IBOUND(J,I,K).LE.0) .AND. (IBOUND(J,I,K+1).LE.0)) GO TO 590
               IF((IBOUND(J,I,K).EQ.0) .OR. (IBOUND(J,I,K+1).EQ.0)) GO TO 590
            ELSE
               IF((IBOUND(J,I,K).EQ.0) .OR. (IBOUND(J,I,K+1).EQ.0)) GO TO 590
            END IF
            HD=HNEW(J,I,K+1)
            TMP=HD
            TOP=dble(BOTM(J,I,LBOTM(K+1)-1))
            !     IF(TMP.LT.TOP) HD=TOP
580         HDIFF=HNEW(J,I,K)-HD
            BUFF(J,I,K)=sngl(HDIFF*CV(J,I,K))
590   CONTINUE
600 CONTINUE
   !
   !5C-----RECORD CONTENTS OF BUFFER AND RETURN.
   IF(IBD.EQ.1)&
   &CALL UBUDSV(KSTP,KPER,TEXT(3),IUPWCB,BUFF,NCOL,NROW,NLAY,IOUT)
   IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT(3),IUPWCB,BUFF,NCOL,NROW,&
   &NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
   RETURN
END
   !
SUBROUTINE GWF2UPWBDS(KSTP,KPER,IGRID)
   !     ******************************************************************
   !     COMPUTE STORAGE BUDGET FLOW TERM.
   !     ******************************************************************
   !
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,ISSFLG,IBOUND,HNEW,HOLD,&
   &BUFF,BOTM,LBOTM,IOUT
   USE GWFBASMODULE,ONLY:MSUM,ICBCFL,VBVL,VBNM,DELT,PERTIM,TOTIM
   USE GWFUPWMODULE,ONLY:IUPWCB,SC1,SC2UPW,So,Sn,IBOUND2,LAYTYPUPW
   USE GWFNWTMODULE,ONLY: Thickfact,Icell
   CHARACTER*16 TEXT
   DOUBLE PRECISION STOIN,STOUT,SSTRG,ZERO,HSING,THICK
   DOUBLE PRECISION STRG, SIN, SOUT
   DOUBLE PRECISION BT, TP, RHO1, RHO2, HLD, BBT, TLED, ONE
   !
   DATA TEXT /'         STORAGE'/
   !     ------------------------------------------------------------------
   !
   CALL SGWF2UPW1PNT(IGRID)
   !
   !1------INITIALIZE BUDGET ACCUMULATORS AND 1/DELT.
   ISS=ISSFLG(KPER)
   ZERO=0.0D0
   STOIN=ZERO
   STOUT=ZERO
   !2------IF STEADY STATE, STORAGE TERM IS ZERO
   IF(ISS.NE.0) GOTO 400
   ONE=1.0D0
   TLED=ONE/DBLE(DELT)
   !
   !3------IF CELL-BY-CELL FLOWS WILL BE SAVED, SET FLAG IBD.
   IBD=0
   IF(IUPWCB.GT.0) IBD=ICBCFL
   !
   !4------CLEAR BUFFER.
   DO 210 K=1,NLAY
      DO 210 I=1,NROW
         DO 210 J=1,NCOL
            BUFF(J,I,K)=ZERO
210 CONTINUE
   !
   !5------LOOP THROUGH EVERY CELL IN THE GRID.
   KT=0
   DO 300 K=1,NLAY
      DO 300 I=1,NROW
         DO 300 J=1,NCOL
   !
   !6------SKIP NO-FLOW AND CONSTANT-HEAD CELLS.
            IF(IBOUND(J,I,K).LE.0) GO TO 300
            HSING=HNEW(J,I,K)
            HLD = DBLE(HOLD(J,I,K))
   !
   !7A----TWO STORAGE CAPACITIES.
            TP=dble(BOTM(J,I,LBOTM(K)-1))
            BT=dble(BOTM(J,I,LBOTM(K)))
            THICK = (TP-BT)
            RHO2 = dble(SC2UPW(J,I,K))*TLED
            ij = Icell(J,I,K)
            RHO1 = dble(SC1(J,I,K))*TLED
            STRG= - THICK*RHO2*(Sn(ij)-So(ij)) - Sn(ij)*THICK*RHO1*(HSING-HLD)
            IF ( LAYTYPUPW(K).GT.0 ) THEN
               IF ( HSING.LE.BT .AND. HLD.LE.BT ) THEN
                  IBOUND2(J,I,K) = 0
               ELSE
                  IBOUND2(J,I,K) = IBOUND(J,I,K)
               END IF
            END IF
            GO TO 288
   !
   !7B----ONE STORAGE CAPACITY.
   !  285 RHO=SC1(J,I,K)*TLED
   !      STRG=RHO*HOLD(J,I,K) - RHO*HSING

   !
   !8-----STORE CELL-BY-CELL FLOW IN BUFFER AND ADD TO ACCUMULATORS.
288         BUFF(J,I,K)=SNGL(STRG)
            SSTRG=STRG
            IF(STRG.LT.ZERO) THEN
               STOUT=STOUT-SSTRG
            ELSE
               STOIN=STOIN+SSTRG
            END IF
   !
300 CONTINUE
   !
   !9-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
   IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,&
   &IUPWCB,BUFF,NCOL,NROW,NLAY,IOUT)
   IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT,IUPWCB,&
   &BUFF,NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND2)
   !
   !10-----ADD TOTAL RATES AND VOLUMES TO VBVL & PUT TITLE IN VBNM.
400 CONTINUE
   SIN=STOIN
   SOUT=STOUT
   VBVL(1,MSUM)=VBVL(1,MSUM)+SNGL(SIN)*DELT
   VBVL(2,MSUM)=VBVL(2,MSUM)+SNGL(SOUT)*DELT
   VBVL(3,MSUM)=SNGL(SIN)
   VBVL(4,MSUM)=SNGL(SOUT)
   VBNM(MSUM)=TEXT
   MSUM=MSUM+1
   !
   !11----RETURN.
   RETURN
END
   !
SUBROUTINE GWF2UPWBDCH(KSTP,KPER,IGRID)
   !     ******************************************************************
   !     COMPUTE FLOW FROM CONSTANT-HEAD CELLS
   !     ******************************************************************
   !
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,CR,CC,CV,&
   &BOTM,LBOTM,IOUT
   USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,DELT,PERTIM,TOTIM,ICBCFL,&
   &ICHFLG
   USE GWFUPWMODULE,ONLY:IUPWCB, Sn, LAYTYPUPW
   USE GWFNWTMODULE,ONLY:Icell, Closezero

   CHARACTER*16 TEXT
   DOUBLE PRECISION HD,CHIN,CHOUT,XX1,XX2,XX3,XX4,XX5,XX6,ZERO
   DOUBLE PRECISION THICK,HDIFF
   DOUBLE PRECISION X1,X2,X3,X4,X5,X6,CIN,COUT
   DOUBLE PRECISION CHCH1,CHCH2,CHCH3,CHCH4,CHCH5,CHCH6
   !
   DATA TEXT /'   CONSTANT HEAD'/
   INTEGER iltyp
   !     ------------------------------------------------------------------
   CALL SGWF2NWT1PNT(IGRID)
   !
   !1------SET IBD TO INDICATE IF CELL-BY-CELL BUDGET VALUES WILL BE SAVED.
   IBD=0
   IF(IUPWCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
   IF(IUPWCB.GT.0) IBD=ICBCFL
   !
   !2------CLEAR BUDGET ACCUMULATORS.
   ZERO=0.0D0
   CHIN=ZERO
   CHOUT=ZERO
   IBDLBL=0
   !
   !3------CLEAR BUFFER.
   DO 5 K=1,NLAY
      DO 5 I=1,NROW
         DO 5 J=1,NCOL
            BUFF(J,I,K)=ZERO
5  CONTINUE
   !
   !3A-----IF SAVING CELL-BY-CELL FLOW IN A LIST, COUNT CONSTANT-HEAD
   !3A-----CELLS AND WRITE HEADER RECORDS.
   IF(IBD.EQ.2) THEN
      NCH=0
      DO 7 K=1,NLAY
         DO 7 I=1,NROW
            DO 7 J=1,NCOL
               IF(IBOUND(J,I,K).LT.0) NCH=NCH+1
7     CONTINUE
      CALL UBDSV2(KSTP,KPER,TEXT,IUPWCB,NCOL,NROW,NLAY,&
      &NCH,IOUT,DELT,PERTIM,TOTIM,IBOUND)
   END IF
   !
   !4------LOOP THROUGH EACH CELL AND CALCULATE FLOW INTO MODEL FROM EACH
   !4------CONSTANT-HEAD CELL.
   DO 200 K=1,NLAY
      iltyp = LAYTYPUPW(K)
      DO 200 I=1,NROW
         DO 200 J=1,NCOL

   !
   !5------IF CELL IS NOT CONSTANT HEAD SKIP IT & GO ON TO NEXT CELL.
            IF (IBOUND(J,I,K).GE.0)GO TO 200
   !
   !6------CLEAR VALUES FOR FLOW RATE THROUGH EACH FACE OF CELL.
            X1=ZERO
            X2=ZERO
            X3=ZERO
            X4=ZERO
            X5=ZERO
            X6=ZERO
            CHCH1=ZERO
            CHCH2=ZERO
            CHCH3=ZERO
            CHCH4=ZERO
            CHCH5=ZERO
            CHCH6=ZERO
   !
   !7------CALCULATE FLOW THROUGH THE LEFT FACE.
   !7------COMMENTS A-C APPEAR ONLY IN THE SECTION HEADED BY COMMENT 7,
   !7------BUT THEY APPLY IN A SIMILAR MANNER TO SECTIONS 8-12.
   !
   !7A-----IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
   !7A-----TO NEXT FACE.  NO FLOW OCCURS AT THE EDGE OF THE GRID, TO AN
   !7A-----ADJACENT NO-FLOW CELL, ??OR TO AN ADJACENT CONSTANT-HEAD CELL??.
            IF(J.EQ.1) GO TO 30
            IF(IBOUND(J-1,I,K).EQ.0) GO TO 30
            IF(IBOUND(J-1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 30
   !
   !7B-----CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
            HDIFF=HNEW(J,I,K)-HNEW(J-1,I,K)
            IF ( HDIFF.GE.-Closezero ) THEN
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J,I,LBOTM(K)-1)) - dble(BOTM(J,I,LBOTM(K)))
                  ij = Icell(J,I,K)
                  CHCH1=HDIFF*CR(J-1,I,K)*THICK*Sn(ij)
               ELSE
                  CHCH1=HDIFF*CR(J-1,I,K)
               END IF
            ELSE
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J-1,I,LBOTM(K)-1)) -&
                  &dble(BOTM(J-1,I,LBOTM(K)))
                  ij = Icell(J-1,I,K)
                  CHCH1=HDIFF*CR(J-1,I,K)*THICK*Sn(ij)
               ELSE
                  CHCH1=HDIFF*CR(J-1,I,K)
               END IF
            END IF
            IF(IBOUND(J-1,I,K).LT.0) GO TO 30
            X1=CHCH1
            XX1=X1
   !
   !7C-----ACCUMULATE POSITIVE AND NEGATIVE FLOW.
            IF(X1.LT.ZERO) THEN
               CHOUT=CHOUT-XX1
            ELSE
               CHIN=CHIN+XX1
            END IF
   !
   !8------CALCULATE FLOW THROUGH THE RIGHT FACE.
30          IF(J.EQ.NCOL) GO TO 60
            IF(IBOUND(J+1,I,K).EQ.0) GO TO 60
            IF(IBOUND(J+1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 60
            HDIFF=HNEW(J,I,K)-HNEW(J+1,I,K)
            IF ( HDIFF.GE.-Closezero ) THEN
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J,I,LBOTM(K)-1)) - dble(BOTM(J,I,LBOTM(K)))
                  ij = Icell(J,I,K)
                  CHCH2=HDIFF*CR(J,I,K)*THICK*Sn(ij)
               ELSE
                  CHCH2=HDIFF*CR(J,I,K)
               END IF
            ELSE
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J+1,I,LBOTM(K)-1)) -&
                  &dble(BOTM(J+1,I,LBOTM(K)))
                  ij = Icell(J+1,I,K)
                  CHCH2=HDIFF*CR(J,I,K)*THICK*Sn(ij)
               ELSE
                  CHCH2=HDIFF*CR(J,I,K)
               END IF
            END IF
            IF(IBOUND(J+1,I,K).LT.0) GO TO 60
            X2=CHCH2
            XX2=X2
            IF(X2.LT.ZERO) THEN
               CHOUT=CHOUT-XX2
            ELSE
               CHIN=CHIN+XX2
            END IF
   !
   !9------CALCULATE FLOW THROUGH THE BACK FACE.
60          IF(I.EQ.1) GO TO 90
            IF (IBOUND(J,I-1,K).EQ.0) GO TO 90
            IF (IBOUND(J,I-1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 90
            HDIFF=HNEW(J,I,K)-HNEW(J,I-1,K)
            IF ( HDIFF.GE.-Closezero ) THEN
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J,I,LBOTM(K)-1)) - dble(BOTM(J,I,LBOTM(K)))
                  ij =  Icell(J,I,K)
                  CHCH3=HDIFF*CC(J,I-1,K)*THICK*Sn(ij)
               ELSE
                  CHCH3=HDIFF*CC(J,I-1,K)
               END IF
            ELSE
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J,I-1,LBOTM(K)-1)) -&
                  &dble(BOTM(J,I-1,LBOTM(K)))
                  ij =  Icell(J,I-1,K)
                  CHCH3=HDIFF*CC(J,I-1,K)*THICK*Sn(ij)
               ELSE
                  CHCH3=HDIFF*CC(J,I-1,K)
               END IF
            END IF
            IF(IBOUND(J,I-1,K).LT.0) GO TO 90
            X3=CHCH3
            XX3=X3
            IF(X3.LT.ZERO) THEN
               CHOUT=CHOUT-XX3
            ELSE
               CHIN=CHIN+XX3
            END IF
   !
   !10-----CALCULATE FLOW THROUGH THE FRONT FACE.
90          IF(I.EQ.NROW) GO TO 120
            IF(IBOUND(J,I+1,K).EQ.0) GO TO 120
            IF(IBOUND(J,I+1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 120
            HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
            IF ( HDIFF.GE.-Closezero ) THEN
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J,I,LBOTM(K)-1)) - dble(BOTM(J,I,LBOTM(K)))
                  ij = Icell(J,I,K)
                  CHCH4=HDIFF*CC(J,I,K)*THICK*Sn(ij)
               ELSE
                  CHCH4=HDIFF*CC(J,I,K)
               END IF
            ELSE
               IF ( iltyp.GT.0 ) THEN
                  THICK = dble(BOTM(J,I+1,LBOTM(K)-1)) -&
                  &dble(BOTM(J,I+1,LBOTM(K)))
                  ij = Icell(J,I+1,K)
                  CHCH4=HDIFF*CC(J,I,K)*THICK*Sn(ij)
               ELSE
                  CHCH4=HDIFF*CC(J,I,K)
               END IF
            END IF
            IF(IBOUND(J,I+1,K).LT.0) GO TO 120
            X4=CHCH4
            XX4=X4
            IF(X4.LT.ZERO) THEN
               CHOUT=CHOUT-XX4
            ELSE
               CHIN=CHIN+XX4
            END IF
   !
   !11-----CALCULATE FLOW THROUGH THE UPPER FACE.
   ! RGN need to correct the flow between layers for case when CV is a function of head.
120         IF(K.EQ.1) GO TO 150
            IF (IBOUND(J,I,K-1).EQ.0) GO TO 150
            IF (IBOUND(J,I,K-1).LT.0 .AND. ICHFLG.EQ.0) GO TO 150
            HD=HNEW(J,I,K)
   !      IF(LAYTYPUPW(K).EQ.0) GO TO 122
   !      TMP=HD
   !      TOP=BOTM(J,I,LBOTM(K)-1)
   !      IF(TMP.LT.TOP) HD=TOP
   !  122 HDIFF=HD-HNEW(J,I,K-1)
            HDIFF=HD-HNEW(J,I,K-1)
            CHCH5=HDIFF*CV(J,I,K-1)
            IF(IBOUND(J,I,K-1).LT.0) GO TO 150
            X5=CHCH5
            XX5=X5
            IF(X5.LT.ZERO) THEN
               CHOUT=CHOUT-XX5
            ELSE
               CHIN=CHIN+XX5
            END IF
   !
   !12-----CALCULATE FLOW THROUGH THE LOWER FACE.
150         IF(K.EQ.NLAY) GO TO 180
            IF(IBOUND(J,I,K+1).EQ.0) GO TO 180
            IF(IBOUND(J,I,K+1).LT.0 .AND. ICHFLG.EQ.0) GO TO 180
            HD=HNEW(J,I,K+1)
   !      IF(LAYTYPUPW(K+1).EQ.0) GO TO 152
   !      TMP=HD
   !      TOP=BOTM(J,I,LBOTM(K+1)-1)
   !      IF(TMP.LT.TOP) HD=TOP
   !  152 HDIFF=HNEW(J,I,K)-HD
            HDIFF=HNEW(J,I,K)-HD
            CHCH6=HDIFF*CV(J,I,K)
            IF(IBOUND(J,I,K+1).LT.0) GO TO 180
            X6=CHCH6
            XX6=X6
            IF(X6.LT.ZERO) THEN
               CHOUT=CHOUT-XX6
            ELSE
               CHIN=CHIN+XX6
            END IF
   !
   !13-----SUM THE FLOWS THROUGH SIX FACES OF CONSTANT HEAD CELL, AND
   !13-----STORE SUM IN BUFFER.
180         RATE=CHCH1+CHCH2+CHCH3+CHCH4+CHCH5+CHCH6
            BUFF(J,I,K)=RATE
   !
   !14-----PRINT THE FLOW FOR THE CELL IF REQUESTED.
            IF(IBD.LT.0) THEN
               IF(IBDLBL.EQ.0) WRITE(IOUT,899) TEXT,KPER,KSTP
899            FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP',I6) !gsf
               WRITE(IOUT,900) K,I,J,RATE
900            FORMAT(1X,'LAYER ',I3,'   ROW ',I5,'   COL ',I5,&
               &'   RATE ',1PG15.6)
               IBDLBL=1
            END IF
   !
   !15-----IF SAVING CELL-BY-CELL FLOW IN LIST, WRITE FLOW FOR CELL.
            IF(IBD.EQ.2) CALL UBDSVA(IUPWCB,NCOL,NROW,J,I,K,RATE,IBOUND,NLAY)
200 CONTINUE
   !
   !16-----IF SAVING CELL-BY-CELL FLOW IN 3-D ARRAY, WRITE THE ARRAY.
   IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,&
   &IUPWCB,BUFF,NCOL,NROW,NLAY,IOUT)
   !
   !17-----SAVE TOTAL CONSTANT HEAD FLOWS AND VOLUMES IN VBVL TABLE
   !17-----FOR INCLUSION IN BUDGET. PUT LABELS IN VBNM TABLE.
   CIN=CHIN
   COUT=CHOUT
   VBVL(1,MSUM)=VBVL(1,MSUM)+sngl(CIN)*DELT
   VBVL(2,MSUM)=VBVL(2,MSUM)+sngl(COUT)*DELT
   VBVL(3,MSUM)=sngl(CIN)
   VBVL(4,MSUM)=sngl(COUT)
   VBNM(MSUM)=TEXT
   MSUM=MSUM+1
   !
   !18-----RETURN.
   RETURN
END
   !
SUBROUTINE SGWF2UPWSC(SC,K,ISPST)
   !     ******************************************************************
   !     COMPUTE STORAGE CAPACITY
   !     ******************************************************************
   !
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,        ONLY:NCOL,NROW,DELR,DELC,BOTM,LBOTM,LAYCBD
   !
   DIMENSION SC(NCOL,NROW)
   !     ------------------------------------------------------------------
   !
   !1------MULTIPLY SPECIFIC STORAGE BY THICKNESS, DELR, AND DELC TO GET
   !1------CONFINED STORAGE CAPACITY.
   IF(ISPST.NE.0) THEN
      DO 80 I=1,NROW
         DO 80 J=1,NCOL
   ! RGN Made this consistent with unconfined storage by not multiplying by thickness.
   !         THICK=BOTM(J,I,LBOTM(K)-1)-BOTM(J,I,LBOTM(K))
   !         SC(J,I)=SC(J,I)*THICK*DELR(J)*DELC(I)
            SC(J,I)=SC(J,I)*DELR(J)*DELC(I)
80    CONTINUE
   ELSE
   !
   !2------MULTIPLY SPECIFIC YIELD BY DELR AND DELC TO GET UNCONFINED
   !2------STORAGEE CAPACITY(SC2).
      DO 85 I=1,NROW
         DO 85 J=1,NCOL
            SC(J,I)=SC(J,I)*DELR(J)*DELC(I)
85    CONTINUE
   END IF
   !
   RETURN
END
   !
SUBROUTINE SGWF2UPWCK(IOUT,NP,PTYP)
   !     ******************************************************************
   !     CHECK THAT JUST-DEFINED PARAMETER OF TYPE 'VK' OR 'VANI' IS USED
   !     CONSISTENTLY WITH LAYVKA ENTRIES FOR LAYERS LISTED IN CLUSTERS FOR
   !     THE PARAMETER
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GWFUPWMODULE,  ONLY:CHANI,LAYVKAUPW
   USE PARAMMODULE
   !
   CHARACTER*4 PTYP
   !     ------------------------------------------------------------------
   !
   !1------LOOP THROUGH THE CLUSTERS FOR THIS PARAMETER.
   DO 10 ICL = IPLOC(1,NP),IPLOC(2,NP)
      LAY = IPCLST(1,ICL)
      LV = LAYVKAUPW(LAY)
      IF (PTYP.EQ.'VK  ' .AND. LV.NE.0) THEN
         WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VK'
590      FORMAT(/,&
         &1X,'LAYVKA entered for layer ',i3,' is: ',i3,'; however,',&
         &' layer ',i3,' is',/,' listed in a cluster for parameter "',a,&
         &'" of type ',a,' and')
         WRITE (IOUT,600)
600      FORMAT(&
         &1X,'parameters of type VK can apply only to layers for which',&
         &/,' LAYVKA is specified as zero -- ',&
         &'STOP EXECUTION (SGWF2UPWCK)')
         CALL USTOP(' ')
      ELSEIF (PTYP.EQ.'VANI' .AND. LV.EQ.0) THEN
         WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VANI'
         WRITE (IOUT,610)
610      FORMAT(&
         &1X,'parameters of type VANI can apply only to layers for which',/,&
         &' LAYVKA is not specified as zero -- STOP EXECUTION',&
         &' (SGWF2UPWCK)')
         CALL USTOP(' ')
      ENDIF
10 CONTINUE
   !
   !2------Return.
   RETURN
END
   !
SUBROUTINE SGWF2UPW1HHARM(K)
   !     ******************************************************************
   !      COMPUTE THE CONSTANT PART OF HORIZONTAL CONDUCTANCE BASED ON THE
   !      HARMONIC AVEARGE K FOR ADJACENT CELLS.
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, iout,&
   &CC, CR
   USE GWFUPWMODULE
   !     ------------------------------------------------------------------
   !
   ZERO=0.
   TWO=2.
   !
   !1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
   !1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
   DO 100 I=1,NROW
      DO 100 J=1,NCOL
   !
   !2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
   !2------TO NEXT CELL.
         IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
            CR(J,I,K)=ZERO
            CC(J,I,K)=ZERO
         ELSE
            T1=HKUPW(J,I,K)
   !3A-----IF THIS IS NOT THE LAST COLUMN (RIGHTMOST), CALCULATE
   !3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
            IF(J.NE.NCOL) THEN
               IF(IBOUND(J+1,I,K).NE.0) THEN
                  T2=HKUPW(J+1,I,K)
                  CR(J,I,K)=TWO*T2*T1*DELC(I)/(T1*DELR(J+1)+T2*DELR(J))
               ELSE
                  CR(J,I,K)=ZERO
               END IF
            ELSE
   !3B-----IF THIS IS THE LAST COLUMN, SET BRANCH CONDUCTANCE=0.
               CR(J,I,K)=ZERO
            END IF
   !
   !3C-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
   !3C-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
            IF(I.NE.NROW) THEN
               IF(IBOUND(J,I+1,K).NE.0) THEN
                  T2=HKUPW(J,I+1,K)
                  IF(CHANI(K).LE.ZERO) THEN
                     KHANI=-CHANI(K)
                     T1=T1*HANI(J,I,KHANI)
                     T2=T2*HANI(J,I+1,KHANI)
                  ELSE
                     T1=T1*CHANI(K)
                     T2=T2*CHANI(K)
                  END IF
                  CC(J,I,K)=TWO*T2*T1*DELR(J)/(T1*DELC(I+1)+T2*DELC(I))
               ELSE
   !3D-----IF THIS IS THE LAST ROW, SET BRANCH CONDUCTANCE=0.
                  CC(J,I,K)=ZERO
               END IF
            ELSE
               CC(J,I,K)=ZERO
            END IF
         END IF
100 CONTINUE
   !
   !4------RETURN
   RETURN
END
   !
SUBROUTINE SGWF2UPW1HLOG(K)
   !     ******************************************************************
   !-----COMPUTE CONSTANT PART OF HORIZONTAL CONDUCTANCE USING LOGARITHMIC
   !-----MEAN HYDRAULIC CONDUCTIVITY -- ACTIVATED BY LAYAVG=1
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, CR, CC
   USE GWFUPWMODULE
   !     ------------------------------------------------------------------
   !
   ZERO=0.
   TWO=2.
   HALF=0.5
   FRAC1=1.005
   FRAC2=0.995
   !
   !1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
   !1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
   DO 100 I=1,NROW
      DO 100 J=1,NCOL
   !
   !2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
   !2------TO NEXT CELL.
         IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
            CR(J,I,K)=ZERO
            CC(J,I,K)=ZERO
         ELSE
   !
            T1=HKUPW(J,I,K)
   !3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
   !3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
            IF(J.NE.NCOL) THEN
               IF(IBOUND(J+1,I,K).NE.0) THEN
   !3A1----LOGARITHMIC MEAN INTERBLOCK TRANSMISSIVITY
                  T2=HKUPW(J+1,I,K)
                  RATIO=T2/T1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     T=(T2-T1)/LOG(RATIO)
                  ELSE
                     T=HALF*(T1+T2)
                  END IF
                  CR(J,I,K)=TWO*DELC(I)*T/(DELR(J+1)+DELR(J))
               ELSE
                  CR(J,I,K)=ZERO
               END IF
            ELSE
               CR(J,I,K)=ZERO
            END IF
   !
   !3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
   !3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
            IF(I.NE.NROW) THEN
               IF(IBOUND(J,I+1,K).NE.0) THEN
                  T2=HKUPW(J,I+1,K)
                  IF(CHANI(K).LE.ZERO) THEN
                     KHANI=-CHANI(K)
                     T1=T1*HANI(J,I,KHANI)
                     T2=T2*HANI(J,I+1,KHANI)
                  ELSE
                     T1=T1*CHANI(K)
                     T2=T2*CHANI(K)
                  END IF
                  RATIO=T2/T1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     T=(T2-T1)/LOG(RATIO)
                  ELSE
                     T=HALF*(T1+T2)
                  END IF
                  CC(J,I,K)=TWO*DELR(J)*T/(DELC(I+1)+DELC(I))
               ELSE
                  CC(J,I,K)=ZERO
               END IF
            ELSE
               CC(J,I,K)=ZERO
            END IF
         END IF
100 CONTINUE
   !
   !4------RETURN
   RETURN
END
   !
SUBROUTINE SGWF2UPW1HUNCNF(K)
   !     ******************************************************************
   !-----COMPUTE CONSTANT PART OF HORIZONTAL CONDUCTANCE USING ARITHMETIC
   !-----MEAN CELL THICKNESS AND LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY.
   !-----THIS IS DIFFERENT FROM SGWF2UPW1HLOG FOR CONFINED LAYERS.
   !-----ACTIVATED BY LAYAVG=2
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,CR,CC,DELR,DELC,BOTM,LBOTM
   USE GWFUPWMODULE,ONLY:HKUPW,CHANI,HANI,LAYTYPUPW
   !     ------------------------------------------------------------------
   !
   ZERO=0.
   HALF=0.5
   FRAC1=1.005
   FRAC2=0.995
   TWO=2.
   !
   !1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
   !1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
   DO 100 I=1,NROW
      DO 100 J=1,NCOL
   !
   !2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
   !2------TO NEXT CELL.
         IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
            CR(J,I,K)=ZERO
            CC(J,I,K)=ZERO
         ELSE
   !
   !3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
            HYC1=HKUPW(J,I,K)
   !3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
   !3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
            IF(J.NE.NCOL) THEN
               IF(IBOUND(J+1,I,K).NE.0) THEN
   !3A1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
                  HYC2=HKUPW(J+1,I,K)
                  RATIO=HYC2/HYC1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     HYC=(HYC2-HYC1)/LOG(RATIO)
                  ELSE
                     HYC=HALF*(HYC1+HYC2)
                  END IF
                  CR(J,I,K)=TWO*DELC(I)*HYC/(DELR(J+1)+DELR(J))
               ELSE
                  CR(J,I,K)=ZERO
               END IF
            ELSE
               CR(J,I,K)=ZERO
            END IF
   !
   !3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
   !3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
            IF(I.NE.NROW) THEN
               IF(IBOUND(J,I+1,K).NE.0) THEN
   !3B1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
                  HYC2=HKUPW(J,I+1,K)
                  IF(CHANI(K).LE.ZERO) THEN
                     KHANI=-CHANI(K)
                     HYC1=HYC1*HANI(J,I,KHANI)
                     HYC2=HYC2*HANI(J,I+1,KHANI)
                  ELSE
                     HYC1=HYC1*CHANI(K)
                     HYC2=HYC2*CHANI(K)
                  END IF
                  RATIO=HYC2/HYC1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     HYC=(HYC2-HYC1)/LOG(RATIO)
                  ELSE
                     HYC=HALF*(HYC1+HYC2)
                  END IF
                  CC(J,I,K)=TWO*DELR(J)*HYC/(DELC(I+1)+DELC(I))
               ELSE
                  CC(J,I,K)=ZERO
               END IF
            ELSE
               CC(J,I,K)=ZERO
            END IF
         END IF
100 CONTINUE
   RETURN
END SUBROUTINE
   !
   !
   !
SUBROUTINE SGWF2UPW1HHARMCON(K)
   !     ******************************************************************
   !      COMPUTE THE HORIZONTAL CONDUCTANCE FOR CONFINED CELLS BASED ON THE
   !      HARMONIC AVEARGE K FOR ADJACENT CELLS.
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, iout,&
   &CC, CR, BOTM, LBOTM
   USE GWFUPWMODULE
   !     ------------------------------------------------------------------
   !
   ZERO=0.
   TWO=2.
   !1A-----Put cell thickness into CC.
   DO 200 I=1,NROW
      DO 200 J=1,NCOL
         TTOP=BOTM(J,I,LBOTM(K)-1)
         BBOT=BOTM(J,I,LBOTM(K))
         THCK=TTOP-BBOT
         CC(J,I,K)=THCK
200 CONTINUE
   !
   !1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
   !1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
   DO 100 I=1,NROW
      DO 100 J=1,NCOL
   !
   !2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
   !2------TO NEXT CELL.
         IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
            CR(J,I,K)=ZERO
            CC(J,I,K)=ZERO
         ELSE
   !
   !3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
            T1=HKUPW(J,I,K)*CC(J,I,K)
   !3A-----IF THIS IS NOT THE LAST COLUMN (RIGHTMOST), CALCULATE
   !3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
            IF(J.NE.NCOL) THEN
               IF(IBOUND(J+1,I,K).NE.0) THEN
                  T2=HKUPW(J+1,I,K)*CC(J+1,I,K)
                  CR(J,I,K)=TWO*T2*T1*DELC(I)/(T1*DELR(J+1)+T2*DELR(J))
               ELSE
                  CR(J,I,K)=ZERO
               END IF
            ELSE
   !3B-----IF THIS IS THE LAST COLUMN, SET BRANCH CONDUCTANCE=0.
               CR(J,I,K)=ZERO
            END IF
   !
   !3C-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
   !3C-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
            IF(I.NE.NROW) THEN
               IF(IBOUND(J,I+1,K).NE.0) THEN
                  T2=HKUPW(J,I+1,K)*CC(J,I+1,K)
                  IF(CHANI(K).LE.ZERO) THEN
                     KHANI=-CHANI(K)
                     T1=T1*HANI(J,I,KHANI)
                     T2=T2*HANI(J,I+1,KHANI)
                  ELSE
                     T1=T1*CHANI(K)
                     T2=T2*CHANI(K)
                  END IF
                  CC(J,I,K)=TWO*T2*T1*DELR(J)/(T1*DELC(I+1)+T2*DELC(I))
               ELSE
   !3D-----IF THIS IS THE LAST ROW, SET BRANCH CONDUCTANCE=0.
                  CC(J,I,K)=ZERO
               END IF
            ELSE
               CC(J,I,K)=ZERO
            END IF
         END IF
100 CONTINUE
   !
   !4------RETURN
   RETURN
END
   !
SUBROUTINE SGWF2UPW1HLOGCON(K)
   !     ******************************************************************
   !-----COMPUTE THE HORIZONTAL CONDUCTANCE FOR CONFINED CELLS BASED ON THE
   !-----LOGARITHMIC  MEAN HYDRAULIC CONDUCTIVITY -- ACTIVATED BY LAYAVG=1
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, CR, CC,&
   &BOTM, LBOTM
   USE GWFUPWMODULE
   !     ------------------------------------------------------------------
   !
   ZERO=0.
   TWO=2.
   HALF=0.5
   FRAC1=1.005
   FRAC2=0.995
   !1A-----Put cell thickness into CC.
   DO 200 I=1,NROW
      DO 200 J=1,NCOL
         TTOP=BOTM(J,I,LBOTM(K)-1)
         BBOT=BOTM(J,I,LBOTM(K))
         THCK=TTOP-BBOT
         CC(J,I,K)=THCK
200 CONTINUE
   !
   !1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
   !1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
   DO 100 I=1,NROW
      DO 100 J=1,NCOL
   !
   !2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
   !2------TO NEXT CELL.
         IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
            CR(J,I,K)=ZERO
            CC(J,I,K)=ZERO
         ELSE
   !
   !3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
            T1=HKUPW(J,I,K)*CC(J,I,K)
   !3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
   !3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
            IF(J.NE.NCOL) THEN
               IF(IBOUND(J+1,I,K).NE.0) THEN
   !3A1----LOGARITHMIC MEAN INTERBLOCK TRANSMISSIVITY
                  T2=HKUPW(J+1,I,K)*CC(J+1,I,K)
                  RATIO=T2/T1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     T=(T2-T1)/LOG(RATIO)
                  ELSE
                     T=HALF*(T1+T2)
                  END IF
                  CR(J,I,K)=TWO*DELC(I)*T/(DELR(J+1)+DELR(J))
               ELSE
                  CR(J,I,K)=ZERO
               END IF
            ELSE
               CR(J,I,K)=ZERO
            END IF
   !
   !3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
   !3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
            IF(I.NE.NROW) THEN
               IF(IBOUND(J,I+1,K).NE.0) THEN
                  T2=HKUPW(J,I+1,K)*CC(J,I+1,K)
                  IF(CHANI(K).LE.ZERO) THEN
                     KHANI=-CHANI(K)
                     T1=T1*HANI(J,I,KHANI)
                     T2=T2*HANI(J,I+1,KHANI)
                  ELSE
                     T1=T1*CHANI(K)
                     T2=T2*CHANI(K)
                  END IF
                  RATIO=T2/T1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     T=(T2-T1)/LOG(RATIO)
                  ELSE
                     T=HALF*(T1+T2)
                  END IF
                  CC(J,I,K)=TWO*DELR(J)*T/(DELC(I+1)+DELC(I))
               ELSE
                  CC(J,I,K)=ZERO
               END IF
            ELSE
               CC(J,I,K)=ZERO
            END IF
         END IF
100 CONTINUE
   !
   !4------RETURN
   RETURN
END
   !
SUBROUTINE SGWF2UPW1HUNCNFCON(K)
   !     ******************************************************************
   !-----COMPUTE HORIZONTAL CONDUCTANCE FOR CONFINED CELLS USING ARITHMETIC
   !-----MEAN CELL THICKNESS AND LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY.
   !-----ACTIVATED BY LAYAVG=2
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,CR,CC,DELR,DELC,BOTM,LBOTM
   USE GWFUPWMODULE,ONLY:HKUPW,CHANI,HANI,LAYTYPUPW
   !     ------------------------------------------------------------------
   !
   ZERO=0.
   HALF=0.5
   FRAC1=1.005
   FRAC2=0.995
   !1A-----Put cell thickness into CC.
   DO 200 I=1,NROW
      DO 200 J=1,NCOL
         TTOP=BOTM(J,I,LBOTM(K)-1)
         BBOT=BOTM(J,I,LBOTM(K))
         THCK=TTOP-BBOT
         CC(J,I,K)=THCK
200 CONTINUE
   !
   !1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
   !1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
   DO 100 I=1,NROW
      DO 100 J=1,NCOL
   !
   !2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
   !2------TO NEXT CELL.
         IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
            CR(J,I,K)=ZERO
            CC(J,I,K)=ZERO
         ELSE
   !
   !3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
            HYC1=HKUPW(J,I,K)
   !3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
   !3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
            IF(J.NE.NCOL) THEN
               IF(IBOUND(J+1,I,K).NE.0) THEN
   !3A1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
                  HYC2=HKUPW(J+1,I,K)
                  RATIO=HYC2/HYC1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     HYC=(HYC2-HYC1)/LOG(RATIO)
                  ELSE
                     HYC=HALF*(HYC1+HYC2)
                  END IF
   !3A2----MULTIPLY LOGARITHMIC K BY ARITMETIC SATURATED THICKNESS.
                  CR(J,I,K)=DELC(I)*HYC*(CC(J,I,K)+CC(J+1,I,K))/&
                  &(DELR(J+1)+DELR(J))
               ELSE
                  CR(J,I,K)=ZERO
               END IF
            ELSE
               CR(J,I,K)=ZERO
            END IF
   !
   !3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
   !3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
            IF(I.NE.NROW) THEN
               IF(IBOUND(J,I+1,K).NE.0) THEN
   !3B1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
                  HYC2=HKUPW(J,I+1,K)
                  IF(CHANI(K).LE.ZERO) THEN
                     KHANI=-CHANI(K)
                     HYC1=HYC1*HANI(J,I,KHANI)
                     HYC2=HYC2*HANI(J,I+1,KHANI)
                  ELSE
                     HYC1=HYC1*CHANI(K)
                     HYC2=HYC2*CHANI(K)
                  END IF
                  RATIO=HYC2/HYC1
                  IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                     HYC=(HYC2-HYC1)/LOG(RATIO)
                  ELSE
                     HYC=HALF*(HYC1+HYC2)
                  END IF
   !3B2----MULTIPLY LOGARITHMIC K BY ARITMETIC SATURATED THICKNESS.
                  CC(J,I,K)=DELR(J)*HYC*(CC(J,I,K)+CC(J,I+1,K))/&
                  &(DELC(I+1)+DELC(I))
               ELSE
                  CC(J,I,K)=ZERO
               END IF
            ELSE
               CC(J,I,K)=ZERO
            END IF
         END IF
100 CONTINUE
   !
   !4------RETURN.
   RETURN
END
   !4------RETURN.

DOUBLE PRECISION FUNCTION SAT_THICK(Hup,Ttop,Bbot,il)
   ! RETURNS SATURATED THICKNESS OF CELL BASED ON SMOOTH FUNCTION
   USE GWFUPWMODULE
   USE GWFNWTMODULE, ONLY: Thickfact
   USE GLOBAL, ONLY: IOUT
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   !     LOCAL VARIABLES
   !     -----------------------------------------------------------------
   INTEGER ic, ir, il, iltyp, METHOD1 ! make METHOD global
   DOUBLE PRECISION hup, bbot, zero, ttop, factor, x, EPS, ACOF
   DOUBLE PRECISION cof1, cof2, s, v, factor1, factor2, Y, Z, EPSQD
   !     -----------------------------------------------------------------
   ! Calculate saturated thickenss
   zero = 0.0D0
   !      v = (ttop-bbot)
   SAT_THICK = 1.0D0
   IF ( LAYTYPUPW(il).LE.0 ) RETURN
   !-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
   EPS = Thickfact
   ACOF = 1.0 / (1.0 - EPS)
   x = (Hup-bbot)/(TTOP-BBOT)
   IF ( x.LT.1.0d-9 )  x = 1.0d-9
   IF(X.LT.EPS)THEN
      Y = ACOF *0.5/EPS * X**2
   ELSEIF(X.LT.1.0-EPS)THEN
      Y = ACOF * X + (1.0-ACOF)*0.5
   ELSEIF(X.LT.1.0)THEN
      X = 1.0 - X
      Y = ACOF *0.5/EPS * X**2
      Y = 1.0-Y
   ELSE
      Y = 1.0
   ENDIF
   factor = Y
   !
   SAT_THICK = factor
END FUNCTION SAT_THICK
   !
   !     -----------------------------------------------------------------
   !
SUBROUTINE SGWF2UPW1VCOND(K)
   !     ******************************************************************
   !     COMPUTE VERTICAL BRANCH CONDUCTANCE BETWEEN A LAYER AND THE NEXT
   !     LOWER LAYER FROM VERTICAL HYDRAULIC CONDUCTIVITY.
   !     ******************************************************************
   !
   !      SPECIFICATIONS:
   !     ------------------------------------------------------------------
   USE GLOBAL,        ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,DELR,DELC,&
   &BOTM,LBOTM,LAYCBD,IOUT,STRT,CV
   USE GWFUPWMODULE
   !
   DOUBLE PRECISION BBOT,TTOP,HHD
   !     ------------------------------------------------------------------
   !
   IF(K.EQ.NLAY) RETURN
   ZERO=0.
   HALF=0.5
   !
   !1------LOOP THROUGH ALL CELLS IN THE LAYER.
   DO 100 I=1,NROW
      DO 100 J=1,NCOL
         CV(J,I,K)=ZERO
         IF(IBOUND(J,I,K).NE.0 .AND. IBOUND(J,I,K+1).NE.0) THEN
   !
   !2------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL.
            IF(LAYVKAUPW(K).EQ.0) THEN
               HYC1=VKAUPW(J,I,K)
            ELSE
               HYC1=HKUPW(J,I,K)/VKAUPW(J,I,K)
            END IF
            IF(HYC1.GT.ZERO) THEN
   !3------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL BELOW.
               IF(LAYVKAUPW(K+1).EQ.0) THEN
                  HYC2=VKAUPW(J,I,K+1)
               ELSE
                  HYC2=(HKUPW(J,I,K+1)/VKAUPW(J,I,K+1))
               END IF
               IF(HYC2.GT.ZERO) THEN
   !
   !4------CALCULATE INVERSE LEAKANCE FOR CELL.  ICONCV FLAG PREVENTS
   !4------CV FROM BEING HEAD DEPENDENT.
                  BBOT=BOTM(J,I,LBOTM(K))
                  TTOP=BOTM(J,I,LBOTM(K)-1)
                  IF(LAYSTRT(K).NE.0) TTOP=STRT(J,I,K)
                  IF(LAYTYPUPW(K).NE.0 .AND. ICONCV.EQ.0) THEN
                     HHD=HNEW(J,I,K)
   !                  IF(HHD.LT.TTOP) TTOP=HHD   !RGN 6/23/09
                  END IF
                  BOVK1=(TTOP-BBOT)*HALF/HYC1
   !
   !5------CALCULATE INVERSE LEAKANCE FOR CELL BELOW.
                  BBOT=BOTM(J,I,LBOTM(K+1))
                  TTOP=BOTM(J,I,LBOTM(K+1)-1)
                  IF(LAYSTRT(K+1).NE.0) TTOP=STRT(J,I,K+1)
                  BBB=(TTOP-BBOT)*HALF
   !
   !5A-----IF CELL BELOW IS NOT SATURATED, DO NOT INCLUDE ITS CONDUCTANCE
   !5A-----IN THE VERTICAL CONDUCTANCE CALULATION, EXCEPT THAT THE NOCVCO
   !5A-----AND ICONCV FLAGS TURN OFF THIS CORRECTION.
                  IF(LAYTYPUPW(K+1).NE.0&
                  &.AND.NOCVCO.EQ.0 .AND. ICONCV.EQ.0) THEN
                     HHD=HNEW(J,I,K+1)
   !                  IF(HHD.LT.TTOP) BBB=ZERO   !RGN 6/23/09
                  END IF
                  BOVK2=BBB/HYC2
   !
   !6------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CONFINING BED.
                  IF(LAYCBD(K).NE.0) THEN
                     IF(VKCB(J,I,LAYCBD(K)).GT.ZERO) THEN
   !
   !7------CALCULATE INVERSE LEAKANCE FOR CONFINING BED.
                        BBB=BOTM(J,I,LBOTM(K))-BOTM(J,I,LBOTM(K)+1)
                        IF(BBB.LT.ZERO) THEN
                           WRITE(IOUT,45) K,I,J
45                         FORMAT(1X,/1X,&
                           &'Negative confining bed thickness below cell (Layer,row,col)',&
                           &I4,',',I5,',',I5)
                           WRITE(IOUT,46) BOTM(J,I,LBOTM(K)),BOTM(J,I,LBOTM(K)+1)
46                         FORMAT(1X,'Top elevation, bottom elevation:',1P,2G13.5)
                           CALL USTOP(' ')
                        END IF
                        CBBOVK=BBB/VKCB(J,I,LAYCBD(K))
                        CV(J,I,K)=DELR(J)*DELC(I)/(BOVK1+CBBOVK+BOVK2)
                     END IF
                  ELSE
                     CV(J,I,K)=DELR(J)*DELC(I)/(BOVK1+BOVK2)
                  END IF
               END IF
            END IF
         END IF
100 CONTINUE
   !
   !8------RETURN.
   RETURN
END
   !
   !     -----------------------------------------------------------------
   !     Updates saturation for previous time step.

SUBROUTINE GWF2UPW1AD(IGRID)
   USE GWFUPWMODULE, ONLY: Sn, So
   USE GWFNWTMODULE, ONLY: Numactive
   IMPLICIT NONE
   !     -----------------------------------------------------------------
   !     ARGUMENTS
   INTEGER IGRID
   !     ------------------------------------------------------------------

   !     LOCAL VARIABLES
   !     -----------------------------------------------------------------
   INTEGER ij
   !     -----------------------------------------------------------------
   !
   !------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2UPW1PNT(IGRID)
   !------Set old saturation to new saturation.
   DO ij = 1, Numactive
      So(ij) = Sn(ij)
   END DO
   RETURN
END SUBROUTINE
   !
   !
   !     SUBROUTINE GWF2UPWUPDATE. UPDATE VALUES AFTER OUTER ITERATION.
SUBROUTINE GWF2UPWUPDATE(Itest, Igrid)
   USE GLOBAL, ONLY: Ncol, Nrow, Nlay, Ibound, HNEW
   USE GWFNWTMODULE, ONLY: A, IA, Numactive, Diag, HITER
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   EXTERNAL SGWF2UPW1PNT
   DOUBLE PRECISION, EXTERNAL:: Sat_thick
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     ------------------------------------------------------------------
   INTEGER Igrid, ij, ic, ir, il, Itest, i, i1, i2
   !     -----------------------------------------------------------------
   !------SET POINTERS FOR THE CURRENT GRID.
   CALL SGWF2UPW1PNT(Igrid)
   DO ij = 1, Numactive
      il = Diag(ij, 1)
      ir = Diag(ij, 2)
      ic = Diag(ij, 3)
      IF (Itest==1) HNEW(ic,ir,il) = HITER(ic,ir,il)
      I1 = IA(ij)
      I2 = IA(ij+1)-1
      DO i = I1, I2
         A(i) = 0.0D0
      END DO
   END DO
   CALL Sn_update()
END SUBROUTINE
   !     -----------------------------------------------------------------
   !     Updates saturation for latest iteration.

SUBROUTINE Sn_update()
   USE GWFUPWMODULE
   USE GWFNWTMODULE, ONLY: Diag, Numactive
   USE GLOBAL,      ONLY: Iout,HNEW,BOTM,LBOTM,NLAY
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     -----------------------------------------------------------------
   DOUBLE PRECISION, EXTERNAL :: Sat_thick
   !     -----------------------------------------------------------------
   !     LOCAL VARIABLES
   !     -----------------------------------------------------------------
   DOUBLE PRECISION HH, TP, BT
   INTEGER ij, ic, ir, il
   !     -----------------------------------------------------------------
   DO ij = 1, Numactive
      il = Diag(ij, 1)
      ir = Diag(ij, 2)
      ic = Diag(ij, 3)
      HH = HNEW(ic,ir,il)
      TP = BOTM(ic,ir,LBOTM(il)-1)
      BT = BOTM(ic,ir,LBOTM(il))
      Sn(ij) = Sat_thick(HH, TP, BT, il)
   END DO
   RETURN
END SUBROUTINE
   !
   !     ------------------------------------------------------------------
   !
   !
   !     -----------------------------------------------------------------
   !
FUNCTION DHORIZUPW(Hup, Ttop, Bbot, il)
   ! RETURNS DERIVATIVE OF HORIZONTAL CONDUCTANCE BASED ON SMOOTH FUNCTION
   ! FUNCTION IS CALCULATED IN UPW PACKAGE IN SUBROUTINE SAT_THICK
   USE GWFNWTMODULE, ONLY: Thickfact
   USE GWFUPWMODULE, ONLY: LAYTYPUPW
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     -----------------------------------------------------------------
   DOUBLE PRECISION Hup, Ttop, Bbot
   !     -----------------------------------------------------------------
   !     LOCAL VARIABLES
   !     -----------------------------------------------------------------
   DOUBLE PRECISION factor, x, s, v, cof1, cof2, EPS, ACOF, Y
   DOUBLE PRECISION EPSQD, z
   DOUBLE PRECISION DHORIZUPW
   INTEGER il
   !     -----------------------------------------------------------------
   !-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
   DHORIZUPW = 0.0D0
   IF ( LAYTYPUPW(il).LE.0 ) RETURN
   EPS = Thickfact
   ACOF = 1.0 / (1.0 - EPS)
   x = (Hup-bbot)/(TTOP-BBOT)
   IF ( x.LT.1.0d-9 )  x = 1.0d-9
   IF(X.LT.EPS)THEN
      Y = ACOF * X / (EPS*(Ttop-Bbot))
   ELSEIF(X.LT.1.0D0-EPS)THEN
      Y = ACOF /(Ttop-Bbot)
   ELSEIF(X.LT.1.0D0)THEN
      X = 1.0 - X
      Y = - ACOF * x / (EPS * (Ttop - Bbot))
   !        Y = 1.0-Y
      Y = -Y     !2-26-16. 1 should go away for derivative
   ELSE
      Y = 0.0
   ENDIF
   factor = Y
   DHORIZUPW = factor
END FUNCTION DHORIZUPW
   !
   !
   !     ------------------------------------------------------------------
   !
SUBROUTINE GWF2UPW1DA(Igrid)
   USE GWFUPWMODULE
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     ------------------------------------------------------------------
   INTEGER Igrid
   !     ------------------------------------------------------------------
   ! Deallocate UPW data.
   DEALLOCATE(Gwfupwdat(IGRID)%Sn)
   DEALLOCATE(Gwfupwdat(IGRID)%So)
   DEALLOCATE(Gwfupwdat(IGRID)%IUPWCB)
   DEALLOCATE(Gwfupwdat(IGRID)%IWDFLG)
   DEALLOCATE(Gwfupwdat(IGRID)%IWETIT)
   DEALLOCATE(Gwfupwdat(IGRID)%IHDWET)
   DEALLOCATE(Gwfupwdat(IGRID)%IPHDRY)
   DEALLOCATE(Gwfupwdat(IGRID)%ISFAC)
   DEALLOCATE(Gwfupwdat(IGRID)%ICONCV)
   DEALLOCATE(Gwfupwdat(IGRID)%ITHFLG)
   DEALLOCATE(Gwfupwdat(IGRID)%NOCVCO)
   DEALLOCATE(Gwfupwdat(IGRID)%NOVFC)
   DEALLOCATE(Gwfupwdat(IGRID)%WETFCT)
   DEALLOCATE(Gwfupwdat(IGRID)%LAYTYPUPW)
   DEALLOCATE(Gwfupwdat(IGRID)%LAYAVG)
   DEALLOCATE(Gwfupwdat(IGRID)%CHANI)
   DEALLOCATE(Gwfupwdat(IGRID)%LAYVKAUPW)
   DEALLOCATE(Gwfupwdat(IGRID)%LAYWET)
   DEALLOCATE(Gwfupwdat(IGRID)%LAYSTRT)
   DEALLOCATE(Gwfupwdat(IGRID)%LAYFLG)
   DEALLOCATE(Gwfupwdat(IGRID)%VKAUPW)
   DEALLOCATE(Gwfupwdat(IGRID)%VKCB)
   DEALLOCATE(Gwfupwdat(IGRID)%SC1)
   DEALLOCATE(Gwfupwdat(IGRID)%SC2UPW)
   DEALLOCATE(Gwfupwdat(IGRID)%HANI)
   DEALLOCATE(Gwfupwdat(IGRID)%WETDRY)
   DEALLOCATE(Gwfupwdat(IGRID)%HKUPW)
   DEALLOCATE(Gwfupwdat(IGRID)%IBOUND2)
END SUBROUTINE GWF2UPW1DA



SUBROUTINE SGWF2UPW1PNT(Igrid)
   USE GWFUPWMODULE
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     ------------------------------------------------------------------
   INTEGER Igrid
   !     ------------------------------------------------------------------
   ! Cell property data
   Sn=>Gwfupwdat(IGRID)%Sn
   So=>Gwfupwdat(IGRID)%So
   IUPWCB=>Gwfupwdat(IGRID)%IUPWCB
   IWDFLG=>Gwfupwdat(IGRID)%IWDFLG
   IWETIT=>Gwfupwdat(IGRID)%IWETIT
   IHDWET=>Gwfupwdat(IGRID)%IHDWET
   IPHDRY=>Gwfupwdat(IGRID)%IPHDRY
   ISFAC=>Gwfupwdat(IGRID)%ISFAC
   ICONCV=>Gwfupwdat(IGRID)%ICONCV
   ITHFLG=>Gwfupwdat(IGRID)%ITHFLG
   NOCVCO=>Gwfupwdat(IGRID)%NOCVCO
   NOVFC=>Gwfupwdat(IGRID)%NOVFC
   WETFCT=>Gwfupwdat(IGRID)%WETFCT
   LAYTYPUPW=>Gwfupwdat(IGRID)%LAYTYPUPW
   LAYAVG=>Gwfupwdat(IGRID)%LAYAVG
   CHANI=>Gwfupwdat(IGRID)%CHANI
   LAYVKAUPW=>Gwfupwdat(IGRID)%LAYVKAUPW
   LAYWET=>Gwfupwdat(IGRID)%LAYWET
   LAYSTRT=>Gwfupwdat(IGRID)%LAYSTRT
   LAYFLG=>Gwfupwdat(IGRID)%LAYFLG
   VKAUPW=>Gwfupwdat(IGRID)%VKAUPW
   VKCB=>Gwfupwdat(IGRID)%VKCB
   SC1=>Gwfupwdat(IGRID)%SC1
   SC2UPW=>Gwfupwdat(IGRID)%SC2UPW
   HANI=>Gwfupwdat(IGRID)%HANI
   WETDRY=>Gwfupwdat(IGRID)%WETDRY
   HKUPW=>Gwfupwdat(IGRID)%HKUPW
   IBOUND2=>Gwfupwdat(IGRID)%IBOUND2
END SUBROUTINE SGWF2UPW1PNT
   !
SUBROUTINE SGWF2UPW1PSV(Igrid)
   USE GWFUPWMODULE
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     ------------------------------------------------------------------
   INTEGER Igrid
   !     ------------------------------------------------------------------
   ! Cell property data
   Gwfupwdat(IGRID)%Sn=>Sn
   Gwfupwdat(IGRID)%So=>So
   Gwfupwdat(IGRID)%IUPWCB=>IUPWCB
   Gwfupwdat(IGRID)%IWDFLG=>IWDFLG
   Gwfupwdat(IGRID)%IWETIT=>IWETIT
   Gwfupwdat(IGRID)%IHDWET=>IHDWET
   Gwfupwdat(IGRID)%IPHDRY=>IPHDRY
   Gwfupwdat(IGRID)%ISFAC=>ISFAC
   Gwfupwdat(IGRID)%ICONCV=>ICONCV
   Gwfupwdat(IGRID)%ITHFLG=>ITHFLG
   Gwfupwdat(IGRID)%NOCVCO=>NOCVCO
   Gwfupwdat(IGRID)%NOVFC=>NOVFC
   Gwfupwdat(IGRID)%WETFCT=>WETFCT
   Gwfupwdat(IGRID)%LAYTYPUPW=>LAYTYPUPW
   Gwfupwdat(IGRID)%LAYAVG=>LAYAVG
   Gwfupwdat(IGRID)%CHANI=>CHANI
   Gwfupwdat(IGRID)%LAYVKAUPW=>LAYVKAUPW
   Gwfupwdat(IGRID)%LAYWET=>LAYWET
   Gwfupwdat(IGRID)%LAYSTRT=>LAYSTRT
   Gwfupwdat(IGRID)%LAYFLG=>LAYFLG
   Gwfupwdat(IGRID)%VKAUPW=>VKAUPW
   Gwfupwdat(IGRID)%VKCB=>VKCB
   Gwfupwdat(IGRID)%SC1=>SC1
   Gwfupwdat(IGRID)%SC2UPW=>SC2UPW
   Gwfupwdat(IGRID)%HANI=>HANI
   Gwfupwdat(IGRID)%WETDRY=>WETDRY
   Gwfupwdat(IGRID)%HKUPW=>HKUPW
   Gwfupwdat(IGRID)%IBOUND2=>IBOUND2
   !
END SUBROUTINE SGWF2UPW1PSV
