   !  This file is provided to facilitate removing the GMG solver from MODFLOW
   !    in case a C language compiler is unavailable.
   !    To remove GMG, replace gmg.f with this file. That is, copy Nogmg.txt
   !    to gmg7.f in the src directory.  Then compile the source files with a
   !    Fortran compiler.
MODULE GMGMODULE
   INTEGER,SAVE,POINTER  ::IITER,IADAMPGMG,ISM,ISC,IOUTGMG
   INTEGER,SAVE,POINTER  ::ISIZ,IPREC,IIOUT
   INTEGER,SAVE,POINTER  ::SITER,TSITER
   INTEGER,SAVE,POINTER  ::GMGID
   INTEGER,SAVE,POINTER  ::IUNITMHC
   REAL   ,SAVE,POINTER  ::HCLOSEGMG,RCLOSEGMG,DAMPGMG
   REAL   ,SAVE,POINTER  ::DUP,DLOW,CHGLIMIT
   REAL   ,SAVE,POINTER,DIMENSION(:,:,:)::HNEWLAST
   DOUBLE PRECISION,SAVE,POINTER :: BIGHEADCHG
   DOUBLE PRECISION,SAVE,POINTER  :: RELAXGMG
END MODULE GMGMODULE
   !
SUBROUTINE GMG7AR(IN,MXITER,IGRID)
   USE GLOBAL, ONLY:IOUT
   WRITE(IOUT,1)
   WRITE(*,1)
1  FORMAT(1X,&
   &'The GMG solver is not included in this version of MODFLOW')
   CALL USTOP(' ')
   RETURN
END
SUBROUTINE GMG7AP(HNEW,RHS,CR,CC,CV,HCOF,HNOFLO,IBOUND,&
&IITER,MXITER,RCLOSE,HCLOSE,&
&KITER,KSTP,KPER,NCOL,NROW,NLAY,&
&ICNVG,SITER,TSITER,DAMP,IADAMP,&
&IOUTGMG,IOUT,GMGID,&
&IUNITMHC,DUP,DLOW,CHGLIMIT,BIGHEADCHG,HNEWLAST)
   !***********************************************************************
   IMPLICIT NONE
   REAL CR(*),CC(*),CV(*)
   REAL HNOFLO,RCLOSE,HCLOSE,DAMP
   DOUBLE PRECISION HNEW(*),RHS(*),HCOF(*)
   INTEGER IBOUND(*)
   INTEGER MXITER,IITER,KITER,KSTP,KPER,ICNVG,IOUTGMG,IOUT
   INTEGER SITER,TSITER
   INTEGER IADAMP
   INTEGER GMGID
   INTEGER NCOL,NROW,NLAY,IUNITMHC
   REAL    DUP,DLOW,CHGLIMIT
   DOUBLE PRECISION   BIGHEADCHG
   REAL    HNEWLAST(*)
   !
   RETURN
END
SUBROUTINE GMG7DA(IGRID)
   RETURN
END
SUBROUTINE GMG7PNT(IGRID)
   RETURN
END
