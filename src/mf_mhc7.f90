   !  The original Ned Banta version of MHC1 has been converted to use data
   !  allocated in the solver rather than having MHC allocate its own data.
SUBROUTINE MHC7IT(NCOL,NROW,NLAY,HNEW,HNEWLAST)
   !     ******************************************************************
   !     Store heads at beginning of iteration
   !     ******************************************************************
   IMPLICIT NONE
   INTEGER,                         INTENT(IN) :: NCOL,NROW,NLAY
   DOUBLE PRECISION, DIMENSION(NCOL,NROW,NLAY),INTENT(IN) ::HNEW
   REAL,             DIMENSION(NCOL,NROW,NLAY),INTENT(OUT)::HNEWLAST
   !  Local variables
   INTEGER :: I, J, K
   !     ------------------------------------------------------------------
   !
   DO K=1,NLAY
      DO I=1,NROW
         DO J=1,NCOL
            HNEWLAST(J,I,K)=HNEW(J,I,K)
         ENDDO
      ENDDO
   ENDDO
   !
   RETURN
END
SUBROUTINE MHC7AP(IUNITMHC,KITER,KSTP,KPER,NCOL,NROW,NLAY,IBOUND,&
&HNEW,HNEWLAST,DDAMP,BIGHEADCHG)
   !     ******************************************************************
   !     Calculate and write head changes
   !     ******************************************************************
   IMPLICIT NONE
   !  Argument-list variables
   INTEGER, INTENT(IN) :: IUNITMHC,KITER,KSTP,KPER,NCOL,NROW,NLAY
   INTEGER, DIMENSION(NCOL,NROW,NLAY),INTENT(IN) ::IBOUND
   DOUBLE PRECISION, DIMENSION(NCOL,NROW,NLAY),INTENT(IN) ::HNEW
   REAL            , DIMENSION(NCOL,NROW,NLAY),INTENT(IN) ::HNEWLAST
   DOUBLE PRECISION, INTENT(IN) :: DDAMP
   DOUBLE PRECISION, INTENT(OUT):: BIGHEADCHG
   !  Local variables
   REAL :: HEADCHG, HEADCHGMAX, HEADCHGMAXNEG, HEADCHGMAXPOS,&
   &HLAST, HLASTNEG, HLASTPOS, HTHIS, HTHISNEG, HTHISPOS
   INTEGER :: I, IMAXHEADCHG, IMAXHEADCHGNEG, IMAXHEADCHGPOS,&
   &J, JMAXHEADCHG, JMAXHEADCHGNEG, JMAXHEADCHGPOS,&
   &K, KMAXHEADCHG, KMAXHEADCHGNEG, KMAXHEADCHGPOS
   !     ------------------------------------------------------------------
   !
   HEADCHG = 0.0
   HEADCHGMAX = 0.0
   HEADCHGMAXNEG = 0.0
   HEADCHGMAXPOS = 0.0
   IMAXHEADCHG = 0
   IMAXHEADCHGNEG = 0
   IMAXHEADCHGPOS = 0
   JMAXHEADCHG = 0
   JMAXHEADCHGNEG = 0
   JMAXHEADCHGPOS = 0
   KMAXHEADCHG = 0
   KMAXHEADCHGNEG = 0
   KMAXHEADCHGPOS = 0
   HLASTNEG=0.0
   HLASTPOS=0.0
   HTHISNEG=0.0
   HTHISPOS=0.0
   !
   !  Find positive and negative max. head changes
   DO K=1,NLAY
      DO I=1,NROW
         DO J=1,NCOL
            IF (IBOUND(J,I,K) > 0) THEN
               HEADCHG = HNEW(J,I,K)-HNEWLAST(J,I,K)
               IF (HEADCHG > 0.0) THEN
                  IF (HEADCHG > HEADCHGMAXPOS) THEN
                     HEADCHGMAXPOS = HEADCHG
                     IMAXHEADCHGPOS = I
                     JMAXHEADCHGPOS = J
                     KMAXHEADCHGPOS = K
                     HLASTPOS = HNEWLAST(J,I,K)
                     HTHISPOS = HNEW(J,I,K)
                  ENDIF
               ELSE
                  IF (HEADCHG < HEADCHGMAXNEG) THEN
                     HEADCHGMAXNEG = HEADCHG
                     IMAXHEADCHGNEG = I
                     JMAXHEADCHGNEG = J
                     KMAXHEADCHGNEG = K
                     HLASTNEG = HNEWLAST(J,I,K)
                     HTHISNEG = HNEW(J,I,K)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   !
   !   Find previous and current head values at cell where
   !   absolute head change is largest
   IF (ABS(HEADCHGMAXPOS) > ABS(HEADCHGMAXNEG)) THEN
      HEADCHGMAX = HEADCHGMAXPOS
      IMAXHEADCHG = IMAXHEADCHGPOS
      JMAXHEADCHG = JMAXHEADCHGPOS
      KMAXHEADCHG = KMAXHEADCHGPOS
      HLAST = HLASTPOS
      HTHIS = HTHISPOS
   ELSE
      HEADCHGMAX = HEADCHGMAXNEG
      IMAXHEADCHG = IMAXHEADCHGNEG
      JMAXHEADCHG = JMAXHEADCHGNEG
      KMAXHEADCHG = KMAXHEADCHGNEG
      HLAST = HLASTNEG
      HTHIS = HTHISNEG
   ENDIF
   !
   !   Store max. head change for use by solver
   BIGHEADCHG=HEADCHGMAX
   !
   IF(IUNITMHC>0) THEN
      IF(KITER .EQ. 1) THEN
         WRITE(IUNITMHC,98)KPER,KSTP
98       FORMAT('"Maximum head changes for Stress Period ',I5,&
         &', Time Step ',I5,'"')
         WRITE(IUNITMHC,99)
99       FORMAT('Iteration,Max_chg,Layer,Row,Column,Damp,',&
         &'Hprev,Hcurr')
      END IF
      WRITE(IUNITMHC,100)KITER,HEADCHGMAX,KMAXHEADCHG,IMAXHEADCHG,&
      &JMAXHEADCHG,DDAMP,HLAST,HTHIS
100   FORMAT(I9,', ',G12.5,', ',I5,', ',I6,', ',I6,', ',&
      &G12.5,', ',G14.7,', ',G14.7)
   ENDIF
   !
   RETURN
END
