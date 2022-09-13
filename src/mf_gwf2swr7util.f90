   !
MODULE GSOLILUTINTERFACE
   INTERFACE
      SUBROUTINE EXPANDILUT( IWK, ISZ, JLU, ALU )
         INTEGER, INTENT(INOUT) :: IWK
         INTEGER, INTENT(IN)    :: ISZ
         INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: JLU
         DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE,&
         &INTENT(INOUT) :: ALU
      END SUBROUTINE EXPANDILUT
   END INTERFACE
END MODULE GSOLILUTINTERFACE
   !
   !-------SOLVERS
   !
   !-------LU DECOMPOSITION USING CROUT'S ALGORITHM
   !       CRS A MATRIX USED TO FILL FULL A MATRIX FOR LU DECOMPOSITION
   !       WITH BACKFILL
SUBROUTINE GSOL_FULL_LUDCAP(ICNVG,ISOLN,NOUTER,NR,NNZ,DAMP,IA,JA,&
&A,AU,X,B,X0,PINDEX,S)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(INOUT) :: ICNVG
   INTEGER, INTENT(INOUT) :: ISOLN
   INTEGER, INTENT(IN) :: NOUTER
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   DOUBLEPRECISION, INTENT(IN) :: DAMP
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)    :: A
   DOUBLEPRECISION, DIMENSION(NR,NR),INTENT(INOUT) :: AU
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)   :: X
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)      :: B
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)   :: X0
   INTEGER, DIMENSION(NR), INTENT(INOUT)           :: PINDEX
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)   :: S
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   DOUBLEPRECISION, PARAMETER :: NEARZERO = 1.0D-20
   INTEGER :: i, j, n, nn
   INTEGER :: npiv
   INTEGER :: ic0, ic1
   DOUBLEPRECISION :: d, t, tm
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   !---------INITIALIZE VARIABLES
   ICNVG = 1
   ISOLN = ISOLN + 1
   DO n = 1, NR
      X0(n) = X(n)
      X(n)  = B(n)
      S(n)  = DZERO
      DO nn = 1, NR
         AU(nn,n) = DZERO
      END DO
   END DO
   !--------FILL AU USING CRS VERSION OF A
   DO n = 1, NR
      ic0 = IA(n)
      ic1 = IA(n+1) - 1
      DO i = ic0, ic1
         AU(n,JA(i)) = A(i)
      END DO
      S(n) = MAXVAL( ABS( A(ic0:ic1) ) )
   END DO
   IF ( ANY(S.EQ.DZERO) ) THEN
      ICNVG = -1
      RETURN
   END IF
   !---------SCALING
   S = DONE / S
   !---------LOOP OVER COLUMNS
   DO n = 1, NR
      DO i = 1, n-1
         t = AU(i,n)
         DO j = 1, i-1
            t = t - AU(i,j) * AU(j,n)
         END DO
         AU(i,n) = t
      END DO
   !-----------SEARCH FOR THE LARGEST PIVOT
      tm = DZERO
      DO i = n, NR
         t = AU(i,n)
         DO j = 1, n-1
            t = t - AU(i,j) * AU(j,n)
         END DO
         AU(i,n) = t
         d = S(i) * ABS(t)
         IF ( d .GE. tm ) THEN
            npiv = i
            tm   = d
         END IF
      END DO
   !-----------INTERCHANGE ROWS
      IF ( n .NE. npiv ) THEN
         DO i = 1, NR
            t = AU(npiv,i)
            AU(npiv,i) = AU(n,i)
            AU(n,i) = t
         END DO
         S(npiv) = S(n)
      END IF
      PINDEX(n) = npiv
   !-----------TEST IF MATRIX IS SINGULAR - PROCEED WITH SMALL PIVOT
      IF ( AU(n,n).EQ.DZERO ) AU(n,n) = NEARZERO
   !-----------DIVIDE BY PIVOT ELEMENT
      IF ( n.NE.NR ) THEN
         t = DONE / AU(n,n)
         DO i = n+1, NR
            AU(i,n) = AU(i,n) * t
         END DO
      END IF
   END DO
   !---------BACKSOLVE FOR X
   nn = 0
   DO n = 1, NR
      j = PINDEX(n)
      t = X(j)
      X(j) = X(n)
      IF ( nn.NE.0 ) THEN
         DO i = nn, n-1
            t = t - AU(n,i) * X(i)
         END DO
      ELSE IF ( t .NE. DZERO ) THEN
         nn = n
      END IF
      X(n) = t
   END DO
   !---------BACKSUBSTITUTION
   DO n = NR, 1, -1
      t = X(n)
      DO i = n+1, NR
         t = t - AU(n,i) * X(i)
      END DO
      X(n) = t / AU(n,n)
   END DO
   !---------APPLY DAMPENING
   IF ( NOUTER.GT.1 ) THEN
      DO n = 1, NR
         X(n) = (DONE-DAMP) * X0(n) + DAMP * X(n)
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_FULL_LUDCAP
   !
   !-------BANDED LU DECOMPOSITION USING CROUT'S METHOD
   !       CRS JACOBIAN USED TO FILL FULL JACOBIAN
SUBROUTINE GSOL_BAND_LUDCAP(ICNVG,ISOLN,NOUTER,NR,NNZ,NBW,NHALFB,&
&DAMP,IA,JA,A,AU,AL,X,B,X0,PINDEX,S)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(INOUT) :: ICNVG
   INTEGER, INTENT(INOUT) :: ISOLN
   INTEGER, INTENT(IN) :: NOUTER
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   INTEGER, INTENT(IN) :: NBW
   INTEGER, INTENT(IN) :: NHALFB
   DOUBLEPRECISION, INTENT(IN) :: DAMP
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)        :: A
   DOUBLEPRECISION, DIMENSION(NR,NBW),INTENT(INOUT)    :: AU
   DOUBLEPRECISION, DIMENSION(NR,NHALFB),INTENT(INOUT) :: AL
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT) :: X
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)    :: B
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT) :: X0
   INTEGER, DIMENSION(NR), INTENT(INOUT)           :: PINDEX
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)   :: S
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   DOUBLEPRECISION, PARAMETER :: NEARZERO = 1.0D-20
   INTEGER :: i, j, l, n, nn
   INTEGER :: ihbp
   INTEGER :: ic0, ic1
   DOUBLEPRECISION :: t, tm
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   !---------INITIALIZE VARIABLES
   ICNVG = 1
   ISOLN = ISOLN + 1
   ihbp  = NHALFB + 1
   DO n = 1, NR
      X0(n) = X(n)
      X(n)  = B(n)
      S(n)  = DZERO
      DO nn = 1, NBW
         AU(n,nn) = DZERO
      END DO
      DO nn = 1, NHALFB
         AL(n,nn) = DZERO
      END DO
   END DO
   !--------FILL JAC%FJACU
   DO n = 1, NR
      t  = DZERO
      tm = DZERO
      ic0 = IA(n)
      ic1 = IA(n+1) - 1
      DO i = ic0, ic1
         t = A(i)
         IF ( ABS(t).GT.tm ) tm = ABS(t)
         j = ihbp + ( JA(i) - n )
         AU(n,j) = t
      END DO
      IF ( tm.EQ.DZERO ) THEN
         ICNVG = -1
         RETURN
      END IF
      S(n) = DONE / tm
   END DO
   !---------REARRANGE BANDED JACOBIAN
   l = NHALFB
   DO n = 1, NHALFB
      DO j = NHALFB+2-n, NBW
         AU(n,j-l) = AU(n,j)
      END DO
      l = l - 1
      DO j = NBW-l, NBW
         AU(n,j) = DZERO
      END DO
   END DO
   l = NHALFB
   !---------LOOP OVER COLUMNS FOR EACH ROW
   DO n = 1, NR
      t = AU(n,1) * S(n)
      i = n
      IF ( l.LT.NR ) l = l + 1
   !-----------FIND THE PIVOT ELEMENT
      DO j = n+1, l
         IF ( ABS(AU(j,1) * S(n)).GT.ABS(t) ) THEN
            t = AU(j,1) * S(n)
            i = j
         END IF
      END DO
      PINDEX(n) = i
   !-----------TEST IF MATRIX IS SINGULAR - PROCEED WITH SMALL PIVOT
      IF ( t.EQ.DZERO ) AU(n,1) = NEARZERO
   !-----------INTERCHANGE ROWS
      IF ( i.NE.n ) THEN
         DO j = 1, NBW
            t = AU(n,j)
            AU(n,j) = AU(i,j)
            AU(i,j) = t
            S(n) = S(j)
         END DO
      END IF
   !-----------ELIMINATION
      DO i = n+1, l
         t = AU(i,1) / AU(n,1)
         AL(n,i-n) = t
         DO j = 2, NBW
            AU(i,j-1) = AU(i,j) - t * AU(n,j)
         END DO
         AU(i,NBW) = DZERO
      END DO
   END DO
   !---------SOLVE BY BACKSUBSTITUTION
   CALL GSOL_BAND_BSOLVAP(NR,NBW,NHALFB,AU,AL,PINDEX,X)
   !---------APPLY DAMPENING
   IF ( NOUTER.GT.1 ) THEN
      DO n = 1, NR
         X(n) = (DONE-DAMP) * X0(n) + DAMP * X(n)
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_BAND_LUDCAP

SUBROUTINE GSOL_BAND_BSOLVAP(NR,NBW,NHALFB,AU,AL,PINDEX,X)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NBW
   INTEGER, INTENT(IN) :: NHALFB
   DOUBLEPRECISION, DIMENSION(NR,NBW), INTENT(IN)    :: AU
   DOUBLEPRECISION, DIMENSION(NR,NHALFB), INTENT(IN) :: AL
   INTEGER, DIMENSION(NR), INTENT(IN)                :: PINDEX
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)     :: X
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: i, l, n
   DOUBLEPRECISION :: t
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   !---------FORWARD SUBSTITUTION
   l = NHALFB
   DO n = 1, NR
      i = PINDEX(n)
      IF ( i.NE.n ) THEN
         t    = X(n)
         X(n) = X(i)
         X(i) = t
      END IF
      IF ( l.LT.NR ) l = l + 1
      DO i = n+1, l
         X(i) = X(i) - AL(n,i-n) * X(n)
      END DO
   END DO
   !---------BACKSUBSTITUTION
   l = 1
   DO n = NR, 1, -1
      t = X(n)
      DO i = 2, l
         t = t - AU(n,i) * X(n+i-1)
      END DO
      X(n) = t / AU(n,1)
      IF ( l.LT.NBW ) l = l + 1
   END DO
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_BAND_BSOLVAP

   !
   !-------CONJUGATE GRADIENT ITERATIVE SOLVER
SUBROUTINE GSOL_CGAP(NCORESM,NCORESV,&
&IPC,ICNVG,ISOLN,NOUTER,NINNER,&
&HCLOSE,RCLOSE,ETA,NR,NNZ,DAMP,IA,JA,A,&
&NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
&X,B,X0,DSCALE1,DSCALE2,&
&D,P,Q,Z,T)
   !        USE, INTRINSIC :: IEEE_ARITHMETIC
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN)    :: NCORESM
   INTEGER, INTENT(IN)    :: NCORESV
   INTEGER, INTENT(IN)    :: IPC
   INTEGER, INTENT(INOUT) :: ICNVG
   INTEGER, INTENT(INOUT) :: ISOLN
   INTEGER, INTENT(IN)    :: NOUTER
   INTEGER, INTENT(IN)    :: NINNER
   DOUBLEPRECISION, INTENT(IN) :: HCLOSE
   DOUBLEPRECISION, INTENT(IN) :: RCLOSE
   DOUBLEPRECISION, INTENT(INOUT) :: ETA
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   DOUBLEPRECISION, INTENT(IN) :: DAMP
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)    :: A
   INTEGER, INTENT(IN) :: NRLU
   INTEGER, INTENT(IN) :: NNZLU
   INTEGER, INTENT(IN) :: NJLU
   INTEGER, INTENT(IN) :: NJW
   INTEGER, INTENT(IN) :: NWLU
   INTEGER, DIMENSION(NJLU), INTENT(IN) :: JLU
   INTEGER, DIMENSION(NR),   INTENT(IN) :: IU
   INTEGER, DIMENSION(NJW),  INTENT(IN) :: JW
   DOUBLEPRECISION, DIMENSION(NWLU),  INTENT(INOUT) :: WLU
   DOUBLEPRECISION, DIMENSION(NNZLU), INTENT(INOUT) :: Mi
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: X
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: B
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: X0
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: DSCALE1
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: DSCALE2
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: D
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: P
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: Q
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: Z
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: T
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   INTEGER :: n
   INTEGER :: its
   INTEGER :: irel
   DOUBLEPRECISION :: sum
   DOUBLEPRECISION :: denom
   DOUBLEPRECISION :: deltax
   DOUBLEPRECISION :: alpha
   DOUBLEPRECISION :: beta
   DOUBLEPRECISION :: rho, rho0
   DOUBLEPRECISION :: rtol, rtolmin
   DOUBLEPRECISION :: dn, dn0
   DOUBLEPRECISION :: dx
   DOUBLEPRECISION :: usdx
   !     + + + FUNCTIONS + + +
   DOUBLEPRECISION :: GSOL_L2NORM
   !     + + + CODE + + +
   ICNVG = 0
   !---------DETERMINE FLOW CONVERGENCE CRITERIA
   irel  = 0
   !        CALL GSOL_RTOLMIN( NR, RCLOSE, rtol )
   CALL GSOL_SRTOLMIN( NR, DSCALE1, RCLOSE, rtol )
   IF ( ETA.GT.DZERO ) THEN
      irel    = 1
      rtolmin = rtol
   END IF
   !---------CALCULATE INITIAL RESIDUAL
   CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,X,T)
   DO n = 1, NR
      X0(n)   = X(n)
      sum     = B(n) - T(n)
      D(n)    = sum
      P(n)    = DZERO
      Q(n)    = DZERO
      Z(n)    = DZERO
   END DO
   !---------SET INITIAL VALUES
   !---------INNER ITERATION
   INNER: DO its = 1, NINNER
      ISOLN = ISOLN + 1
   !           APPLY PRECONDITIONER TO UPDATE Z
      CALL GSOL_PCA(IPC,NR,NNZ,IA,JA,&
      &NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
      &D,Z)
      rho = DOT_PRODUCT(D, Z)
      IF ( rho.eq.DZERO ) THEN
         ICNVG = -1
         EXIT INNER
      END IF
   !-----------COMPUTE DIRECTIONAL VECTORS
      IF ( its.EQ.1 ) THEN
         DO n = 1, NR
            P(n) = Z(n)
         END DO
   !-------------CALCULATE INITIAL FORCING FUNCTION FOR CURRENT INNER ITERATION
         dn  = GSOL_L2NORM( NR, D )
         dn0 = dn
         IF ( irel.GT.0 ) THEN
            rtol = ETA * dn
         END IF
      ELSE
         beta = rho / rho0
         CALL GSOL_AXPY(NCORESV, NR, Z,   beta, P, P)
      END IF
   !-----------COMPUTE ITERATES
   !           UPDATE Q WITH A AND P
      CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,P,Q)
      denom = DOT_PRODUCT(P, Q)
      IF ( denom.EQ.DZERO .OR. (ISNAN(denom).EQV. .TRUE.) )  THEN
         ICNVG = -1
         EXIT INNER
      END IF
      alpha = rho / denom
   !-----------UPDATE X
      CALL GSOL_SETX(NR, DZERO, T)
      CALL GSOL_AXPY(NCORESV, NR, T, alpha,    P, T)
      CALL GSOL_AXPY(NCORESV, NR, X,  DONE,    T, X)
   !-----------FIND MAXIMUM CHANGE IN X
      deltax = DZERO
      DO n = 1, NR
         dx     = T(n)
         usdx   = ABS( dx * DSCALE2(n) )
         IF ( usdx.GT.deltax ) THEN
            deltax = usdx
         END IF
      END DO
   !-----------UPDATE RESIDUAL
      CALL GSOL_AXPY(NCORESV, NR, D, -alpha, Q, D)
   !-----------TEST FOR OVERSOLVING
      dn = GSOL_L2NORM( NR, D )
      IF ( dn.LT.rtol .AND. deltax.LT.HCLOSE ) THEN
         ICNVG = 1
         EXIT INNER
      END IF
   !-----------SAVE CURRENT INNER ITERATES
      rho0   = rho
   END DO INNER
   !---------UPDATE ETA
   IF ( irel.GT.0 ) THEN
      CALL GSOL_ETAC(ETA,dn0,dn)
   END IF
   !---------APPLY DAMPENING
   IF ( NOUTER.GT.1 .AND. DAMP.NE.DONE ) THEN
      DO n = 1, NR
         X(n) = (DONE-DAMP) * X0(n) + DAMP * X(n)
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_CGAP
   !
   !-------BICONJUGATE GRADIENT STABILIZED ITERATIVE SOLVER
SUBROUTINE GSOL_BCGSAP(NCORESM,NCORESV,&
&IPC,ICNVG,ISOLN,NOUTER,NINNER,&
&HCLOSE,RCLOSE,ETA,NR,NNZ,DAMP,IA,JA,A,&
&NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
&X,B,X0,DSCALE1,DSCALE2,&
&D,DHAT,P,PHAT,S,SHAT,V,T)
   !        USE, INTRINSIC :: IEEE_ARITHMETIC
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN)    :: NCORESM
   INTEGER, INTENT(IN)    :: NCORESV
   INTEGER, INTENT(IN)    :: IPC
   INTEGER, INTENT(INOUT) :: ICNVG
   INTEGER, INTENT(INOUT) :: ISOLN
   INTEGER, INTENT(IN)    :: NOUTER
   INTEGER, INTENT(IN)    :: NINNER
   DOUBLEPRECISION, INTENT(IN) :: HCLOSE
   DOUBLEPRECISION, INTENT(IN) :: RCLOSE
   DOUBLEPRECISION, INTENT(INOUT) :: ETA
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   DOUBLEPRECISION, INTENT(IN) :: DAMP
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)    :: A
   INTEGER, INTENT(IN) :: NRLU
   INTEGER, INTENT(IN) :: NNZLU
   INTEGER, INTENT(IN) :: NJLU
   INTEGER, INTENT(IN) :: NJW
   INTEGER, INTENT(IN) :: NWLU
   INTEGER, DIMENSION(NJLU), INTENT(IN) :: JLU
   INTEGER, DIMENSION(NR),   INTENT(IN) :: IU
   INTEGER, DIMENSION(NJW),  INTENT(IN) :: JW
   DOUBLEPRECISION, DIMENSION(NWLU),  INTENT(INOUT) :: WLU
   DOUBLEPRECISION, DIMENSION(NNZLU), INTENT(INOUT) :: Mi
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: X
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: B
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: X0
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: DSCALE1
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: DSCALE2
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: D
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: DHAT
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: P
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: PHAT
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: S
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: SHAT
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: V
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: T
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   INTEGER :: n
   INTEGER :: its
   INTEGER :: irel
   DOUBLEPRECISION :: sum
   DOUBLEPRECISION :: denom
   DOUBLEPRECISION :: deltax
   DOUBLEPRECISION :: alpha, alpha0
   DOUBLEPRECISION :: beta
   DOUBLEPRECISION :: omega, omega0
   DOUBLEPRECISION :: rho, rho0
   DOUBLEPRECISION :: rtol, rtolmin
   DOUBLEPRECISION :: dn, dn0
   DOUBLEPRECISION :: sn
   DOUBLEPRECISION :: dx
   DOUBLEPRECISION :: usdx
   !     + + + FUNCTIONS + + +
   DOUBLEPRECISION :: GSOL_L2NORM
   !     + + + CODE + + +
   ICNVG = 0
   !---------DETERMINE FLOW CONVERGENCE CRITERIA
   irel  = 0
   !        CALL GSOL_RTOLMIN( NR, RCLOSE, rtol )
   CALL GSOL_SRTOLMIN( NR, DSCALE1, RCLOSE, rtol )
   IF ( ETA.GT.DZERO ) THEN
      irel    = 1
      rtolmin = rtol
   END IF
   !---------CALCULATE INITIAL RESIDUAL
   CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,X,T)
   DO n = 1, NR
      X0(n)   = X(n)
      sum     = B(n) - T(n)
      D(n)    = sum
      DHAT(n) = sum
   END DO
   !---------SET INITIAL VALUES
   !---------INNER ITERATION
   INNER: DO its = 1, NINNER
      ISOLN = ISOLN + 1
      rho = DOT_PRODUCT(DHAT, D)
      IF ( rho.eq.DZERO ) THEN
         ICNVG = -1
         EXIT INNER
      END IF
   !-----------COMPUTE DIRECTIONAL VECTORS
      IF (its.EQ.1) THEN
         DO n = 1, NR
            P(n) = D(n)
         END DO
   !-------------CALCULATE INITIAL FORCING FUNCTION FOR CURRENT INNER ITERATION
         dn  = GSOL_L2NORM( NR, D )
         dn0 = dn
         IF ( irel.GT.0 ) THEN
            rtol = ETA * dn
         END IF
      ELSE
         beta = (rho / rho0) * (alpha0 / omega0)
         !DO n = 1, NR
         !  P(n) = D(n) + beta * (P(n) - omega * V(n))
         !END DO
         CALL GSOL_AXPY(NCORESV, NR, P, -omega, V, P)
         CALL GSOL_AXPY(NCORESV, NR, D,   beta, P, P)
      END IF
   !-----------COMPUTE ITERATES
   !           APPLY PRECONDITIONER TO UPDATE PHAT
      CALL GSOL_PCA(IPC,NR,NNZ,IA,JA,&
      &NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
      &P,PHAT)
   !           UPDATE V WITH A AND PHAT
      CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,PHAT,V)
      denom = DOT_PRODUCT(DHAT, V)
      IF ( denom.EQ.DZERO .OR. (ISNAN(denom).EQV. .TRUE.) )  THEN
         ICNVG = -1
         EXIT INNER
      END IF
      alpha = rho / denom
      !DO n = 1, NR
      !  S(n) = D(n) - alpha * V(n)
      !END DO
      CALL GSOL_AXPY(NCORESV, NR, D, -alpha, V, S)
      sn    = GSOL_L2NORM( NR, S )
      IF ( sn.LT.rtol ) THEN
         !DO n = 1, NR
         !  dx    = alpha * PHAT(n)
         !  X(n)  = X(n) + dx
         !END DO
         CALL GSOL_AXPY(NCORESV, NR, X, alpha, PHAT, X)
         dn    = sn
         ICNVG = 1
         EXIT INNER
      END IF
   !           APPLY PRECONDITIONER TO UPDATE SHAT
      CALL GSOL_PCA(IPC,NR,NNZ,IA,JA,&
      &NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
      &S,SHAT)
   !           UPDATE T WITH A AND SHAT
      CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,SHAT,T)
      denom = DOT_PRODUCT(T, T)
      IF ( denom.EQ.DZERO .OR. (ISNAN(denom).EQV. .TRUE.) )  THEN
         ICNVG = -1
         EXIT INNER
      END IF
      omega = DOT_PRODUCT(T, S) / denom
   !-----------UPDATE X AND RESIDUAL
      CALL GSOL_SETX(NR, DZERO, D)
      CALL GSOL_AXPY(NCORESV, NR, D, omega, SHAT, D)
      CALL GSOL_AXPY(NCORESV, NR, D, alpha, PHAT, D)
      CALL GSOL_AXPY(NCORESV, NR, X,  DONE,    D, X)
      deltax = DZERO
   !          DO n = 1, NR
   !            dx     = alpha * PHAT(n) + omega * SHAT(n)
   !            X(n)   = X(n) + dx
   !!            deltax = MAX( ABS( dx ), deltax )
   !            usdx   = ABS( dx * DSCALE2(n) )
   !            IF ( usdx.GT.deltax ) THEN
   !              deltax = usdx
   !            END IF
   !            D(n)   = S(n) - omega * T(n)
   !          END DO
      DO n = 1, NR
         dx     = D(n)
         usdx   = ABS( dx * DSCALE2(n) )
         IF ( usdx.GT.deltax ) THEN
            deltax = usdx
         END IF
      END DO
      CALL GSOL_AXPY(NCORESV, NR, S, -omega, T, D)
   !-----------TEST FOR OVERSOLVING
      dn = GSOL_L2NORM( NR, D )
      IF ( dn.LT.rtol .AND. deltax.LT.HCLOSE ) THEN
         ICNVG = 1
         EXIT INNER
      END IF
   !-----------SAVE CURRENT INNER ITERATES
      rho0   = rho
      alpha0 = alpha
      omega0 = omega
   END DO INNER
   !---------UPDATE ETA
   IF ( irel.GT.0 ) THEN
   !          CALL GSOL_ETAC(ETA,dn0,dn,rtolmin)
      CALL GSOL_ETAC(ETA,dn0,dn)
   END IF
   !---------APPLY DAMPENING
   IF ( NOUTER.GT.1 .AND. DAMP.NE.DONE ) THEN
      DO n = 1, NR
         X(n) = (DONE-DAMP) * X0(n) + DAMP * X(n)
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_BCGSAP
   !
   !-------RESTARTED GENERALIZED MINIMUM RESIDUAL ITERATIVE SOLVER
SUBROUTINE GSOL_GMRMAP(NCORESM,NCORESV,&
&IPC,ICNVG,ISOLN,NOUTER,NINNER,&
&RCLOSE,ETA,NR,NNZ,DAMP,&
&IA,JA,A,&
&NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
&X,B,X0,DSCALE1,DSCALE2,&
&R,Z,D,T,CS,SN,S,Y,H,V)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN)    :: NCORESM
   INTEGER, INTENT(IN)    :: NCORESV
   INTEGER, INTENT(IN)    :: IPC
   INTEGER, INTENT(INOUT) :: ICNVG
   INTEGER, INTENT(INOUT) :: ISOLN
   INTEGER, INTENT(IN)    :: NOUTER
   INTEGER, INTENT(IN)    :: NINNER
   !        DOUBLEPRECISION, INTENT(IN) :: HCLOSE
   DOUBLEPRECISION, INTENT(IN) :: RCLOSE
   DOUBLEPRECISION, INTENT(INOUT) :: ETA
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   DOUBLEPRECISION, INTENT(IN) :: DAMP
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)    :: A
   INTEGER, INTENT(IN) :: NRLU
   INTEGER, INTENT(IN) :: NNZLU
   INTEGER, INTENT(IN) :: NJLU
   INTEGER, INTENT(IN) :: NJW
   INTEGER, INTENT(IN) :: NWLU
   INTEGER, DIMENSION(NJLU), INTENT(IN) :: JLU
   INTEGER, DIMENSION(NR),   INTENT(IN) :: IU
   INTEGER, DIMENSION(NJW),  INTENT(IN) :: JW
   DOUBLEPRECISION, DIMENSION(NWLU),  INTENT(INOUT) :: WLU
   DOUBLEPRECISION, DIMENSION(NNZLU), INTENT(INOUT) :: Mi
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: X
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: B
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)    :: X0
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: DSCALE1
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)       :: DSCALE2
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)              :: R
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)              :: Z
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)              :: D
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT)              :: T
   DOUBLEPRECISION, DIMENSION(NINNER), INTENT(INOUT)          :: CS
   DOUBLEPRECISION, DIMENSION(NINNER), INTENT(INOUT)          :: SN
   DOUBLEPRECISION, DIMENSION(NINNER+1), INTENT(INOUT)        :: S
   DOUBLEPRECISION, DIMENSION(NINNER+1), INTENT(INOUT)        :: Y
   DOUBLEPRECISION, DIMENSION(NINNER+1,NINNER), INTENT(INOUT) :: H
   DOUBLEPRECISION, DIMENSION(NR,NINNER+1), INTENT(INOUT)     :: V
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   DOUBLEPRECISION, PARAMETER :: delta = 1.0D-03
   !        INTEGER :: i, n, ii, nn
   !        INTEGER :: ots, igmr, irst
   INTEGER :: i, n, ii
   INTEGER :: igmr, irst
   INTEGER :: irel
   DOUBLEPRECISION :: sum
   DOUBLEPRECISION :: deltax
   DOUBLEPRECISION :: alpha, beta
   DOUBLEPRECISION :: ht, mu
   DOUBLEPRECISION :: dn, dn0
   DOUBLEPRECISION :: rtol, rtolmin
   DOUBLEPRECISION :: dx
   DOUBLEPRECISION :: usdx
   !     + + + FUNCTIONS + + +
   DOUBLEPRECISION :: GSOL_L2NORM
   !     + + + CODE + + +
   !---------DETERMINE FLOW CONVERGENCE CRITERIA
   irel  = 0
   !        CALL GSOL_RTOLMIN( NR, RCLOSE, rtol )
   CALL GSOL_SRTOLMIN( NR, DSCALE1, RCLOSE, rtol )
   IF ( ETA.GT.DZERO ) THEN
      irel    = 1
      rtolmin = rtol
   END IF
   !---------SET INITIAL VALUES
   ICNVG = 0
   !---------CALCULATE INITIAL RESIDUAL
   CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,X,R)
   DO n = 1, NR
      X0(n)   = X(n)
      sum     = B(n) - R(n)
      R(n)    = sum
   END DO
   !
   !---------RESTART LOOP FOR LINEAR GMRES
   irst = 1
   GMRST: DO
   !-----------APPLY PRECONDITIONER TO UPDATE z
      CALL GSOL_PCA(IPC,NR,NNZ,IA,JA,&
      &NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
      &R,Z)
   !-----------INITIALIZE GMRES VARIABLES PRIOR TO EACH RESTART (igmr) ITERATION
   !          d               = -z
      D               = Z
      H               = DZERO
      V               = DZERO
      CS              = DZERO
      SN              = DZERO
      beta            = GSOL_L2NORM( NR, Z )
      S(1)            = beta
      S(2:NINNER+1)   = DZERO
      V(1:NR,1)       = D / beta
      GMITER: DO igmr = 1, NINNER
         ISOLN = ISOLN + 1
         ii    = igmr
   !-------------CALCULATE INITIAL FORCING FUNCTION FOR CURRENT INNER ITERATION
         IF ( igmr.EQ.1 ) THEN
            dn   = GSOL_L2NORM( NR, D )
            dn0  = dn
            IF ( irel.GT.0 ) THEN
               rtol = ETA * dn
            END IF
         END IF
   !             UPDATE t WITH A AND v(:,igmr)
         CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,V(:,igmr),T)
   !             APPLY PRECONDITIONER TO UPDATE v(:,igmr+1)
         CALL GSOL_PCA(IPC,NR,NNZ,IA,JA,&
         &NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
         &T,V(:,igmr+1))
         alpha = GSOL_L2NORM( NR, V(:,igmr+1) )
         DO i = 1, igmr
            H(i,igmr)   = DOT_PRODUCT( V(:,igmr+1), V(:,i) )
            V(:,igmr+1) = V(:,igmr+1) - H(i,igmr) * V(:,i)
         END DO
         H(igmr+1,igmr) = GSOL_L2NORM( NR, V(:,igmr+1) )
         IF ( (alpha + delta * H(igmr+1,igmr) ).EQ.alpha ) THEN
            DO i = 1, igmr
               ht          = DOT_PRODUCT( V(:,igmr+1), V(:,i) )
               H(i,igmr)   = H(i,igmr) + ht
               v(:,igmr+1) = v(:,igmr+1) - ht * v(:,i)
            END DO
            H(igmr+1,igmr) = GSOL_L2NORM( NR, V(:,igmr+1) )
         END IF
         IF ( H(igmr+1,igmr).NE.DZERO ) THEN
            V(:,igmr+1) = V(:,igmr+1) / H(igmr+1,igmr)
         END IF
         IF ( igmr.GT.1 ) THEN
            Y(1:igmr+1) = H(1:igmr+1,igmr)
            DO i = 1, igmr-1
               CALL GSOL_MGVNS( i, CS(i), SN(i), Y(1:igmr+1) )
            END DO
            H(1:igmr+1,igmr) = Y(1:igmr+1)
         END IF
         mu             = SQRT(H(igmr,igmr)**2 + H(igmr+1,igmr)**2)
         CS(igmr)       = H(igmr,igmr) / mu
         SN(igmr)       = -H(igmr+1,igmr) / mu
         H(igmr,igmr)   = CS(igmr) * H(igmr,igmr) -&
         &SN(igmr) * H(igmr+1,igmr)
         H(igmr+1,igmr) = DZERO
         CALL GSOL_MGVNS( igmr, CS(igmr), SN(igmr), S(1:igmr+1) )
         dn  = ABS( S(igmr+1) )
   !-------------TEST FOR OVERSOLVING
         IF ( dn.LT.rtol .OR. igmr.EQ.NINNER ) THEN
            ii = ii - 1
            Y(ii+1) = S(ii+1) / H(ii+1,ii+1)
            DO i = ii, 1, -1
               ht   = DOT_PRODUCT( H(i,i+1:ii+1), Y(i+1:ii+1) )
               Y(i) = (S(i) - ht) / H(i,i)
            END DO
   !---------------UPDATE X
            deltax = DZERO
            DO n = 1, NR
               dx   = DOT_PRODUCT( V(n,1:ii+1), Y(1:ii+1) )
               X(n) = X(n) + dx
   !                deltax = MAX( ABS( dx ), deltax )
               usdx   = ABS( dx * DSCALE2(n) )
               IF ( usdx.GT.deltax ) THEN
                  deltax = usdx
               END IF
            END DO
   !---------------TERMINATE RESTART
   !              IF ( dn.LT.rtol .AND. deltax.LT.HCLOSE ) THEN
            IF ( dn.LT.rtol ) THEN
               ICNVG = 1
               EXIT GMRST
            END IF
            IF ( irst.EQ.NINNER ) THEN
               ICNVG = 0
               EXIT GMRST
            END IF
         END IF
      END DO GMITER
   !-----------UPDATE RESIDUAL
   !           UPDATE t WITH A AND X
      CALL GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,X,T)
      DO n = 1, NR
         sum     = B(n) - T(n)
         R(n)    = sum
      END DO
   !-----------TERMINATE RESTARTS IF NUMBER OF RESTARTS GREATER THAN NINNER
      irst = irst + 1
   END DO GMRST
   !---------UPDATE ETA
   IF ( irel.GT.0 ) THEN
   !          CALL GSOL_ETAC(ETA,dn0,dn,rtolmin)
      CALL GSOL_ETAC(ETA,dn0,dn)
   END IF
   !---------APPLY DAMPENING
   IF ( NOUTER.GT.1 ) THEN
      DO n = 1, NR
         X(n) = (DONE-DAMP) * X0(n) + DAMP * X(n)
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_GMRMAP

SUBROUTINE GSOL_MGVNS(N, Cs, Sn, G)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN)                            :: N
   DOUBLEPRECISION, INTENT(IN)                    :: Cs
   DOUBLEPRECISION, INTENT(IN)                    :: Sn
   DOUBLEPRECISION, DIMENSION(N+1), INTENT(INOUT) :: G
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION :: g1, g2
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   g1     = Cs * G(N) - Sn * G(N+1)
   g2     = Sn * G(N) + Cs * G(N+1)
   G(N)   = g1
   G(N+1) = g2
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_MGVNS

   !
   !-------ROUTINE TO SCALE THE COEFFICIENT MATRIX (A),
   !       THE FORCING VECTOR (B), AND THE ESTIMATE OF THE SOLUTION VECTOR (X)
SUBROUTINE GSOL_SCALE(IOPT,ISCL,NR,NNZ,IA,JA,A,X,B,&
&DSCALE,DSCALE2)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: IOPT
   INTEGER, INTENT(IN) :: ISCL
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   DOUBLEPRECISION, DIMENSION(NNZ), INTENT(INOUT) :: A
   DOUBLEPRECISION, DIMENSION(NR),  INTENT(INOUT) :: X
   DOUBLEPRECISION, DIMENSION(NR),  INTENT(INOUT) :: B
   DOUBLEPRECISION, DIMENSION(NR),  INTENT(INOUT) :: DSCALE
   DOUBLEPRECISION, DIMENSION(NR),  INTENT(INOUT) :: DSCALE2
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: i, n
   INTEGER :: id, jc
   INTEGER :: i0, i1
   DOUBLEPRECISION :: v, c1, c2
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   !
   !---------SCALE SCALE A, X, AND B
   IF ( IOPT.EQ.0 ) THEN
   !-----------SYMMETRIC SCALING
      SELECT CASE ( ISCL )
       CASE ( 1 )
         DO n = 1, NR
            id   = IA(n)
            v    = A(id)
            c1   = 1.0D0 / SQRT( ABS( v ) )
            DSCALE(n)  = c1
            DSCALE2(n) = c1
         END DO
   !               SCALE A -- A = DSCALE(row) * A(i) * DSCALE2(col)
         DO n = 1, NR
            c1 = DSCALE(n)
            i0 = IA(n)
            i1 = IA(n+1) - 1
            DO i = i0, i1
               jc = JA(i)
               c2 = DSCALE2(jc)
               A(i) = c1 * A(i) * c2
            END DO
         END DO
   !-----------L-2 NORM SCALING
       CASE ( 2 )
   !               SCALE EACH ROW SO THAT THE L-2 NORM IS 1
         DO n = 1, NR
            c1 = 0.0D0
            i0 = IA(n)
            i1 = IA(n+1) - 1
            DO i = i0, i1
               c1 = c1 + A(i) * A(i)
            END DO
            c1 = SQRT( c1 )
            IF ( c1.EQ.0.0D0 ) THEN
               c1 = 1.0D0
            ELSE
               c1 = 1.0D0 / c1
            END IF
            DSCALE(n) = c1
   !                 INITIAL SCALING OF A -- A = DSCALE(row) * A(i)
            DO i = i0, i1
               A(i) = c1 * A(i)
            END DO
         END DO
   !               SCALE EACH COLUMN SO THAT THE L-2 NORM IS 1
         DO n = 1, NR
            DSCALE2(n) = 0.0D0
         END DO
         c2 = 0.0D0
         DO n = 1, NR
            i0 = IA(n)
            i1 = IA(n+1) - 1
            DO i = i0, i1
               jc = JA(i)
               c2 = A(i)
               DSCALE2(jc) = DSCALE2(jc) + c2 * c2
            END DO
         END DO
         DO n = 1, NR
            c2 = DSCALE2(n)
            IF ( c2.EQ.0.0D0 ) THEN
               c2 = 1.0D0
            ELSE
               c2 = 1.0D0 / SQRT( c2 )
            END IF
            DSCALE2(n) = c2
         END DO
   !               FINAL SCALING OF A -- A = DSCALE2(col) * A(i)
         DO n = 1, NR
            i0 = IA(n)
            i1 = IA(n+1) - 1
            DO i = i0, i1
               jc = JA(i)
               c2 = DSCALE2(jc)
               A(i) = c2 * A(i)
            END DO
         END DO
      END SELECT
   !-----------SCALE X AND B
      DO n = 1, NR
         c1    = DSCALE(n)
         c2    = DSCALE2(n)
         X(n)  = X(n) / c2
         B(n)  = B(n) * c1
      END DO
   !---------UNSCALE SCALE A, X, AND B
   ELSE
      DO n = 1, NR
         c1 = DSCALE(n)
         i0 = IA(n)
         i1 = IA(n+1) - 1
   !             UNSCALE A
         DO i = i0, i1
            jc = JA(i)
            c2 = DSCALE2(jc)
            A(i) = ( 1.0D0 / c1 ) * A(i) * ( 1.0D0 / c2 )
         END DO
   !             UNSCALE X AND B
         c2   = DSCALE2(n)
         X(n) = X(n) * c2
         B(n) = B(n) / c1
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_SCALE
   !
   !-------UPDATE PRECONDITIONER
   !
   !-------ROUTINE TO UPDATE THE PRECONDITIONER
SUBROUTINE GSOL_PCU(IOUT,IPC,NR,NNZ,IA,JA,A,&
&NRLU,NNZLU,NJLU,NJW,NWLU,&
&JLU,IU,JW,WLU,Mi,&
&NLEVELS,DROPTOL)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: IOUT
   INTEGER, INTENT(IN) :: IPC
   !         COEFFICIENT MATRIX
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   DOUBLEPRECISION, DIMENSION(NNZ),   INTENT(IN)    :: A
   !         PRECONDITIONER MATRIX
   INTEGER, INTENT(IN) :: NRLU
   INTEGER, INTENT(IN) :: NNZLU
   INTEGER, INTENT(IN) :: NJLU
   INTEGER, INTENT(IN) :: NJW
   INTEGER, INTENT(IN) :: NWLU
   INTEGER, DIMENSION(NJLU), INTENT(INOUT) :: JLU
   INTEGER, DIMENSION(NR),   INTENT(INOUT) :: IU
   INTEGER, DIMENSION(NJW),  INTENT(INOUT) :: JW
   DOUBLEPRECISION, DIMENSION(NWLU),  INTENT(INOUT) :: WLU
   DOUBLEPRECISION, DIMENSION(NNZLU), INTENT(INOUT) :: Mi
   INTEGER, INTENT(IN) :: NLEVELS
   DOUBLEPRECISION, INTENT(IN) :: DROPTOL
   !     + + + LOCAL DEFINITIONS + + +
   CHARACTER (LEN=72) :: ciluterr
   INTEGER :: ierr
   INTEGER :: izero
   DOUBLEPRECISION :: delta
   !     + + + FUNCTIONS + + +
   !     + + + FORMATS + + +
2000 FORMAT (/,' MATRIX IS SEVERELY NON-DIAGONALLY DOMINANT.  CHECK',&
   &' INPUT FILES.',/,' -- STOP EXECUTION (GSOL_PCU)')
   !     + + + CODE + + +
   izero = 0
   delta = 0.0D0
   SELECT CASE(IPC)
   !           NO PRE-CONDITIONER
    CASE (0)
      Mi = A
   !           JACOBI PRE-CONDITIONER
    CASE (1)
      CALL GSOL_PC_JACO(NR,NNZ,IA,A,Mi)
   !           ILU0 AND MILU0
    CASE (2,3)
      LUPC: DO
         CALL GSOL_PC_ILU0(IPC,NR,NNZ,IA,JA,IU,JW,&
         &A,Mi,delta,izero)
         IF ( izero.NE.1 ) THEN
            EXIT LUPC
         END IF
         delta = 1.5D0 * delta + 0.001
         IF ( delta.GT.0.5D0 ) THEN
            WRITE (IOUT,2000)
            CALL USTOP('MATRIX IS SEVERELY NON-DIAGONALLY DOMINANT')
         END IF
      END DO LUPC
   !           ILU WITH DUAL TRUCATION STRATEGY AND LEVEL FILL
    CASE (4)
      ierr = 0
      CALL ilut(NR,A,JA,IA,NLEVELS,DROPTOL,&
      &Mi,JLU,IU,NNZLU,WLU,JW,ierr)
      IF ( ierr.NE.0 ) THEN
         WRITE (ciluterr,'(A,1X,I10)') 'ILUT ERROR: ', ierr
         WRITE (IOUT,'(A)') ciluterr
         CALL USTOP(ciluterr)
      END IF
   !           ADDITIONAL PRECONDITIONERS - ILU, etc.
   END SELECT
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_PCU
   !
   !-------JACOBI PRECONDITIONER - INVERSE OF DIAGONAL
   !      SUBROUTINE GSOL_PC_JACO(NR,NNZ,IA,JA,A,Mi)
SUBROUTINE GSOL_PC_JACO(NR,NNZ,IA,A,Mi)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)     :: A
   DOUBLEPRECISION, DIMENSION(NR),   INTENT(INOUT)  :: Mi
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   INTEGER :: idiag
   !        INTEGER :: i
   INTEGER :: n
   !        INTEGER :: ic0, ic1
   DOUBLEPRECISION :: sum
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   DO n = 1, NR
      Mi(n) = DZERO
   END DO
   DO n = 1, NR
      idiag = IA(n)
      sum = DZERO
      IF (ABS(A(idiag)).GT.DZERO) sum = DONE / A(idiag)
      Mi(n) = sum
   END DO
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_PC_JACO

SUBROUTINE GSOL_PC_ILU0(IPC,NR,NNZ,IA,JA,IU,IW,A,Mi,DELTA,IZERO)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: IPC
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   INTEGER, DIMENSION(NR+1),  INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),   INTENT(IN) :: JA
   INTEGER, DIMENSION(NR),    INTENT(IN) :: IU
   INTEGER, DIMENSION(NR), INTENT(INOUT) :: IW !RENAMED FROM JW TO IW WITHIN THIS SUBROUTINE
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)     :: A
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(INOUT)  :: Mi
   DOUBLEPRECISION, INTENT(IN) :: DELTA
   INTEGER, INTENT(INOUT) :: IZERO
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   INTEGER :: ic0, ic1, id0, iu1
   INTEGER :: j, n
   INTEGER :: jj, jj0, jj1
   !        INTEGER :: jpos, jrow, jw
   INTEGER :: jrow, jw
   !        INTEGER, DIMENSION(NR)  :: id, iw
   DOUBLEPRECISION :: milumult
   DOUBLEPRECISION :: rs
   DOUBLEPRECISION :: tl
   DOUBLEPRECISION :: d
   DOUBLEPRECISION :: sd1
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   milumult = DZERO
   IZERO    = 0
   IF ( IPC.EQ.3 ) milumult = DONE
   !
   !---------FILL PRECONDITIONED MATRIX (Mi) WITH A
   DO n = 1, NNZ
      Mi(n) = A(n)
   END DO
   !---------INITIALIZE POINTERS TO DIAGONAL AND NON-ZERO COLUMN ENTRIES
   DO n = 1, NR
   !          id(n) = IA(n)
      IW(n) = 0
   END DO
   MAIN: DO n = 1, NR
      ic0 = IA(n)
      ic1 = IA(n+1) - 1
      iu1 = IU(n) - 1
      DO j = ic0, ic1
         iw(JA(j)) = j
      END DO
      rs   = DZERO
      LOWER: DO j = ic0 + 1, iu1
         jrow = JA(j)
   !            tl = Mi(j) * Mi( id(jrow) )
         tl = Mi(j) * Mi( IA(jrow) )
         Mi(j) = tl
         jj0 = IU(jrow)
         jj1 = IA(jrow+1) - 1
         DO jj = jj0, jj1
            jw = iw( JA(jj) )
            IF ( jw.NE.0 ) THEN
               Mi(jw) = Mi(jw) - tl * Mi(jj)
            ELSE
               rs = rs + tl * Mi(jj)
            END IF
         END DO
      END DO LOWER
   !          id0 = id(n)
      id0 = IA(n)
      d   = Mi(id0)
      tl  = ( DONE + DELTA ) * d - rs * milumult
   !-----------CALCULATE INVERSE OF JACOBIAN DIAGONAL FOR SOLUTION
      sd1 = SIGN(d,tl)
      IF ( sd1.NE.d ) THEN
         IZERO = 1
         EXIT MAIN
      END IF
      IF ( ABS(tl).EQ.DZERO ) THEN
         IZERO = 1
         EXIT MAIN
      END IF
      Mi(id0) = DONE / tl
      !IF ( ABS(tl).GT.DZERO ) THEN
      !  Mi(id0) = DONE / tl
      !ELSE
      !  Izero = 1
      !  EXIT MAIN
      !END IF
   !-----------RESET POINTER FOR IW TO ZERO
      DO j = ic0, ic1
         jj = JA(j)
         IW(jj) = 0
      END DO
   END DO MAIN
   !C---------REVERT TO A IF ZERO ON DIAGONAL ENCOUNTERED
   !        IF ( IPC.NE.3 ) THEN
   !          IF ( Izero.GT.0 ) THEN
   !            DO n = 1, NNZ
   !              Mi(n) = A(n)
   !            END DO
   !            Izero = 0
   !          END IF
   !        END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_PC_ILU0

   !
   !
   !-------APPLY PRECONDITIONER
   !
SUBROUTINE GSOL_PCA(IPC,NR,NNZ,IA,JA,&
&NRLU,NNZLU,NJLU,NJW,NWLU,JLU,IU,JW,WLU,Mi,&
&R,D)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: IPC
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ),  INTENT(IN) :: JA
   INTEGER, INTENT(IN) :: NRLU
   INTEGER, INTENT(IN) :: NNZLU
   INTEGER, INTENT(IN) :: NJLU
   INTEGER, INTENT(IN) :: NJW
   INTEGER, INTENT(IN) :: NWLU
   INTEGER, DIMENSION(NJLU), INTENT(IN) :: JLU
   INTEGER, DIMENSION(NR),   INTENT(IN) :: IU
   INTEGER, DIMENSION(NJW),  INTENT(IN) :: JW
   DOUBLEPRECISION, DIMENSION(NWLU),  INTENT(INOUT) :: WLU
   DOUBLEPRECISION, DIMENSION(NNZLU), INTENT(INOUT) :: Mi
   DOUBLEPRECISION, DIMENSION(NR),    INTENT(IN)    :: R
   DOUBLEPRECISION, DIMENSION(NR),    INTENT(INOUT) :: D
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: ic0, ic1
   INTEGER :: j, n
   DOUBLEPRECISION :: sum
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   SELECT CASE (IPC)
   !           NO PRECONDITIONER - Mi IS A - THEREFORE SET D TO R
    CASE (0)
      DO n = 1, NR
         D(n) = R(n)
      END DO
   !           JACOBI PRECONDITIONER
    CASE (1)
      DO n = 1, NR
         D(n) = R(n) * Mi(n)
      END DO
   !           INCOMPLETE LU FACTORIZATION
    CASE (2,3)
   !               FORWARD SOLVE - Mi * D = R
      DO n = 1, NR
         D(n) = R(n)
         ic0 = IA(n) + 1
         ic1 = IU(n) - 1
         sum = D(n)
         DO j = ic0, ic1
            sum = sum - Mi(j) * D(JA(j))
         END DO
         D(n) = sum
      END DO
   !               BACKWARD SOLVE - D = D / U
      DO n = NR, 1, -1
         ic0 = IA(n)
         ic1 = IA(n+1) - 1
         sum = D(n)
         DO j = IU(n), ic1
            sum = sum - Mi(j) * D(JA(j))
         END DO
   !                 COMPUTE D FOR DIAGONAL - D = D / U
         D(n) = Mi(ic0) * sum
      END DO
    CASE (4)
   !	subroutine lusol(n, y, x, alu, jlu, ju)
      CALL LUSOL(NR, R, D, Mi, JLU, IU)

   END SELECT
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_PCA

   !-------MATRIX AND VECTOR ROUTINES
   !       MATRIX VECTOR PRODUCT
SUBROUTINE GSOL_CMATVEC(NCORESM,NR,NNZ,IA,JA,A,D1,D2)
   !USE GWFSWRMODULE, ONLY: DZERO, NR, JAC
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: NCORESM
   INTEGER, INTENT(IN) :: NR
   INTEGER, INTENT(IN) :: NNZ
   INTEGER, DIMENSION(NR+1), INTENT(IN) :: IA
   INTEGER, DIMENSION(NNZ), INTENT(IN)  :: JA
   DOUBLEPRECISION, DIMENSION(NNZ),  INTENT(IN)    :: A
   DOUBLEPRECISION, DIMENSION(NR),   INTENT(IN)    :: D1
   DOUBLEPRECISION, DIMENSION(NR),   INTENT(INOUT) :: D2
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
   DOUBLEPRECISION :: sum
   INTEGER :: ic0, ic1
   INTEGER :: icol
   INTEGER :: m, n
   INTEGER :: n0, iblksize
   INTEGER :: i
   INTEGER :: istart, iend
   INTEGER :: jstart, jend
   INTEGER :: jlen
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   IF ( NCORESM.NE.1 ) THEN
      n0 = ( ( NR - 1 ) / NCORESM ) + 1
   !$    OMP  PARALLEL
   !$    OMP& SHARED(NR,A,D1,D2)
   !$    OMP& PRIVATE(iblksize,istart,iend,jstart,jend,jlen)
   !$    OMP& NUM_THREADS(NCORESM)
   !$    OMP  DO SCHEDULE(STATIC)
      DO i = 1, NCORESM
         iblksize = MIN( n0, NR - ( i - 1 ) * n0 )
         istart   = ( i - 1 ) * n0 + 1
         iend     = istart + iblksize
         jstart   = IA(istart)
         jend     = IA(iend)
         jlen     = jend - jstart + 1
         IF ( iblksize.GT.0 ) THEN
            CALL GSOL_GEMV(iblksize,NR,jlen,jstart,&
            &IA(istart),JA(jstart),&
            &A(jstart),D1,D2(istart))
         END IF
      END DO
   !$    OMP  END DO NOWAIT
   !$    OMP  END PARALLEL
   ELSE
      DO n = 1, NR
         sum = DZERO
   !             ADD DIAGONAL AND OFF-DIAGONAL TERMS
         ic0 = IA(n)
         ic1 = IA(n+1) - 1
         DO m = ic0, ic1
            icol = JA(m)
            sum = sum + A(m) * D1(icol)
         END DO
         D2(n) = sum
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_CMATVEC

SUBROUTINE GSOL_GEMV(IBLKSIZE,NIAC,JLEN,JSTART,IA,JA,A,D1,D2)
   IMPLICIT NONE
   !         + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: IBLKSIZE
   INTEGER, INTENT(IN) :: NIAC
   INTEGER, INTENT(IN) :: JLEN
   INTEGER, INTENT(IN) :: JSTART
   INTEGER, DIMENSION(IBLKSIZE+1),       INTENT(IN)    :: IA
   INTEGER, DIMENSION(JLEN),             INTENT(IN)    :: JA
   DOUBLEPRECISION, DIMENSION(JLEN),     INTENT(IN)    :: A
   DOUBLEPRECISION, DIMENSION(NIAC),     INTENT(IN)    :: D1
   DOUBLEPRECISION, DIMENSION(IBLKSIZE), INTENT(INOUT) :: D2
   !         + + + LOCAL DEFINITIONS + + +
   INTEGER :: i, j
   INTEGER :: j0, j1, jcol
   DOUBLEPRECISION, PARAMETER :: dzero = 0.0d0
   !         + + + FUNCTIONS + + +
   !         + + + CODE + + +
   DO i = 1, IBLKSIZE
      j0 = IA(i)   - JSTART + 1
      j1 = IA(i+1) - JSTART
      D2(i) = dzero
      DO j = j0, j1
         jcol = JA(j)
         D2(i) = D2(i) + A(j) * D1(jcol)
      END DO
   END DO
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_GEMV

   !       ADD PRODUCT OF CONSTANT AND VECTOR TO VECTOR
SUBROUTINE GSOL_AXPY(NCORESV, NR, V1, C, V2, R)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: NCORESV
   INTEGER, INTENT(IN) :: NR
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)    :: V1
   DOUBLEPRECISION, INTENT(IN) :: C
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN)    :: V2
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT) :: R
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: n
   INTEGER :: n0
   INTEGER :: i
   INTEGER :: iblksize, istart
   DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   IF ( NCORESV.NE.1 ) THEN
      n0 = ( ( NR - 1 ) / NCORESV ) + 1
   !$    OMP  PARALLEL
   !$    OMP& SHARED(NR,V1,C,V2,R)
   !$    OMP& PRIVATE(iblksize,istart)
   !$    OMP& NUM_THREADS(NCORESV)
   !$    OMP  DO SCHEDULE(STATIC)
      DO i = 1, NCORESV
         iblksize = MIN( n0, NR - ( i - 1 ) * n0 )
         istart   = ( i - 1 ) * n0 + 1
         !iend     = istart + iblksize
         IF ( iblksize.GT.0 ) THEN
            CALL GSOL_MCAXPY(iblksize,V1(istart),&
            &C,V2(istart),R(istart))
         END IF
      END DO
   !$    OMP  END DO NOWAIT
   !$    OMP  END PARALLEL
   ELSE
      DO n = 1, NR
         R(n) = V1(n) + C * V2(n)
      END DO
   END IF
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_AXPY

SUBROUTINE GSOL_MCAXPY(IBLKSIZE, D1, DC, D2, DR)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: IBLKSIZE
   DOUBLEPRECISION, DIMENSION(IBLKSIZE), INTENT(IN)    :: D1
   DOUBLEPRECISION, INTENT(IN) :: DC
   DOUBLEPRECISION, DIMENSION(IBLKSIZE), INTENT(IN)    :: D2
   DOUBLEPRECISION, DIMENSION(IBLKSIZE), INTENT(INOUT) :: DR
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: n
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   DO n = 1, IBLKSIZE
      DR(n) = D1(n) + DC * D2(n)
   END DO
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_MCAXPY

   !       SET VECTOR TO A CONSTANT
SUBROUTINE GSOL_SETX(NR, C, R)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: NR
   DOUBLEPRECISION, INTENT(IN) :: C
   DOUBLEPRECISION, DIMENSION(NR), INTENT(INOUT) :: R
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: n
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   DO n = 1, NR
      R(n) = C
   END DO
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_SETX

   !-------CALCULATE THE L2 NORM OF A VECTOR
DOUBLEPRECISION FUNCTION GSOL_L2NORM(NR,V) RESULT(value)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN) :: NR
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN) :: V
   !     + + + LOCAL DEFINITIONS + + +
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   value = SQRT(DOT_PRODUCT(V,V))
   !---------RETURN
   RETURN
END FUNCTION GSOL_L2NORM
   !
   !-------CALCULATE FORCING FUNCTION TO PREVENT OVERSOLVING IN
   !       KRYLOV SOLVERS
   !      SUBROUTINE GSOL_ETAC(Eta,Dn0,Dn,Dmin)
SUBROUTINE GSOL_ETAC(Eta,Dn0,Dn)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   DOUBLEPRECISION, INTENT(INOUT) :: Eta
   DOUBLEPRECISION, INTENT(IN)    :: Dn0
   DOUBLEPRECISION, INTENT(IN)    :: Dn
   !        DOUBLEPRECISION, INTENT(IN)    :: Dmin
   !     + + + LOCAL DEFINITIONS + + +
   DOUBLEPRECISION, PARAMETER :: DONE = 1.0D0
   DOUBLEPRECISION :: gamma
   DOUBLEPRECISION :: etamin
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   gamma  = 0.9D+00
   etamin = gamma * Eta**2.0D0
   Eta    = gamma * ( Dn / Dn0 )**2.0D0
   IF ( Eta.LT.etamin ) Eta = etamin
   IF ( Eta.GT.DONE   ) Eta = DONE
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_ETAC

   !
   !-------CALCULATE MINIMUM L2NORM BASED ON USER-SPECIFIED
   !       FLOW CONVERGENCE CRITERIA (RCLOSE) - USED TO CONSTRAIN
   !       Eta WHEN INEXACT NEWTON METHODS ARE BEING USED
SUBROUTINE GSOL_RTOLMIN(NR,RCLOSE,Rtol)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN)            :: NR
   DOUBLEPRECISION, INTENT(IN)    :: RCLOSE
   DOUBLEPRECISION, INTENT(INOUT) :: Rtol
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: n
   DOUBLEPRECISION :: sum
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   sum = 0.0D0
   DO n = 1, NR
      sum = sum + RCLOSE * RCLOSE
   END DO
   Rtol = SQRT( sum )
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_RTOLMIN
   !
   !-------CALCULATE SCALED MINIMUM L2NORM BASED ON USER-SPECIFIED
   !       FLOW CONVERGENCE CRITERIA (RCLOSE) - USED TO CONSTRAIN
   !       Eta WHEN INEXACT NEWTON METHODS ARE BEING USED
SUBROUTINE GSOL_SRTOLMIN(NR,DSCALE1,RCLOSE,Rtol)
   IMPLICIT NONE
   !     + + + DUMMY ARGUMENTS + + +
   INTEGER, INTENT(IN)            :: NR
   DOUBLEPRECISION, DIMENSION(NR), INTENT(IN) :: DSCALE1
   DOUBLEPRECISION, INTENT(IN)    :: RCLOSE
   DOUBLEPRECISION, INTENT(INOUT) :: Rtol
   !     + + + LOCAL DEFINITIONS + + +
   INTEGER :: n
   DOUBLEPRECISION :: t
   DOUBLEPRECISION :: sum
   !     + + + FUNCTIONS + + +
   !     + + + CODE + + +
   sum = 0.0D0
   DO n = 1, NR
      t   = RCLOSE * DSCALE1(n)
      sum = sum + ( t * t )
   END DO
   Rtol = SQRT( sum )
   !---------RETURN
   RETURN
END SUBROUTINE GSOL_SRTOLMIN
   !
   !-------END OF GSOL ROUTINES
   !-------BEGINNING OF SUBROUTINES FROM OTHER LIBRARIES
   !
   !      SUBROUTINES ILUTSIZE AND EXPANDILUT ARE DERIVATIVE OF SPARSKIT
   !      VERSION 2 ILUT SUBROUTINE - USED TO CALCULATE THE MAXIMUM SIZE OF
   !      STORAGE NEEDED FOR THE ILUT PRECONDITIONER GIVEN A USER SPECIFIED
   !      LEVEL OF FILL (lfil) FOR THE ILUT PRECONDITIONER
   !
   !
SUBROUTINE ILUTSIZE(N,NNZ,JA,IA,LFIL,IWK)
   !-----------------------------------------------------------------------
   USE GSOLILUTINTERFACE
   IMPLICIT NONE
   INTEGER, INTENT(IN)    :: N
   INTEGER, INTENT(IN)    :: NNZ
   INTEGER, INTENT(IN)    :: JA(NNZ)
   INTEGER, INTENT(IN)    :: IA(N+1)
   INTEGER, INTENT(IN)    :: LFIL
   INTEGER, INTENT(INOUT) :: IWK
   !----------------------------------------------------------------------
   !     LOCALS
   INTEGER :: ju0,k,j1,j2,j,ii,i,lenl,lend,lenu,jj,jrow,jpos,len
   INTEGER :: lennonzero
   INTEGER :: iv
   DOUBLEPRECISION :: t, tnorm
   DOUBLEPRECISION :: s, fact
   DOUBLEPRECISION :: droptol
   DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: a
   INTEGER, DIMENSION (:), ALLOCATABLE :: ju
   INTEGER, DIMENSION (:), ALLOCATABLE :: jw
   INTEGER, DIMENSION (:), ALLOCATABLE :: jlu
   DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: alu
   DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: w
   !
   !
   !     ALLOCATE LOCAL STORAGE
   ALLOCATE( a(NNZ)   )
   ALLOCATE( ju(N)    )
   ALLOCATE( jw(2*N)  )
   ALLOCATE( jlu(IWK) )
   ALLOCATE( alu(IWK) )
   ALLOCATE( w(N+1)   )
   !     INITIALIZE LOCAL STORAGE
   droptol = 0.0D-6
   DO j = 1, NNZ
      a(j) = 1.0D0
   END DO
   DO j = 1, N
      ju(j) = 0
   END DO
   DO j = 1, 2*N
      jw(j) = 0
   END DO
   DO j = 1, N+1
      w(j)  = 0.0D0
   END DO
   DO j = 1, IWK
      jlu(j) = 0
      alu(j) = 0.0D0
   END DO
   !-----------------------------------------------------------------------
   !     INITIALIZE JU0 (POINTS TO NEXT ELEMENT TO BE ADDED TO JLU)
   !     AND POINTER ARRAY.
   !-----------------------------------------------------------------------
   ju0 = N+2
   jlu(1) = ju0
   !
   !     INITIALIZE NONZERO INDICATOR ARRAY.
   !
   DO j = 1, N
      jw(N+j)  = 0
   END DO
   !-----------------------------------------------------------------------
   !     BEGINNING OF MAIN LOOP.
   !-----------------------------------------------------------------------
   EACHROWL: DO ii = 1, N
      j1 = IA(ii)
      j2 = IA(ii+1) - 1

      tnorm = 0.0D0
      DO k = j1, j2
         tnorm = tnorm + ABS( a(k) )
      END DO
      tnorm = tnorm / REAL( (j2-j1+1), 8 )
   !
   !     UNPACK L-PART AND U-PART OF ROW OF A IN ARRAYS W
   !
      lenu       = 1
      lenl       = 0
      lend       = 0
      lennonzero = 0
      w(ii)      = 0.0D0
      jw(ii)     = ii
      jw(N+ii)   = ii
   !
      DO j = j1, j2
         k = ja(j)
         t = a(j)
   !             LOWER
         IF ( k.LT.ii ) THEN
            lenl     = lenl + 1
            jw(lenl) = k
            w(lenl)  = t
            jw(N+k)  = lenl
   !             DIAGONAL
         ELSE IF ( k.EQ.ii ) THEN
            lend     = lend + 1
            w(ii)    = t
   !            UPPER
         ELSE
            lenu     = lenu + 1
            jpos     = ii + lenu - 1
            jw(jpos) = k
            w(jpos)  = t
            jw(N+k)  = jpos
         ENDIF
      END DO
      jj  = 0
      len = 0
   !
   !     ELIMINATE PREVIOUS ROWS
   !
150   jj = jj+1
      IF ( jj.GT.lenl ) GOTO 160
   !-----------------------------------------------------------------------
   !     IN ORDER TO DO THE ELIMINATION IN THE CORRECT ORDER WE MUST SELECT
   !     THE SMALLEST COLUMN INDEX AMONG JW(K), K=JJ+1, ..., LENL.
   !-----------------------------------------------------------------------
      jrow = jw(jj)
      k    = jj
   !
   !     DETERMINE SMALLEST COLUMN INDEX
   !
      DO j = jj+1, lenl
         if ( jw(j).LT.jrow ) then
            jrow = jw(j)
            k    = j
         ENDIF
      END DO
   !
      IF ( k.NE.jj ) THEN
   !     EXCHANGE IN JW
         j          = jw(jj)
         jw(jj)     = jw(k)
         jw(k)      = j
   !     EXCHANGE IN JR
         jw(N+jrow) = jj
         jw(N+j)    = k
   !     EXCHANGE IN W
         s          = w(jj)
         w(jj)      = w(k)
         w(k)       = s
      ENDIF
   !
   !     ZERO OUT ELEMENT IN ROW BY SETTING JW(N+JROW) TO ZERO.
   !
      jw(N+jrow) = 0
   !
   !     GET THE MULTIPLIER FOR ROW TO BE ELIMINATED (JROW).
   !
      fact = w(jj) * alu(jrow)
      IF ( ABS(fact).LE.droptol ) GOTO 150
   !
   !     COMBINE CURRENT ROW AND ROW JROW
   !
      COMBINEL: DO k = ju(jrow), jlu(jrow+1)-1
         s    = fact * alu(k)
         j    = jlu(k)
         jpos = jw(N+j)
         IF ( j.GE.ii ) THEN
   !
   !     DEALING WITH UPPER PART.
   !
            IF ( jpos.EQ.0 ) THEN
   !
   !     THIS IS A FILL-IN ELEMENT
   !
               lenu    = lenu + 1
   !                  IF ( lenu.GT.N ) GOTO 995
               i       = ii + lenu - 1
               jw(i)   = j
               jw(N+j) =  i
               w(i)    = -s
            ELSE
   !
   !     THIS IS NOT A FILL-IN ELEMENT
   !
               lennonzero = lennonzero + 1
               w(jpos)    = w(jpos) - s

            ENDIF
         ELSE
   !
   !     DEALING  WITH LOWER PART.
   !
            IF ( jpos.EQ.0 ) THEN
   !
   !     THIS IS A FILL-IN ELEMENT
   !
               lenl     = lenl + 1
   !                  IF ( lenl.GT.N ) GOTO 995
               jw(lenl) = j
               jw(N+j)  = lenl
               w(lenl)  = -s
            ELSE
   !
   !     THIS IS NOT A FILL-IN ELEMENT
   !
               lennonzero = lennonzero + 1
               w(jpos)    = w(jpos) - s
            ENDIF
         ENDIF
      END DO COMBINEL
   !
   !     STORE THIS PIVOT ELEMENT -- (FROM LEFT TO RIGHT -- NO DANGER OF
   !     OVERLAP WITH THE WORKING ELEMENTS IN L (PIVOTS).
   !
      len      = len+1
      w(len)   = fact
      jw(len)  = jrow
      GOTO 150

160   CONTINUE
   !
   !     RESET DOUBLE-POINTER TO ZERO (U-PART)
   !
      DO k = 1, lenu
         jw(N+jw(ii+k-1)) = 0
      END DO
   !
   !     UPDATE L-MATRIX
   !
      lenl = len
      len  = MIN0(lenl,LFIL)
   !
   !     SORT BY QUICK-SPLIT
   !
      CALL qsplit (w,jw,lenl,len)
   !
   !     STORE L-PART
   !
      DO k = 1, len
         !IF ( ju0.GT.iwk ) GOTO 996
         IF ( ju0.GT.iwk) THEN
            CALL EXPANDILUT( IWK, ju0, jlu, alu )
         END IF
         alu(ju0) = w(k)
         jlu(ju0) = jw(k)
         ju0      = ju0 + 1
      END DO
   !
   !     SAVE POINTER TO BEGINNING OF ROW II OF U
   !
      ju(ii) = ju0
   !
   !     UPDATE U-MATRIX
   !
      len = 0
      DO k = 1, lenu-1
         len        = len+1
         w(ii+len)  = w(ii+k)
         jw(ii+len) = jw(ii+k)
      ENDDO
      lenu = len + 1
      len  = MIN0(lenu,LFIL)
   !
      CALL qsplit (w(ii+1), jw(ii+1), lenu-1,len)
   !
   !     COPY
   !
      t = ABS( w(ii) )
   !         IF ( len + ju0.GT.IWK ) GOTO 997
      IF ( len + ju0.GT.IWK) THEN
         CALL EXPANDILUT( IWK, (len + ju0), jlu, alu )
      END IF
      DO k = ii+1, ii+len-1
         jlu(ju0) = jw(k)
         alu(ju0) = w(k)
         t        = t + ABS( w(k) )
         ju0      = ju0 + 1
      END DO
   !
   !     STORE INVERSE OF DIAGONAL ELEMENT OF U
   !
      IF ( w(ii).EQ.0.0D0) w(ii) = (0.0001 + droptol)*tnorm
   !
      alu(ii) = 1.0D0 / w(ii)
   !
   !     UPDATE POINTER TO BEGINNING OF NEXT ROW OF U.
   !
      jlu(ii+1) = ju0
   !-----------------------------------------------------------------------
   !     END MAIN LOOP
   !-----------------------------------------------------------------------
   END DO EACHROWL

   !     DETERMINE THE MAXIMUM NUMBER OF NON-ZERO ELEMENTS OF JLU
   iv = 0
   DO ii = 1, IWK
      if ( jlu(ii).NE.0 ) iv = iv + 1
   END DO
   !      iwk = iv * 2
   iwk = iv + N
   !
   !     CLEAN UP TEMPORARY STORAGE
   DEALLOCATE( a   )
   DEALLOCATE( ju  )
   DEALLOCATE( jw  )
   DEALLOCATE( jlu )
   DEALLOCATE( alu )
   DEALLOCATE( w   )
   !
   !     RETURN
   RETURN
   !----------------END-OF-ILUTSIZE----------------------------------------
   !-----------------------------------------------------------------------
END SUBROUTINE ILUTSIZE

SUBROUTINE EXPANDILUT( IWK, ISZ, JLU, ALU )
   IMPLICIT NONE
   !
   INTEGER, INTENT(INOUT) :: IWK
   INTEGER, INTENT(IN)    :: ISZ
   INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: JLU
   DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ALU
   !         LOCALS
   INTEGER :: i
   INTEGER, DIMENSION(:), ALLOCATABLE :: j
   DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: a
   !         CODE
   !         ALLOCATE LOCAL STORAGE
   ALLOCATE ( j(ISZ) )
   ALLOCATE ( a(ISZ) )
   !         FILL LOCAL STORAGE WITH JLU
   DO i = 1, ISZ
      IF ( i.GT.IWK ) THEN
         j(i) = 0
         a(i) = 0.0D0
      ELSE
         j(i) = JLU(i)
         a(i) = ALU(i)
      END IF
   END DO
   !         DEALLOCATE AND ALLOCATE JLU TO NEW SIZE
   DEALLOCATE( JLU )
   DEALLOCATE( ALU )
   ALLOCATE( JLU(ISZ) )
   ALLOCATE( ALU(ISZ) )
   !         FILL REDIMENSIONED JLU WITH DATA FROM LOCAL STORAGE
   DO i = 1, ISZ
      JLU(i) = j(i)
      ALU(i) = a(i)
   END DO
   !         SET IWK TO NEW SIZE OF IWK
   IWK = ISZ
   !        CLEAN UP LOCAL STORAGE
   DEALLOCATE( j )
   DEALLOCATE( a )
   !         RETURN
   RETURN
END SUBROUTINE EXPANDILUT


   !       SUBSET OF SPARSKIT VERSION 2 SOURCE CODE
   !
   !  SPARSKIT VERSION 2 SUBROUTINES INCLUDED INCLUDE:
   !
   !    1 - ilut
   !    2 - lusol
   !    3 - qsplit
   !
   !-----------------------------------------------------------------------
   !                   S P A R S K I T   V E R S I O N  2.
   !-----------------------------------------------------------------------
   !
   !Latest update : Tue Mar  8 11:01:12 CST 2005
   !
   !-----------------------------------------------------------------------
   !
   !Welcome  to SPARSKIT  VERSION  2.  SPARSKIT is  a  package of  FORTRAN
   !subroutines  for working  with  sparse matrices.  It includes  general
   !sparse  matrix  manipulation  routines  as  well as  a  few  iterative
   !solvers, see detailed description of contents below.
   !
   ! Copyright (C) 2005, the Regents of the University of Minnesota
   !
   !SPARSKIT is  free software; you  can redistribute it and/or  modify it
   !under the terms of the  GNU Lesser General Public License as published
   !by the  Free Software Foundation [version  2.1 of the  License, or any
   !later version.]
   !
   !A copy of  the licencing agreement is attached in  the file LGPL.  For
   !additional information  contact the Free Software  Foundation Inc., 59
   !Temple Place - Suite 330, Boston, MA 02111, USA or visit the web-site
   !
   ! http://www.gnu.org/copyleft/lesser.html
   !
   !
   !DISCLAIMER
   !----------
   !
   !SPARSKIT  is distributed  in  the hope  that  it will  be useful,  but
   !WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
   !MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
   !Lesser General Public License for more details.
   !
   !For more information contact saad@cs.umn.edu
   !
   !

subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
   !-----------------------------------------------------------------------
   implicit none
   integer n
   doubleprecision a(*),alu(*),w(n+1),droptol
   integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr
   !----------------------------------------------------------------------*
   !                      *** ILUT preconditioner ***                     *
   !      incomplete LU factorization with dual truncation mechanism      *
   !----------------------------------------------------------------------*
   !     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
   !----------------------------------------------------------------------*
   ! PARAMETERS
   !-----------
   !
   ! on entry:
   !==========
   ! n       = integer. The row dimension of the matrix A. The matrix
   !
   ! a,ja,ia = matrix stored in Compressed Sparse Row format.
   !
   ! lfil    = integer. The fill-in parameter. Each row of L and each row
   !           of U will have a maximum of lfil elements (excluding the
   !           diagonal element). lfil must be .ge. 0.
   !           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
   !           EARLIER VERSIONS.
   !
   ! droptol = doubleprecision. Sets the threshold for dropping small terms in the
   !           factorization. See below for details on dropping strategy.
   !
   !
   ! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
   !           are not big enough to store the ILU factorizations, ilut
   !           will stop with an error message.
   !
   ! On return:
   !===========
   !
   ! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
   !           the L and U factors together. The diagonal (stored in
   !           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
   !           contains the i-th row of L (excluding the diagonal entry=1)
   !           followed by the i-th row of U.
   !
   ! ju      = integer array of length n containing the pointers to
   !           the beginning of each row of U in the matrix alu,jlu.
   !
   ! ierr    = integer. Error message with the following meaning.
   !           ierr  = 0    --> successful return.
   !           ierr .gt. 0  --> zero pivot encountered at step number ierr.
   !           ierr  = -1   --> Error. input matrix may be wrong.
   !                            (The elimination process has generated a
   !                            row in L or U whose length is .gt.  n.)
   !           ierr  = -2   --> The matrix L overflows the array al.
   !           ierr  = -3   --> The matrix U overflows the array alu.
   !           ierr  = -4   --> Illegal value for lfil.
   !           ierr  = -5   --> zero row encountered.
   !
   ! work arrays:
   !=============
   ! jw      = integer work array of length 2*n.
   ! w       = real work array of length n+1.
   !
   !----------------------------------------------------------------------
   ! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
   ! jw(n+1:2n)  stores nonzero indicators
   !
   ! Notes:
   ! ------
   ! The diagonal elements of the input matrix must be  nonzero (at least
   ! 'structurally').
   !
   !----------------------------------------------------------------------*
   !---- Dual drop strategy works as follows.                             *
   !                                                                      *
   !     1) Theresholding in L and U as set by droptol. Any element whose *
   !        magnitude is less than some tolerance (relative to the abs    *
   !        value of diagonal element in u) is dropped.                   *
   !                                                                      *
   !     2) Keeping only the largest lfil elements in the i-th row of L   *
   !        and the largest lfil elements in the i-th row of U (excluding *
   !        diagonal elements).                                           *
   !                                                                      *
   ! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
   ! keeping  the largest  elements in  each row  of L  and U.   Taking   *
   ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
   ! (however, fill-in is then mpredictible).                             *
   !----------------------------------------------------------------------*
   !     locals
   integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len
   doubleprecision tnorm, t, abs, s, fact
   if (lfil .lt. 0) goto 998
   !-----------------------------------------------------------------------
   !     initialize ju0 (points to next element to be added to alu,jlu)
   !     and pointer array.
   !-----------------------------------------------------------------------
   ju0 = n+2
   jlu(1) = ju0
   !
   !     initialize nonzero indicator array.
   !
   do 1 j=1,n
      jw(n+j)  = 0
1  continue
   !-----------------------------------------------------------------------
   !     beginning of main loop.
   !-----------------------------------------------------------------------
   do 500 ii = 1, n
      j1 = ia(ii)
      j2 = ia(ii+1) - 1
      tnorm = 0.0d0
      do 501 k=j1,j2
         tnorm = tnorm+abs(a(k))
501   continue
      if (tnorm .eq. 0.0) goto 999
      tnorm = tnorm/real(j2-j1+1)
   !
   !     unpack L-part and U-part of row of A in arrays w
   !
      lenu = 1
      lenl = 0
      jw(ii) = ii
      w(ii) = 0.0
      jw(n+ii) = ii
   !
      do 170  j = j1, j2
         k = ja(j)
         t = a(j)
         if (k .lt. ii) then
            lenl = lenl+1
            jw(lenl) = k
            w(lenl) = t
            jw(n+k) = lenl
         else if (k .eq. ii) then
            w(ii) = t
         else
            lenu = lenu+1
            jpos = ii+lenu-1
            jw(jpos) = k
            w(jpos) = t
            jw(n+k) = jpos
         endif
170   continue
      jj = 0
      len = 0
   !
   !     eliminate previous rows
   !
150   jj = jj+1
      if (jj .gt. lenl) goto 160
   !-----------------------------------------------------------------------
   !     in order to do the elimination in the correct order we must select
   !     the smallest column index among jw(k), k=jj+1, ..., lenl.
   !-----------------------------------------------------------------------
      jrow = jw(jj)
      k = jj
   !
   !     determine smallest column index
   !
      do 151 j=jj+1,lenl
         if (jw(j) .lt. jrow) then
            jrow = jw(j)
            k = j
         endif
151   continue
   !
      if (k .ne. jj) then
   !     exchange in jw
         j = jw(jj)
         jw(jj) = jw(k)
         jw(k) = j
   !     exchange in jr
         jw(n+jrow) = jj
         jw(n+j) = k
   !     exchange in w
         s = w(jj)
         w(jj) = w(k)
         w(k) = s
      endif
   !
   !     zero out element in row by setting jw(n+jrow) to zero.
   !
      jw(n+jrow) = 0
   !
   !     get the multiplier for row to be eliminated (jrow).
   !
      fact = w(jj)*alu(jrow)
      if (abs(fact) .le. droptol) goto 150
   !
   !     combine current row and row jrow
   !
      do 203 k = ju(jrow), jlu(jrow+1)-1
         s = fact*alu(k)
         j = jlu(k)
         jpos = jw(n+j)
         if (j .ge. ii) then
   !
   !     dealing with upper part.
   !
            if (jpos .eq. 0) then
   !
   !     this is a fill-in element
   !
               lenu = lenu+1
               if (lenu .gt. n) goto 995
               i = ii+lenu-1
               jw(i) = j
               jw(n+j) = i
               w(i) = - s
            else
   !
   !     this is not a fill-in element
   !
               w(jpos) = w(jpos) - s

            endif
         else
   !
   !     dealing  with lower part.
   !
            if (jpos .eq. 0) then
   !
   !     this is a fill-in element
   !
               lenl = lenl+1
               if (lenl .gt. n) goto 995
               jw(lenl) = j
               jw(n+j) = lenl
               w(lenl) = - s
            else
   !
   !     this is not a fill-in element
   !
               w(jpos) = w(jpos) - s
            endif
         endif
203   continue
   !
   !     store this pivot element -- (from left to right -- no danger of
   !     overlap with the working elements in L (pivots).
   !
      len = len+1
      w(len) = fact
      jw(len)  = jrow
      goto 150
160   continue
   !
   !     reset double-pointer to zero (U-part)
   !
      do 308 k=1, lenu
         jw(n+jw(ii+k-1)) = 0
308   continue
   !
   !     update L-matrix
   !
      lenl = len
      len = min0(lenl,lfil)
   !
   !     sort by quick-split
   !
      call qsplit (w,jw,lenl,len)
   !
   !     store L-part
   !
      do 204 k=1, len
   !            if (ju0 .gt. iwk) goto 996
         if (ju0 .gt. iwk) then
            write (*,'(//1x,2i10)') ju0, iwk
            goto 996
         end if
         alu(ju0) =  w(k)
         jlu(ju0) =  jw(k)
         ju0 = ju0+1
204   continue
   !
   !     save pointer to beginning of row ii of U
   !
      ju(ii) = ju0
   !
   !     update U-matrix -- first apply dropping strategy
   !
      len = 0
      do k=1, lenu-1
         if (abs(w(ii+k)) .gt. droptol*tnorm) then
            len = len+1
            w(ii+len) = w(ii+k)
            jw(ii+len) = jw(ii+k)
         endif
      enddo
      lenu = len+1
      len = min0(lenu,lfil)
   !
      call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
   !
   !     copy
   !
      t = abs(w(ii))
   !         if (len + ju0 .gt. iwk) goto 997
      if (len + ju0 .gt. iwk) then
         write (*,'(//1x,2i10)') (len + ju0), iwk
         goto 997
      end if
      do 302 k=ii+1,ii+len-1
         jlu(ju0) = jw(k)
         alu(ju0) = w(k)
         t = t + abs(w(k) )
         ju0 = ju0+1
302   continue
   !
   !     store inverse of diagonal element of u
   !
      if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
   !
      alu(ii) = 1.0d0/ w(ii)
   !
   !     update pointer to beginning of next row of U.
   !
      jlu(ii+1) = ju0
   !-----------------------------------------------------------------------
   !     end main loop
   !-----------------------------------------------------------------------
500 continue
   ierr = 0
   return
   !
   !     incomprehensible error. Matrix must be wrong.
   !
995 ierr = -1
   return
   !
   !     insufficient storage in L.
   !
996 ierr = -2
   return
   !
   !     insufficient storage in U.
   !
997 ierr = -3
   return
   !
   !     illegal lfil entered.
   !
998 ierr = -4
   return
   !
   !     zero row encountered
   !
999 ierr = -5
   return
   !----------------end-of-ilut--------------------------------------------
   !-----------------------------------------------------------------------
end
   !-----------------------------------------------------------------------
subroutine lusol(n, y, x, alu, jlu, ju)
   doubleprecision x(n), y(n), alu(*)
   integer n, jlu(*), ju(*)
   !-----------------------------------------------------------------------
   !
   ! This routine solves the system (LU) x = y,
   ! given an LU decomposition of a matrix stored in (alu, jlu, ju)
   ! modified sparse row format
   !
   !-----------------------------------------------------------------------
   ! on entry:
   ! n   = dimension of system
   ! y   = the right-hand-side vector
   ! alu, jlu, ju
   !     = the LU matrix as provided from the ILU routines.
   !
   ! on return
   ! x   = solution of LU x = y.
   !-----------------------------------------------------------------------
   !
   ! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)
   !       will solve the system with rhs x and overwrite the result on x .
   !
   !-----------------------------------------------------------------------
   ! local variables
   !
   integer i,k
   !
   ! forward solve
   !
   do 40 i = 1, n
      x(i) = y(i)
      do 41 k=jlu(i),ju(i)-1
         x(i) = x(i) - alu(k)* x(jlu(k))
41    continue
40 continue
   !
   !     backward solve.
   !
   do 90 i = n, 1, -1
      do 91 k=ju(i),jlu(i+1)-1
         x(i) = x(i) - alu(k)*x(jlu(k))
91    continue
      x(i) = alu(i)*x(i)
90 continue
   !
   return
   !----------------end of lusol ------------------------------------------
   !-----------------------------------------------------------------------
end
   !-----------------------------------------------------------------------
subroutine qsplit(a,ind,n,ncut)
   doubleprecision a(n)
   integer ind(n), n, ncut
   !-----------------------------------------------------------------------
   !     does a quick-sort split of a real array.
   !     on input a(1:n). is a real array
   !     on output a(1:n) is permuted such that its elements satisfy:
   !
   !     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
   !     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
   !
   !     ind(1:n) is an integer array which permuted in the same way as a(*).
   !-----------------------------------------------------------------------
   doubleprecision tmp, abskey
   integer itmp, first, last
   !-----
   first = 1
   last = n
   if (ncut .lt. first .or. ncut .gt. last) return
   !
   !     outer loop -- while mid .ne. ncut do
   !
1  mid = first
   abskey = abs(a(mid))
   do 2 j=first+1, last
      if (abs(a(j)) .gt. abskey) then
         mid = mid+1
   !     interchange
         tmp = a(mid)
         itmp = ind(mid)
         a(mid) = a(j)
         ind(mid) = ind(j)
         a(j)  = tmp
         ind(j) = itmp
      endif
2  continue
   !
   !     interchange
   !
   tmp = a(mid)
   a(mid) = a(first)
   a(first)  = tmp
   !
   itmp = ind(mid)
   ind(mid) = ind(first)
   ind(first) = itmp
   !
   !     test for while loop
   !
   if (mid .eq. ncut) return
   if (mid .gt. ncut) then
      last = mid-1
   else
      first = mid+1
   endif
   goto 1
   !----------------end-of-qsplit------------------------------------------
   !-----------------------------------------------------------------------
end
   !
   !------REORDERING ROUTINES
   !      SUBSET OF THE REVERSE CUTHILL MCKEE REORDERING LIBRARY
   !      MAINTAINED BY JOHN BURKARDT AT FLORIDA STATE UNIVERSITY
SUBROUTINE GSOL_GENRCM ( node_num, adj_num, adj_row, adj, perm )

   !*****************************************************************************80
   !
   !! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
   !
   !  Discussion:
   !
   !    For each connected component in the graph, the routine obtains
   !    an ordering by calling RCM.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    04 January 2003
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Alan George, Joseph Liu.
   !    FORTRAN90 version by John Burkardt
   !
   !  Reference:
   !
   !    Alan George, Joseph Liu,
   !    Computer Solution of Large Sparse Positive Definite Systems,
   !    Prentice Hall, 1981.
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
   !
   !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
   !
   !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about
   !    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
   !
   !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
   !    For each row, it contains the column indices of the nonzero entries.
   !
   !    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
   !
   !  Local Parameters:
   !
   !    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
   !    structure.  The level structure is stored in the currently unused
   !    spaces in the permutation vector PERM.
   !
   !    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
   !
   implicit none

   integer ( kind = 4 ) adj_num
   integer ( kind = 4 ) node_num

   integer ( kind = 4 ) adj(adj_num)
   integer ( kind = 4 ) adj_row(node_num+1)
   integer ( kind = 4 ) i
   integer ( kind = 4 ) iccsze
   integer ( kind = 4 ) mask(node_num)
   integer ( kind = 4 ) level_num
   integer ( kind = 4 ) level_row(node_num+1)
   integer ( kind = 4 ) num
   integer ( kind = 4 ) perm(node_num)
   integer ( kind = 4 ) root

   mask(1:node_num) = 1

   num = 1

   do i = 1, node_num
      !
      !  For each masked connected component...
      !
      if ( mask(i) /= 0 ) then

         root = i
         !
         !  Find a pseudo-peripheral node ROOT.  The level structure found by
         !  ROOT_FIND is stored starting at PERM(NUM).
         !
         call root_find ( root, adj_num, adj_row, adj, mask,&
         &level_num, level_row, perm(num), node_num )
         !
         !  RCM orders the component using ROOT as the starting node.
         !
         call GSOL_RCM( root, adj_num, adj_row, adj, mask, perm(num),&
         &iccsze, node_num )

         num = num + iccsze
         !
         !  We can stop once every node is in one of the connected components.
         !
         if ( node_num < num ) then
            return
         end if

      end if

   end do

   return
END SUBROUTINE GSOL_GENRCM

SUBROUTINE ROOT_FIND ( root, adj_num, adj_row, adj, mask,&
&level_num, level_row, level, node_num )

   !*****************************************************************************80
   !
   !! ROOT_FIND finds a pseudo-peripheral node.
   !
   !  Discussion:
   !
   !    The diameter of a graph is the maximum distance (number of edges)
   !    between any two nodes of the graph.
   !
   !    The eccentricity of a node is the maximum distance between that
   !    node and any other node of the graph.
   !
   !    A peripheral node is a node whose eccentricity equals the
   !    diameter of the graph.
   !
   !    A pseudo-peripheral node is an approximation to a peripheral node;
   !    it may be a peripheral node, but all we know is that we tried our
   !    best.
   !
   !    The routine is given a graph, and seeks pseudo-peripheral nodes,
   !    using a modified version of the scheme of Gibbs, Poole and
   !    Stockmeyer.  It determines such a node for the section subgraph
   !    specified by MASK and ROOT.
   !
   !    The routine also determines the level structure associated with
   !    the given pseudo-peripheral node; that is, how far each node
   !    is from the pseudo-peripheral node.  The level structure is
   !    returned as a list of nodes LS, and pointers to the beginning
   !    of the list of nodes that are at a distance of 0, 1, 2, ...,
   !    NODE_NUM-1 from the pseudo-peripheral node.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    28 October 2003
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Alan George, Joseph Liu.
   !    FORTRAN90 version by John Burkardt
   !
   !  Reference:
   !
   !    Alan George, Joseph Liu,
   !    Computer Solution of Large Sparse Positive Definite Systems,
   !    Prentice Hall, 1981.
   !
   !    Norman Gibbs, William Poole, Paul Stockmeyer,
   !    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
   !    SIAM Journal on Numerical Analysis,
   !    Volume 13, pages 236-250, 1976.
   !
   !    Norman Gibbs,
   !    Algorithm 509: A Hybrid Profile Reduction Algorithm,
   !    ACM Transactions on Mathematical Software,
   !    Volume 2, pages 378-387, 1976.
   !
   !  Parameters:
   !
   !    Input/output, integer ( kind = 4 ) ROOT.  On input, ROOT is a node in the
   !    the component of the graph for which a pseudo-peripheral node is
   !    sought.  On output, ROOT is the pseudo-peripheral node obtained.
   !
   !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
   !
   !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about
   !    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
   !
   !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
   !    For each row, it contains the column indices of the nonzero entries.
   !
   !    Input, integer ( kind = 4 ) MASK(NODE_NUM), specifies a section subgraph.
   !    Nodes for which MASK is zero are ignored by FNROOT.
   !
   !    Output, integer ( kind = 4 ) LEVEL_NUM, is the number of levels in the
   !    level structure rooted at the node ROOT.
   !
   !    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
   !    level structure array pair containing the level structure found.
   !
   !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
   !
   implicit none

   integer ( kind = 4 ) adj_num
   integer ( kind = 4 ) node_num

   integer ( kind = 4 ) adj(adj_num)
   integer ( kind = 4 ) adj_row(node_num+1)
   integer ( kind = 4 ) iccsze
   integer ( kind = 4 ) j
   integer ( kind = 4 ) jstrt
   integer ( kind = 4 ) k
   integer ( kind = 4 ) kstop
   integer ( kind = 4 ) kstrt
   integer ( kind = 4 ) level(node_num)
   integer ( kind = 4 ) level_num
   integer ( kind = 4 ) level_num2
   integer ( kind = 4 ) level_row(node_num+1)
   integer ( kind = 4 ) mask(node_num)
   integer ( kind = 4 ) mindeg
   integer ( kind = 4 ) nabor
   integer ( kind = 4 ) ndeg
   integer ( kind = 4 ) node
   integer ( kind = 4 ) root
   !
   !  Determine the level structure rooted at ROOT.
   !
   call level_set ( root, adj_num, adj_row, adj, mask, level_num,&
   &level_row, level, node_num )
   !
   !  Count the number of nodes in this level structure.
   !
   iccsze = level_row(level_num+1) - 1
   !
   !  Extreme case:
   !    A complete graph has a level set of only a single level.
   !    Every node is equally good (or bad).
   !
   if ( level_num == 1 ) then
      return
   end if
   !
   !  Extreme case:
   !    A "line graph" 0--0--0--0--0 has every node in its only level.
   !    By chance, we've stumbled on the ideal root.
   !
   if ( level_num == iccsze ) then
      return
   end if
   !
   !  Pick any node from the last level that has minimum degree
   !  as the starting point to generate a new level set.
   !
   do

      mindeg = iccsze

      jstrt = level_row(level_num)
      root = level(jstrt)

      if ( jstrt < iccsze ) then

         do j = jstrt, iccsze

            node = level(j)
            ndeg = 0
            kstrt = adj_row(node)
            kstop = adj_row(node+1) - 1

            do k = kstrt, kstop
               nabor = adj(k)
               if ( 0 < mask(nabor) ) then
                  ndeg = ndeg + 1
               end if
            end do

            if ( ndeg < mindeg ) then
               root = node
               mindeg = ndeg
            end if

         end do

      end if
      !
      !  Generate the rooted level structure associated with this node.
      !
      call level_set ( root, adj_num, adj_row, adj, mask,&
      &level_num2, level_row, level, node_num )
      !
      !  If the number of levels did not increase, accept the new ROOT.
      !
      if ( level_num2 <= level_num ) then
         exit
      end if

      level_num = level_num2
      !
      !  In the unlikely case that ROOT is one endpoint of a line graph,
      !  we can exit now.
      !
      if ( iccsze <= level_num ) then
         exit
      end if

   end do

   return

END SUBROUTINE ROOT_FIND

SUBROUTINE LEVEL_SET ( root, adj_num, adj_row, adj, mask,&
&level_num, level_row, level, node_num )

   !*****************************************************************************80
   !
   !! LEVEL_SET generates the connected level structure rooted at a given node.
   !
   !  Discussion:
   !
   !    Only nodes for which MASK is nonzero will be considered.
   !
   !    The root node chosen by the user is assigned level 1, and masked.
   !    All (unmasked) nodes reachable from a node in level 1 are
   !    assigned level 2 and masked.  The process continues until there
   !    are no unmasked nodes adjacent to any node in the current level.
   !    The number of levels may vary between 2 and NODE_NUM.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    28 October 2003
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Alan George, Joseph Liu.
   !    FORTRAN90 version by John Burkardt
   !
   !  Reference:
   !
   !    Alan George, Joseph Liu,
   !    Computer Solution of Large Sparse Positive Definite Systems,
   !    Prentice Hall, 1981.
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) ROOT, the node at which the level structure
   !    is to be rooted.
   !
   !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
   !
   !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about
   !    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
   !
   !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
   !    For each row, it contains the column indices of the nonzero entries.
   !
   !    Input/output, integer ( kind = 4 ) MASK(NODE_NUM).  On input, only nodes
   !    with nonzero MASK are to be processed.  On output, those nodes which were
   !    included in the level set have MASK set to 1.
   !
   !    Output, integer ( kind = 4 ) LEVEL_NUM, the number of levels in the level
   !    structure.  ROOT is in level 1.  The neighbors of ROOT
   !    are in level 2, and so on.
   !
   !    Output, integer ( kind = 4 ) LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM),
   !    the rooted level structure.
   !
   !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
   !
   implicit none

   integer ( kind = 4 ) adj_num
   integer ( kind = 4 ) node_num

   integer ( kind = 4 ) adj(adj_num)
   integer ( kind = 4 ) adj_row(node_num+1)
   integer ( kind = 4 ) i
   integer ( kind = 4 ) iccsze
   integer ( kind = 4 ) j
   integer ( kind = 4 ) jstop
   integer ( kind = 4 ) jstrt
   integer ( kind = 4 ) lbegin
   integer ( kind = 4 ) level_num
   integer ( kind = 4 ) level_row(node_num+1)
   integer ( kind = 4 ) level(node_num)
   integer ( kind = 4 ) lvlend
   integer ( kind = 4 ) lvsize
   integer ( kind = 4 ) mask(node_num)
   integer ( kind = 4 ) nbr
   integer ( kind = 4 ) node
   integer ( kind = 4 ) root

   mask(root) = 0
   level(1) = root
   level_num = 0
   lvlend = 0
   iccsze = 1
   !
   !  LBEGIN is the pointer to the beginning of the current level, and
   !  LVLEND points to the end of this level.
   !
   do

      lbegin = lvlend + 1
      lvlend = iccsze
      level_num = level_num + 1
      level_row(level_num) = lbegin
      !
      !  Generate the next level by finding all the masked neighbors of nodes
      !  in the current level.
      !
      do i = lbegin, lvlend

         node = level(i)
         jstrt = adj_row(node)
         jstop = adj_row(node+1) - 1

         do j = jstrt, jstop

            nbr = adj(j)

            if ( mask(nbr) /= 0 ) then
               iccsze = iccsze + 1
               level(iccsze) = nbr
               mask(nbr) = 0
            end if

         end do

      end do
      !
      !  Compute the current level width (the number of nodes encountered.)
      !  If it is positive, generate the next level.
      !
      lvsize = iccsze - lvlend

      if ( lvsize <= 0 ) then
         exit
      end if

   end do

   level_row(level_num+1) = lvlend + 1
   !
   !  Reset MASK to 1 for the nodes in the level structure.
   !
   mask(level(1:iccsze)) = 1

   return
END SUBROUTINE LEVEL_SET

SUBROUTINE GSOL_RCM( root, adj_num, adj_row, adj, mask, perm,&
&iccsze, node_num )

   !*****************************************************************************80
   !
   !! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
   !
   !  Discussion:
   !
   !    The connected component is specified by a node ROOT and a mask.
   !    The numbering starts at the root node.
   !
   !    An outline of the algorithm is as follows:
   !
   !    X(1) = ROOT.
   !
   !    for ( I = 1 to N-1)
   !      Find all unlabeled neighbors of X(I),
   !      assign them the next available labels, in order of increasing degree.
   !
   !    When done, reverse the ordering.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    02 January 2007
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Alan George, Joseph Liu.
   !    FORTRAN90 version by John Burkardt
   !
   !  Reference:
   !
   !    Alan George, Joseph Liu,
   !    Computer Solution of Large Sparse Positive Definite Systems,
   !    Prentice Hall, 1981.
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) ROOT, the node that defines the connected
   !    component.  It is used as the starting point for the RCM ordering.
   !
   !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
   !
   !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about
   !    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
   !
   !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
   !    For each row, it contains the column indices of the nonzero entries.
   !
   !    Input/output, integer ( kind = 4 ) MASK(NODE_NUM), a mask for the nodes.
   !    Only those nodes with nonzero input mask values are considered by the
   !    routine.  The nodes numbered by RCM will have their mask values
   !    set to zero.
   !
   !    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
   !
   !    Output, integer ( kind = 4 ) ICCSZE, the size of the connected component
   !    that has been numbered.
   !
   !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
   !
   !  Local Parameters:
   !
   !    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold
   !    the degree of the nodes in the section graph specified by mask and root.
   !
   implicit none

   integer ( kind = 4 ) adj_num
   integer ( kind = 4 ) node_num

   integer ( kind = 4 ) adj(adj_num)
   integer ( kind = 4 ) adj_row(node_num+1)
   integer ( kind = 4 ) deg(node_num)
   integer ( kind = 4 ) fnbr
   integer ( kind = 4 ) i
   integer ( kind = 4 ) iccsze
   integer ( kind = 4 ) j
   integer ( kind = 4 ) jstop
   integer ( kind = 4 ) jstrt
   integer ( kind = 4 ) k
   integer ( kind = 4 ) l
   integer ( kind = 4 ) lbegin
   integer ( kind = 4 ) lnbr
   integer ( kind = 4 ) lperm
   integer ( kind = 4 ) lvlend
   integer ( kind = 4 ) mask(node_num)
   integer ( kind = 4 ) nbr
   integer ( kind = 4 ) node
   integer ( kind = 4 ) perm(node_num)
   integer ( kind = 4 ) root
   !
   !  Find the degrees of the nodes in the component specified by MASK and ROOT.
   !
   call gsol_degree ( root, adj_num, adj_row, adj, mask, deg,&
   &iccsze, perm, node_num )

   mask(root) = 0

   if ( iccsze <= 1 ) then
      return
   end if

   lvlend = 0
   lnbr = 1
   !
   !  LBEGIN and LVLEND point to the beginning and
   !  the end of the current level respectively.
   !
   do while ( lvlend < lnbr )

      lbegin = lvlend + 1
      lvlend = lnbr

      do i = lbegin, lvlend
         !
         !  For each node in the current level...
         !
         node = perm(i)
         jstrt = adj_row(node)
         jstop = adj_row(node+1) - 1
         !
         !  Find the unnumbered neighbors of NODE.
         !
         !  FNBR and LNBR point to the first and last neighbors
         !  of the current node in PERM.
         !
         fnbr = lnbr + 1

         do j = jstrt, jstop

            nbr = adj(j)

            if ( mask(nbr) /= 0 ) then
               lnbr = lnbr + 1
               mask(nbr) = 0
               perm(lnbr) = nbr
            end if

         end do
         !
         !  If no neighbors, skip to next node in this level.
         !
         if ( lnbr <= fnbr ) then
            cycle
         end if
         !
         !  Sort the neighbors of NODE in increasing order by degree.
         !  Linear insertion is used.
         !
         k = fnbr

         do while ( k < lnbr )

            l = k
            k = k + 1
            nbr = perm(k)

            do while ( fnbr < l )

               lperm = perm(l)

               if ( deg(lperm) <= deg(nbr) ) then
                  exit
               end if

               perm(l+1) = lperm
               l = l - 1

            end do

            perm(l+1) = nbr

         end do

      end do

   end do
   !
   !  We now have the Cuthill-McKee ordering.  Reverse it.
   !
   call i4vec_reverse ( iccsze, perm )

   return
END SUBROUTINE GSOL_RCM

SUBROUTINE GSOL_DEGREE ( root, adj_num, adj_row, adj, mask, deg,&
&iccsze, ls, node_num )

   !*****************************************************************************80
   !
   !! DEGREE computes the degrees of the nodes in the connected component.
   !
   !  Discussion:
   !
   !    The connected component is specified by MASK and ROOT.
   !    Nodes for which MASK is zero are ignored.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    05 January 2003
   !
   !  Author:
   !
   !    Original FORTRAN77 version by Alan George, Joseph Liu.
   !    FORTRAN90 version by John Burkardt
   !
   !  Reference:
   !
   !    Alan George, Joseph Liu,
   !    Computer Solution of Large Sparse Positive Definite Systems,
   !    Prentice Hall, 1981.
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) ROOT, the node that defines the connected
   !    component.
   !
   !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
   !
   !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about
   !    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
   !
   !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
   !    For each row, it contains the column indices of the nonzero entries.
   !
   !    Input, integer ( kind = 4 ) MASK(NODE_NUM), is nonzero for those nodes
   !    which are to be considered.
   !
   !    Output, integer ( kind = 4 ) DEG(NODE_NUM), contains, for each  node in
   !    the connected component, its degree.
   !
   !    Output, integer ( kind = 4 ) ICCSIZE, the number of nodes in the
   !    connected component.
   !
   !    Output, integer ( kind = 4 ) LS(NODE_NUM), stores in entries 1 through
   !    ICCSIZE the nodes in the connected component, starting with ROOT, and
   !    proceeding by levels.
   !
   !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
   !
   implicit none

   integer ( kind = 4 ) adj_num
   integer ( kind = 4 ) node_num

   integer ( kind = 4 ) adj(adj_num)
   integer ( kind = 4 ) adj_row(node_num+1)
   integer ( kind = 4 ) deg(node_num)
   integer ( kind = 4 ) i
   integer ( kind = 4 ) iccsze
   integer ( kind = 4 ) ideg
   integer ( kind = 4 ) j
   integer ( kind = 4 ) jstop
   integer ( kind = 4 ) jstrt
   integer ( kind = 4 ) lbegin
   integer ( kind = 4 ) ls(node_num)
   integer ( kind = 4 ) lvlend
   integer ( kind = 4 ) lvsize
   integer ( kind = 4 ) mask(node_num)
   integer ( kind = 4 ) nbr
   integer ( kind = 4 ) node
   integer ( kind = 4 ) root
   !
   !  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
   !
   ls(1) = root
   adj_row(root) = -adj_row(root)
   lvlend = 0
   iccsze = 1
   !
   !  LBEGIN is the pointer to the beginning of the current level, and
   !  LVLEND points to the end of this level.
   !
   do

      lbegin = lvlend + 1
      lvlend = iccsze
      !
      !  Find the degrees of nodes in the current level,
      !  and at the same time, generate the next level.
      !
      do i = lbegin, lvlend

         node = ls(i)
         jstrt = -adj_row(node)
         jstop = abs ( adj_row(node+1) ) - 1
         ideg = 0

         do j = jstrt, jstop

            nbr = adj(j)

            if ( mask(nbr) /= 0 ) then

               ideg = ideg + 1

               if ( 0 <= adj_row(nbr) ) then
                  adj_row(nbr) = -adj_row(nbr)
                  iccsze = iccsze + 1
                  ls(iccsze) = nbr
               end if

            end if

         end do

         deg(node) = ideg

      end do
      !
      !  Compute the current level width.
      !
      lvsize = iccsze - lvlend
      !
      !  If the current level width is nonzero, generate another level.
      !
      if ( lvsize == 0 ) then
         exit
      end if

   end do
   !
   !  Reset ADJ_ROW to its correct sign and return.
   !
   do i = 1, iccsze
      node = ls(i)
      adj_row(node) = -adj_row(node)
   end do

   return
END SUBROUTINE GSOL_DEGREE

SUBROUTINE I4VEC_REVERSE ( n, a )

   !*****************************************************************************80
   !
   !! I4VEC_REVERSE reverses the elements of an I4VEC.
   !
   !  Example:
   !
   !    Input:
   !
   !      N = 5,
   !      A = ( 11, 12, 13, 14, 15 ).
   !
   !    Output:
   !
   !      A = ( 15, 14, 13, 12, 11 ).
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    26 July 1999
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) N, the number of entries in the array.
   !
   !    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
   !
   implicit none

   integer ( kind = 4 ) n

   integer ( kind = 4 ) a(n)
   integer ( kind = 4 ) i

   do i = 1, n/2
      call i4_swap ( a(i), a(n+1-i) )
   end do

   return
END SUBROUTINE I4VEC_REVERSE

SUBROUTINE I4_SWAP ( i, j )

   !*****************************************************************************80
   !
   !! I4_SWAP swaps two I4's.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    30 November 1998
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Parameters:
   !
   !    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
   !    J have been interchanged.
   !
   implicit none

   integer ( kind = 4 ) i
   integer ( kind = 4 ) j
   integer ( kind = 4 ) k

   k = i
   i = j
   j = k

   return
END SUBROUTINE I4_SWAP

SUBROUTINE GSOL_PERM_INVERSE3 ( n, perm, perm_inv )

   !*****************************************************************************80
   !
   !! PERM_INVERSE3 produces the inverse of a given permutation.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    28 October 2003
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) N, the number of items permuted.
   !
   !    Input, integer ( kind = 4 ) PERM(N), a permutation.
   !
   !    Output, integer ( kind = 4 ) PERM_INV(N), the inverse permutation.
   !
   implicit none

   integer ( kind = 4 ) n

   integer ( kind = 4 ) i
   integer ( kind = 4 ) perm(n)
   integer ( kind = 4 ) perm_inv(n)

   do i = 1, n
      perm_inv(perm(i)) = i
   end do

   return
END SUBROUTINE GSOL_PERM_INVERSE3

FUNCTION GSOL_ADJ_PERM_BANDWIDTH( node_num, adj_num, adj_row, adj,&
&perm, perm_inv )

   !*****************************************************************************80
   !
   !! ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
   !
   !  Discussion:
   !
   !    The matrix is defined by the adjacency information and a permutation.
   !
   !    The routine also computes the bandwidth and the size of the envelope.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    11 March 2005
   !
   !  Author:
   !
   !    John Burkardt
   !
   !  Reference:
   !
   !    Alan George, Joseph Liu,
   !    Computer Solution of Large Sparse Positive Definite Systems,
   !    Prentice Hall, 1981.
   !
   !  Parameters:
   !
   !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
   !
   !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
   !
   !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about
   !    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
   !
   !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
   !    For each row, it contains the column indices of the nonzero entries.
   !
   !    Input, integer ( kind = 4 ) PERM(NODE_NUM), PERM_INV(NODE_NUM), the
   !    permutation and inverse permutation.
   !
   !    Output, integer ( kind = 4 ) ADJ_PERM_BANDWIDTH, the bandwidth of the
   !    permuted adjacency matrix.
   !
   implicit none

   integer ( kind = 4 ) adj_num
   integer ( kind = 4 ) node_num

   integer ( kind = 4 ) adj(adj_num)
   integer ( kind = 4 ) gsol_adj_perm_bandwidth
   integer ( kind = 4 ) adj_row(node_num+1)
   integer ( kind = 4 ) band_hi
   integer ( kind = 4 ) band_lo
   integer ( kind = 4 ) col
   integer ( kind = 4 ) i
   integer ( kind = 4 ) j
   integer ( kind = 4 ) perm(node_num)
   integer ( kind = 4 ) perm_inv(node_num)

   band_lo = 0
   band_hi = 0

   do i = 1, node_num

      do j = adj_row(perm(i)), adj_row(perm(i)+1) - 1
         col = perm_inv(adj(j))
         band_lo = max ( band_lo, i - col )
         band_hi = max ( band_hi, col - i )
      end do

   end do

   gsol_adj_perm_bandwidth = band_lo + 1 + band_hi

   return
END FUNCTION GSOL_ADJ_PERM_BANDWIDTH

   !-------END OF SUBROUTINES FROM OTHER LIBRARIES
