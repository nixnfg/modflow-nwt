
   !      for XMD package for NWT solver (NWT version 1.3; 07-01-2022)
   !
   !                     M. Ibaraki
   !
   !             Thu Oct 29 10:03:25 EDT 2009
   !
   !     variable definitions:
   !
   !      iacl               choice of acceleration method
   !                         = 0; conjugate gradient
   !                         = 1; ORTHOMIN
   !                         = 2; CGSTAB
   !      n                  number of unknowns
   !      norder             = 0; original ordering
   !                         = 1; RCM ordering
   !                         = 2; Minimum Degree ordering
   !      nja                size of ja, a, arrays
   !      njaf               size of af, jaf arrays
   !      level              level of ILU
   !      itmax              number of maximum allowable iterations
   !      north              number of orthogonalization
   !      liwrk              size of integer work array
   !      lrwrk              size of real work array
   !
   !      ia(n+1),ja(nja)    usual ia, ja arrays for coefficient matrix
   !      lorder(n)          ordering vector: lorder( new_order ) = old_order
   !      iwork(liwrk)      temporary work array
   !
   !      dptol              flag for the drop tolerance
   !                         =.true. perform the drop tolerance
   !                         =.false. do NOT perform the drop tolerance
   !
   !      epsrn              drop tolerance
   !      ctol               absolute convergence criteria
   !      rrctol             residual reduction convergence criteria
   !
   !      rwork(lrwrk)       temporary work array
   !      a(nja)             matrix stored as linear array
   !      af(njaf)           factored matrix (each row of af contains a row L\U)
   !                         where A = LU
   !      b(n)               right hand side vector
   !      x(n)               solution vector
   !
   !
   !      nx,ny,nz           graph of matrix is regular rectangular grid
   !                         of size nx * ny * nz
   !
   !mi
   !      MODULE XMDMODULE
   !      IMPLICIT NONE
   !      LOGICAL, SAVE, POINTER ::  REDSYS,LDCOMB
   !      DOUBLE PRECISION, SAVE, POINTER ::  EPSRN,RRCTOL,HCLOSEXMD
   !      INTEGER, SAVE, POINTER ::  MXITERXMD
   !      DOUBLE PRECISION, SAVE, DIMENSION(:), POINTER ::  RWORK,AF,DGSCAL
   !      INTEGER, SAVE, DIMENSION(:), POINTER ::  LORDER,IWORK,MSINDX
   !      INTEGER, SAVE, POINTER :: IACL,NORDER,NJAF,LEVEL,NORTH,LIWRK,
   !     *  LRWRK,IDROPTOL,NBLACK,IERR,IDSCALE,MBLACK
   !      END MODULE XMDMODULE
   !mi


MODULE XMDMODULE
   IMPLICIT NONE
   LOGICAL, SAVE, POINTER ::  REDSYS,LDCOMB
   DOUBLE PRECISION, SAVE, POINTER ::  EPSRN,RRCTOL,HCLOSEXMD
   INTEGER, SAVE, POINTER ::  MXITERXMD
   INTEGER, SAVE, POINTER :: IACL,NORDER,LEVEL,NORTH,IDROPTOL,&
   &IERR
END MODULE XMDMODULE
   !------------------------------------------------------------------
SUBROUTINE XMD7AR(IN)

   USE GLOBAL, ONLY: IOUT
   USE XMDMODULE
   USE GWFNWTMODULE, ONLY: IPRNWT,NUMACTIVE,IA,JA,NJA,IFDPARAM
   !mi
   use xmdcmn
   !mi
   IMPLICIT NONE
   !     ------------------------------------------------------------------
   !     SPECIFICATIONS:
   !     ------------------------------------------------------------------
   INTRINSIC INT
   EXTERNAL URDCOM, URWORD
   !mi
   !     include 'xmdcmn.com'
   !mi
   !     ------------------------------------------------------------------
   !     ARGUMENTS
   !     ------------------------------------------------------------------
   INTEGER IN,NODES
   !     ------------------------------------------------------------------
   !     LOCAL VARIABLES
   !     ------------------------------------------------------------------
   INTEGER lloc, istart, istop, i, IREDSYS
   CHARACTER(LEN=200) line
   REAL R,RRCTOLS,EPSRNS,HCLOSEXMDDUM
   !     LOCAL VARIABLES FOR XMD SOLVER
   !     ------------------------------------------------------------------
   !
   !1------IDENTIFY PACKAGE AND INITIALIZE.
   WRITE(IOUT,1)IN
1  FORMAT(1X,'XMD -- LINEAR SOLUTION BY XMD PACKAGE VERSION',&
   &1X,'1.30',&
   &/1X,'    BY MOTOMU IBARAKI, OHIO STATE UNIVERSITY, COLOMBUS, OH',&
   &/1X,'                INPUT READ FROM UNIT',I3)
   !
   !mi
   !      ALLOCATE (IACL,NORDER,NJAF,LEVEL,NORTH,LIWRK,LRWRK,IDROPTOL,
   !     *  NBLACK,IERR,IDSCALE,MBLACK,HCLOSEXMD,MXITERXMD)
   !      ALLOCATE (EPSRN,RRCTOL)
   !      ALLOCATE (REDSYS,LDCOMB)
   !mi
   ALLOCATE (IACL,NORDER,LEVEL,NORTH,IDROPTOL,&
   &IERR,HCLOSEXMD,MXITERXMD)
   ALLOCATE (EPSRN,RRCTOL)
   ALLOCATE (REDSYS,LDCOMB)
   !-----XMD INPUT
   IF ( IFDPARAM.EQ.4 )CALL URDCOM(In, Iout, line)
   lloc = 1
   i = 1
   IF ( IFDPARAM.EQ.4 ) THEN
      CALL URWORD(line, lloc, istart, istop, 2, IACL, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, NORDER, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, LEVEL, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, NORTH, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, IREDSYS, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, I, RRCTOLS, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, IDROPTOL, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, I, EPSRNS, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, I,HCLOSEXMDDUM,Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, MXITERXMD,r,Iout,In)  !---added this 12/29/10
   ELSEIF ( IFDPARAM.EQ.1 ) THEN
      IACL = 1
      NORDER = 0
      LEVEL = 1
      NORTH = 5
      IREDSYS = 1
      RRCTOLS = 0.0
      IDROPTOL = 1
      EPSRNS = 5.0e-3
      HCLOSEXMDDUM = 1.0e-3
      MXITERXMD = 50
   ELSEIF ( IFDPARAM.EQ.2 ) THEN
      IACL = 1               !3/25/16
      NORDER = 0
      LEVEL = 1             !3/25/16
      NORTH = 10             !3/25/16
      IREDSYS = 1
      RRCTOLS = 0.0
      IDROPTOL = 1
      EPSRNS = 5.0e-3        !3/25/16
      HCLOSEXMDDUM = 1.0e-3  !3/25/16
      MXITERXMD = 100        !3/25/16
   ELSEIF ( IFDPARAM.EQ.3 ) THEN
      IACL = 2
      NORDER = 0  !3/25/16
      LEVEL = 15  !3/25/16
      NORTH = 1  !3/25/16
      IREDSYS = 1
      RRCTOLS = 0.0
      IDROPTOL = 1
      EPSRNS = 1.0e-3        !3/25/16
      HCLOSEXMDDUM = 1.0e-4  !3/25/16
      MXITERXMD = 200        !3/25/16
   END IF
   HCLOSEXMD = dble(HCLOSEXMDDUM)
   !
   !
   !
   !      read (*,*) IACL, NORDER, LEVEL, NORTH, IREDSYS, IDROPTOL, Icomb,
   !     [            EPSRNS
   !
   !      write (*,1192) IACL, NORDER, LEVEL, NORTH, IREDSYS, IDROPTOL,
   !     [               Icomb, EPSRNS
   !
   ! 1192 format(7i8, 1pe10.3)
   !
   IF ( HCLOSEXMD.LT.1.0e-8 ) HCLOSEXMD = 1.0e-8
   RRCTOL = RRCTOLS
   EPSRN = EPSRNS
   MIUNIT = IOUT
   MIOUT = IPRNWT - 2
   if(north.eq.0)north = 7   !!!!
   IF(EPSRN.LT.1.0E-20) EPSRN = 1.0e-3
   REDSYS = .FALSE.
   IF(IREDSYS.EQ.1) THEN
      REDSYS = .TRUE.
   ENDIF
   !
   WRITE(IOUT,23) IACL,NORDER,LEVEL,NORTH,IREDSYS,RRCTOL,&
   &IDROPTOL, EPSRN, Hclosexmd, MXITERXMD
23 FORMAT(1X,'ACCELERATION METHOD                    (IACL) = ',I9/&
   &1X,'NODE ORDERING FLAG                   (NORDER) = ',I9/&
   &1X,'LEVEL OF FILL                         (LEVEL) = ',I9/&
   &1X,'MAXIMUM NUMBER OF ORTHOGONALIZATIONS  (NORTH) = ',I9/&
   &1X,'INDEX FOR USING REDUCED SYSTEM      (IREDSYS) = ',I9/&
   &1X,'RESID. REDUCTION CONVERGE CRITERION  (RRCTOL) = ',E13.6/&
   &1X,'INDEX FOR USING DROP TOLERANCE     (IDROPTOL) = ',I9/&
   &1X,'DROP TOLERANCE VALUE                  (EPSRN) = ',E13.6/&
   &1X,'CONVERGENCE CRITERIA OF           (HCLOSEXMD) = ',E13.6/&
   &1X,'MAX. NUMBER OF LINEAR ITERATIONS  (MXITERXMD) = ',I9/)
   !
   !4-----ALLOCATE SPACE USED BY SOLVER
   NODES = NUMACTIVE
   !  -----------------
   !     preprocessor
   !  -----------------

   call xmdprpc(ia, ja, nja, numactive, norder, ierr, redsys)

   !  ------------------------------------
   !     check array sizes and structure
   !  ------------------------------------
   call xmdcheck(ia, ja, numactive, nja, ierr)

   !  ---------------------------------------------------
   !     PERFORM SYMBOLIC FACTORIZATION FOR LEVEL BASED PRECONDITIONING
   !  ---------------------------------------------------
   IF(IDROPTOL.EQ.0)THEN
   !  --------------------------------
   !     level based preconditioning
   !  --------------------------------
      call xmdprecl(ia, ja, level, nja, nodes, ierr)
   ENDIF
   !
   RETURN
END
   !-----------------------------------------------------------------------------------
SUBROUTINE XMD7DA(IGRID)
   !  DEALLOCATE GLOBAL DATA
   use xmdcmn
   use xmdmatrix
   INTEGER ALLOC_ERR, IGRID
   !
   DEALLOCATE(icolour, STAT = ALLOC_ERR)
   DEALLOCATE(RBorder, STAT = ALLOC_ERR)
   DEALLOCATE(iblackend, STAT = ALLOC_ERR)
   DEALLOCATE(lorder, STAT = ALLOC_ERR)
   DEALLOCATE(iaf, STAT = ALLOC_ERR)
   DEALLOCATE(idiagf, STAT = ALLOC_ERR)
   RETURN
END



