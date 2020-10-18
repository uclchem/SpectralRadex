MODULE Slatec
  USE types
  IMPLICIT INTEGER (I-N)
  IMPLICIT REAL(dp) (A-H,O-Z)
  !integer, parameter :: dp2 = selected_real_kind(15,307)

CONTAINS
  !LUDCMP can be skipped entirely in SLATEC (SGEIR does both)

  SUBROUTINE ludcmp(a,n,np,indx,d)
    INTEGER :: n,np,indx(n)
    REAL(dp) :: d,a(np,np)
    CONTINUE
    RETURN
  END

  SUBROUTINE lubksb(a,n,np,indx,b)
    INTEGER :: n,np,indx(np)
    REAL(dp) :: a(np,np),b(np)

    REAL(dp) :: ra(n-1,n-1),rb(n-1)
    INTEGER :: itask,iwork(n-1),ind,i,j
    REAL(dp) :: work((n-1)*(n-1+1))

    itask=1

    !Requires "smaller" matrix: nlev*nlev only, with last equation
    !replaced by sum_of_pops=1
    !This is a bit of a waste of computing time...

    do j=1,n-1
      do i=1,n-2
        ra(i,j)=a(i,j)
      enddo
      ra(n-1,j)=1.
    enddo
    do i=1,n-2
      rb(i)=0.
    enddo
    rb(n-1)=1.
    CALL sgeir(ra,n-1,n-1,rb,itask,ind,work,iwork)

    do i=1,n-1
      b(i)=rb(i)
    enddo

    RETURN
  END SUBROUTINE lubksb

  ! IMPLICIT INTEGER :: (I-N)
  ! IMPLICIT REAL(dp) :: (A-H,O-Z)
  !***BEGIN PROLOGUE  SGEIR
  !***PURPOSE  Solve a general system of linear equations.  Iterative
  !  refinement is used to obtain an error estimate.
  !***LIBRARY   SLATEC
  !***CATEGORY  D2A1
  !***TYPE      SINGLE PRECISION (SGEIR-S, CGEIR-C)
  !***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
  !   GENERAL SYSTEM OF LINEAR EQUATIONS
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
  !
  !  Subroutine SGEIR solves a general NxN system of single
  !  precision linear equations using LINPACK subroutines SGEFA and
  !  SGESL.  One pass of iterative refinement is used only to obtain
  !  an estimate of the accuracy.  That is, if A is an NxN real
  !  matrix and if X and B are real N-vectors, then SGEIR solves
  !  the equation
  !
  !         A*X=B.
  !
  !  The matrix A is first factored into upper and lower tri-
  !  angular matrices U and L using partial pivoting.  These
  !  factors and the pivoting information are used to calculate
  !  the solution, X.  Then the residual vector is found and
  !  used to calculate an estimate of the relative error, IND.
  !  IND estimates the accuracy of the solution only when the
  !  input matrix and the right hand side are represented
  !  exactly in the computer and does not take into account
  !  any errors in the input data.
  !
  !  If the equation A*X=B is to be solved for more than one vector
  !  B, the factoring of A does not need to be performed again and
  !  the option to solve only (ITASK .GT. 1) will be faster for
  !  the succeeding solutions.  In this case, the contents of A,
  !  LDA, N, WORK, and IWORK must not have been altered by the
  !  user following factorization (ITASK=1).  IND will not be
  !  changed by SGEIR in this case.
  !
  !  Argument Description ***
  !
  !  A      DOUBLE PRECISION(LDA,N)
  !   the doubly subscripted array with dimension (LDA,N)
  !   which contains the coefficient matrix.  A is not
  !   altered by the routine.
  !  LDA    INTEGER
  !   the leading dimension of the array A.  LDA must be great-
  !   er than or equal to N.  (terminal error message IND=-1)
  !  N      INTEGER
  !   the order of the matrix A.  The first N elements of
  !   the array A are the elements of the first column of
  !   matrix A.  N must be greater than or equal to 1.
  !   (terminal error message IND=-2)
  !  V      DOUBLE PRECISION(N)
  !   on entry, the singly subscripted array(vector) of di-
  !     mension N which contains the right hand side B of a
  !     system of simultaneous linear equations A*X=B.
  !   on return, V contains the solution vector, X .
  !  ITASK  INTEGER
  !   If ITASK=1, the matrix A is factored and then the
  !     linear equation is solved.
  !   If ITASK .GT. 1, the equation is solved using the existing
  !     factored matrix A (stored in WORK).
  !   If ITASK .LT. 1, then terminal error message IND=-3 is
  !     printed.
  !  IND    INTEGER
  !   GT. 0  IND is a rough estimate of the number of digits
  !    of accuracy in the solution, X.  IND=75 means
  !    that the solution vector X is zero.
  !   LT. 0  see error message corresponding to IND below.
  !  WORK   DOUBLE PRECISION(N*(N+1))
  !   a singly subscripted array of dimension at least N*(N+1).
  !  IWORK  INTEGER(N)
  !   a singly subscripted array of dimension at least N.
  !
  !  Error Messages Printed ***
  !
  !  IND=-1  terminal   N is greater than LDA.
  !  IND=-2  terminal   N is less than one.
  !  IND=-3  terminal   ITASK is less than one.
  !  IND=-4  terminal   The matrix A is computationally singular.
  !        A solution has not been computed.
  !  IND=-10 warning    The solution has no apparent significance.
  !        The solution may be inaccurate or the matrix
  !        A may be poorly scaled.
  !
  !     Note-  The above terminal(*fatal*) error messages are
  !     designed to be handled by XERMSG in which
  !     LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
  !     for warning error messages from XERMSG.  Unless
  !     the user provides otherwise, an error message
  !     will be printed followed by an abort.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  R1MACH, SASUM, SCOPY, SDSDOT, SGEFA, SGESL, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  ! 800430  DATE WRITTEN
  ! 890531  Changed all specific intrinsics to generic.  (WRB)
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  ! 900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SGEIR
  !
  SUBROUTINE SGEIR (A, LDA, N, V, ITASK, IND, WORK, IWORK)

    INTEGER :: LDA,N,ITASK,IND,IWORK(*),INFO,J
    REAL(dp) :: A(LDA,*),V(*),WORK(N,*),XNORM,DNORM
    CHARACTER*8 XERN1, XERN2
    !***FIRST EXECUTABLE STATEMENT  SGEIR
    IF (LDA.LT.N) THEN
       IND = -1
       WRITE (XERN1, '(I8)') LDA
       WRITE (XERN2, '(I8)') N
       CALL XERMSG ('SLATEC', 'SGEIR', 'LDA = ' // XERN1 //&
          &' IS LESS THAN N = ' // XERN2, -1, 1)
       RETURN
    ENDIF
    !
    IF (N.LE.0) THEN
       IND = -2
       WRITE (XERN1, '(I8)') N
       CALL XERMSG ('SLATEC', 'SGEIR', 'N = ' // XERN1 //&
          &' IS LESS THAN 1', -2, 1)
       RETURN
    ENDIF
    !
    IF (ITASK.LT.1) THEN
       IND = -3
       WRITE (XERN1, '(I8)') ITASK
       CALL XERMSG ('SLATEC', 'SGEIR', 'ITASK = ' // XERN1 //&
         & ' IS LESS THAN 1', -3, 1)
       RETURN
    ENDIF
    !
    IF (ITASK.EQ.1) THEN
      !
      ! MOVE MATRIX A TO WORK
      !
      DO J=1,N
        CALL SCOPY(N,A(1,J),1,WORK(1,J),1)
      END DO
      !
      ! FACTOR MATRIX A INTO LU
      !
      CALL SGEFA(WORK,N,N,IWORK,INFO)
      !
      ! CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
      !
      IF (INFO.NE.0) THEN
        IND = -4
        CALL XERMSG ('SLATEC', 'SGEIR',&
           &'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
        RETURN
      ENDIF
    ENDIF
    !
    !SOLVE WHEN FACTORING COMPLETE
    !MOVE VECTOR B TO WORK
    !
    CALL SCOPY(N,V(1),1,WORK(1,N+1),1)
    CALL SGESL(WORK,N,N,IWORK,V,0)
    !
    !FORM NORM OF X0
    !
    XNORM=SASUM(N,V(1),1)
    IF (XNORM.EQ.0.0) THEN
       IND = 75
       RETURN
    ENDIF
    !
    !COMPUTE  RESIDUAL
    !
    DO J=1,N
       WORK(J,N+1) = SDSDOT(N,-WORK(J,N+1),A(J,1),LDA,V,1)
    END DO
    !
    !SOLVE A*DELTA=R
    !
    CALL SGESL(WORK,N,N,IWORK,WORK(1,N+1),0)
    !
    !FORM NORM OF DELTA
    !
    DNORM = SASUM(N,WORK(1,N+1),1)
    !
    !COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
    !AND CHECK FOR IND GREATER THAN ZERO
    !
    IND = -LOG10(MAX(R1MACH(4),DNORM/XNORM))
    IF (IND.LE.0) THEN
       IND = -10
    !cc         CALL XERMSG ('SLATEC', 'SGEIR',
    !cc     *      'SOLUTION MAY HAVE NO SIGNIFICANCE ...', -10, 0)
    !cc   Error message disabled by mrh 03jan02. In solution of stateq this
    !message is triggered, but the solution does seem to have
    !significance. No need to print to screen, then.
    !Checked on 2-level system and standard example files: OK.
    ENDIF
    RETURN
  END SUBROUTINE SGEIR

  !***BEGIN PROLOGUE  R1MACH
  !***PURPOSE  Return floating point machine dependent constants.
  !***LIBRARY   SLATEC
  !***CATEGORY  R1
  !***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
  !***KEYWORDS  MACHINE CONSTANTS
  !***AUTHOR  Fox, P. A., (Bell Labs)
  ! Hall, A. D., (Bell Labs)
  ! Schryer, N. L., (Bell Labs)
  !***DESCRIPTION
  !
  ! R1MACH can be used to obtain machine-dependent parameters for the
  ! local machine environment.  It is a function subprogram with one
  ! (input) argument, and can be referenced as follows:
  !
  ! A = R1MACH(I)
  !
  ! where I=1,...,5.  The (output) value of A above is determined by
  ! the (input) value of I.  The results for various values of I are
  ! discussed below.
  !
  ! R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
  ! R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
  ! R1MACH(3) = B**(-T), the smallest relative spacing.
  ! R1MACH(4) = B**(1-T), the largest relative spacing.
  ! R1MACH(5) = LOG10(B)
  !
  ! Assume single precision numbers are represented in the T-digit,
  ! base-B form
  !
  !    sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
  !
  ! where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
  ! EMIN .LE. E .LE. EMAX.
  !
  ! The values of B, T, EMIN and EMAX are provided in I1MACH as
  ! follows:
  ! I1MACH(10) = B, the base.
  ! I1MACH(11) = T, the number of base-B digits.
  ! I1MACH(12) = EMIN, the smallest exponent E.
  ! I1MACH(13) = EMAX, the largest exponent E.
  !
  ! To alter this function for a particular environment, the desired
  ! set of DATA statements should be activated by removing the C from
  ! column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
  ! checked for consistency with the local operating system.
  !
  !***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
  !a portable library, ACM Transactions on Mathematical
  !Software 4, 2 (June 1978), pp. 177-188.
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  ! 790101  DATE WRITTEN
  ! 890213  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  ! 900618  Added DEC RISC constants.  (WRB)
  ! 900723  Added IBM RS 6000 constants.  (WRB)
  ! 910710  Added HP 730 constants.  (SMR)
  ! 911114  Added Convex IEEE constants.  (WRB)
  ! 920121  Added SUN -r8 compiler option constants.  (WRB)
  ! 920229  Added Touchstone Delta i860 constants.  (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  ! 920625  Added CONVEX -p8 and -pd8 compiler option constants.
  ! (BKS, WRB)
  ! 930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
  !***END PROLOGUE  R1MACH
  REAL(dp) FUNCTION R1MACH (I)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)
    
    !
    INTEGER :: SMALL(2)
    INTEGER :: LARGE(2)
    INTEGER :: RIGHT(2)
    INTEGER :: DIVER(2)
    INTEGER :: LOG10(2)
    !
    REAL(dp) :: RMACH(5)
    SAVE RMACH
    !
    EQUIVALENCE (RMACH(1),SMALL(1))
    EQUIVALENCE (RMACH(2),LARGE(1))
    EQUIVALENCE (RMACH(3),RIGHT(1))
    EQUIVALENCE (RMACH(4),DIVER(1))
    EQUIVALENCE (RMACH(5),LOG10(1))
    !
    !MACHINE CONSTANTS FOR THE AMIGA
    !ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
    !
    !DATA SMALL(1) / Z'00800000' /
    !DATA LARGE(1) / Z'7F7FFFFF' /
    !DATA RIGHT(1) / Z'33800000' /
    !DATA DIVER(1) / Z'34000000' /
    !DATA LOG10(1) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE AMIGA
    !ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
    !
    !DATA SMALL(1) / Z'00800000' /
    !DATA LARGE(1) / Z'7EFFFFFF' /
    !DATA RIGHT(1) / Z'33800000' /
    !DATA DIVER(1) / Z'34000000' /
    !DATA LOG10(1) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE APOLLO
    !
    !DATA SMALL(1) / 16#00800000 /
    !DATA LARGE(1) / 16#7FFFFFFF /
    !DATA RIGHT(1) / 16#33800000 /
    !DATA DIVER(1) / 16#34000000 /
    !DATA LOG10(1) / 16#3E9A209B /
    !
    !MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
    !
    !DATA RMACH(1) / Z400800000 /
    !DATA RMACH(2) / Z5FFFFFFFF /
    !DATA RMACH(3) / Z4E9800000 /
    !DATA RMACH(4) / Z4EA800000 /
    !DATA RMACH(5) / Z500E730E8 /
    !
    !MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
    !
    !DATA RMACH(1) / O1771000000000000 /
    !DATA RMACH(2) / O0777777777777777 /
    !DATA RMACH(3) / O1311000000000000 /
    !DATA RMACH(4) / O1301000000000000 /
    !DATA RMACH(5) / O1157163034761675 /
    !
    !MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
    !
    !DATA RMACH(1) / Z"3001800000000000" /
    !DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
    !DATA RMACH(3) / Z"3FD2800000000000" /
    !DATA RMACH(4) / Z"3FD3800000000000" /
    !DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
    !
    !MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
    !
    !DATA RMACH(1) / 00564000000000000000B /
    !DATA RMACH(2) / 37767777777777777776B /
    !DATA RMACH(3) / 16414000000000000000B /
    !DATA RMACH(4) / 16424000000000000000B /
    !DATA RMACH(5) / 17164642023241175720B /
    !
    !MACHINE CONSTANTS FOR THE CELERITY C1260
    !
    !DATA SMALL(1) / Z'00800000' /
    !DATA LARGE(1) / Z'7F7FFFFF' /
    !DATA RIGHT(1) / Z'33800000' /
    !DATA DIVER(1) / Z'34000000' /
    !DATA LOG10(1) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE CONVEX
    !USING THE -fn COMPILER OPTION
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7FFFFFFF' /
    !DATA RMACH(3) / Z'34800000' /
    !DATA RMACH(4) / Z'35000000' /
    !DATA RMACH(5) / Z'3F9A209B' /
    !
    !MACHINE CONSTANTS FOR THE CONVEX
    !USING THE -fi COMPILER OPTION
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7F7FFFFF' /
    !DATA RMACH(3) / Z'33800000' /
    !DATA RMACH(4) / Z'34000000' /
    !DATA RMACH(5) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE CONVEX
    !USING THE -p8 OR -pd8 COMPILER OPTION
    !
    !DATA RMACH(1) / Z'0010000000000000' /
    !DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
    !DATA RMACH(3) / Z'3CC0000000000000' /
    !DATA RMACH(4) / Z'3CD0000000000000' /
    !DATA RMACH(5) / Z'3FF34413509F79FF' /
    !
    !MACHINE CONSTANTS FOR THE CRAY
    !
    !DATA RMACH(1) / 200034000000000000000B /
    !DATA RMACH(2) / 577767777777777777776B /
    !DATA RMACH(3) / 377224000000000000000B /
    !DATA RMACH(4) / 377234000000000000000B /
    !DATA RMACH(5) / 377774642023241175720B /
    !
    !MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
    !NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
    !STATIC RMACH(5)
    !
    !DATA SMALL /    20K,       0 /
    !DATA LARGE / 77777K, 177777K /
    !DATA RIGHT / 35420K,       0 /
    !DATA DIVER / 36020K,       0 /
    !DATA LOG10 / 40423K,  42023K /
    !
    !MACHINE CONSTANTS FOR THE DEC ALPHA
    !USING G_FLOAT
    !
    !DATA RMACH(1) / '00000080'X /
    !DATA RMACH(2) / 'FFFF7FFF'X /
    !DATA RMACH(3) / '00003480'X /
    !DATA RMACH(4) / '00003500'X /
    !DATA RMACH(5) / '209B3F9A'X /
    !
    !MACHINE CONSTANTS FOR THE DEC ALPHA
    !USING IEEE_FLOAT
    !
    !DATA RMACH(1) / '00800000'X /
    !DATA RMACH(2) / '7F7FFFFF'X /
    !DATA RMACH(3) / '33800000'X /
    !DATA RMACH(4) / '34000000'X /
    !DATA RMACH(5) / '3E9A209B'X /
    !
    !MACHINE CONSTANTS FOR THE DEC RISC
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7F7FFFFF' /
    !DATA RMACH(3) / Z'33800000' /
    !DATA RMACH(4) / Z'34000000' /
    !DATA RMACH(5) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE DEC VAX
    !(EXPRESSED IN INTEGER :: AND HEXADECIMAL)
    !THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
    !THE INTEGER :: FORMAT SHOULD BE OK FOR UNIX SYSTEMS
    !
    !DATA SMALL(1) /       128 /
    !DATA LARGE(1) /    -32769 /
    !DATA RIGHT(1) /     13440 /
    !DATA DIVER(1) /     13568 /
    !DATA LOG10(1) / 547045274 /
    !
    !DATA SMALL(1) / Z00000080 /
    !DATA LARGE(1) / ZFFFF7FFF /
    !DATA RIGHT(1) / Z00003480 /
    !DATA DIVER(1) / Z00003500 /
    !DATA LOG10(1) / Z209B3F9A /
    !
    !MACHINE CONSTANTS FOR THE ELXSI 6400
    !(ASSUMING DOUBLE PRECISION*4 IS THE DEFAULT DOUBLE PRECISION)
    !
    !DATA SMALL(1) / '00800000'X /
    !DATA LARGE(1) / '7F7FFFFF'X /
    !DATA RIGHT(1) / '33800000'X /
    !DATA DIVER(1) / '34000000'X /
    !DATA LOG10(1) / '3E9A209B'X /
    !
    !MACHINE CONSTANTS FOR THE HARRIS 220
    !
    !DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
    !DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
    !DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
    !DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
    !DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
    !
    !MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
    !
    !DATA RMACH(1) / O402400000000 /
    !DATA RMACH(2) / O376777777777 /
    !DATA RMACH(3) / O714400000000 /
    !DATA RMACH(4) / O716400000000 /
    !DATA RMACH(5) / O776464202324 /
    !
    !MACHINE CONSTANTS FOR THE HP 730
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7F7FFFFF' /
    !DATA RMACH(3) / Z'33800000' /
    !DATA RMACH(4) / Z'34000000' /
    !DATA RMACH(5) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE HP 2100
    !3 WORD REAL(dp) :: WITH FTN4
    !
    !DATA SMALL(1), SMALL(2) / 40000B,       1 /
    !DATA LARGE(1), LARGE(2) / 77777B, 177776B /
    !DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
    !DATA DIVER(1), DIVER(2) / 40000B,    327B /
    !DATA LOG10(1), LOG10(2) / 46420B,  46777B /
    !
    !MACHINE CONSTANTS FOR THE HP 2100
    !4 WORD REAL(dp) :: WITH FTN4
    !
    !DATA SMALL(1), SMALL(2) / 40000B,       1 /
    !DATA LARGE(1), LARGE(2) / 77777B, 177776B /
    !DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
    !DATA DIVER(1), DIVER(2) / 40000B,    327B /
    !DATA LOG10(1), LOG10(2) / 46420B,  46777B /
    !
    !MACHINE CONSTANTS FOR THE HP 9000
    !
    !DATA SMALL(1) / 00004000000B /
    !DATA LARGE(1) / 17677777777B /
    !DATA RIGHT(1) / 06340000000B /
    !DATA DIVER(1) / 06400000000B /
    !DATA LOG10(1) / 07646420233B /
    !
    !MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
    !THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
    !THE PERKIN ELMER (INTERDATA) 7/32.
    !
    !DATA RMACH(1) / Z00100000 /
    !DATA RMACH(2) / Z7FFFFFFF /
    !DATA RMACH(3) / Z3B100000 /
    !DATA RMACH(4) / Z3C100000 /
    !DATA RMACH(5) / Z41134413 /
    !
    !MACHINE CONSTANTS FOR THE IBM PC
    !
    !DATA SMALL(1) / 1.18E-38      /
    !DATA LARGE(1) / 3.40E+38      /
    !DATA RIGHT(1) / 0.595E-07     /
    !DATA DIVER(1) / 1.19E-07      /
    !DATA LOG10(1) / 0.30102999566 /
    !
    !MACHINE CONSTANTS FOR THE IBM RS 6000
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7F7FFFFF' /
    !DATA RMACH(3) / Z'33800000' /
    !DATA RMACH(4) / Z'34000000' /
    !DATA RMACH(5) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE INTEL i860
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7F7FFFFF' /
    !DATA RMACH(3) / Z'33800000' /
    !DATA RMACH(4) / Z'34000000' /
    !DATA RMACH(5) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
    !
    !DATA RMACH(1) / "000400000000 /
    !DATA RMACH(2) / "377777777777 /
    !DATA RMACH(3) / "146400000000 /
    !DATA RMACH(4) / "147400000000 /
    !DATA RMACH(5) / "177464202324 /
    !
    !MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
    !32-BIT INTEGERS (EXPRESSED IN INTEGER :: AND OCTAL).
    !
    !DATA SMALL(1) /    8388608 /
    !DATA LARGE(1) / 2147483647 /
    !DATA RIGHT(1) /  880803840 /
    !DATA DIVER(1) /  889192448 /
    !DATA LOG10(1) / 1067065499 /
    !
    !DATA RMACH(1) / O00040000000 /
    !DATA RMACH(2) / O17777777777 /
    !DATA RMACH(3) / O06440000000 /
    !DATA RMACH(4) / O06500000000 /
    !DATA RMACH(5) / O07746420233 /
    !
    !MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
    !16-BIT INTEGERS  (EXPRESSED IN INTEGER :: AND OCTAL).
    !
    !DATA SMALL(1), SMALL(2) /   128,     0 /
    !DATA LARGE(1), LARGE(2) / 32767,    -1 /
    !DATA RIGHT(1), RIGHT(2) / 13440,     0 /
    !DATA DIVER(1), DIVER(2) / 13568,     0 /
    !DATA LOG10(1), LOG10(2) / 16282,  8347 /
    !
    !DATA SMALL(1), SMALL(2) / O000200, O000000 /
    !DATA LARGE(1), LARGE(2) / O077777, O177777 /
    !DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
    !DATA DIVER(1), DIVER(2) / O032400, O000000 /
    !DATA LOG10(1), LOG10(2) / O037632, O020233 /
    !
    !MACHINE CONSTANTS FOR THE SILICON GRAPHICS
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7F7FFFFF' /
    !DATA RMACH(3) / Z'33800000' /
    !DATA RMACH(4) / Z'34000000' /
    !DATA RMACH(5) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE SUN
    !
    !DATA RMACH(1) / Z'00800000' /
    !DATA RMACH(2) / Z'7F7FFFFF' /
    !DATA RMACH(3) / Z'33800000' /
    !DATA RMACH(4) / Z'34000000' /
    !DATA RMACH(5) / Z'3E9A209B' /
    !
    !MACHINE CONSTANTS FOR THE SUN
    !USING THE -r8 COMPILER OPTION
    !
    !DATA RMACH(1) / Z'0010000000000000' /
    !DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
    !DATA RMACH(3) / Z'3CA0000000000000' /
    !DATA RMACH(4) / Z'3CB0000000000000' /
    !DATA RMACH(5) / Z'3FD34413509F79FF' /
    !
    !MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
    !
    !DATA RMACH(1) / O000400000000 /
    !DATA RMACH(2) / O377777777777 /
    !DATA RMACH(3) / O146400000000 /
    !DATA RMACH(4) / O147400000000 /
    !DATA RMACH(5) / O177464202324 /
    !
    !MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
    !
    !DATA SMALL(1), SMALL(2) /     0,    256/
    !DATA LARGE(1), LARGE(2) /    -1,   -129/
    !DATA RIGHT(1), RIGHT(2) /     0,  26880/
    !DATA DIVER(1), DIVER(2) /     0,  27136/
    !DATA LOG10(1), LOG10(2) /  8347,  32538/
    !
    !***FIRST EXECUTABLE STATEMENT  R1MACH
    IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'R1MACH',&
     &'I OUT OF BOUNDS', 1, 2)
!
      R1MACH = RMACH(I)
      RETURN
!
  END FUNCTION R1MACH

  !***BEGIN PROLOGUE  SASUM
  !***PURPOSE  Compute the sum of the magnitudes of the elements of a
  !  vector.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A3A
  !***TYPE      SINGLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  ! Hanson, R. J., (SNLA)
  ! Kincaid, D. R., (U. of Texas)
  ! Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !      B L A S  Subprogram
  !  Description of Parameters
  !
  !--Input--
  ! N  number of elements in input vector(S)
  !  SX  single precision vector with N elements
  !INCX  storage spacing between elements of SX
  !
  !--Output--
  !  SASUM  single precision result (zero if N .LE. 0)
  !
  !Returns sum of magnitudes of single precision SX.
  !SASUM = sum from 0 to N-1 of ABS(SX(IX+I*INCX)),
  !where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !Krogh, Basic linear algebra subprograms for Fortran
  !usage, Algorithm No. 539, Transactions on Mathematical
  !Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 791001  DATE WRITTEN
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900821  Modified to correct problem with a negative increment.
  ! (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SASUM
  REAL(dp) FUNCTION SASUM (N, SX, INCX)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    REAL(dp) :: SX(*)
    INTEGER :: I, INCX, IX, M, MP1, N

    SASUM = 0.0D0
    IF (N .GT. 0) THEN
      !
      IF (INCX .NE. 1) THEN
        !
        !Code for increment not equal to 1.
        !
        IX = 1
        IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
        DO I = 1,N
          SASUM = SASUM + DABS(SX(IX))
          IX = IX + INCX
        END DO
      ELSE
        !
        !Code for increment equal to 1.
        !
        !Clean-up loop so remaining vector length is a multiple of 6.
        !
        M = MOD(N,6)
        IF (M .NE. 0) THEN
          DO I = 1,M
            SASUM = SASUM + DABS(SX(I))
          END DO
          IF (N .LT. 6) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1,N,6
          SASUM = SASUM + DABS(SX(I)) + DABS(SX(I+1)) + DABS(SX(I+2)) +&
                  &DABS(SX(I+3)) + DABS(SX(I+4)) + DABS(SX(I+5))
        END DO
      END IF
    END IF
    RETURN
  END FUNCTION SASUM

  !***BEGIN PROLOGUE  SCOPY
  !***PURPOSE  Copy a vector.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A5
  !***TYPE      SINGLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
  !***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  ! Hanson, R. J., (SNLA)
  ! Kincaid, D. R., (U. of Texas)
  ! Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !      B L A S  Subprogram
  !  Description of Parameters
  !
  !--Input--
  ! N  number of elements in input vector(s)
  !  SX  single precision vector with N elements
  !INCX  storage spacing between elements of SX
  !  SY  single precision vector with N elements
  !INCY  storage spacing between elements of SY
  !
  !--Output--
  !  SY  copy of vector SX (unchanged if N .LE. 0)
  !
  !Copy single precision SX to single precision SY.
  !For I = 0 to N-1, copy  SX(LX+I*INCX) to SY(LY+I*INCY),
  !where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !defined in a similar way using INCY.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !Krogh, Basic linear algebra subprograms for Fortran
  !usage, Algorithm No. 539, Transactions on Mathematical
  !Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 791001  DATE WRITTEN
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SCOPY
  SUBROUTINE SCOPY (N, SX, INCX, SY, INCY)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    REAL(dp) :: SX(*), SY(*)
    IF (N .GT. 0) THEN
      IF ((INCX .NE. INCY) .or. ((INCX-1) .lt. 0.0)) THEN
        !
        !Code for unequal or nonpositive increments.
        !
        IX = 1
        IY = 1
        IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
        IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
        DO I = 1,N
          SY(IY) = SX(IX)
          IX = IX + INCX
          IY = IY + INCY
        END DO
      ELSE IF ((INCX-1) .eq. 0.0) THEN
        !
        !Code for both increments equal to 1.
        !
        !Clean-up loop so remaining vector length is a multiple of 7.
        !
        M = MOD(N,7)
        IF (M .NE. 0) THEN
          DO I = 1,M
            SY(I) = SX(I)
          END DO
          IF (N .LT. 7) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1,N,7
          SY(I) = SX(I)
          SY(I+1) = SX(I+1)
          SY(I+2) = SX(I+2)
          SY(I+3) = SX(I+3)
          SY(I+4) = SX(I+4)
          SY(I+5) = SX(I+5)
          SY(I+6) = SX(I+6)
        END DO
      ELSE 
        !
        !Code for equal, positive, non-unit increments.
        !
        NS = N*INCX
        DO I = 1,NS,INCX
          SY(I) = SX(I)
        END DO
      END IF
    END IF
    RETURN
  END SUBROUTINE SCOPY


  !***BEGIN PROLOGUE  SDSDOT
  !***PURPOSE  Compute the inner product of two vectors with extended
  !  precision accumulation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A4
  !***TYPE      SINGLE PRECISION (SDSDOT-S, CDCDOT-C)
  !***KEYWORDS  BLAS, DOT PRODUCT, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  ! Hanson, R. J., (SNLA)
  ! Kincaid, D. R., (U. of Texas)
  ! Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !      B L A S  Subprogram
  !  Description of Parameters
  !
  !--Input--
  ! N  number of elements in input vector(s)
  !  SB  single precision scalar to be added to inner product
  !  SX  single precision vector with N elements
  !INCX  storage spacing between elements of SX
  !  SY  single precision vector with N elements
  !INCY  storage spacing between elements of SY
  !
  !--Output--
  ! SDSDOT  single precision dot product (SB if N .LE. 0)
  !
  !Returns S.P. result with dot product accumulated in D.P.
  !SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
  !where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !defined in a similar way using INCY.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !Krogh, Basic linear algebra subprograms for Fortran
  !usage, Algorithm No. 539, Transactions on Mathematical
  !Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 791001  DATE WRITTEN
  ! 890531  Changed all specific intrinsics to generic.  (WRB)
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SDSDOT
  REAL(dp) FUNCTION SDSDOT (N, SB, SX, INCX, SY, INCY)
    !IMPLICIT INTEGER :: (I-N)
    !IMPLICIT REAL(dp) :: (A-H,O-Z)

    REAL(dp) :: SX(*), SY(*), SB
    REAL(dp) :: DSDOT

    DSDOT = SB
    IF (N .GT. 0) THEN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
        !
        !Code for equal and positive increments.
        !
        NS = N*INCX
        DO I = 1,NS,INCX
           DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
        END DO
        SDSDOT = DSDOT
      ELSE
        !
        !Code for unequal or nonpositive increments.
        !
        KX = 1
        KY = 1
        IF (INCX .LT. 0) KX = 1+(1-N)*INCX
        IF (INCY .LT. 0) KY = 1+(1-N)*INCY
        DO 10 I = 1,N
          DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
          KX = KX + INCX
          KY = KY + INCY
        10 CONTINUE
      END IF
    END IF
    SDSDOT = DSDOT
  END FUNCTION SDSDOT




  !***BEGIN PROLOGUE  SGEFA
  !***PURPOSE  Factor a matrix using Gaussian elimination.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2A1
  !***TYPE      SINGLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
  !***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
  !   MATRIX FACTORIZATION
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !SGEFA factors a real matrix by Gaussian elimination.
  !
  !SGEFA is usually called by SGECO, but it can be called
  !directly with a saving in time if  RCOND  is not needed.
  !(Time for SGECO) = (1 + 9/N)*(Time for SGEFA) .
  !
  !On Entry
  !
  ! A       DOUBLE PRECISION(LDA, N)
  !      the matrix to be factored.
  !
  ! LDA     INTEGER
  !      the leading dimension of the array  A .
  !
  ! N       INTEGER
  !      the order of the matrix  A .
  !
  !On Return
  !
  ! A       an upper triangular matrix and the multipliers
  !      which were used to obtain it.
  !      The factorization can be written  A = L*U , where
  !      L  is a product of permutation and unit lower
  !      triangular matrices and  U  is upper triangular.
  !
  ! IPVT    INTEGER(N)
  !      an INTEGER :: vector of pivot indices.
  !
  ! INFO    INTEGER
  !      = 0  normal value.
  !      = K  if  U(K,K) .EQ. 0.0 .  This is not an error
  !    condition for this subroutine, but it does
  !    indicate that SGESL or SGEDI will divide by zero
  !    if called.  Use  RCOND  in SGECO for a reliable
  !    indication of singularity.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  ISAMAX, SAXPY, SSCAL
  !***REVISION HISTORY  (YYMMDD)
  ! 780814  DATE WRITTEN
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900326  Removed duplicate information from DESCRIPTION section.
  ! (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SGEFA
  SUBROUTINE SGEFA (A, LDA, N, IPVT, INFO)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    INTEGER :: LDA,N,IPVT(*),INFO
    REAL(dp) :: A(LDA,*)

    REAL(dp) :: T
    INTEGER :: J,K,KP1,L,NM1

    !GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
    !
    !***FIRST EXECUTABLE STATEMENT  SGEFA
    INFO = 0
    NM1 = N - 1
    IF (NM1 .LT. 1) THEN
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
    ELSE 
      DO K = 1, NM1
        KP1 = K + 1
  
        ! FIND L = PIVOT INDEX
        L = ISAMAX(N-K+1,A(K,K),1) + K - 1
        IPVT(K) = L
      
        ! ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
        IF (A(L,K) .EQ. 0.0D0) THEN
          INFO=K
        ELSE

          ! INTERCHANGE IF NECESSARY
          IF (L .NE. K) THEN
            T = A(L,K)
            A(L,K) = A(K,K)
            A(K,K) = T
          END IF
            
          ! COMPUTE MULTIPLIERS
          T = -1.0D0/A(K,K)
          CALL SSCAL(N-K,T,A(K+1,K),1)

          ! ROW ELIMINATION WITH COLUMN INDEXING
          DO J = KP1, N
            T = A(L,J)
            IF (L .NE. K) THEN
              A(L,J) = A(K,J)
              A(K,J) = T
            END IF
            CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
          END DO
        END IF
      END DO
    END IF
  END SUBROUTINE SGEFA



  !***BEGIN PROLOGUE  ISAMAX
  !***PURPOSE  Find the smallest index of that component of a vector
  !  having the maximum magnitude.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A2
  !***TYPE      SINGLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  ! Hanson, R. J., (SNLA)
  ! Kincaid, D. R., (U. of Texas)
  ! Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !      B L A S  Subprogram
  !  Description of Parameters
  !
  !--Input--
  ! N  number of elements in input vector(s)
  !  SX  single precision vector with N elements
  !INCX  storage spacing between elements of SX
  !
  !--Output--
  ! ISAMAX  smallest index (zero if N .LE. 0)
  !
  !Find smallest index of maximum magnitude of single precision SX.
  !ISAMAX = first I, I = 1 to N, to maximize  ABS(SX(IX+(I-1)*INCX)),
  !where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !Krogh, Basic linear algebra subprograms for Fortran
  !usage, Algorithm No. 539, Transactions on Mathematical
  !Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 791001  DATE WRITTEN
  ! 861211  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900821  Modified to correct problem with a negative increment.
  ! (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  ! 920618  Slight restructuring of code.  (RWC, WRB)
  !***END PROLOGUE  ISAMAX
  INTEGER FUNCTION ISAMAX (N, SX, INCX)
    !IMPLICIT INTEGER :: (I-N)
    !IMPLICIT REAL(dp) :: (A-H,O-Z)
    REAL(dp) :: SX(*), SMAX, XMAG
    INTEGER :: I, INCX, IX, N

    ISAMAX = 0
    IF (N .LE. 0) RETURN
    ISAMAX = 1
    IF (N .EQ. 1) RETURN

    IF (INCX .NE. 1) THEN
    !
    !Code for increment not equal to 1.
    !
    IX = 1
    IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      SMAX = DABS(SX(IX))
      IX = IX + INCX
      DO I = 2,N
        XMAG = DABS(SX(IX))
        IF (XMAG .GT. SMAX) THEN
          ISAMAX = I
          SMAX = XMAG
        ENDIF
        IX = IX + INCX
      END DO
    ELSE
      !
      !Code for increments equal to 1.
      !
      SMAX = DABS(SX(1))
      DO I = 2,N
        XMAG = DABS(SX(I))
        IF (XMAG .GT. SMAX) THEN
          ISAMAX = I
          SMAX = XMAG
        ENDIF
      END DO
    END IF
  END FUNCTION ISAMAX

  !***BEGIN PROLOGUE  SAXPY
  !***PURPOSE  Compute a constant times a vector plus a vector.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A7
  !***TYPE      SINGLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  ! Hanson, R. J., (SNLA)
  ! Kincaid, D. R., (U. of Texas)
  ! Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !      B L A S  Subprogram
  !  Description of Parameters
  !
  !--Input--
  ! N  number of elements in input vector(s)
  !  SA  single precision scalar multiplier
  !  SX  single precision vector with N elements
  !INCX  storage spacing between elements of SX
  !  SY  single precision vector with N elements
  !INCY  storage spacing between elements of SY
  !
  !--Output--
  !  SY  single precision result (unchanged if N .LE. 0)
  !
  !Overwrite single precision SY with single precision SA*SX +SY.
  !For I = 0 to N-1, replace  SY(LY+I*INCY) with SA*SX(LX+I*INCX) +
  !  SY(LY+I*INCY),
  !where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !defined in a similar way using INCY.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !Krogh, Basic linear algebra subprograms for Fortran
  !usage, Algorithm No. 539, Transactions on Mathematical
  !Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 791001  DATE WRITTEN
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SAXPY
  SUBROUTINE SAXPY (N, SA, SX, INCX, SY, INCY)
    !IMPLICIT INTEGER :: (I-N)
    !IMPLICIT REAL(dp) :: (A-H,O-Z)
    
    REAL(dp) :: SX(*), SY(*), SA
    !***FIRST EXECUTABLE STATEMENT  SAXPY
    IF (N.LE.0 .OR. SA.EQ.0.0D0) RETURN
    IF (((INCX-1) .lt. 0.0) .or. (INCX .NE. INCY))  THEN
      !
      !Code for unequal or nonpositive increments.
      !
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
      END DO
    ELSE IF ((INCX-1).eq. 0.0) THEN
    !
    !Code for both increments equal to 1.
    !
    !Clean-up loop so remaining vector length is a multiple of 4.
    !
      M = MOD(N,4)
      IF (M .NE. 0) THEN
        DO I = 1,M
          SY(I) = SY(I) + SA*SX(I)
        END DO
        IF (N .LT. 4) RETURN
      END IF
      MP1 = M + 1
        DO I = MP1,N,4
          SY(I) = SY(I) + SA*SX(I)
          SY(I+1) = SY(I+1) + SA*SX(I+1)
          SY(I+2) = SY(I+2) + SA*SX(I+2)
          SY(I+3) = SY(I+3) + SA*SX(I+3)
        END DO
    ELSE
      !
      !Code for equal, positive, non-unit increments.
      !
      NS = N*INCX
      DO I = 1,NS,INCX
        SY(I) = SA*SX(I) + SY(I)
      END DO
    END IF
  END SUBROUTINE SAXPY

  !***BEGIN PROLOGUE  SSCAL
  !***PURPOSE  Multiply a vector by a constant.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A6
  !***TYPE      SINGLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  ! Hanson, R. J., (SNLA)
  ! Kincaid, D. R., (U. of Texas)
  ! Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !      B L A S  Subprogram
  !  Description of Parameters
  !
  !--Input--
  ! N  number of elements in input vector(s)
  !  SA  single precision scale factor
  !  SX  single precision vector with N elements
  !INCX  storage spacing between elements of SX
  !
  !--Output--
  !  SX  single precision result (unchanged if N .LE. 0)
  !
  !Replace single precision SX by single precision SA*SX.
  !For I = 0 to N-1, replace SX(IX+I*INCX) with  SA * SX(IX+I*INCX),
  !where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !Krogh, Basic linear algebra subprograms for Fortran
  !usage, Algorithm No. 539, Transactions on Mathematical
  !Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 791001  DATE WRITTEN
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900821  Modified to correct problem with a negative increment.
  ! (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SSCAL
  SUBROUTINE SSCAL (N, SA, SX, INCX)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    REAL(dp) :: SA, SX(*)
    INTEGER :: I, INCX, IX, M, MP1, N

    IF (N .GT. 0) THEN
      IF (INCX .NE. 1) THEN
        !Code for increment not equal to 1.
        IX = 1
        IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
        DO I = 1,N
          SX(IX) = SA*SX(IX)
          IX = IX + INCX
        END DO
      !Code for increment equal to 1.
      !Clean-up loop so remaining vector length is a multiple of 5.
      ELSE
        M = MOD(N,5)
        IF (M .NE. 0) THEN
          DO I = 1,M
            SX(I) = SA*SX(I)
          END DO
          IF (N .LT. 5) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1,N,5
          SX(I) = SA*SX(I)
          SX(I+1) = SA*SX(I+1)
          SX(I+2) = SA*SX(I+2)
          SX(I+3) = SA*SX(I+3)
          SX(I+4) = SA*SX(I+4)
        END DO
      END IF
    END IF
    RETURN
  END SUBROUTINE SSCAL


  !***BEGIN PROLOGUE  SGESL
  !***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
  !  factors of SGECO or SGEFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2A1
  !***TYPE      SINGLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !SGESL solves the real system
  !A * X = B  or  TRANS(A) * X = B
  !using the factors computed by SGECO or SGEFA.
  !
  !On Entry
  !
  ! A       DOUBLE PRECISION(LDA, N)
  !      the output from SGECO or SGEFA.
  !
  ! LDA     INTEGER
  !      the leading dimension of the array  A .
  !
  ! N       INTEGER
  !      the order of the matrix  A .
  !
  ! IPVT    INTEGER(N)
  !      the pivot vector from SGECO or SGEFA.
  !
  ! B       DOUBLE PRECISION(N)
  !      the right hand side vector.
  !
  ! JOB     INTEGER
  !      = 0         to solve  A*X = B ,
  !      = nonzero   to solve  TRANS(A)*X = B  where
  !           TRANS(A)  is the transpose.
  !
  !On Return
  !
  ! B       the solution vector  X .
  !
  !Error Condition
  !
  ! A division by zero will occur if the input factor contains a
  ! zero on the diagonal.  Technically, this indicates singularity,
  ! but it is often caused by improper arguments or improper
  ! setting of LDA .  It will not occur if the subroutines are
  ! called correctly and if SGECO has set RCOND .GT. 0.0
  ! or SGEFA has set INFO .EQ. 0 .
  !
  !To compute  INVERSE(A) * C  where  C  is a matrix
  !with  P  columns
  ! CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
  ! IF (RCOND is too small) GO TO ...
  ! DO 10 J = 1, P
  !    CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
  ! 10 CONTINUE
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  SAXPY, SDOT
  !***REVISION HISTORY  (YYMMDD)
  ! 780814  DATE WRITTEN
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900326  Removed duplicate information from DESCRIPTION section.
  ! (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SGESL
  SUBROUTINE SGESL (A, LDA, N, IPVT, B, JOB)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    INTEGER :: LDA,N,IPVT(*),JOB
    REAL(dp) :: A(LDA,*),B(*)
    REAL(dp) :: T
    INTEGER :: K,KB,L,NM1

    NM1 = N - 1
    IF (JOB .EQ. 0) THEN
      ! JOB = 0 , SOLVE  A * X = B
      ! FIRST SOLVE  L*Y = B
      IF (NM1 .GE. 1) THEN
        DO K = 1, NM1
          L = IPVT(K)
          T = B(L)
          IF (L .NE. K) THEN
            B(L) = B(K)
            B(K) = T
          END IF
          CALL SAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
        END DO
      END IF

      ! NOW SOLVE  U*X = Y
      DO KB = 1, N
        K = N + 1 - KB
        B(K) = B(K)/A(K,K)
        T = -B(K)
        CALL SAXPY(K-1,T,A(1,K),1,B(1),1)
      END DO
    ELSE
      ! JOB = NONZERO, SOLVE  TRANS(A) * X = B
      ! FIRST SOLVE  TRANS(U)*Y = B
      DO K = 1, N
        T = SDOT(K-1,A(1,K),1,B(1),1)
        B(K) = (B(K) - T)/A(K,K)
      END DO

      ! NOW SOLVE TRANS(L)*X = Y
      IF (NM1 .GE. 1) THEN
        DO KB = 1, NM1
          K = N - KB
          B(K) = B(K) + SDOT(N-K,A(K+1,K),1,B(K+1),1)
          L = IPVT(K)
          IF (L .NE. K) THEN
            T = B(L)
            B(L) = B(K)
            B(K) = T
          END IF
        END DO
      END IF
    END IF
    RETURN
  END SUBROUTINE SGESL

  !***BEGIN PROLOGUE  SDOT
  !***PURPOSE  Compute the inner product of two vectors.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A4
  !***TYPE      SINGLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
  !***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  ! Hanson, R. J., (SNLA)
  ! Kincaid, D. R., (U. of Texas)
  ! Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !      B L A S  Subprogram
  !  Description of Parameters
  !
  !--Input--
  ! N  number of elements in input vector(s)
  !  SX  single precision vector with N elements
  !INCX  storage spacing between elements of SX
  !  SY  single precision vector with N elements
  !INCY  storage spacing between elements of SY
  !
  !--Output--
  !SDOT  single precision dot product (zero if N .LE. 0)
  !
  !Returns the dot product of single precision SX and SY.
  !SDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
  !where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !defined in a similar way using INCY.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !Krogh, Basic linear algebra subprograms for Fortran
  !usage, Algorithm No. 539, Transactions on Mathematical
  !Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 791001  DATE WRITTEN
  ! 890831  Modified array declarations.  (WRB)
  ! 890831  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SDOT
  REAL(dp)  FUNCTION SDOT (N, SX, INCX, SY, INCY)
  ! IMPLICIT INTEGER :: (I-N)
  ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    REAL(dp) :: SX(*), SY(*)
    SDOT = 0.0D0
    IF (N .GT. 0 ) THEN
      IF ((INCX .NE. INCY) .or. ((INCX-1) .lt. 0)) THEN
        !
        !Code for unequal or nonpositive increments.
        !
        IX = 1
        IY = 1
        IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
        IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
        DO I = 1,N
            SDOT = SDOT + SX(IX)*SY(IY)
            IX = IX + INCX
            IY = IY + INCY
        END DO
      ELSE IF ((INCX-1) .eq. 0) THEN
        !Code for both increments equal to 1.
        !Clean-up loop so remaining vector length is a multiple of 5.
        M = MOD(N,5)
        IF (M .NE. 0) THEN
          DO I = 1,M
            SDOT = SDOT + SX(I)*SY(I)
          END DO
          IF (N .LT. 5) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1,N,5
          SDOT = SDOT + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2) +&
                        &SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
        END DO
      ELSE
        !Code for equal, positive, non-unit increments.
        NS = N*INCX
        DO I = 1,NS,INCX
          SDOT = SDOT + SX(I)*SY(I)
        END DO
      END IF
    END IF
  END FUNCTION SDOT

  !***BEGIN PROLOGUE  XERMSG
  !***PURPOSE  Process error messages for SLATEC and other libraries.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERMSG-A)
  !***KEYWORDS  ERROR MESSAGE, XERROR
  !***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
  !***DESCRIPTION
  !
  ! XERMSG processes a diagnostic message in a manner determined by the
  ! value of LEVEL and the current value of the library error control
  ! flag, KONTRL.  See subroutine XSETF for details.
  !
  !  LIBRAR   A character constant (or character variable) with the name
  !   of the library.  This will be 'SLATEC' for the SLATEC
  !   Common Math Library.  The error handling package is
  !   general enough to be used by many libraries
  !   simultaneously, so it is desirable for the routine that
  !   detects and reports an error to identify the library name
  !   as well as the routine name.
  !
  !  SUBROU   A character constant (or character variable) with the name
  !   of the routine that detected the error.  Usually it is the
  !   name of the routine that is calling XERMSG.  There are
  !   some instances where a user callable library routine calls
  !   lower level subsidiary routines where the error is
  !   detected.  In such cases it may be more informative to
  !   supply the name of the routine the user called rather than
  !   the name of the subsidiary routine that detected the
  !   error.
  !
  !  MESSG    A character constant (or character variable) with the text
  !   of the error or warning message.  In the example below,
  !   the message is a character constant that contains a
  !   generic message.
  !
  !  CALL XERMSG ('SLATEC', 'MMPY',
  ! *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
  ! *3, 1)
  !
  !   It is possible (and is sometimes desirable) to generate a
  !   specific message--e.g., one that contains actual numeric
  !   values.  Specific numeric values can be converted into
  !   character strings using formatted WRITE statements into
  !   character variables.  This is called standard Fortran
  !   internal file I/O and is exemplified in the first three
  !   lines of the following example.  You can also catenate
  !   substrings of characters to construct the error message.
  !   Here is an example showing the use of both writing to
  !   an internal file and catenating character strings.
  !
  !  CHARACTER*5 CHARN, CHARL
  !  WRITE (CHARN,10) N
  !  WRITE (CHARL,10) LDA
  !      10 FORMAT(I5)
  !  CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
  ! *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
  ! *   CHARL, 3, 1)
  !
  !   There are two subtleties worth mentioning.  One is that
  !   the // for character catenation is used to construct the
  !   error message so that no single character constant is
  !   continued to the next line.  This avoids confusion as to
  !   whether there are trailing blanks at the end of the line.
  !   The second is that by catenating the parts of the message
  !   as an actual argument rather than encoding the entire
  !   message into one large character variable, we avoid
  !   having to know how long the message will be in order to
  !   declare an adequate length for that large character
  !   variable.  XERMSG calls XERPRN to print the message using
  !   multiple lines if necessary.  If the message is very long,
  !   XERPRN will break it into pieces of 72 characters (as
  !   requested by XERMSG) for printing on multiple lines.
  !   Also, XERMSG asks XERPRN to prefix each line with ' *  '
  !   so that the total line length could be 76 characters.
  !   Note also that XERPRN scans the error message backwards
  !   to ignore trailing blanks.  Another feature is that
  !   the substring '$$' is treated as a new line sentinel
  !   by XERPRN.  If you want to construct a multiline
  !   message without having to count out multiples of 72
  !   characters, just use '$$' as a separator.  '$$'
  !   obviously must occur within 72 characters of the
  !   start of each line to have its intended effect since
  !   XERPRN is asked to wrap around at 72 characters in
  !   addition to looking for '$$'.
  !
  !  NERR     An INTEGER :: value that is chosen by the library routine's
  !   author.  It must be in the range -99 to 999 (three
  !   printable digits).  Each distinct error should have its
  !   own error number.  These error numbers should be described
  !   in the machine readable documentation for the routine.
  !   The error numbers need be unique only within each routine,
  !   so it is reasonable for each routine to start enumerating
  !   errors from 1 and proceeding to the next integer.
  !
  !  LEVEL    An INTEGER :: value in the range 0 to 2 that indicates the
  !   level (severity) of the error.  Their meanings are
  !
  !  -1  A warning message.  This is used if it is not clear
  !      that there really is an error, but the user's attention
  !      may be needed.  An attempt is made to only print this
  !      message once.
  !
  !   0  A warning message.  This is used if it is not clear
  !      that there really is an error, but the user's attention
  !      may be needed.
  !
  !   1  A recoverable error.  This is used even if the error is
  !      so serious that the routine cannot return any useful
  !      answer.  If the user has told the error package to
  !      return after recoverable errors, then XERMSG will
  !      return to the Library routine which can then return to
  !      the user's routine.  The user may also permit the error
  !      package to terminate the program upon encountering a
  !      recoverable error.
  !
  !   2  A fatal error.  XERMSG will not return to its caller
  !      after it receives a fatal error.  This level should
  !      hardly ever be used; it is much better to allow the
  !      user a chance to recover.  An example of one of the few
  !      cases in which it is permissible to declare a level 2
  !      error is a reverse communication Library routine that
  !      is likely to be called repeatedly until it integrates
  !      across some interval.  If there is a serious error in
  !      the input such that another step cannot be taken and
  !      the Library routine is called again without the input
  !      error having been corrected by the caller, the Library
  !      routine will probably be called forever with improper
  !      input.  In this case, it is reasonable to declare the
  !      error to be fatal.
  !
  !  Each of the arguments to XERMSG is input; none will be modified by
  !  XERMSG.  A routine may make multiple calls to XERMSG with warning
  !  level messages; however, after a call to XERMSG with a recoverable
  !  error, the routine should return to the user.  Do not try to call
  !  XERMSG with a second recoverable error after the first recoverable
  !  error because the error package saves the error number.  The user
  !  can retrieve this error number by calling another entry point in
  !  the error handling package and then clear the error number when
  !  recovering from the error.  Calling XERMSG in succession causes the
  !  old error number to be overwritten by the latest error number.
  !  This is considered harmless for error numbers associated with
  !  warning messages but must not be done for error numbers of serious
  !  errors.  After a call to XERMSG with a recoverable error, the user
  !  must be given a chance to call NUMXER or XERCLR to retrieve or
  !  clear the error number.
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !Error-handling Package, SAND82-0800, Sandia
  !Laboratories, 1982.
  !***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
  !***REVISION HISTORY  (YYMMDD)
  ! 880101  DATE WRITTEN
  ! 880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
  ! THERE ARE TWO BASIC CHANGES.
  ! 1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
  !     PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
  !     INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
  !     ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
  !     ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
  !     ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
  !     72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
  !     LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
  ! 2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
  !     FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
  !     OF LOWER CASE.
  ! 880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
  ! THE PRINCIPAL CHANGES ARE
  ! 1.  CLARIFY COMMENTS IN THE PROLOGUES
  ! 2.  RENAME XRPRNT TO XERPRN
  ! 3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
  !     SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
  !     CHARACTER FOR NEW RECORDS.
  ! 890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
  ! CLEAN UP THE CODING.
  ! 890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
  ! PREFIX.
  ! 891013  REVISED TO CORRECT COMMENTS.
  ! 891214  Prologue converted to Version 4.0 format.  (WRB)
  ! 900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
  ! NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
  ! LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
  ! XERCTL to XERCNT.  (RWC)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERMSG

  SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
  ! IMPLICIT INTEGER :: (I-N)
  ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    CHARACTER*(*) LIBRAR, SUBROU, MESSG
    CHARACTER*8 XLIBR, XSUBR
    CHARACTER*72  TEMP
    CHARACTER*20  LFIRST
    LKNTRL = J4SAVE (2, 0, .FALSE.)
    MAXMES = J4SAVE (4, 0, .FALSE.)
    !  LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
    !  MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
    !   SHOULD BE PRINTED.
    !
    !  WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
    !   CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
    !   AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.

      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.&
     &   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //&
           &'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//&
          & 'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
!
!  RECORD THE MESSAGE.
!
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!  HANDLE PRINT-ONCE WARNING MESSAGES.
!
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
!
!  ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
!
!  SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!  ZERO AND THE ERROR IS NOT FATAL.
!
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
!
!  ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!  MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!  AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!  IS NOT ZERO.
!
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!  IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!  PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!  FROM EACH OF THE FOLLOWING THREE OPTIONS.
!  1.  LEVEL OF THE MESSAGE
!    'INFORMATIVE MESSAGE'
!    'POTENTIALLY RECOVERABLE ERROR'
!    'FATAL ERROR'
!  2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!    'PROG CONTINUES'
!    'PROG ABORTED'
!  3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
! MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
! WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!    'TRACEBACK REQUESTED'
!    'TRACEBACK NOT REQUESTED'
!  NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!  EXCEED 74 CHARACTERS.
!  WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
!
      IF (LKNTRL .GT. 0) THEN
!
!  THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
!
!  THEN WHETHER THE PROGRAM WILL CONTINUE.
!
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.&
     &       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
!
!  FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!  NOW SEND OUT THE MESSAGE.
!
      CALL XERPRN (' *  ', -1, MESSG, 72)
!
!  IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!   TRACEBACK.
!
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
!
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
!
!  IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
!
!  IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!  CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
!
!  THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!  FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!  SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN&
     &         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
  END SUBROUTINE XERMSG

  !***BEGIN PROLOGUE  J4SAVE
  !***SUBSIDIARY
  !***PURPOSE  Save or recall global variables needed by error
  !  handling routines.
  !***LIBRARY   SLATEC (XERROR)
  !***TYPE      INTEGER :: (J4SAVE-I)
  !***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !Abstract
  ! J4SAVE saves and recalls several global variables needed
  ! by the library error handling routines.
  !
  !Description of Parameters
  ! --Input--
  ! IWHICH - Index of item desired.
  !      = 1 Refers to current error number.
  !      = 2 Refers to current error control flag.
  !      = 3 Refers to current unit number to which error
  !   messages are to be sent.  (0 means use standard.)
  !      = 4 Refers to the maximum number of times any
  !    message is to be printed (as set by XERMAX).
  !      = 5 Refers to the total number of units to which
  !    each error message is to be written.
  !      = 6 Refers to the 2nd unit for error messages
  !      = 7 Refers to the 3rd unit for error messages
  !      = 8 Refers to the 4th unit for error messages
  !      = 9 Refers to the 5th unit for error messages
  ! IVALUE - The value to be set for the IWHICH-th parameter,
  !if ISET is .TRUE. .
  ! ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
  !given the value, IVALUE.  If ISET=.FALSE., the
  !IWHICH-th parameter will be unchanged, and IVALUE
  !is a dummy parameter.
  ! --Output--
  ! The (old) value of the IWHICH-th parameter will be returned
  ! in the function value, J4SAVE.
  !
  !***SEE ALSO  XERMSG
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !Error-handling Package, SAND82-0800, Sandia
  !Laboratories, 1982.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 790801  DATE WRITTEN
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900205  Minor modifications to prologue.  (WRB)
  ! 900402  Added TYPE section.  (WRB)
  ! 910411  Added KEYWORDS section.  (WRB)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  J4SAVE
  FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)

      LOGICAL ISET
      INTEGER :: IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
  END FUNCTION J4SAVE


  !***BEGIN PROLOGUE  XERCNT
  !***SUBSIDIARY
  !***PURPOSE  Allow user control over handling of errors.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERCNT-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !Abstract
  ! Allows user control over handling of individual errors.
  ! Just after each message is recorded, but before it is
  ! processed any further (i.e., before it is printed or
  ! a decision to abort is made), a call is made to XERCNT.
  ! If the user has provided his own version of XERCNT, he
  ! can then override the value of KONTROL used in processing
  ! this message by redefining its value.
  ! KONTRL may be set to any value from -2 to 2.
  ! The meanings for KONTRL are the same as in XSETF, except
  ! that the value of KONTRL changes only for this message.
  ! If KONTRL is set to a value outside the range from -2 to 2,
  ! it will be moved back into that range.
  !
  !Description of Parameters
  !
  ! --Input--
  ! LIBRAR - the library that the routine is in.
  ! SUBROU - the subroutine that XERMSG is being called from
  ! MESSG  - the first 20 characters of the error message.
  ! NERR   - same as in the call to XERMSG.
  ! LEVEL  - same as in the call to XERMSG.
  ! KONTRL - the current value of the control flag as set
  !by a call to XSETF.
  !
  ! --Output--
  ! KONTRL - the new value of KONTRL.  If KONTRL is not
  !defined, it will remain at its original value.
  !This changed value of control affects only
  !the current occurrence of the current message.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !Error-handling Package, SAND82-0800, Sandia
  !Laboratories, 1982.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 790801  DATE WRITTEN
  ! 861211  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900206  Routine changed from user-callable to subsidiary.  (WRB)
  ! 900510  Changed calling sequence to include LIBRARY and SUBROUTINE
  ! names, changed routine name from XERCTL to XERCNT.  (RWC)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERCNT
  SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)
    CHARACTER*(*) LIBRAR, SUBROU, MESSG
    RETURN
  END SUBROUTINE XERCNT

  !***BEGIN PROLOGUE  XERHLT
  !***SUBSIDIARY
  !***PURPOSE  Abort program execution and print error message.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERHLT-A)
  !***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !Abstract
  ! ***Note*** machine dependent routine
  ! XERHLT aborts the execution of the program.
  ! The error message causing the abort is given in the calling
  ! sequence, in case one needs it for printing on a dayfile,
  ! for example.
  !
  !Description of Parameters
  ! MESSG is as in XERMSG.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !Error-handling Package, SAND82-0800, Sandia
  !Laboratories, 1982.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 790801  DATE WRITTEN
  ! 861211  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! 900206  Routine changed from user-callable to subsidiary.  (WRB)
  ! 900510  Changed calling sequence to delete length of character
  ! and changed routine name from XERABT to XERHLT.  (RWC)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERHLT
  SUBROUTINE XERHLT (MESSG)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)
    CHARACTER*(*) MESSG
    STOP
  END SUBROUTINE XERHLT

  !***BEGIN PROLOGUE  XERPRN
  !***SUBSIDIARY
  !***PURPOSE  Print error messages processed by XERMSG.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERPRN-A)
  !***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
  !***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
  !***DESCRIPTION
  !
  ! This routine sends one or more lines to each of the (up to five)
  ! logical units to which error messages are to be sent.  This routine
  ! is called several times by XERMSG, sometimes with a single line to
  ! print and sometimes with a (potentially very long) message that may
  ! wrap around into multiple lines.
  !
  ! PREFIX  Input argument of type CHARACTER.  This argument contains
  !  characters to be put at the beginning of each line before
  !  the body of the message.  No more than 16 characters of
  !  PREFIX will be used.
  !
  ! NPREF   Input argument of type INTEGER.  This argument is the number
  !  of characters to use from PREFIX.  If it is negative, the
  !  intrinsic function LEN is used to determine its length.  If
  !  it is zero, PREFIX is not used.  If it exceeds 16 or if
  !  LEN(PREFIX) exceeds 16, only the first 16 characters will be
  !  used.  If NPREF is positive and the length of PREFIX is less
  !  than NPREF, a copy of PREFIX extended with blanks to length
  !  NPREF will be used.
  !
  ! MESSG   Input argument of type CHARACTER.  This is the text of a
  !  message to be printed.  If it is a long message, it will be
  !  broken into pieces for printing on multiple lines.  Each line
  !  will start with the appropriate prefix and be followed by a
  !  piece of the message.  NWRAP is the number of characters per
  !  piece; that is, after each NWRAP characters, we break and
  !  start a new line.  In addition the characters '$$' embedded
  !  in MESSG are a sentinel for a new line.  The counting of
  !  characters up to NWRAP starts over for each new line.  The
  !  value of NWRAP typically used by XERMSG is 72 since many
  !  older error messages in the SLATEC Library are laid out to
  !  rely on wrap-around every 72 characters.
  !
  ! NWRAP   Input argument of type INTEGER.  This gives the maximum size
  !  piece into which to break MESSG for printing on multiple
  !  lines.  An embedded '$$' ends a line, and the count restarts
  !  at the following character.  If a line break does not occur
  !  on a blank (it would split a word) that word is moved to the
  !  next line.  Values of NWRAP less than 16 will be treated as
  !  16.  Values of NWRAP greater than 132 will be treated as 132.
  !  The actual line length will be NPREF + NWRAP after NPREF has
  !  been adjusted to fall between 0 and 16 and NWRAP has been
  !  adjusted to fall between 16 and 132.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !Error-handling Package, SAND82-0800, Sandia
  !Laboratories, 1982.
  !***ROUTINES CALLED  I1MACH, XGETUA
  !***REVISION HISTORY  (YYMMDD)
  ! 880621  DATE WRITTEN
  ! 880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
  ! JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
  ! THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
  ! SLASH CHARACTER IN FORMAT STATEMENTS.
  ! 890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
  ! STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
  ! LINES TO BE PRINTED.
  ! 890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
  ! CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
  ! 891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
  ! 891214  Prologue converted to Version 4.0 format.  (WRB)
  ! 900510  Added code to break messages between words.  (RWC)
  ! 920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERPRN
  SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
    ! IMPLICIT INTEGER :: (I-N)
    ! IMPLICIT REAL(dp) :: (A-H,O-Z)

    CHARACTER(*) :: PREFIX, MESSG
    INTEGER :: NPREF, NWRAP
    CHARACTER(148) :: CBUFF
    INTEGER :: IU(5), NUNIT
    CHARACTER(2), PARAMETER :: NEWLIN ='$$'

    CALL XGETUA(IU,NUNIT)
    
    !  A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
    !  ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
    !  ERROR MESSAGE UNIT.
    N = I1MACH(4)
    DO I=1,NUNIT
       IF (IU(I) .EQ. 0) IU(I) = N
    END DO

    !  LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
    !  BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
    !  THE REST OF THIS ROUTINE.
    IF ( NPREF .LT. 0 ) THEN
       LPREF = LEN(PREFIX)
    ELSE
       LPREF = NPREF
    ENDIF
    LPREF = MIN(16, LPREF)
    IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX

    !LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
    !TIME FROM MESSG TO PRINT ON ONE LINE.
    LWRAP = MAX(16, MIN(132, NWRAP))

    !  SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
    LENMSG = LEN(MESSG)
    N = LENMSG
    DO I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') EXIT
         LENMSG = LENMSG - 1
    END DO
   

    !IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
    IF (LENMSG .EQ. 0) THEN
      CBUFF(LPREF+1:LPREF+1) = ' '
      DO I=1,NUNIT
        WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
      END DO
      RETURN
    ENDIF

    !  SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
    !  STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
    !  WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
    !  WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
    !
    !  WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
    !  INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
    !  OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
    !  OF THE SECOND ARGUMENT.
    !
    !  THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
    !  FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
    !  OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
    !  POSITION NEXTC.
    !
    !  LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
    !      REMAINDER OF THE CHARACTER STRING.  LPIECE
    !      SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
    !      WHICHEVER IS LESS.
    !
    !  LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
    !      NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
    !      PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
    !      BLANK LINES.  THIS TAKES CARE OF THE SITUATION
    !      WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
    !      EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
    !      SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
    !      SHOULD BE INCREMENTED BY 2.
    !
    !  LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
    !
    !  ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
    !      RESET LPIECE = LPIECE-1.  NOTE THAT THIS
    !      PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
    !      LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
    !      AT THE END OF A LINE.
    
    NEXTC = 1
    DO WHILE (NEXTC .LE. LENMSG)
      LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
        !
        !  THERE WAS NO NEW LINE SENTINEL FOUND.
        !
        IDELTA = 0
        LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
        IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
          DO I=LPIECE+1,2,-1
             IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                LPIECE = I-1
                IDELTA = 1
                EXIT
             ENDIF
          END DO
        ENDIF
        CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
        NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
        !
        !  WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
        !  DON'T PRINT A BLANK LINE.
        !
        NEXTC = NEXTC + 2
        CYCLE
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
      !
      !  LPIECE SHOULD BE SET DOWN TO LWRAP.
      !
        IDELTA = 0
        LPIECE = LWRAP
        DO I=LPIECE+1,2,-1
          IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
             LPIECE = I-1
             IDELTA = 1
             EXIT
          ENDIF
        END DO
        CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
        NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
        !
        !  IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
        !  WE SHOULD DECREMENT LPIECE BY ONE.
        !
        LPIECE = LPIECE - 1
        CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
        NEXTC  = NEXTC + LPIECE + 2
      ENDIF
      !
      !  PRINT
      !
      DO I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
      END DO
      !
    END DO
    RETURN
  END SUBROUTINE XERPRN


! !***BEGIN PROLOGUE  I1MACH
! !***PURPOSE  Return INTEGER :: machine dependent constants.
! !***LIBRARY   SLATEC
! !***CATEGORY  R1
! !***TYPE      INTEGER :: (I1MACH-I)
! !***KEYWORDS  MACHINE CONSTANTS
! !***AUTHOR  Fox, P. A., (Bell Labs)
! ! Hall, A. D., (Bell Labs)
! ! Schryer, N. L., (Bell Labs)
! !***DESCRIPTION
! !
! ! I1MACH can be used to obtain machine-dependent parameters for the
! ! local machine environment.  It is a function subprogram with one
! ! (input) argument and can be referenced as follows:
! !
! ! K = I1MACH(I)
! !
! ! where I=1,...,16.  The (output) value of K above is determined by
! ! the (input) value of I.  The results for various values of I are
! ! discussed below.
! !
! ! I/O unit numbers:
! !I1MACH( 1) = the standard input unit.
! !I1MACH( 2) = the standard output unit.
! !I1MACH( 3) = the standard punch unit.
! !I1MACH( 4) = the standard error message unit.
! !
! ! Words:
! !I1MACH( 5) = the number of bits per INTEGER :: storage unit.
! !I1MACH( 6) = the number of characters per INTEGER :: storage unit.
! !
! ! Integers:
! !assume integers are represented in the S-digit, base-A form
! !
! !      sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
! !
! !      where 0 .LE. X(I) .LT. A for I=0,...,S-1.
! !I1MACH( 7) = A, the base.
! !I1MACH( 8) = S, the number of base-A digits.
! !I1MACH( 9) = A**S - 1, the largest magnitude.
! !
! ! Floating-Point Numbers:
! !Assume floating-point numbers are represented in the T-digit,
! !base-B form
! !      sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
! !
! !      where 0 .LE. X(I) .LT. B for I=1,...,T,
! !      0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
! !I1MACH(10) = B, the base.
! !
! ! Single-Precision:
! !I1MACH(11) = T, the number of base-B digits.
! !I1MACH(12) = EMIN, the smallest exponent E.
! !I1MACH(13) = EMAX, the largest exponent E.
! !
! ! Double-Precision:
! !I1MACH(14) = T, the number of base-B digits.
! !I1MACH(15) = EMIN, the smallest exponent E.
! !I1MACH(16) = EMAX, the largest exponent E.
! !
! ! To alter this function for a particular environment, the desired
! ! set of DATA statements should be activated by removing the C from
! ! column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
! ! checked for consistency with the local operating system.
! !
! !***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
! !a portable library, ACM Transactions on Mathematical
! !Software 4, 2 (June 1978), pp. 177-188.
! !***ROUTINES CALLED  (NONE)
! !***REVISION HISTORY  (YYMMDD)
! ! 750101  DATE WRITTEN
! ! 891012  Added VAX G-floating constants.  (WRB)
! ! 891012  REVISION DATE from Version 3.2
! ! 891214  Prologue converted to Version 4.0 format.  (BAB)
! ! 900618  Added DEC RISC constants.  (WRB)
! ! 900723  Added IBM RS 6000 constants.  (WRB)
! ! 901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
! ! (RWC)
! ! 910710  Added HP 730 constants.  (SMR)
! ! 911114  Added Convex IEEE constants.  (WRB)
! ! 920121  Added SUN -r8 compiler option constants.  (WRB)
! ! 920229  Added Touchstone Delta i860 constants.  (WRB)
! ! 920501  Reformatted the REFERENCES section.  (WRB)
! ! 920625  Added Convex -p8 and -pd8 compiler option constants.
! ! (BKS, WRB)
! ! 930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
! ! 930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
! ! options.  (DWL, RWC and WRB).
! !***END PROLOGUE  I1MACH
! !
! !
! !MACHINE CONSTANTS FOR THE AMIGA
! !ABSOFT COMPILER
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -126 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1022 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE APOLLO
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        129 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1025 /
! !
! !MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
! !
! !DATA IMACH( 1) /          7 /
! !DATA IMACH( 2) /          2 /
! !DATA IMACH( 3) /          2 /
! !DATA IMACH( 4) /          2 /
! !DATA IMACH( 5) /         36 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         33 /
! !DATA IMACH( 9) / Z1FFFFFFFF /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -256 /
! !DATA IMACH(13) /        255 /
! !DATA IMACH(14) /         60 /
! !DATA IMACH(15) /       -256 /
! !DATA IMACH(16) /        255 /
! !
! !MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         48 /
! !DATA IMACH( 6) /          6 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         39 /
! !DATA IMACH( 9) / O0007777777777777 /
! !DATA IMACH(10) /          8 /
! !DATA IMACH(11) /         13 /
! !DATA IMACH(12) /        -50 /
! !DATA IMACH(13) /         76 /
! !DATA IMACH(14) /         26 /
! !DATA IMACH(15) /        -50 /
! !DATA IMACH(16) /         76 /
! !
! !MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         48 /
! !DATA IMACH( 6) /          6 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         39 /
! !DATA IMACH( 9) / O0007777777777777 /
! !DATA IMACH(10) /          8 /
! !DATA IMACH(11) /         13 /
! !DATA IMACH(12) /        -50 /
! !DATA IMACH(13) /         76 /
! !DATA IMACH(14) /         26 /
! !DATA IMACH(15) /     -32754 /
! !DATA IMACH(16) /      32780 /
! !
! !MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         64 /
! !DATA IMACH( 6) /          8 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         63 /
! !DATA IMACH( 9) / 9223372036854775807 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         47 /
! !DATA IMACH(12) /      -4095 /
! !DATA IMACH(13) /       4094 /
! !DATA IMACH(14) /         94 /
! !DATA IMACH(15) /      -4095 /
! !DATA IMACH(16) /       4094 /
! !
! !MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /    6LOUTPUT/
! !DATA IMACH( 5) /         60 /
! !DATA IMACH( 6) /         10 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         48 /
! !DATA IMACH( 9) / 00007777777777777777B /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         47 /
! !DATA IMACH(12) /       -929 /
! !DATA IMACH(13) /       1070 /
! !DATA IMACH(14) /         94 /
! !DATA IMACH(15) /       -929 /
! !DATA IMACH(16) /       1069 /
! !
! !MACHINE CONSTANTS FOR THE CELERITY C1260
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          0 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / Z'7FFFFFFF' /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -126 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1022 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE CONVEX
! !USING THE -fn COMPILER OPTION
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1023 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE CONVEX
! !USING THE -fi COMPILER OPTION
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE CONVEX
! !USING THE -p8 COMPILER OPTION
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         64 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         63 /
! !DATA IMACH( 9) / 9223372036854775807 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         53 /
! !DATA IMACH(12) /      -1023 /
! !DATA IMACH(13) /       1023 /
! !DATA IMACH(14) /        113 /
! !DATA IMACH(15) /     -16383 /
! !DATA IMACH(16) /      16383 /
! !
! !MACHINE CONSTANTS FOR THE CONVEX
! !USING THE -pd8 COMPILER OPTION
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         64 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         63 /
! !DATA IMACH( 9) / 9223372036854775807 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         53 /
! !DATA IMACH(12) /      -1023 /
! !DATA IMACH(13) /       1023 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1023 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE CRAY
! !USING THE 46 BIT INTEGER :: COMPILER OPTION
! !
! !DATA IMACH( 1) /        100 /
! !DATA IMACH( 2) /        101 /
! !DATA IMACH( 3) /        102 /
! !DATA IMACH( 4) /        101 /
! !DATA IMACH( 5) /         64 /
! !DATA IMACH( 6) /          8 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         46 /
! !DATA IMACH( 9) / 1777777777777777B /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         47 /
! !DATA IMACH(12) /      -8189 /
! !DATA IMACH(13) /       8190 /
! !DATA IMACH(14) /         94 /
! !DATA IMACH(15) /      -8099 /
! !DATA IMACH(16) /       8190 /
! !
! !MACHINE CONSTANTS FOR THE CRAY
! !USING THE 64 BIT INTEGER :: COMPILER OPTION
! !
! !DATA IMACH( 1) /        100 /
! !DATA IMACH( 2) /        101 /
! !DATA IMACH( 3) /        102 /
! !DATA IMACH( 4) /        101 /
! !DATA IMACH( 5) /         64 /
! !DATA IMACH( 6) /          8 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         63 /
! !DATA IMACH( 9) / 777777777777777777777B /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         47 /
! !DATA IMACH(12) /      -8189 /
! !DATA IMACH(13) /       8190 /
! !DATA IMACH(14) /         94 /
! !DATA IMACH(15) /      -8099 /
! !DATA IMACH(16) /       8190 /
! !
! !MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
! !
! !DATA IMACH( 1) /         11 /
! !DATA IMACH( 2) /         12 /
! !DATA IMACH( 3) /          8 /
! !DATA IMACH( 4) /         10 /
! !DATA IMACH( 5) /         16 /
! !DATA IMACH( 6) /          2 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         15 /
! !DATA IMACH( 9) /      32767 /
! !DATA IMACH(10) /         16 /
! !DATA IMACH(11) /          6 /
! !DATA IMACH(12) /        -64 /
! !DATA IMACH(13) /         63 /
! !DATA IMACH(14) /         14 /
! !DATA IMACH(15) /        -64 /
! !DATA IMACH(16) /         63 /
! !
! !MACHINE CONSTANTS FOR THE DEC ALPHA
! !USING G_FLOAT
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1023 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE DEC ALPHA
! !USING IEEE_FLOAT
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE DEC RISC
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE DEC VAX
! !USING D_FLOATING
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         56 /
! !DATA IMACH(15) /       -127 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR THE DEC VAX
! !USING G_FLOATING
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1023 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE ELXSI 6400
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         32 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -126 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1022 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE HARRIS 220
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          0 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         24 /
! !DATA IMACH( 6) /          3 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         23 /
! !DATA IMACH( 9) /    8388607 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         23 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         38 /
! !DATA IMACH(15) /       -127 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /         43 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         36 /
! !DATA IMACH( 6) /          6 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         35 /
! !DATA IMACH( 9) / O377777777777 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         27 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         63 /
! !DATA IMACH(15) /       -127 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR THE HP 730
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE HP 2100
! !3 WORD REAL(dp) :: OPTION WITH FTN4
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          4 /
! !DATA IMACH( 4) /          1 /
! !DATA IMACH( 5) /         16 /
! !DATA IMACH( 6) /          2 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         15 /
! !DATA IMACH( 9) /      32767 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         23 /
! !DATA IMACH(12) /       -128 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         39 /
! !DATA IMACH(15) /       -128 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR THE HP 2100
! !4 WORD REAL(dp) :: OPTION WITH FTN4
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          4 /
! !DATA IMACH( 4) /          1 /
! !DATA IMACH( 5) /         16 /
! !DATA IMACH( 6) /          2 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         15 /
! !DATA IMACH( 9) /      32767 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         23 /
! !DATA IMACH(12) /       -128 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         55 /
! !DATA IMACH(15) /       -128 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR THE HP 9000
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          7 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         32 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -126 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1015 /
! !DATA IMACH(16) /       1017 /
! !
! !MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
! !THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
! !THE PERKIN ELMER (INTERDATA) 7/32.
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          7 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) /  Z7FFFFFFF /
! !DATA IMACH(10) /         16 /
! !DATA IMACH(11) /          6 /
! !DATA IMACH(12) /        -64 /
! !DATA IMACH(13) /         63 /
! !DATA IMACH(14) /         14 /
! !DATA IMACH(15) /        -64 /
! !DATA IMACH(16) /         63 /
! !
! !MACHINE CONSTANTS FOR THE IBM PC
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          0 /
! !DATA IMACH( 4) /          0 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE IBM RS 6000
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          0 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE INTEL i860
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         36 /
! !DATA IMACH( 6) /          5 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         35 /
! !DATA IMACH( 9) / "377777777777 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         27 /
! !DATA IMACH(12) /       -128 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         54 /
! !DATA IMACH(15) /       -101 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         36 /
! !DATA IMACH( 6) /          5 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         35 /
! !DATA IMACH( 9) / "377777777777 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         27 /
! !DATA IMACH(12) /       -128 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         62 /
! !DATA IMACH(15) /       -128 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
! !32-BIT INTEGER :: ARITHMETIC.
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         56 /
! !DATA IMACH(15) /       -127 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
! !16-BIT INTEGER :: ARITHMETIC.
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          5 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         16 /
! !DATA IMACH( 6) /          2 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         15 /
! !DATA IMACH( 9) /      32767 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         56 /
! !DATA IMACH(15) /       -127 /
! !DATA IMACH(16) /        127 /
! !
! !MACHINE CONSTANTS FOR THE SILICON GRAPHICS
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE SUN
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -125 /
! !DATA IMACH(13) /        128 /
! !DATA IMACH(14) /         53 /
! !DATA IMACH(15) /      -1021 /
! !DATA IMACH(16) /       1024 /
! !
! !MACHINE CONSTANTS FOR THE SUN
! !USING THE -r8 COMPILER OPTION
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          6 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         32 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         31 /
! !DATA IMACH( 9) / 2147483647 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         53 /
! !DATA IMACH(12) /      -1021 /
! !DATA IMACH(13) /       1024 /
! !DATA IMACH(14) /        113 /
! !DATA IMACH(15) /     -16381 /
! !DATA IMACH(16) /      16384 /
! !
! !MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
! !
! !DATA IMACH( 1) /          5 /
! !DATA IMACH( 2) /          6 /
! !DATA IMACH( 3) /          1 /
! !DATA IMACH( 4) /          6 /
! !DATA IMACH( 5) /         36 /
! !DATA IMACH( 6) /          4 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         35 /
! !DATA IMACH( 9) / O377777777777 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         27 /
! !DATA IMACH(12) /       -128 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         60 /
! !DATA IMACH(15) /      -1024 /
! !DATA IMACH(16) /       1023 /
! !
! !MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
! !
! !DATA IMACH( 1) /          1 /
! !DATA IMACH( 2) /          1 /
! !DATA IMACH( 3) /          0 /
! !DATA IMACH( 4) /          1 /
! !DATA IMACH( 5) /         16 /
! !DATA IMACH( 6) /          2 /
! !DATA IMACH( 7) /          2 /
! !DATA IMACH( 8) /         15 /
! !DATA IMACH( 9) /      32767 /
! !DATA IMACH(10) /          2 /
! !DATA IMACH(11) /         24 /
! !DATA IMACH(12) /       -127 /
! !DATA IMACH(13) /        127 /
! !DATA IMACH(14) /         56 /
! !DATA IMACH(15) /       -127 /
! !DATA IMACH(16) /        127 /
! !
  INTEGER FUNCTION I1MACH (I)
    !IMPLICIT INTEGER :: (I-N)
    !IMPLICIT REAL(dp) :: (A-H,O-Z)
    INTEGER :: IMACH(16),OUTPUT
    SAVE IMACH
    EQUIVALENCE (IMACH(4),OUTPUT)
    !***FIRST EXECUTABLE STATEMENT  I1MACH
    IF (I .GE. 1  .OR.  I .LE. 16) THEN
      I1MACH = IMACH(I)
      RETURN
    END IF
    WRITE (*, FMT = 9000)
    9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
    !
    !CALL FDUMP
    !
    STOP
  END FUNCTION I1MACH


  ! !***BEGIN PROLOGUE  XGETUA
  ! !***PURPOSE  Return unit number(s) to which error messages are being
  ! !  sent.
  ! !***LIBRARY   SLATEC (XERROR)
  ! !***CATEGORY  R3C
  ! !***TYPE      ALL (XGETUA-A)
  ! !***KEYWORDS  ERROR, XERROR
  ! !***AUTHOR  Jones, R. E., (SNLA)
  ! !***DESCRIPTION
  ! !
  ! !Abstract
  ! ! XGETUA may be called to determine the unit number or numbers
  ! ! to which error messages are being sent.
  ! ! These unit numbers may have been set by a call to XSETUN,
  ! ! or a call to XSETUA, or may be a default value.
  ! !
  ! !Description of Parameters
  ! ! --Output--
  ! ! IUNIT - an array of one to five unit numbers, depending
  ! !      on the value of N.  A value of zero refers to the
  ! !      default unit, as defined by the I1MACH machine
  ! !      constant routine.  Only IUNIT(1),...,IUNIT(N) are
  ! !      defined by XGETUA.  The values of IUNIT(N+1),...,
  ! !      IUNIT(5) are not defined (for N .LT. 5) or altered
  ! !      in any way by XGETUA.
  ! ! N     - the number of units to which copies of the
  ! !      error messages are being sent.  N will be in the
  ! !      range from 1 to 5.
  ! !
  ! !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  ! !Error-handling Package, SAND82-0800, Sandia
  ! !Laboratories, 1982.
  ! !***ROUTINES CALLED  J4SAVE
  ! !***REVISION HISTORY  (YYMMDD)
  ! ! 790801  DATE WRITTEN
  ! ! 861211  REVISION DATE from Version 3.2
  ! ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! ! 920501  Reformatted the REFERENCES section.  (WRB)
  ! !***END PROLOGUE  XGETUA
  SUBROUTINE XGETUA (IUNITA, N)
    DIMENSION IUNITA(5)
    N = J4SAVE(5,0,.FALSE.)
    DO I=1,N
      INDEX = I+4
      IF (I.EQ.1) INDEX = 3
      IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
    END DO
    RETURN
  END SUBROUTINE XGETUA


  ! !***BEGIN PROLOGUE  XERSVE
  ! !***SUBSIDIARY
  ! !***PURPOSE  Record that an error has occurred.
  ! !***LIBRARY   SLATEC (XERROR)
  ! !***CATEGORY  R3
  ! !***TYPE      ALL (XERSVE-A)
  ! !***KEYWORDS  ERROR, XERROR
  ! !***AUTHOR  Jones, R. E., (SNLA)
  ! !***DESCRIPTION
  ! !
  ! ! *Usage:
  ! !
  ! ! INTEGER ::  KFLAG, NERR, LEVEL, ICOUNT
  ! ! CHARACTER * (len) LIBRAR, SUBROU, MESSG
  ! !
  ! ! CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
  ! !
  ! ! *Arguments:
  ! !
  ! ! LIBRAR :IN    is the library that the message is from.
  ! ! SUBROU :IN    is the subroutine that the message is from.
  ! ! MESSG  :IN    is the message to be saved.
  ! ! KFLAG  :IN    indicates the action to be performed.
  ! !     when KFLAG > 0, the message in MESSG is saved.
  ! !     when KFLAG=0 the tables will be dumped and
  ! !     cleared.
  ! !     when KFLAG < 0, the tables will be dumped and
  ! !     not cleared.
  ! ! NERR   :IN    is the error number.
  ! ! LEVEL  :IN    is the error severity.
  ! ! ICOUNT :OUT   the number of times this message has been seen,
  ! !     or zero if the table has overflowed and does not
  ! !     contain this message specifically.  When KFLAG=0,
  ! !     ICOUNT will not be altered.
  ! !
  ! ! *Description:
  ! !
  ! ! Record that this error occurred and possibly dump and clear the
  ! ! tables.
  ! !
  ! !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  ! !Error-handling Package, SAND82-0800, Sandia
  ! !Laboratories, 1982.
  ! !***ROUTINES CALLED  I1MACH, XGETUA
  ! !***REVISION HISTORY  (YYMMDD)
  ! ! 800319  DATE WRITTEN
  ! ! 861211  REVISION DATE from Version 3.2
  ! ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  ! ! 900413  Routine modified to remove reference to KFLAG.  (WRB)
  ! ! 900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
  ! ! sequence, use IF-THEN-ELSE, make number of saved entries
  ! ! easily changeable, changed routine name from XERSAV to
  ! ! XERSVE.  (RWC)
  ! ! 910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
  ! ! 920501  Reformatted the REFERENCES section.  (WRB)
  ! !***END PROLOGUE  XERSVE
  SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,ICOUNT)
    INTEGER, PARAMETER :: LENTAB=10
    INTEGER :: LUN(5)
    CHARACTER(*) :: LIBRAR, SUBROU, MESSG
    CHARACTER(8) :: LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
    CHARACTER(20) :: MESTAB(LENTAB), MES
    DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
    SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
    DATA KOUNTX/0/, NMSG/0/
    !***FIRST EXECUTABLE STATEMENT  XERSVE
    !
    IF (KFLAG.LE.0) THEN
      ! Dump the table.
      IF (NMSG.EQ.0) RETURN
   
      ! Print to each unit.
      CALL XGETUA (LUN, NUNIT)
      DO KUNIT = 1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
        ! Print the table header.
        WRITE (IUNIT,9000)
        
        ! Print body of table.
        DO I = 1,NMSG
          WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),&
                &NERTAB(I),LEVTAB(I),KOUNT(I)
        END DO

        ! Print number of other errors.
        IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
        WRITE (IUNIT,9030)
      END DO

      ! Clear the error tables.
      IF (KFLAG.EQ.0) THEN
         NMSG = 0
         KOUNTX = 0
      ENDIF
    ELSE
      ! PROCESS A MESSAGE...
      ! SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
      ! OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      LIB = LIBRAR
      SUB = SUBROU
      MES = MESSG
      DO I = 1,NMSG
        IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.&
         &   MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.&
         &   LEVEL.EQ.LEVTAB(I)) THEN
           KOUNT(I) = KOUNT(I) + 1
           ICOUNT = KOUNT(I)
           RETURN
        ENDIF
      END DO
    
      IF (NMSG.LT.LENTAB) THEN
      ! Empty slot found for new message.
        NMSG = NMSG + 1
        LIBTAB(I) = LIB
        SUBTAB(I) = SUB
        MESTAB(I) = MES
        NERTAB(I) = NERR
        LEVTAB(I) = LEVEL
        KOUNT (I) = 1
        ICOUNT    = 1
      ELSE

        ! Table is full.
        KOUNTX = KOUNTX+1
        ICOUNT = 0
      ENDIF
    ENDIF
    RETURN
    !
    !Formats.
    !
     9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /&
          &  ' LIBRARY    SUBROUTINE MESSAGE START             NERR',&
         &   '     LEVEL     COUNT')
     9010 FORMAT (1X,A,3X,A,3X,A,3I10)
     9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
     9030 FORMAT (1X)
  END SUBROUTINE XERSVE


  SUBROUTINE FDUMP
  !***BEGIN PROLOGUE  FDUMP
  !***PURPOSE  Symbolic dump (should be locally written).
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3
  !***TYPE      ALL (FDUMP-A)
  !***KEYWORDS  ERROR, XERMSG
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  ! ***Note*** Machine Dependent Routine
  ! FDUMP is intended to be replaced by a locally written
  ! version which produces a symbolic dump.  Failing this,
  ! it should be replaced by a version which prints the
  ! subprogram nesting list.  Note that this dump must be
  ! printed on each of up to five files, as indicated by the
  ! XGETUA routine.  See XSETUA and XGETUA for details.
  !
  !Written by Ron Jones, with SLATEC Common Math Library Subcommittee
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  ! 790801  DATE WRITTEN
  ! 861211  REVISION DATE from Version 3.2
  ! 891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  FDUMP
  !***FIRST EXECUTABLE STATEMENT  FDUMP
  RETURN
  END SUBROUTINE FDUMP
END MODULE SLATEC