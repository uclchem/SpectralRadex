!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	RADEX
!	Original code by Van der Tak et al. 2007
!	All publications using this code should reference their release paper:
!	A&A 468, 627 (2007)
!
!   This amended version by Jon Holdship has been updated to Modern Fortran 
!   with a view to removing common blocks and making an F2PY module for python.
!   
!	This code has been tested against the original but the authors do not
!	guarantee there are no errors. We advise users to check results they intend
!	to publish with another version of RADEX or another code.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

SUBROUTINE RADEX(inputDictionary,success_flag,nlines,Qup,Qlow,lineOutputs)
USE IO
USE Solver
USE Background
IMPLICIT NONE
    !Main program: controls program flow and drives subroutines
    CHARACTER(*) :: inputDictionary
    INTEGER :: niter,nlines,success_flag,iline   ! iteration counters
    DOUBLE PRECISION :: lineOutputs(10,500)
    CHARACTER(6) :: Qup(3000),Qlow(3000)
    LOGICAL :: conv    ! are we converged?
    !f2py intent(in) inputDictionary
    !f2py intent(out) success_flag,nlines, Qup,Qlow,lineOutputs
    success_flag=1
    !     Get input parameters
    IF (DEBUG) write(*,*) 'calling getinputs'

    CALL parseInputDictionary(inputDictionary,success_flag)
    !     Read data file
    IF (success_flag .ne. 1) RETURN
    IF (DEBUG) write(*,*) 'calling readdata'
    CALL ReadData

    !     Calculate background radiation field
    IF (DEBUG) write(*,*) 'calling backrad'
    CALL backrad

    niter = 0
    conv  = .false.

    !Set up rate matrix, splitting it in radiative and collisional parts
    !Invert rate matrix to get `thin' starting condition
    IF (DEBUG) write(*,*) 'calling matrix'
    CALL matrix(niter,conv)

    !Start iterating
    DO niter=1,maxiter
        !Invert rate matrix using escape probability for line photons
        CALL matrix(niter,conv)
        IF (conv) THEN
            EXIT
        END IF
    END DO

    IF (.NOT. conv) write(*,*) '   Warning: Calculation did not converge in ',maxiter&
        &,' iterations.'

    !Prepare output by calculating final quantities
    IF (DEBUG) write(*,*) 'calculating output summary variables'
    CALL CalcOutputArrays(nlines)
    niter=1
    DO iline=1,nlines
      !Check if line within output freq range
      IF (spfreq(iline).lt.fmax.and.spfreq(iline).gt.fmin) THEN
        lineOutputs(:,niter)=(/eup(iline),spfreq(iline),wavelength(iline),&
            &tex(iline),taul(iline),antennaTemp(iline),upperPops(iline),lowerPops(iline),&
            &intensityKkms(iline),intensityErgs(iline)/)
        Qup(niter)=upperQNum(iline)
        Qlow(niter)=lowQNum(iline)
        niter=niter+1
      END IF
    END DO
END SUBROUTINE RADEX
