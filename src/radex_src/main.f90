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

PROGRAM Radex
USE IO
USE Solver
USE Background
IMPLICIT NONE
    !Main program: controls program flow and drives subroutines

    INTEGER :: niter,nlines   ! iteration counter
    INTEGER :: imore=1   ! are we running again?
    LOGICAL :: conv    ! are we converged?

    ! Begin executable statements
    write(*,*)
    write(*,*)'   Welcome to Radex, Moden Fortran Edition'
    write(*,*)

    method=1

    DO WHILE (imore.eq.1)
        !     Get input parameters
        IF (DEBUG) write(*,*) 'calling getinputs'

        CALL getinputs

        !     Read data file
        IF (DEBUG) write(*,*) 'calling readdata'
        CALL ReadData

        write(*,*) "Beginning Calculation"

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
                write(*,*) 'Finished in ',niter,' iterations.'
                EXIT
            END IF
        END DO

        IF (.NOT. conv) write(*,*) '   Warning: Calculation did not converge in ',maxiter&
            &,' iterations.'

        !Prepare output by calculating final quantities
        IF (DEBUG) write(*,*) 'calculating output summary variables'
        CALL CalcOutputArrays(nlines)

        !Write output
        IF (DEBUG) write(*,*) 'calling output'
        CALL output(niter)

        !     See if user wants more, else call it a day
        51   format(A,$)
        write(*,51) '  Another calculation [0/1] ? '
        read(*,*) imore
        write(13,52) imore
        52   format(i2)
    END DO
    write(*,*) '   Have a nice day.'
    !     Done! Now close log file ...
    close(13)
    !     ...and output file.
    close(8)
END PROGRAM Radex
