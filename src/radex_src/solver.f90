MODULE Solver
USE CommonData
USE Slatec
USE types
IMPLICIT NONE
Contains
  SUBROUTINE Matrix(niter,conv,success_flag)
    !     Set up rate matrix
    INTEGER :: niter            ! iteration counter
    INTEGER :: ilev,jlev,klev   ! to loop over energy levels
    INTEGER :: nplus            ! to solve statistical equilibrium
    INTEGER :: iline            ! to loop over lines
    INTEGER :: m,n              ! line upper/lower levels
    INTEGER :: nthick           ! counts optically thick lines
    INTEGER :: nfat             ! counts highly optically thick lines
    INTEGER :: nreduce          ! size of reduced rate matrix
    INTEGER :: indx(maxlev),dsign       ! needed for NumRep equation solver
    INTEGER :: terminate,success_flag

    REAL(dp) :: rhs(maxlev)           ! RHS of rate equation
    REAL(dp) :: yrate(maxlev,maxlev)  ! rate matrix

    REAL(dp) :: etr,exr               ! to calculate radiative rates
    REAL(dp) :: xt                    ! frequency cubed
    REAL(dp) :: hnu                   ! photon energy
    REAL(dp) :: bnutex                ! line source function
    REAL(dp) :: cddv                  ! N(mol) / delta V
    REAL(dp) :: beta          ! escape probability
    REAL(dp) :: bnu                   ! Planck function
    REAL(dp) :: uarray(maxlev,maxlev) ! reduced rate matrix
    REAL(dp) :: redcrit               ! reduction criterion
    REAL(dp) :: sumx                  ! summed radiative rate
    REAL(dp) :: total                 ! to normalize populations
    REAL(dp) :: tsum,thistex          ! to check convergence
    REAL(dp) :: lr                    ! learning rate to solve

    LOGICAL :: conv                 ! are we converged?

    !	keep old population for underrelaxation procedure -- sb/fvdt 30nov2011
    REAL(dp) :: xpopold(maxlev)		 

    !     db
    LOGICAL :: reduce
    reduce=.false.

    !     Executable statements begin here

    IF (debug) write(*,*) 'niter = ',niter

    !     Clear array of level populations.
    DO ilev=1,nlev
      rhs(ilev) = 0.0
      DO jlev=1,nlev
        yrate(ilev,jlev) = 0.0
      END DO
    END DO

    !     Initialize rate matrix
    nplus = nlev +1
    DO ilev=1,nlev
    !$$$	ctot(ilev)=0.0
    DO jlev=1,nlev
      yrate(ilev,jlev) = -1.0d-30*totdens
    END DO
    !     Add conservation equation
    yrate(nplus,ilev) = 1.0
    rhs(ilev) = 1.0e-30*totdens
    yrate(ilev,nplus) = 1.0d-30*totdens
    END DO
    rhs(nplus) = 1.0e-30*totdens

    !$$$c     Total collision rate
    !$$$      DO ilev=1,nlev
    !$$$	DO jlev=1,nlev
    !$$$	  ctot(ilev)=colld(ilev,jlev)+ctot(ilev)
    !$$$	END DO
    !$$$      END DO

    !     Contribution of radiative processes to the rate matrix.

    !     First iteration: use background intensity
    IF (niter.eq.0) THEN
       DO 51 iline=1,nline
          m   = iupp(iline)
          n   = ilow(iline)
          etr = fk*xnu(iline)/trj(iline)
          if(etr.ge.160.0d0) THEN         
             exr = 0.0d0
          else
             exr = 1.0/(dexp(etr)-1.0d0)
          END IF
    yrate(m,m) = yrate(m,m) + aeinst(iline)*(1.0+exr)
    yrate(n,n) = yrate(n,n) + aeinst(iline)*gstat(m)*exr/gstat(n)
    yrate(m,n) = yrate(m,n) - aeinst(iline)*(gstat(m)/gstat(n))*exr
    yrate(n,m) = yrate(n,m) - aeinst(iline)*(1.0+exr)
    51      END DO
    else
    !     Subsequent iterations: use escape probability.
       cddv = cdmol / deltav
    !     Count optically thick lines
       nthick = 0
       nfat   = 0

       DO 52 iline=1,nline
          xt  = xnu(iline)**3.0
          m   = iupp(iline)
          n   = ilow(iline)

    !     Calculate source function
          hnu = fk * xnu(iline) / tex(iline)
          if(hnu.ge.160.0) THEN
             bnutex = 0.0d0
          else
          bnutex = thc*xt/(dexp(fk*xnu(iline)/tex(iline))-1.0)
          END IF

    !     Calculate line optical depth.
          taul(iline) = cddv*(xpop(n)*gstat(m)/gstat(n)-xpop(m))&
               &/(fgaus*xt/aeinst(iline))
          if(taul(iline).gt.1.d-2) nthick = nthick+1
          if(taul(iline).gt.1.d05) nfat   = nfat+1

    !     Use escape probability approximation for internal intensity.
          beta = EscProb(taul(iline))
    !     Split off local contribution to radiation field  sb/fvdt 30nov2011
          bnu  = totalb(iline)*beta
    !c up to 30nov2011:     bnu  = totalb(iline)*beta+(1.0d0-beta)*bnutex
          exr  = bnu/(thc*xt)

    !     Radiative contribution to the rate matrix
    	yrate(m,m) = yrate(m,m)+aeinst(iline)*(beta+exr)
    !c up to 30nov2011:		yrate(m,m) = yrate(m,m)+aeinst(iline)*(1.0+exr)
        	yrate(n,n) = yrate(n,n)+aeinst(iline)*(gstat(m)*exr/gstat(n))
        	yrate(m,n) = yrate(m,n)-aeinst(iline)*(gstat(m)/gstat(n))*exr
        	yrate(n,m) = yrate(n,m)-aeinst(iline)*(beta+exr)
    !c up to 30nov2011:    	yrate(n,m) = yrate(n,m)-aeinst(iline)*(1.0+exr)
    52      END DO
    END IF

    !     Warn user if convergence problems expected
    IF ((niter.eq.1).and.(nfat.gt.0)) WRITE(*,*)&
         &"*** Warning: Some lines have very high optical depth"

    IF (debug) THEN
       WRITE(*,*) yrate(1,1),yrate(1,2),yrate(1,3),yrate(1,4)
       WRITE(*,*) yrate(2,1),yrate(2,2),yrate(2,3),yrate(2,4)
       WRITE(*,*) yrate(3,1),yrate(3,2),yrate(3,3),yrate(3,4)
       WRITE(*,*) yrate(4,1),yrate(4,2),yrate(4,3),yrate(4,4)
    END IF

    !     Contribution of collisional processes to the rate matrix.
    DO ilev=1,nlev
    yrate(ilev,ilev) = yrate(ilev,ilev) + ctot(ilev)
    DO jlev=1,nlev
         if(ilev.ne.jlev) yrate(ilev,jlev) = yrate(ilev,jlev) &
              &- crate(jlev,ilev)
    END DO
    END DO

    IF (debug) THEN
       WRITE(*,*) yrate(1,1),yrate(1,2),yrate(1,3),yrate(1,4)
       WRITE(*,*) yrate(2,1),yrate(2,2),yrate(2,3),yrate(2,4)
       WRITE(*,*) yrate(3,1),yrate(3,2),yrate(3,3),yrate(3,4)
       WRITE(*,*) yrate(4,1),yrate(4,2),yrate(4,3),yrate(4,4)
    END IF

    !     db
    IF (reduce) THEN

      !     An auxiliary array is passed to the linear equation solver after
      !     renormalization. The array Y retains the original matrix elements.

      IF (debug) WRITE(*,*) 'reducing matrix...' 

      DO ilev=1,nlev
         DO jlev=1,nlev
            uarray(ilev,jlev) = yrate(ilev,jlev)
         END DO
      END DO

      !     Now test whether the matrix should be reduced
      !     to exclude the radiatively coupled levels.

      redcrit = 10.0*tkin/fk
      nreduce = 0
      DO ilev=1,nlev
         if(eterm(ilev).le.redcrit) nreduce = nreduce+1
      END DO

      IF (debug) WRITE(*,*) 'nreduce=',nreduce

      !     We now separate the collisionally coupled levels from those that
      !     are coupled mainly by radiative processes, compute an effective
      !     cascade matrix for rates of transfer from one low-lying level
      !     to another and then solve this reduced system of equations
      !     explicitly for the low-lying levels only.

      DO jlev=1,nreduce
         DO ilev=1,nreduce
            DO klev=nreduce+1,nlev
               uarray(ilev,jlev) = abs(yrate(klev,jlev)*yrate(ilev,klev)&
                   & /yrate(klev,klev)) + uarray(ilev,jlev) 
            END DO
         END DO
      END DO

      !     Invert this reduced matrix

      IF (debug) WRITE(*,*) 'inverting reduced matrix...'

      !call ludcmp(uarray,nreduce+1,maxlev,indx,dsign)
      call lubksb(uarray,nreduce+1,maxlev,indx,rhs,success_flag)
      IF (success_flag .eq.0) RETURN

      IF (debug) WRITE(*,*) 'computing cascade...'

      !     Compute the populations of the highly excited states
      if(nlev.gt.nreduce) THEN
         DO klev=nreduce+1,nlev
            sumx = 0.0
            DO jlev=1,nreduce
               sumx = rhs(jlev)*yrate(klev,jlev) + sumx
            END DO
            rhs(klev) = abs(sumx/yrate(klev,klev))
         END DO
      END IF

    ELSE  !if we don't want to reduce
       IF (debug) WRITE(*,*) 'inverting non-reduced matrix...'
       !call ludcmp(yrate,nplus,maxlev,indx,dsign)
       call lubksb(yrate,nplus,maxlev,indx,rhs,success_flag)
       IF (success_flag .eq.0) RETURN
    END IF

    !     Level populations are the normalized RHS components
    total = 0.0d0
    DO ilev=1,nlev
       total = rhs(ilev)+total   
    END DO
    IF (debug) WRITE(*,*) 'total rhs=',total
    IF (debug) WRITE(*,*) 'rhs=',(rhs(ilev),ilev=1,nlev)

    !     Limit population to minpop
    DO ilev=1,nlev

    !sb301111 Store old population 
       xpopold(ilev)=xpop(ilev)

       xpop(ilev) = dmax1(minpop,rhs(ilev)/total)

    !sb301111 first iteration: there is no old population
        if(niter.eq.0) xpopold(ilev)=xpop(ilev)

    END DO

    IF (debug) WRITE(*,*) 'computing T_ex...'

    !     Compute excitation temperatures of the lines
    tsum = 0.0
    DO iline=1,nline
    m  = iupp(iline)
    n  = ilow(iline)
      xt = xnu(iline)**3.d0
      IF (niter.eq.0) THEN
         IF ((xpop(n).le.minpop).or.(xpop(m).le.minpop)) THEN
            tex(iline) = totalb(iline)
         else
            tex(iline) = fk*xnu(iline)/&
                 &(dlog(xpop(n)*gstat(m)/(xpop(m)*gstat(n))))
         END IF
      else
         IF ((xpop(n).le.minpop).or.(xpop(m).le.minpop)) THEN
            thistex = tex(iline)
         else
            thistex = fk*xnu(iline)/&
                 &(dlog(xpop(n)*gstat(m)/(xpop(m)*gstat(n))))
         END IF

    !     Only thick lines count for convergence
    if(taul(iline).gt.0.01) tsum = tsum&
              &+ abs((thistex-tex(iline))/thistex)

    !     Update excitation temperature & optical depth
    tex(iline)  = 0.5*(thistex+tex(iline))
    taul(iline) = cddv*(xpop(n)*gstat(m)/gstat(n)-xpop(m))&
            & /(fgaus*xt/aeinst(iline))
      END IF
    END DO

    IF (debug) write(*,*) niter,(tex(iline),iline=1,nline),&
         &(xpop(ilev),ilev=1,nlev),(taul(iline),iline=1,nline)

    !     Introduce a minimum number of iterations
    if(niter.ge.miniter) THEN
       if(nthick.eq.0) conv = .true.
       if(tsum/nthick.lt.ccrit) conv = .true.
    IF (debug) WRITE(*,*) niter,nthick,tsum,tsum/nthick,conv
    END IF

    !sb301111 now DO the underrelaxation!
    IF(niter.le.20) THEN
       lr=0.1+0.01*niter
    else
       lr=5.0/niter+0.05
    END IF
    DO ilev=1,nlev
       xpop(ilev)=lr*xpop(ilev)+(1.0-lr)*xpopold(ilev)
    end do

    return
  END SUBROUTINE Matrix

!     ------------------------------------------------------------

  FUNCTION EscProb(TAU)

    REAL(dp) :: EscProb,beta,tau
    REAL(dp) :: taur  !optical radius


    SELECT CASE (method)
      CASE(1)
        !Uniform sphere formula from Osterbrock (Astrophysics of
        !Gaseous Nebulae and Active Galactic Nuclei) Appendix 2
        !with power law approximations for large and small tau
        taur=tau/2.0
        IF (taur.lt.0) THEN
          beta=1.0-tau
        ELSE IF (taur.lt.0.1) THEN
          beta = 1.d0-0.75d0*taur+(taur**2.)/2.5d0&
            &-(taur**3.)/6.d0+(taur**4.)/17.5d0
        ELSE IF(taur.gt.5.d1) THEN
            beta = 0.75d0/taur
        ELSE
          beta = 0.75d0/taur*(1.d0-1.d0/(2.d0*(taur**2.))+&
              &(1.d0/taur+1.d0/(2.d0*(taur**2.)))*dexp(-2.*taur))
        END IF


      CASE(2)
        !Expanding sphere = Large Velocity Gradient (LVG) or Sobolev case.
        !Formula from De Jong, Boland and Dalgarno (1980, A&A 91, 68)
        !corrected by factor 2 in order to match EscProb(TAU=0)=1
        taur=2.34*tau
        IF (taur.lt.0) THEN
          beta=1.0-tau     
        ELSE IF (taur.lt.0.01) THEN
           beta = 1.0-taur/2.0+(taur**2.)/6.d0
        ELSE IF(taur.gt.5.d1) THEN
           beta = 1.0/taur
        ELSE
           beta = (1.0-dexp(-taur))/taur
        ENDIF

      CASE (3)
        !Slab geometry (e.g., shocks): de Jong, Dalgarno & Chu 1975, 
        !ApJ 199, 69 (again with power law approximations)
        taur=3.0*tau
        IF (taur.lt.0) THEN
           beta=1.0-tau     
        ELSE IF (taur.lt.0.01) THEN
           beta = 1.0-taur/2.0+(taur**2.)/6.d0
        ELSE IF(taur.gt.5.d1) THEN
           beta = 1.0/taur
        ELSE
           beta = (1.0-dexp(-taur))/taur
        ENDIF
      CASE DEFAULT 
         WRITE(*,*) 'Error: Escape probability method undefined'
         STOP
    END SELECT
    EscProb = beta
    RETURN
  END FUNCTION EscProb

  SUBROUTINE CalcOutputArrays(nlines)
    INTEGER, INTENT(INOUT) :: nlines !number of output lines
    INTEGER :: iline    ! to loop over lines
    INTEGER :: m,n      ! upper & lower level of the line

    REAL(dp) ::  xt        ! frequency cubed
    REAL(dp) ::  hnu       ! photon energy
    REAL(dp) ::  bnutex    ! line source function
    REAL(dp) ::  ftau      ! exp(-tau)
    REAL(dp) ::  toti      ! background intensity
    REAL(dp) ::  tbl       ! black body temperature
    REAL(dp) ::  wh        ! Planck correction
    REAL(dp) ::  tback     ! background temperature
    REAL(dp) ::  ta        ! line antenna temperature
    REAL(dp) ::  tr        ! line radiation temperature
    REAL(dp) ::  beta ! escape probability
    REAL(dp) ::  bnu       ! Planck function
    REAL(dp) ::  kkms      ! line integrated intensity (K km/s)
    REAL(dp) ::  ergs      ! line flux (erg / s / cm^2)
    REAL(dp) ::  wavel     ! line wavelength (micron)
    nlines=0
    DO iline=1,nline
      m  = iupp(iline)
      n  = ilow(iline)
      xt = xnu(iline)**3.
      !Calculate source function
      hnu = fk*xnu(iline)/tex(iline)
      IF (hnu.ge.160.0d0) THEN
        bnutex = 0.0d0
      else
        bnutex = thc*xt/(dexp(fk*xnu(iline)/tex(iline))-1.d0)
      END IF
      !Calculate line brightness in excess of background
      ftau = 0.0d0
      IF (abs(taul(iline)).le.3.d2) ftau = dexp(-taul(iline))
      toti = backi(iline)*ftau+bnutex*(1.d0-ftau)
      IF (toti.eq.0.0d0) THEN
        tbl = 0.0d0
      else
        wh = thc*xt/toti+1.d0
        IF (wh.le.0.d0) THEN
          tbl = toti/(thc*xnu(iline)*xnu(iline)/fk)
        else
          tbl = fk*xnu(iline)/dlog(wh)
        END IF
      END IF
      IF (backi(iline).eq.0.0d0) THEN
        tback = 0.0d0
      else
        tback = fk*xnu(iline)/dlog(thc*xt/backi(iline)+1.d0)
      END IF
      !Calculate antenna temperature
      tbl = tbl-tback
            hnu = fk*xnu(iline)
            IF (abs(tback/hnu).le.0.02) THEN
        ta = toti
      else
              ta = toti-backi(iline)
            END IF
            ta = ta/(thc*xnu(iline)*xnu(iline)/fk)
      !Calculate radiation temperature
      beta = escprob(taul(iline))
      bnu  = totalb(iline)*beta+(1.d0-beta)*bnutex
      IF (bnu.eq.0.0d0) THEN
        tr = totalb(iline)
      else
        wh = thc*xt/bnu+1.0
        IF (wh.le.0.0) THEN
          tr = bnu/(thc*xnu(iline)*xnu(iline)/fk)
        else
          tr = fk*xnu(iline)/dlog(wh)
        END IF
      END IF

      wavel = clight / spfreq(iline) / 1.0d5 ! unit =  micron
      kkms  = 1.0645*deltav*ta
      ergs  = fgaus*kboltz*deltav*ta*(xnu(iline)**3.)
      antennaTemp(iline)=ta
      wavelength(iline)=wavel
      lowerPops(iline)=(xpop(n))
      upperPops(iline)=(xpop(m))
      intensityKkms(iline)=(kkms/1.0d5)
      intensityErgs(iline)=ergs
      lowQNum(iline)=qnum(n)
      upperQNum(iline)=qnum(m)
      nlines=nlines+1
    END DO

  END SUBROUTINE CalcOutputArrays
END MODULE Solver      
