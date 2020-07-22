MODULE Background
USE CommonData
IMPLICIT NONE
CONTAINS

  SUBROUTINE BACKRAD
    !----- This routine returns the intensity of continuum radiation that is
    !----- felt by the radiating molecules.  Intensity is computed at the
    !----- line frequencies only.

    !----- OPTIONS:
    !-----  1 - Single blackbody; default is the cosmic microwave background
    !-----      at T_CMB=2.725 K.  For values different from the default, the
    !-----      corresponding redshift is determined according to 
    !-----      T=T_CMB(1+z)
    !-----  2 - The mean Galactic background radiation field plus CMB. This 
    !-----      is a slightly revised version of Black, J. H. 1994, in 
    !-----      The First Symposium on the Infrared Cirrus and Diffuse 
    !-----      Interstellar Clouds. ASP Conference Series, Vol. 58, 
    !-----      R.M. Cutri and W.B. Latter, Eds., p.355. Details are posted
    !-----      at http://www.oso.chalmers.se/~jblack/RESEARCH/isrf.html
    !-----      This spectrum is NOT adjustable by a scale factor because
    !-----      it consists of several components that are not expected to
    !-----      scale linearly with respect to each other.
    !-----  3 - A user-defined radiation field that is specified by NRAD
    !-----      values of frequency [cm-1], intensity [Jy nsr-1], and dilution
    !-----      factor [dimensionless]. Spline interpolation is applied to 
    !-----      this table. The intensity need not be specified at all 
    !-----      frequencies of the line list, but a warning message will
    !-----      appear if extrapolation (rather than interpolation) is required.


    CHARACTER(80) :: bgfile,title
    INTEGER :: iline,irad,nrad
    REAL(dp) :: xnubr(maxline),spinbr(maxline),dilbr(maxline),hnu,tbb3
    REAL(dp) :: logfreq(maxline),logflux(maxline),fpp(maxline)
    REAL(dp), PARAMETER :: huge=1.0e38 ! highest allowed by f90 (fdvt 28apr06)
    REAL(dp) :: aa,bb,fout
    REAL(dp) :: xnumin,xnumax,cbi,xln

    !bgfile:  file with user's background
    !title:   one-liner describing user's background
    !nrad:    number of wavelength points in user's background
    !irad:    to loop over background wavelength
    !iline:   to loop over lines
    !xnubr:   frequencies of user's background [cm^-1]
    !spinbr:  intensities of user's background [Jy/nsr]
    !dilbr:   dilution factors of user's background
    !hnu:     helps to calculate intensity
    !tbb3:    cmb addition to user's background
    !cbi:     cmb intensity

    !logfreq: base 10 log of frequency
    !logflux: base 10 log of intensity
    !fpp:     helps to calculate splines
    !aa,bb:   help to calculate Planck function
    !fout:    interpolated intensity
    !xnumin,xnumax: min/max line frequencies
    !xln:     log of line frequency


    IF (tbg.gt.0.0) THEN
    !option 1: Single black body
      DO iline = 1,nline
        hnu = fk*xnu(iline)/tbg
        if (debug) write(*,*) iline,hnu,xnu(iline)
        if(hnu.ge.160.0d0) THEN
          backi(iline) = eps
          !	sb/fvdt 30nov2011 Do not set backi to zero: may propagate into T(ex) 
          !	and if line is thick, Tsum = NaN and code never converges ...
        else
          backi(iline) = thc*(xnu(iline)**3.0)&
                   &/(dexp(fk*xnu(iline)/tbg)-1.d0)
        END IF
        trj(iline)    = tbg
        totalb(iline) = backi(iline)
      END DO


    ELSE IF (tbg.EQ.0.0) THEN
      !option 2:  mean Galactic background radiation
      CALL galbr

    ELSE IF (tbg.lt.0.0) THEN
      !option 3:  user-specified radiation field with spline interpolation
      21      format(A,$)
      write(*,21) 'File with observed background intensity? '
      read*, bgfile
      open(unit=31,file=bgfile,status='old',err=666)

      tbb3 = -1.0*tbg

      read(31,310) title
      read(31,*) nrad

      !Read and dilute the input intensity
      DO irad = 1,nrad
        read(31,*) xnubr(irad),spinbr(irad),dilbr(irad)
        spinbr(irad) = spinbr(irad)*dilbr(irad)
      END DO

      close(31)
      310   format(a)

      ! In most cases it is safest to interpolate within the logarithms
      ! of the frequencies and intensities
      DO irad = 1,nrad
      logfreq(irad) = dlog(xnubr(irad))
      logflux(irad) = dlog(spinbr(irad))
      END DO

      !Set the spline coefficients
      call splcoeff(logfreq,logflux,nrad,huge,huge,fpp)

      ! Interpolate continuum brightness at frequencies in the line list
      !---------------------------------------------------
      !------ NOTE:  Interpolation is done in the input --
      !------ spectrum after dilution factors have      --
      !------ been applied but before the CMB has been  --
      !------ added. Units converted from               --
      !------ [Jy nsr-1] to [erg s-1 cm-2 Hz-1 sr-1]    --
      !---------------------------------------------------

      xnumin = xnu(1)
      xnumax = 0.0d0

      DO iline = 1,nline
        xln = dlog(xnu(iline))
        call splintrp(logfreq,logflux,fpp,NRAD,xln,fout)
        aa  = thc*(xnu(iline)**3.0d0)
        bb  = fk*xnu(iline)
        !Add CMB if applicable
        if(tbb3.gt.0.0d0) THEN
           cbi = aa/(dexp(bb/tbb3)-1.0d0)
        else
           cbi = 0.0d0
        endif        
        if(xnu(iline).ge.xnubr(1).and.xnu(iline).le.xnubr(nrad)) THEN
           backi(iline) = 1.0d-14*dexp(fout) + cbi
        else
           backi(iline) = cbi
        endif
        trj(iline) = bb/dlog(aa/backi(iline) + 1.0d0)
        if(xnu(iline).lt.xnumin) xnumin = xnu(iline)
        if(xnu(iline).gt.xnumax) xnumax = xnu(iline)
        !added 24aug2011 (thanks Simon Bruderer)
        totalb(iline) = backi(iline)
      END DO

      IF ((xnumin.lt.xnubr(1)).or.(xnumax.gt.xnubr(nrad))) WRITE(*,*)& 
        & 'Warning: the line list requires extrapolation of the background'
    END IF

    RETURN
    666  WRITE(*,*)'Error opening background file'
    STOP
  END SUBROUTINE BACKRAD


  SUBROUTINE GALBR
    !
    !.....Computes the mean background radiation near the Sun's location
    !.....in the Galaxy:  the cosmic microwave background radiation is a
    !.....blackbody with T_CMB=2.725 (Fixsen & Mather 2002, ApJ, 581, 817)       
    !.....The far-IR/submm component is based upon an earlier analysis of 
    !.....COBE data and is described by a single-temperature fit (Wright 
    !.....et al. 1991, ApJ, 381, 200).  At frequencies below 10 cm-1 
    !.....(29.9 GHz), there is a background contribution from non-thermal 
    !.....radiation. The UV/Visible/near-IR part of the spectrum is based 
    !.....on the model of average Galactic starlight in the solar 
    !.....neighborhood of Mathis, Mezger, and Panagia (1983, Astron. 
    !.....Astrophys., 128, 212).
    !
    INTEGER :: iline
    REAL(dp) :: aa,hnuk,cbi,cmi,cmib,yy,xla,ylg


    !aa,hnuk: help to calculate Planck function
    !tcmb:    CMB temperature
    !cbi:     CMB intensity
    !cmi:     synchrotron radiation intensity
    !cmib:    dust radiation intensity
    !yy,xla,ylg: to calculate stellar radiation field

    DO iline = 1,nline
      aa   = thc*(xnu(iline)**3.0d0)
      hnuk = fk*xnu(iline)/tcmb

      IF (xnu(iline).le.10.0d0) THEN
        cbi = aa/(dexp(hnuk) - 1.0d0) 
        cmi = 0.3d0*1.767d-19/(xnu(iline)**0.75d0)     ! synchrotron component

      ELSE IF (xnu(iline).le.104.98d0) THEN
        cbi  = aa/(dexp(hnuk) - 1.0d0)
        cmib = aa/(dexp(fk*xnu(iline)/23.3d0) - 1.0d0)
        cmi  = 0.3d0*5.846d-7*(xnu(iline)**1.65d0)*cmib  ! COBE single-T dust

      ELSE IF (xnu(iline).le.1113.126d0) THEN
        cmi = 1.3853d-12*(xnu(iline)**(-1.8381d0))
        cbi = 0.0d0                                      

      ELSE IF (xnu(iline).le.4461.40d0) THEN
        cbi = 0.0d0
        cmi = 1.0d-18*(18.213601d0 - 0.023017717d0*xnu(iline)&
           &+ 1.1029705d-5*(xnu(iline)**2.0d0) &
           &- 2.1887383d-9*(xnu(iline)**3.0d0)&
           &+ 1.5728533d-13*(xnu(iline)**4.0d0))   

      ELSE IF (xnu(iline).le.8333.33d0) THEN
        cbi = 0.0d0
        cmi = 1.d-18*(-2.4304726d0 + 0.0020261152d0*xnu(iline)&
           &- 2.0830715d-7*(xnu(iline)**2.0d0)&
           &+ 6.1703393d-12*(xnu(iline)**3.0d0)) 

      ELSE IF (xnu(iline).le.14286.d0) THEN
        yy  = -17.092474d0 - 4.2153656d-5*xnu(iline)
        cbi = 0.0d0
        cmi = 10.d0**yy

      ELSE IF (xnu(iline).le.40000.d0) THEN
        xla = 1.0d8/xnu(iline)
        ylg = -1.7506877d-14*(xla**4.d0) &
              &+ 3.9030189d-10*(xla**3.d0)&
              &+ 3.1282174d-7*(xla*xla) &
              &- 3.0189024d-3*xla &
              &+ 2.0845155d0
        cmi = 1.581d-24*xnu(iline)*ylg
        cbi = 0.0d0

      ELSE IF (xnu(iline).le.55556.d0) THEN
        xla = 1.0d8/xnu(iline)
        ylg = -0.56020085d0 + 9.806303d-4*xla
        cmi = 1.581d-24*xnu(iline)*ylg
        cbi = 0.0d0

      ELSE IF (xnu(iline).le.90909.d0) THEN
        xla = 1.0d8/xnu(iline)
        ylg = -21.822255d0 + 3.2072800d-2*xla - 7.3408518d-6*xla*xla
        cmi = 1.581d-25*xnu(iline)*ylg
        cbi = 0.0d0

      ELSE IF (xnu(iline).le.109678.76d0) THEN
        xla = 1.0d8/xnu(iline)
        ylg = 30.955076d0 - 7.3393509d-2*xla + 4.4906433d-5*xla*xla
        cmi = 1.581d-25*xnu(iline)*ylg
        cbi = 0.0d0

      ELSE
        !  radiation field extends to Lyman limit of H
        write(*,202) xnu(iline)                                                     
         202 FORMAT(' ** XNU = ',1PE13.6,' IS OUTSIDE THE RANGE OF THE FITTING'&       
          & //'FITTINGUNCTION AND BEYOND LYMAN LIMIT')
        cbi=0.0d0
        cmi=0.0d0

      END IF

      backi(iline) = cbi+cmi                                                        
      trj(iline)   = fk*xnu(iline)/dlog(1.0d0+aa/backi(iline))    ! brightness temperature
      !	  added 24aug2011 (thanks Simon Bruderer)
      totalb(iline) = backi(iline)

    END DO

    RETURN                                                                  
  END SUBROUTINE GALBR                                                                     

!
!.....This package contains several routines for applications of 
!.....cubic-spline interpolation.  A further routine, splinteg, is
!.....available on request (John.Black@chalmers.se) and does
!.....numerical quadrature through use of the spline coefficients
!.....determined for a set of points in routine splcoeff.  These
!.....routines have been written by J. H. Black.  The basic spline
!.....interpolation routines are based upon the routines of Numerical
!.....Recipes (Ch. 3.3), but are not identical to them.
!
  SUBROUTINE SPLCOEFF(x,f,N,fp1,fpn,fpp)
    INTEGER :: N,NMAX
    INTEGER :: I,K
    PARAMETER (NMAX=2500)
    REAL(dp) :: fp1,fpn,x(nmax),f(nmax),fpp(nmax)
    REAL(dp) :: p,qn,sig,un,u(NMAX)
    !.....N values of a function f(x) are tabulated at points x(i), in 
    !.....order of increasing x(i).  fpp(x) is the evaluated second derivative
    !.....of f(x).  fp1 and fpn are the values of the first derivatives at
    !.....the beginning and ending points, i=1 and i=N.  For natural splines,
    !.....the second derivative is set to zero at a boundary.  Set fp1 and/or
    !.....fpn > 1.0D60 for natural splines.  Otherwise, the derivatives can
    !.....be specified and matched to those of extrapolating functions.
    !
    !.....Lower boundary condition (use 1e60 for f77/g77, 1e38 for f90):
    if (fp1.gt.0.99d38) THEN
      fpp(1)=0.d0
      u(1)=0.d0
    else
      fpp(1)=-0.5d0
      u(1)=(3.d0/(x(2)-x(1)))*((f(2)-f(1))/(x(2)-x(1))-fp1)
    endif
    !
    do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*fpp(i-1)+2.d0
      fpp(i)=(sig-1.d0)/p
      u(i)=(6.d0*((f(i+1)-f(i))/(x(i+1)-x(i))-(f(i)-f(i-1))/&
       &(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    !
    !.....Upper boundary condition (see above):
    if (fpn.gt.0.99d38) THEN
      qn=0.d0
      un=0.d0
    else
      qn=0.5d0
      un=(3.d0/(x(n)-x(n-1)))*(fpn-(f(n)-f(n-1))/(x(n)-x(n-1)))
    endif
    !
    fpp(n)=(un-qn*u(n-1))/(qn*fpp(n-1)+1.)
    do k=n-1,1,-1
      fpp(k)=fpp(k)*fpp(k+1)+u(k)
    enddo
    !
    RETURN
  END SUBROUTINE SPLCOEFF
 
  SUBROUTINE splintrp(xin,fin,fppin,N,x,fout)
    !.....For N tabulated values xin(i) and fin(x) of a function, and
    !.....the array fppin(x), which is the 2nd derivative of fin(x) delivered
    !.....by SPLCOEFF above, an interpolated value of fout is delivered
    !.....at x.  The routine can also, if desired, return the 
    !.....values of the first and second derivatives of the fitted 
    !.....function, foutp and foutpp.
    INTEGER :: N
    INTEGER :: k,khi,klo
    REAL(dp) :: x,fout,xin(n),fppin(n),fin(n)
    REAL(dp) :: a,b,h,foutp,foutpp
    klo=1
    khi=N
    1     if (khi-klo.gt.1) THEN
      k=(khi+klo)/2
      if(xin(k).gt.x)THEN
        khi=k
      else
        klo=k
      endif
    goto 1
    endif
    h=xin(khi)-xin(klo)
    if (h.eq.0.d0) WRITE(*,*) 'Warning: bad xin input in splintrp '
    a=(xin(khi)-x)/h
    b=(x-xin(klo))/h
    fout=a*fin(klo)+b*fin(khi)+((a**3.d0-a)*fppin(klo)+&
       &(b**3.d0-b)*fppin(khi))*(h**2.d0)/6.d0
    !
    !.....first derivative
    !
    foutp=(fin(khi)-fin(klo))/h - (3.d0*a*a-1.d0)*h*fppin(klo)/6.d0 +&
         &(3.d0*b*b-1.d0)*h*fppin(khi)/6.d0
    !
    !.....second derivative
    !
    foutpp=a*fppin(klo)+b*fppin(khi)
    !
    RETURN
  END SUBROUTINE splintrp
END MODULE Background
      





