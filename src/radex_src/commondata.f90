MODULE CommonData
USE types
IMPLICIT NONE

    !file for input and output
    character(200) :: outfile,molfile,molfileold,specref
    character(*), PARAMETER :: radat='./data/'
    character(*), PARAMETER :: version = '30nov2011'
    character(*), PARAMETER :: logfile = './radex.log'

    !Escape probability method (uncomment your choice)
    INTEGER :: method=1 !1=uniform sphere,2=LVG, 3=slab


    !---------------------------------------------------------
    !Physical and astronomical constants (CODATA 2002)
    REAL(dp), PARAMETER :: clight  = 2.99792458d10  ! speed of light     (cm/s)
    REAL(dp), PARAMETER :: hplanck = 6.6260963d-27  ! Planck constant    (erg/Hz)
    REAL(dp), PARAMETER :: kboltz  = 1.3806505d-16  ! Boltzmann constant (erg/K)
    REAL(dp), PARAMETER :: pi      = 3.14159265d0  ! pi
    REAL(dp), PARAMETER :: amu     = 1.67262171d-24 ! atomic mass unit   (g)
    REAL(dp), PARAMETER :: tcmb = 2.725 ! CMB background temperature (K)
    !---------------------------------------------------------
      
    !Array sizes
    INTEGER, PARAMETER :: maxpart = 9     ! maximum no. of collision partners (seven defined)
    INTEGER, PARAMETER :: maxtemp = 99    ! maximum no. of collision temperatures
    INTEGER, PARAMETER :: maxlev  = 2999  ! maximum no. of energy levels
    INTEGER, PARAMETER :: maxline = 99999  ! maximum no. of radiative transitions
    INTEGER, PARAMETER :: maxcoll = 99999 ! maximum no. of collisional transitions

    
    !---------------------------------------------------------
    !   Molecular data
    !---------------------------------------------------------
    INTEGER :: nlev,nline,ncoll,npart,ntemp,iupp(maxline),ilow(maxline)
    !nlev:  actual number of levels
    !nline: actual number of lines
    !ncoll: actual number of transitions
    !npart: actual number of partners
    !ntemp: actual number of collision temperatures

    !iupp(i): upper level of line i
    !ilow(i): lower level of line i

    REAL(dp) :: amass,eterm(maxlev),gstat(maxlev),aeinst(maxline)&
         &,eup(maxline)
    ! amass:  molecular mass              (amu)
    ! eterm:  energy levels               (1/cm)
    ! gstat:  statistical weights
    ! aeinst: Einstein A coefficients     (1/s)
    ! eup:    line upper level energy     (K)
    ! colld:  downward rate coefficients  (cm^3 /s)
    ! xpop:   level populations

 
    !--------------------------------------------------------- 
    !Physical conditions
    !---------------------------------------------------------
    REAL(dp) :: density(maxpart),tkin,tbg,cdmol,deltav,totdens

    ! density:  number densities of collision partners  (cm^-3)
    ! totdens:  total number density of all partners    (cm^-3)
    ! tkin:     kinetic temperature                     (K)
    ! tbg:      temperature of background radiation     (K)
    ! cdmol:    molecular column density                (cm^-2)
    ! deltav:   FWHM line width                         (cm/s)


    !---------------------------------------------------------      
    !  Numerical parameters
    !---------------------------------------------------------
    INTEGER, PARAMETER :: miniter=10  ! minimum number of iterations
    INTEGER, PARAMETER :: maxiter=9999 ! maximum number of iterations

    REAL(dp) :: fmin,fmax !minimum/maximum output frequency
    REAL(dp), PARAMETER :: ccrit=1.0e-6   ! relative tolerance on solution
    REAL(dp), PARAMETER :: eps=1.0d-30    ! round-off error
    REAL(dp), PARAMETER :: minpop=1.0d-20 ! minimum level population

    !---------------------------------------------------------
    ! Radiative quantities
    !---------------------------------------------------------
    REAL(dp) :: taul(maxline),tex(maxline),backi(maxline),xnu(maxline)
    REAL(dp) :: trj(maxline),totalb(maxline),spfreq(maxline),antennaTemp(maxline)
    REAL(dp) :: upperPops(maxline),lowerPops(maxline),wavelength(maxline)
    REAL(dp) :: intensityKkms(maxline),intensityErgs(maxline)

    CHARACTER(6) :: qnum(maxlev),lowQNum(maxlev),upperQNum(maxlev)
    REAL(dp), PARAMETER :: fk    = hplanck*clight/kboltz
    REAL(dp), PARAMETER :: thc   = 2.d0*hplanck*clight
    REAL(dp), PARAMETER :: fgaus = 1.0645*8.0*pi

    ! xnu:    line frequency (cm^-1)
    ! taul:   line optical depth
    ! tex:    line excitation temperature

    ! trj:    background brightness (RJ)
    ! backi:  background intensity [erg s-1 cm-2 Hz-1 sr-1]
    ! totalb: background temperature (BB)

    ! fk,thc: help to calculate intensities
    ! fgaus:  accounts for Gaussian line shape

    ! spfreq: spectroscopic line frequency (GHz), not used in
    !       calculation but only to print output
    ! qnum:   quantum numbers of levels

    !---------------------------------------------------------
    !   Collisional quantities
    !---------------------------------------------------------
    REAL(dp) :: ctot(maxlev),crate(maxlev,maxlev),xpop(maxlev)

    ! crate: collision rate matrix (density * rate coefficient)
    ! ctot:  total collision rate 
    ! xpop:  level populations

    INTEGER :: id(maxpart)      ! to identify collision partners
    INTEGER :: lcu(maxpart,maxcoll),lcl(maxpart,maxcoll)
    REAL(dp) :: coll(maxpart,maxcoll,maxtemp)
    REAL(dp) :: temp(maxpart,maxtemp) ! collision temperatures
    character*120 collref ! text about source of collisional data

    !For development / maintenance purposes:
    LOGICAL, PARAMETER :: debug =.False.
END MODULE CommonData
