! Integrating d2s/dx/dy over x, y, z and Enu (4d)
! For tau neutrino regeneration
! A. Vincent 8/11/2016
! Input: 1) output file name
! Output is

	module pars
	implicit none
	double precision, parameter :: pi=3.141592653, mN=0.9385,s2t=0.23,g2w = 0.42,NAVO=6.022d23
    double precision, parameter :: Lfu=0.3467, Lfd=-0.4233, Rfu=-0.1533,Rfd=0.0767
    double precision, parameter :: MZ = 91.18d0, MW = 80.385d0,mtau = 1.776
    double precision, parameter :: hbarc=1.97d-14
	end module pars
	

    
    MODULE Integrand
    IMPLICIT NONE
    !PUBLIC :: F
	!all units are GeV

    double precision :: Enu,EnuOut,flavor !these are set in the main program
    !PRIVATE
    CONTAINS

    FUNCTION HEAVISIDE(xin)
    double precision, INTENT(IN) :: xin
    double precision :: heaviside
    if (xin .gt. 0.d0) then
    HEAVISIDE =  1.d0
    else if (xin .lt. 0.d0) then
    HEAVISIDE = 0.d0
    else
    HEAVISIDE =  0.50d0
    END IF
    END FUNCTION HEAVISIDE

    !B13 of Palomares
    FUNCTION Dtau(Etau)
    double precision,intent(IN) :: Etau
    double precision :: Dtau, p1, q1, q2
    p1 = 0.833
    q1 = 1.66
    q2 = 1.15
    Dtau = (1.d0+p1*Etau/10.d6)/(1.+q1*Etau/10.d6 + q2*(Etau/10.d6)**2)
    END FUNCTION Dtau


    !secondary production of nue, numu
    function dnnottaudE(z,pol)
    double precision, intent(in) :: z
    integer pol
    double precision dnnottaudE
    dnnottaudE = 0.18*(4.0-12.*z+12.*z**2-4.*z**3)

    end function dnnottaudE
    !tau decay energy distribution from https://arxiv.org/pdf/hep-ph/0005310v2.pdf appendix A
    function dntaudE(z,channel,pol)
    integer channel,pol !tau- (from nutau) --> pol = -1; nutaubar --> pol = +1
    double precision,intent(IN) :: z
    double precision :: dntaudE,g0, g1,r
    if (channel == 0) then  !electron channel
     g0 = 5./3. - 3.d0*z**2 + 4.d0/3.d0*z**3
     g1 = 1./3. - 3.d0*z**2 + 8.d0/3.d0*z**3
     dntaudE =  .18*(g0+pol*g1)
    else if (channel == 1) then !hadron channel
     dntaudE = 0.
!     1)pions
     r = .1395**2/1.776**2
     if (1. .gt. r + z) then
        g0 = 1./(1.-r) !*heaviside(1.d0-r-z)
        g1 = -(2.*z-1.+r)/(1.-r)**2!*heaviside(1.d0-r-z)
        dntaudE = dntaudE +  .12*(g0+pol*g1)
!    2) rho
        r = .77**2/1.776**2
        if (1. .gt. r+z) then
            g0 = 1./(1.-r) !*heaviside(1.d0-r-z)
            g1 = -(2.*z-1.+r)/(1.-r)*(1.-2*r)/(1+2*r) !*heaviside(1.d0-r-z)
            dntaudE = dntaudE  + .26*(g0+pol*g1)
!    3) a1 (yeah this is a thing!)
            r = 1.26**2/1.776**2
            if (1. .gt. r+z) then
                g0 = 1./(1.-r) !*heaviside(1.d0-r-z)
                g1 = -(2.*z-1.+r)/(1.-r)*(1.-2*r)/(1+2*r) !*heaviside(1.d0-r-z)
                dntaudE = dntaudE + .13*(g0+pol*g1)
            end if
        end if
      end if
!   4) X (everything else hadronic)
    if (z .gt. 0.3) then
    g0 = 1./.30!*heaviside(0.3-z)
     dntaudE = dntaudE + .13*g0
    end if
    end if

    end function dntaudE

    !bit to account for tau mass from https://arxiv.org/pdf/hep-ph/0407371v2.pdf
    FUNCTION DSDXDY_nubar(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nubar,deltatau
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))
    PDFarray = getPDF(x,Q)
    deltatau = mtau**2/(s(Enu)-mN**2)

    if (x .lt. deltatau/(1.-mtau/Enu)) then
        dsdxdy_nubar = 0.0
    else
	dsdxdy_nubar = mN*Enu/(16.*pi)*g2w**2&
    /(Q**2+MW**2)**2*((PDFarray(2)+PDFarray(4))*(1.d0-deltatau/x-y)*(1.-y) &
    +(1.-deltatau/x)*(PDFarray(-1) + PDFarray(-3)+PDFarray(-5)))
 
    !for testing
!     Q = 100.0d0 ! Q! 
!     x=1.0d-04 ! x    
!     PDFarray = getPDF(x,Q)
!     write(*,*) x, PDFarray(-2) 
    end if
    end function dsdxdy_nubar
    
    FUNCTION DSDXDY_nu(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nu, deltatau
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))
    PDFarray = getPDF(x,Q)
    deltatau = mtau**2/(s(Enu)-mN**2)
    if (x .lt. deltatau/(1.-mtau/Enu)) then
    dsdxdy_nu = 0.0
    else
    !overall x factor is inside PDFs
    dsdxdy_nu =  mN*Enu/(16.*pi)*g2w**2&
    /(Q**2+MW**2)**2*((1.-deltatau/x)*(PDFarray(1)+PDFarray(3)+PDFarray(5)) &
    +(1.d0-deltatau/x-y)*(1-y)*(PDFarray(-2) + PDFarray(-4)))
    end if


    end function dsdxdy_nu
    
    !! END CROSS SECTION, BEGIN KINEMATICS
    
    !mandelstam s
	FUNCTION s(Enu)
	use pars
	double precision :: Enu, s
    s= 2.d0*mN*Enu
    END FUNCTION s
    
    
    !Q^2
	FUNCTION Q2(x,y,Enu)
	use pars
	double precision :: x, y, Enu, Q2
    Q2= 2.d0*x*y*mN*Enu
    END FUNCTION Q2
    

    FUNCTION getPDF(x,Q) ! Gets all flavors: 0: glue, 1-> d 2-> u, 3-> s, 4-> c, 5-> b and 6 -> t
    double precision x, Q, getPDF(-6:6)
    double precision xpdfCV(-6:6)
    call InitPDF(0)     
    call evolvePDF(x,Q,xpdfCV)
    getPDF = xpdfCV
    END FUNCTION getPDF


    
    !!!!!!    !!!!!!    !!!!!!    !!!!!!    !!!!!!    !!!!!!    !!!!!!
    !Function to be integrated below
    
    FUNCTION F(NUMFUN,TS) RESULT(Value)
    USE Precision_Model
    use pars

    
    double precision u,x,y,z,xs
    double precision dndz
    double precision bork1, bork2
    
    INTEGER, INTENT(IN) :: NUMFUN
    REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: TS
    REAL(kind=stnd), DIMENSION(1:NUMFUN) :: Value
    x = TS(1)
    y = TS(2)

    x = exp(x)
    y = exp(y)
    z = EnuOut/Enu/(1-y)
    if (z .gt. 1) then !kinematically forbidden
    u = 0.0
    else

    if (flavor**2 .ne. 9) then !electron or muon regeneration
        dndz = dnnottaudE(z,1)
        if (flavor .gt. 0) xs = DSDXDY_nu(x,y,Enu)
        if (flavor .lt. 0) xs = DSDXDY_nubar(x,y,Enu)
    else if (flavor .gt. 0) then !produce tau-
        dndz = dntaudE(z,0,-1) + dntaudE(z,1,-1)
        xs = DSDXDY_nu(x,y,Enu)
    else if (flavor .lt. 0) then !tau+
        dndz = dntaudE(z,0,+1) + dntaudE(z,1,+1)
        xs = DSDXDY_nubar(x,y,Enu)

    end if

		u = x*y*z/EnuOut*xs*dndz*hbarc**2

    end if

    Value(1) = u
    
    RETURN
    END FUNCTION F
    

    END MODULE Integrand
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                   Begin program here                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PROGRAM Example2D


    USE Precision_Model
    USE CUI                      ! Cubpack User Interface
    USE Integrand
    USE pars
    implicit none

    character*50 :: leminin,lemaxin, flavorIn
    INTEGER, PARAMETER :: n=2, & ! the dimension
                          m=1, & ! the number of simple regions
                          l=1, & ! the length of the integrand vector
    Cube=2

    INTEGER, DIMENSION(1:m) :: RgType
    INTEGER :: NEval, j,i,k
    integer, parameter :: nybins = 200
    REAL(kind=stnd), DIMENSION(1:n,0:n,1:m) :: Vertices
    REAL(kind=stnd), DIMENSION(1:l) ::  IntegralValue, AbsErr
    REAL(kind=stnd) ::  EpsRel, EpsAbs
    real :: x0,y0,x1,y1
    real :: starttime, finishtime
    LOGICAL :: Restart
    double precision dNdy(nybins,nybins),dNdy_p(nybins,nybins),dNdy_n(nybins,nybins)
    double precision Evec (nybins)

    !    integer  nevarray(res,res)
    double precision :: z1,z2 !, logxpy,logxmy,logymx
    double precision :: exposure,logEmin,logEmax,dE
    double precision :: bork1, bork2
    character(len=1024) filename
    character(len=8) frmt
    character (len=1024) muout

    REAL(kind=stnd), DIMENSION(2):: TS

    !    NNPDF tools
    integer NumberReplicas




    call get_command_argument(1,leminin)
    call get_command_argument(2,lemaxin)
    call get_command_argument(3,flavorIn)
    call get_command_argument(4,filename)

    read( leminin,*) logEmin
    read( lemaxin,*) logEmax
    read( flavorIn,*) flavor

    if (flavor .eq. -3) print*, "Producing nutaubar regeneration cross sections. This will take forever"
    if (flavor .eq. -2) print*, "Producing numubar (and nuebar) regeneration cross sections"
    if (flavor .eq. -1) print*, "Producing numubar (and nuebar) regeneration cross sections"
    if (flavor .eq. 3) print*, "Producing nutau regeneration cross sections. This will take forever"
    if (flavor .eq. 2) print*, "Producing numu (and nue) regeneration cross sections"
    if (flavor .eq. 1) print*, "Producing numu (and nue) regeneration cross sections"

    dE = (logEmax-logEmin)/nybins



    epsrel = 1.d-3
    epsabs = 1.d-100
    Restart= .false.

    RgType(1) = Cube
    x0 = log(1.d-7)
    y0 = log(1.d-7)
    x1 = 0.
    y1 = -1e-7 !avoids division by zero


    !		   x      log10(Enu)
    Vertices(1:n,0,1) = (/ x0, y0 /)
    Vertices(1:n,1,1) = (/ x1 ,y0 /)
    Vertices(1:n,2,1) = (/ x0 , y1  /)

            ! 1) Protons:
    call InitPDFsetByName("CT14qed_inc_proton")
    call numberPDF(NumberReplicas)



    call cpu_time(starttime)

    do j=1,nybins
        Enu = 10.d0**(logEmin+dE*dble(j-1))
        call cpu_time(bork1)

        do k = 1,nybins
            EnuOut = 10.d0**(logEmin+dE*dble(k-1))

            if  (Enuout .ge. Enu) then !kinematically not allowed
            z1 = 0
            else

            CALL CUBATR(n,l,F,m,Vertices(:,:,1:m),RgType,IntegralValue,AbsErr, &
            NEval=NEval,EpsRel=epsrel,EpsAbs=EpsAbs,Restart=Restart,MaxPts=50000000)
            call cubatr()

            z1 = IntegralValue(1)

            evec(j) = Enu
            end if
            dNdy_p(j,k) = z1

        end do

        call cpu_time(bork2)
        print*, j, bork2-bork1
    end do

    ! 2) Neutrons
    call InitPDFsetByName("CT14qed_inc_neutron")
    call numberPDF(NumberReplicas)
    do j=1,nybins
        call cpu_time(bork1)
        Enu = 10.d0**(logEmin+dE*dble(j-1))
        do k = 1,nybins
            EnuOut = 10.d0**(logEmin+dE*dble(k-1))
!            y = 1.-EnuOut/Enu
            if (Enuout .ge. Enu) then
            z2 = 0
            else

            CALL CUBATR(n,l,F,m,Vertices(:,:,1:m),RgType,IntegralValue,AbsErr, &
            NEval=NEval,EpsRel=epsrel,EpsAbs=EpsAbs,Restart=Restart,MaxPts=50000000)

            call cubatr()

            z2 = IntegralValue(1)
            end if
            dNdy_n(j,k) = z2
        end do
        call cpu_time(bork2)
        print*, j, bork2-bork1
    end do
    dNdy = (dNdy_n+dNdy_p)/2


    call cpu_time(finishtime)

    print*,"Time = ",finishtime-starttime," seconds."
    frmt = '(I3.3)'
    open(unit=55,file=filename,action='write')
    do i=1,nybins
    write(55,*) (dNdy(i,j), j=1,nybins) !, dNdy(i)
    end do
    close(55)
    print*, "Mission accomplished"
    STOP
    END PROGRAM Example2D
