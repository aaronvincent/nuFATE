! Integrating d2s/dx/dy over x and y
!input:
!logemin (GeV)
!log E max (GeV)
!flavor (1,2,3) = e, mu tau, (-1,-2,-3) = antineutrinos
!output file name

!output: vector of total cross sections (cm^2)
!total cross section
! A. Vincent 03/2017




	module pars
	implicit none
	double precision, parameter :: pi=3.141592653, mN=0.9385,s2t=0.23,g2w = 0.42,NAVO=6.022d23
    double precision, parameter :: Lfu=0.3467, Lfd=-0.4233, Rfu=-0.1533,Rfd=0.0767,GW = 2.085
    double precision, parameter :: MZ = 91.18d0, MW = 80.385d0,GF = 1.16d-5,mmu=.106d0,me=511.d-6
    double precision, parameter :: hbarc=1.97d-14, mtau = 1.776
	end module pars
	

    
    MODULE Integrand
    IMPLICIT NONE
    !PUBLIC :: F
	!all units are GeV

    double precision :: Enu, flavor !these are set in the main program
    !PRIVATE
    CONTAINS


!!!! CROSS SECTION DEFINITIONS
!!!! CHARGED CURRENT


!    CC nubar-nucleus (e and mu)
    FUNCTION DSDXDY_nubar_CC(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nubar_CC
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))
    PDFarray = getPDF(x,Q)
	dsdxdy_nubar_CC = mN*Enu/(16.*pi)*g2w**2&
    /(Q**2+MW**2)**2*((PDFarray(2)+PDFarray(4))*(1.d0-y)**2 &
    +(PDFarray(-1) + PDFarray(-3)+PDFarray(-5)))
    end function dsdxdy_nubar_CC

!   CC nu-nucleus (e and mu)
    FUNCTION DSDXDY_nu_CC(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nu_CC
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))
    PDFarray = getPDF(x,Q)
    !overall x factor is inside PDFs
    dsdxdy_nu_CC =  mN*Enu/(16.*pi)*g2w**2&
    /(Q**2+MW**2)**2*((PDFarray(1)+PDFarray(3)+PDFarray(5)) &
    +(1.d0-y)**2*(PDFarray(-2) + PDFarray(-4)))
    end function dsdxdy_nu_CC

!    CC nubar-electron (no regeneration)
    FUNCTION sigeh(Enu)
    use pars
    double precision :: Enu, sigeh
    sigeh = 1.d0/3.d0*GF**2*selectron(Enu)/pi*(1.d0-(mmu**2-me**2)/selectron(Enu))**2 &
    /((1.d0-selectron(Enu)/MW**2)**2+GW**2/MW**2)*.676/.1057
    end function sigeh

    !    CC nubar-nucleus (tau)
    !bit to account for tau mass from https://arxiv.org/pdf/hep-ph/0407371v2.pdf
    FUNCTION DSDXDY_nubar_CCtau(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nubar_CCtau,deltatau
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))
    PDFarray = getPDF(x,Q)
    deltatau = 1.76**2/(s(Enu)-mN**2)
    if (x .lt. deltatau/(1-mtau/Enu)) then
        dsdxdy_nubar_CCtau = 0.0
    else
        dsdxdy_nubar_CCtau = mN*Enu/(16.*pi)*g2w**2&
        /(Q**2+MW**2)**2*((PDFarray(2)+PDFarray(4))*(1.d0-deltatau/x-y)*(1.-y) &
        +(1.-deltatau/x)*(PDFarray(-1) + PDFarray(-3)+PDFarray(-5)))
    end if
    end function dsdxdy_nubar_CCtau

    !    CC nu-nucleus (tau)
    FUNCTION DSDXDY_nu_CCtau(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nu_CCtau, deltatau
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))
    PDFarray = getPDF(x,Q)
    deltatau = 1.76**2/(s(Enu)-mN**2)
    if (x .lt. deltatau/(1-mtau/Enu)) then
        dsdxdy_nu_CCtau = 0.0
    else
    !overall x factor is inside PDFs
        dsdxdy_nu_CCtau =  mN*Enu/(16.*pi)*g2w**2&
        /(Q**2+MW**2)**2*((1.-deltatau/x)*(PDFarray(1)+PDFarray(3)+PDFarray(5)) &
        +(1.d0-deltatau/x-y)*(1-y)*(PDFarray(-2) + PDFarray(-4)))
    end if
    end function dsdxdy_nu_CCtau


!!!! NEUTRAL CURRENT
    !nubar-nucleus
    FUNCTION DSDXDY_nubar_NC(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nubar_NC
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))

    PDFarray = getPDF(x,Q)

    dsdxdy_nubar_NC = mN*Enu/(16.*pi)*((afu(x,y,Enu)**2 + bfu(x,y,Enu)**2*(1.-y)**2) &
    *(PDFarray(2)/2.d0+pdfarray(4)) &
    +(afsb(x,y,Enu)**2 + bfsb(x,y,Enu)**2*(1.-y)**2)*(PDFarray(3) + PDFarray(5))  &
    +(bfu(x,y,Enu)**2 + afu(x,y,Enu)**2*(1.-y)**2)*(PDFarray(-2)/2.d0+PDFarray(-4))&
    +(bfsb(x,y,Enu)**2 + afsb(x,y,Enu)**2*(1-y)**2)*(PDFarray(-3) + PDFarray(-5))&
    +(bfd(x,y,Enu)**2 + afd(x,y,Enu)**2*(1-y)**2)*PDFarray(1)/2. &
    +(afd(x,y,Enu)**2 + bfd(x,y,Enu)**2*(1-y)**2)*PDFarray(-1)/2.d0)

    !for testing
    !     Q = 100.0d0 ! Q!
    !     x=1.0d-04 ! x
    !     PDFarray = getPDF(x,Q)
    !     write(*,*) x, PDFarray(-2)

    end function dsdxdy_nubar_NC

    !nu-nucleus, neutral current
    FUNCTION DSDXDY_nu_NC(x,y,Enu)
    use pars
    double precision :: x,y,Enu,Q,DSDXDY_nu_NC
    double precision pdfarray(-6:6)
    Q = sqrt(Q2(x,y,Enu))

    PDFarray = getPDF(x,Q)

    dsdxdy_nu_NC =  mN*Enu/(16.*pi)*((bfu(x,y,Enu)**2 + afu(x,y,Enu)**2*(1.-y)**2) &
    *(PDFarray(2)/2.d0+pdfarray(4)) &
    +(bfsb(x,y,Enu)**2 + afsb(x,y,Enu)**2*(1.-y)**2)*(PDFarray(3) + PDFarray(5))  &
    +(afu(x,y,Enu)**2 + bfu(x,y,Enu)**2*(1.-y)**2)*(PDFarray(-2)/2.d0+PDFarray(-4))&
    +(afsb(x,y,Enu)**2 + bfsb(x,y,Enu)**2*(1-y)**2)*(PDFarray(-3) + PDFarray(-5))&
    +(afd(x,y,Enu)**2 + bfd(x,y,Enu)**2*(1-y)**2)*PDFarray(1)/2. &
    +(bfd(x,y,Enu)**2 + afd(x,y,Enu)**2*(1-y)**2)*PDFarray(-1)/2.d0)
    !for testing
    !     Q = 100.0d0 ! Q!
    !     x=1.0d-04 ! x
    !     PDFarray = getPDF(x,Q)
    !     write(*,*) x, PDFarray(-2)

    end function dsdxdy_nu_NC


    !aux fct bits: these are like this for historical reasons
    FUNCTION afu(x,y,Enu)
    use pars
    double precision :: x, y, Enu, afu
    afu = g2w/(1-s2t)*Lfu/(Q2(x,y,Enu) + MZ**2)
    END FUNCTION afu

    FUNCTION bfu(x,y,Enu)
    use pars
    double precision :: x, y, Enu, bfu
    bfu = g2w/(1-s2t)*Rfu/(Q2(x,y,Enu) + MZ**2)
    END FUNCTION bfu

    FUNCTION afsb(x,y,Enu)
    use pars
    double precision :: x, y, Enu, afsb
    afsb = g2w/(1-s2t)*Lfd/(Q2(x,y,Enu) + MZ**2)
    END FUNCTION afsb

    FUNCTION bfsb(x,y,Enu)
    use pars
    double precision :: x, y, Enu, bfsb
    bfsb = g2w/(1-s2t)*Rfd/(Q2(x,y,Enu) + MZ**2)
    END FUNCTION bfsb

    FUNCTION afd(x,y,Enu)
    use pars
    double precision :: x, y, Enu, afd
    afd = g2w/(1-s2t)*Lfd/(Q2(x,y,Enu) + MZ**2)
    END FUNCTION afd

    FUNCTION bfd(x,y,Enu)
    use pars
    double precision :: x,y,Enu
    double precision :: bfd
    bfd = g2w/(1-s2t)*Rfd/(Q2(x,y,Enu) + MZ**2)
    END FUNCTION bfd





    !! END CROSS SECTION, BEGIN KINEMATICS

    !mandelstam s
	FUNCTION s(Enu)
	use pars
	double precision :: Enu, s
    s= 2.d0*mN*Enu
    END FUNCTION s

    FUNCTION selectron(Enu)
    use pars
    double precision :: Enu, selectron
    selectron= 2.d0*me*Enu
    END FUNCTION selectron


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
!    use heavimod  
    
    double precision s,t,u,x,y,xs,xsatt,totattnu,totattnubar
    
    INTEGER, INTENT(IN) :: NUMFUN
    REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: TS
    REAL(kind=stnd), DIMENSION(1:NUMFUN) :: Value
    x = TS(1)
    y = TS(2)

    x = exp(x)
    y = exp(y)
    !equal nu and nutau
    !NC bit is the same for everyone
    if (flavor .gt. 0 )then
        xs = dsdxdy_nu_NC(x,y,Enu)
    else
        xs = dsdxdy_nubar_NC(x,y,Enu)
    end if

    if (flavor .eq. 1) then !nue
        xs = xs + dsdxdy_nu_CC(x,y,Enu)
    else if (flavor .eq. -1) then !nuebar
        xs = xs + dsdxdy_nubar_CC(x,y,Enu)
    else if (flavor .eq. 2) then
        xs = xs + dsdxdy_nu_CC(x,y,Enu)
    else if (flavor .eq. -2) then
        xs = xs + dsdxdy_nubar_CC(x,y,Enu)
    else if (flavor .eq. 3) then
        xs = xs + dsdxdy_nu_CCtau(x,y,Enu)
    else if (flavor .eq. -3) then
        xs = xs + dsdxdy_nubar_CCtau(x,y,Enu)
    end if
!    print*, Enu,x,y,xs
    u = x*y*xs*hbarc**2

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
    INTEGER :: NEval, j,i
    integer, parameter :: nybins = 200
    REAL(kind=stnd), DIMENSION(1:n,0:n,1:m) :: Vertices
    REAL(kind=stnd), DIMENSION(1:l) ::  IntegralValue, AbsErr
    REAL(kind=stnd) ::  EpsRel, EpsAbs
    real :: x0,y0,x1,y1
    real :: starttime, finishtime
    LOGICAL :: Restart
    double precision dNdy(nybins),dNdy_p(nybins),dNdy_n(nybins)
    double precision Evec (nybins)

!    integer  nevarray(res,res)
    double precision :: z1,z2 !, logxpy,logxmy,logymx
    double precision :: exposure,logEmin,logEmax,dE
            character(len=1024) filename
            character(len=8) frmt
            character (len=1024) muout
!    double precision x,y,mu,sigma
! 	double precision y
!     double precision muplus
    REAL(kind=stnd), DIMENSION(2):: TS
    
!    NNPDF tools
    integer NumberReplicas

!    Parse commands


    call get_command_argument(1,leminin)
    call get_command_argument(2,lemaxin)
    call get_command_argument(3,flavorIn)
    call get_command_argument(4,filename)

    read( leminin,*) logEmin
    read( lemaxin,*) logEmax
    read( flavorIn,*) flavor


    dE = (logEmax-logEmin)/nybins

!    print*, "MLQ = ",MLQ, " GeV"
!    print*, "x1j = ",x1j, " "





    epsrel = 1.d-3
    epsabs = 1.d-100
     Restart= .false.

!model and experimental parameters



    RgType(1) = Cube
    x0 = log(1.d-7)
    y0 = log(1.d-7)
    x1 = 0.d0
    y1 = 0.d0

    			!		   x      log10(Enu)
    Vertices(1:n,0,1) = (/ x0, y0 /)
    Vertices(1:n,1,1) = (/ x1 ,y0 /)
    Vertices(1:n,2,1) = (/ x0 , y1  /)


! 1) Protons:

    call InitPDFsetByName("CT14qed_inc_proton")
    call numberPDF(NumberReplicas)


    do j=1,nybins
    call cpu_time(starttime)

    Enu = 10.d0**(logEmin+dE*dble(j-1))

      





      CALL CUBATR(n,l,F,m,Vertices(:,:,1:m),RgType,IntegralValue,AbsErr, &
           NEval=NEval,EpsRel=epsrel,EpsAbs=EpsAbs,Restart=Restart,MaxPts=50000000) !00000)
        !Print *,"-> Integral approximation = ",IntegralValue
        !Print *,"-> with estimated error < ",AbsErr
        !Print*," -> and relative error < " , EpsRel
        !Print *,"-> The number of integrand evaluations used = ",NEval
       ! Print *,""
        call cubatr()
            !TS(1) = 10.d0
            !TS(2) = 10.d0 
            !IntegralValue = F(1,TS)
    
            z1 = IntegralValue(1)
!             cin(i,j) = z1*((2.d0*muplus)**4)*x/y
!             nevarray(i,j) = neval
		evec(j) = Enu
		dNdy_p(j) = z1
    end do
    call InitPDFsetByName("CT14qed_inc_neutron")
    call numberPDF(NumberReplicas)
        do j=1,nybins
        call cpu_time(starttime)

        Enu = 10.d0**(logEmin+dE*dble(j-1))
! 2) Neutrons

  
      CALL CUBATR(n,l,F,m,Vertices(:,:,1:m),RgType,IntegralValue,AbsErr, &
           NEval=NEval,EpsRel=epsrel,EpsAbs=EpsAbs,Restart=Restart,MaxPts=50000000)
        !Print *,"-> Integral approximation = ",IntegralValue
        !Print *,"-> with estimated error < ",AbsErr
        !Print*," -> and relative error < " , EpsRel
        !Print *,"-> The number of integrand evaluations used = ",NEval
       ! Print *,""
        call cubatr()
            !TS(1) = 10.d0
            !TS(2) = 10.d0 
            !IntegralValue = F(1,TS)
    
            z2 = IntegralValue(1)
!             cin(i,j) = z1*((2.d0*muplus)**4)*x/y
!             nevarray(i,j) = neval
		evec(j) = Enu
		dNdy_n(j) = z2		

        if (flavor .eq. -1) then
            dNdy_n(j) = dNdy_n(j) + sigeh(Enu)*hbarc**2     !doesn't matter where we put this
        end if
		
            !print*, "x: ", x
        print*, Enu, z1+z2
!             print*, "z1: ", z1
!     		print*,Edep,flux(Edep,-2.6d0)
		call cpu_time(finishtime)
		print*,"Time = ",finishtime-starttime," seconds."
    end do

    dNdy = (dNdy_n + dNdy_p)/2 !average protons & neutrons (electrons get divided by 2, which is correct)

!    dNdy = (10.*dNdy_p + 8.*dNdy_n)/18.
    frmt = '(I3.3)'
            open(unit=55,file=filename,action='write')
            do i=1,nybins
            write(55,*) dNdy(i) !, dNdy(i)
            end do
            close(55)
    print*, "Which is cool"
    STOP
    END PROGRAM Example2D
