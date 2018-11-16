clc
clear all
%%photonuclear, ALLM parameterization
%parameters for ALLM
global R cp1 cp2 cp3 ap1 ap2 ap3 bp1 bp2 bp3 cr1 cr2 cr3 ar1 ar2 ar3
global br1 br2 br3 m02 mp2 mr2 Q02 L2 alpha2 M mtau mpi
R=0;
cp1=0.28067;
cp2=0.22291;
cp3=2.1979;
ap1=-0.0808;
ap2=-0.44812;
ap3=1.1709;
bp1=0.36292;
bp2=1.8917;
bp3=1.8439;
cr1=0.80107;
cr2=0.97307;
cr3=3.4942;
ar1=0.58400;
ar2=0.37888;
ar3=2.6063;
br1=0.01147;
br2=3.7582;
br3=0.49338;
m02=0.31985;%in GeV^2
mp2=49.457;%in GeV^2
mr2=0.15052;%in GeV^2
Q02=0.52544;%in GeV^2
L2=0.06527;%in GeV^2
%other parameters
hbar=6.582e-25;%in GeV*s
c=3e10;%in cm/s
alpha2=(1./137)^2;
mp=0.9382;
mn=0.9396;
M=(mp+mn)/2;
mtau=1.777;%tau mass
mpi=0.1396;%pion mass
A=18;
Z=8;

%precision of solution
abserr=1e-10;
relerr=1e-6;


NumE=200;
Etau=logspace(3,10,NumE);%incoming tau
Etauout=logspace(3,10,NumE);%outgoing tau
dsigmadE=zeros(NumE);
for i=1:NumE
    E=Etau(i);
    ymin=(mpi+mpi^2/(2*M))/E;
    ymax=1-M/(2*E)*(1+mtau^2/M^2);
    for j=1:NumE
        Eout=Etauout(j);
        y=(E-Eout)/E;
        if y>ymin && y<ymax
            Q2min=mtau^2*(E*y)^2/(E*Eout)-mtau^4/(2*E*Eout);
            Q2max=2*M*(E*y-mpi)-mpi^2;
            dsig=@(Q2) getsigma_pn(Q2,y,E,Z,A)/E;%convert from dsigmady to dsigmadE
            dsigmadE(i,j)=integral(dsig,Q2min,Q2max,'AbsTol',abserr,'RelTol',relerr);
        end
    end
end
dsigmadE=dsigmadE*(hbar*c)^2;%convert from GeV^-3 to cm^2*GeV^-1
