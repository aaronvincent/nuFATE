clc
clear all
%%electron pair production
%parameters
global mtau me alpha2 re2 gamma1 gamma2
gamma1=1.95e-5;
gamma2=5.3e-5;
mtau=1.777;%tau mass
me=0.511e-3;%electron mass
alpha2=(1./137)^2;
re2=(2.8179e-13)^2;%in cm
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
    ymin=4*me/E;
    ymax=1-3*sqrt(exp(1))/4*mtau/E*Z^(1/3);
    for j=1:NumE
        Eout=Etauout(j);
        y=(E-Eout)/E;
        if y>ymin && y<ymax
            rhomax=sqrt(1-4*me/E/y)*(1-6*mtau^2/E^2/(1-y));
            dsig=@(rho) getsigma_epp(rho,y,E,Z)/E;%convert from dsigmady to dsigmadE
            dsigmadE(i,j)=2*integral(dsig,0,rhomax,'AbsTol',abserr,'RelTol',relerr);
        end
    end
end

dsigmadE;

% Z=1:100;
% for i=1:length(Z)
%     B(i)=getB(Z(i));
% end
% plot(Z,B);
% xlabel('Z');
% ylabel('B');
