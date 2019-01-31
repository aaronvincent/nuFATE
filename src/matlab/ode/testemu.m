clear all;
close all;

N=200;
lEmin=3;
lEmax=10;
dlg=(lEmax-lEmin)/(N-1.);
Ei=logspace(lEmin-dlg/2,lEmax+dlg/2,N+1);

%% get the spectrum for a single direction
g = 0; %spectral index
% defaults:
logemin = 3; %Min log energy (GeV) (do not touch unless you have recomputed cross sections)
logemax = 10; % Max log E

NA = 6.02e+23; %Avogadro's constant
energy_nodes=logspace(lEmin,lEmax,N);

c = 3.e+10; %speed of light in cm/s
mtau = 1.777; %tau mass in GeV
ttau = 2.906e-13; %tau lifetime in rest frame in s
rhobar=0.917;


%tau energy loss in propogation through ice
xs_epp_O=load('epp/O.txt');
xs_epp_H=load('epp/H.txt');
xs_pn_O=load('pn/O.txt');
xs_pn_H=load('pn/H.txt');
dsigma_crust=16/18*xs_epp_O/16+2/18*xs_epp_H/1+16/18*xs_pn_O/16+2/18*xs_pn_H;
%uncomment this line to remove the diagonal entries
dsigma_crust=dsigma_crust-diag(diag(dsigma_crust));
RHSMatrix_crust=get_matrices(Ei,energy_nodes,dsigma_crust);

%total loss cross section
diffE=diff(Ei);
xsdiag_crust=sum(bsxfun(@times,diffE,dsigma_crust),2);
RHSMatrix_crust=-diag(xsdiag_crust)+RHSMatrix_crust;

RHSMatrix_tau=RHSMatrix_crust;
%uncomment this line to include tau decay
RHSMatrix_tau=RHSMatrix_tau-diag(mtau./(NA*c*ttau*rhobar*energy_nodes)); %tau decay

phi_0 = energy_nodes.^(2-g)';
%start with E=1e10 GeV
phi_0(1:end-1)=0;
%10km
x=10*1e5;
t = rhobar*x*NA; %optical depth
M=RHSMatrix_tau;
odefun=@(t,phi) M*phi;
[ts,phis]=ode15s(odefun,[0,t],phi_0);
phi_1=phis(end,:);

loglog(energy_nodes,phi_1./energy_nodes.^2);
hold on

%%from proposal
figure
%data=load('~/PROPOSAL/src/output/log_10km_1.1582e7_1e4.txt');
%data=load('~/PROPOSAL/src/output/log_10km_2.0998e8_1e4.txt');
%data=load('~/PROPOSAL/src/output/log_10km_1.4491e9_1e4.txt');
%data=load('~/PROPOSAL/src/output/log_10km_ALLM_NOHARD_CONT_1e5.txt');
data=load('log_10km_1e5.txt');
%data=load('~/PROPOSAL/src/output/log_1km_1e8_1e4.txt');
%data=load('~/PROPOSAL/src/output/log_10km_rand.txt');

Erebin=Ei;
%data_eff=data(:,3)./1e3;
data_eff=data(:,2)./1e3;
count=zeros(1,length(Erebin)-1);
for i=1:length(data_eff)
    for j=1:length(Erebin)-1
        if data_eff(i)>=Erebin(j) && data_eff(i)<Erebin(j+1)
            count(j)=count(j)+1;
        end
    end
end

diffE_nuFate=diff(Ei);
norm_nuFate=diffE_nuFate(:,end)*1;
%norm_nuFate=sum(diffE_nuFate);
phi_nuFate=phi_1./energy_nodes.^2/norm_nuFate;
loglog(energy_nodes,phi_nuFate);
hold on;

diffE=diff(Erebin);
difflogE=diff(log10(Erebin));
count_norm=count./diffE;
%count_norm=count./difflogE;
count_norm=count_norm/length(data_eff);
%count_norm=count_norm*max(phi_nuFate)/max(count_norm);
loglog(energy_nodes,count_norm,'-r');


%nuFate
N1=100;
N2=10;
lEmin=3;
lEmid=6;
lEmax=10;
dlg1=(lEmid-lEmin)/(N1-1.);
a=lEmid+dlg1/2;
b=(2*lEmax*N2-a)/(2*N2-1);
dlg2=(b-a)/N2;
Ei1=logspace(lEmin-dlg1/2,lEmid+dlg1/2,N1+1);
Ei2=logspace(a,b,N2+1);
Ei=[Ei1(1:end-1) Ei2];
Erec1=log10(Ei1)+dlg1/2;
Erec1(end)=[];
Erec2=log10(Ei2)-dlg2/2;
Erec2(1)=[];
Erec=[Erec1 Erec2];
energy_tau=10.^Erec;

dsigma_crust=load('epp/O_110.txt')+load('pn/O_110.txt');
dsigma_crust=dsigma_crust/16;
dsigma_crust=dsigma_crust-diag(diag(dsigma_crust));
RHSMatrix_crust=get_matrices(Ei,energy_tau,dsigma_crust);

%total loss cross section
diffE=diff(Ei);
xsdiag_crust=sum(bsxfun(@times,diffE,dsigma_crust),2);
RHSMatrix_crust=-diag(xsdiag_crust)+RHSMatrix_crust;

RHSMatrix_tau=RHSMatrix_crust;
RHSMatrix_tau=RHSMatrix_tau-diag(mtau./(NA*c*ttau*rhobar*energy_tau)); %tau decay

phi_0 = energy_tau.^(2-g)';
phi_0(1:end-1)=0;
x=10*1e5;
M=RHSMatrix_tau;
[v,w] = eig(M);
ci = v^-1*phi_0;
w = diag(w);
t = rhobar*x*NA; %optical depth
phi_110 = v*(ci.*exp(w*t));
phi_out=phi_110./energy_tau'.^2/diffE(end);
loglog(energy_tau,phi_out,'k');
legend('ODE','PROPOSAL','nuFATE','location','northwest');

function  RHSMatrix = get_matrices(energy_edges,energy_nodes,dsigmady)
NumNodes = length(dsigmady);
DeltaE = diff(log(energy_edges));
RHSMatrix = zeros(NumNodes);
for i = 1:NumNodes
    for j = i+1:NumNodes
        RHSMatrix(i,j) = DeltaE(j)*dsigmady(j,i)*energy_nodes(j).^-1*energy_nodes(i).^2;
        %RHSMatrix(i,j) = DeltaE(j-1)*dsigmady(j,i)*energy_nodes(j);
    end
end
end
