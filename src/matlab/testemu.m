clear all;
close all;

%% get the spectrum for a single direction
g = 2; %spectral index
zenith = 91; %choose a zenith.
% defaults:
logemin = 3; %Min log energy (GeV) (do not touch unless you have recomputed cross sections)
logemax = 10; % Max log E
NumNodes = 200;
NA = 6.02e+23; %Avogadro's constant
energy_nodes = logspace(logemin,logemax,NumNodes);

NA = 6.02e+23; %Avogadro's constant
c = 3.e+10; %speed of light in cm/s
mtau = 1.777; %tau mass in GeV
ttau = 2.906e-13; %tau lifetime in rest frame in s
rhobar=2;

elossfname='../../resources/eloss.h5';
elements={'O','Si','Al','Fe','Ca','Na','K','Mg','Ni','S'};
frac_crust=[0.467,0.277,0.08,0.05,0.03,0.03,0.03,0.02,0,0];
frac_mantle=[0.448,0.215,0.02,0.06,0.02,0,0,0.228,0,0];
frac_core=[0,0,0,0.89,0,0,0,0,0.06,0.05];

xsecfname = '../../resources/NuFATECrossSections.h5';
dxs_fname = '/differential_cross_sections/dxsnubar';
sigma_fname={'/total_cross_sections/nuebarxs', ...
             '/total_cross_sections/numubarxs', ...
             '/total_cross_sections/nutaubarxs'};
dsigmady = h5read(xsecfname, dxs_fname)';
RHSMatrix_NC = get_matrices(energy_nodes,dsigmady);
%sigma_array = trapz(energy_nodes,dsigmady');
sigma_array = h5read(xsecfname,sigma_fname{2})';
RHSMatrix = -diag(sigma_array) + RHSMatrix_NC;

%tau energy loss in propogation
dsigma_crust=zeros(NumNodes);
dsigma_mantle=zeros(NumNodes);
dsigma_core=zeros(NumNodes);
for i=1:length(elements)
    dsigma_crust=dsigma_crust+h5read(elossfname,['/epp/',elements{i}])'*frac_crust(i);
    dsigma_crust=dsigma_crust+h5read(elossfname,['/pn/',elements{i}])'*frac_crust(i);
    dsigma_mantle=dsigma_mantle+h5read(elossfname,['/epp/',elements{i}])'*frac_mantle(i);
    dsigma_mantle=dsigma_mantle+h5read(elossfname,['/pn/',elements{i}])'*frac_mantle(i);
    dsigma_core=dsigma_core+h5read(elossfname,['/epp/',elements{i}])'*frac_core(i);
    dsigma_core=dsigma_core+h5read(elossfname,['/pn/',elements{i}])'*frac_core(i);
end
%uncomment this line to remove the diagonal entries
%dsigma_crust=dsigma_crust-diag(diag(dsigma_crust));
RHSMatrix_crust=get_matrices(energy_nodes,dsigma_crust);
RHSMatrix_mantle=get_matrices(energy_nodes,dsigma_mantle);
RHSMatrix_core=get_matrices(energy_nodes,dsigma_core);

%total loss cross section
RHSMatrix_crust=-diag(trapz(energy_nodes,dsigma_crust'))+RHSMatrix_crust;
RHSMatrix_mantle=-diag(trapz(energy_nodes,dsigma_mantle'))+RHSMatrix_mantle;
RHSMatrix_core=-diag(trapz(energy_nodes,dsigma_core'))+RHSMatrix_core;
RHSMatrix_tau=RHSMatrix_crust;
%uncomment this line to include tau decay
%RHSMatrix_tau=RHSMatrix_tau-diag(mtau./(NA*c*ttau*rhobar*energy_nodes)); %tau decay

y=(energy_nodes-energy_nodes')./energy_nodes;
beta=trapz(energy_nodes,y'.*dsigma_crust,2);
sig=trapz(energy_nodes,dsigma_crust');

figure
loglog(energy_nodes,beta/20,'-r',energy_nodes,sig,'-b');hold on;
loglog(energy_nodes,mtau./(NA*c*ttau*rhobar*energy_nodes),'-k');
grid on;

%plot the matrix as a function of in-going and out-going particles
% figure
% [X,Y]=meshgrid(energy_nodes,energy_nodes);
% surf(X,Y,dsigma_crust);
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% set(gca,'zscale','log');

%solve for tau flux with energy loss only
phi_0 = energy_nodes.^(2-g)';
x=0.0001*1e5;%cm
M=RHSMatrix_tau;
[v,w] = eig(M);
ci = v\phi_0;
w = diag(w);
t = rhobar*x*NA; %optical depth
phi_1 = v*(ci.*exp(w*t));

figure
loglog(energy_nodes,phi_0./energy_nodes'.^2,'-r',energy_nodes,phi_1./energy_nodes'.^2,'-b');
hold on

%excluding gain
M=-diag(trapz(energy_nodes,dsigma_crust'));
[v,w] = eig(M);
ci = (v^-1)*phi_0;
w = diag(w);
t = rhobar*x*NA; %optical depth
phi_2 = v*(ci.*exp(w*t));

loglog(energy_nodes,phi_2./energy_nodes'.^2,'--g');



function  RHSMatrix = get_matrices(energy_nodes,dsigmady)
NumNodes = length(energy_nodes);
DeltaE = diff(log(energy_nodes));
RHSMatrix = zeros(NumNodes);
for i = 1:NumNodes
    for j = i+1:NumNodes
        RHSMatrix(i,j) = DeltaE(j-1)*dsigmady(j,i)*energy_nodes(j).^-1*energy_nodes(i).^2;
    end
end
end