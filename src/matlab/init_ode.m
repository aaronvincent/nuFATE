function [RHSMatrix,energy_nodes,energy_tau] = init_ode(flavor,varargin)
%%precompute the matrix element to speed up
%%input: flavor (only specify neutrino or antineutrino),
%%cross section file (optional)
%%output: RHS matri elements RHSMatrix, energy nodes-energy_nodes

if flavor==0
    error(['You must specify a nonzero number for the flavor, a positive'...
    ' number for neutrinos or a negtive number for antineutrinos']);
end
% defaults:
logemin = 3; %Min log energy (GeV) (do not touch unless you have recomputed cross sections)
logemax = 10; % Max log E
NumNodes = 200;
%NumTau = 20; %number of energy bins for tau propogation


NumTau=200;
dlg=(logemax-logemin)/(NumTau-1);
Ei=logspace(logemin-dlg/2,logemax+dlg/2,NumTau+1);

if nargin >= 2
    xsecfname = varargin{2};
else
    xsecfname = '../../resources/NuFATECrossSections.h5';
end
%get cross section locations
if flavor>0 %neutrinos
    sigma_fname={'/total_cross_sections/nuexs', ...
                 '/total_cross_sections/numuxs', ...
                 '/total_cross_sections/nutauxs'};
    dxs_fname = '/differential_cross_sections/dxsnu';
else %antineutrinos
    sigma_fname={'/total_cross_sections/nuebarxs', ...
                 '/total_cross_sections/numubarxs', ...
                 '/total_cross_sections/nutaubarxs'};
    dxs_fname = '/differential_cross_sections/dxsnubar';
end
elossfname='../../resources/eloss.h5';
elements={'O','Si','Al','Fe','Ca','Na','K','Mg','Ni','S'};
%atomic number
A=[16,28,27,56,40,23,39,24,58.7,32];
frac_crust=[0.467,0.277,0.08,0.05,0.03,0.03,0.03,0.02,0,0];
frac_mantle=[0.448,0.215,0.02,0.06,0.02,0,0,0.228,0,0];
frac_core=[0,0,0,0.89,0,0,0,0,0.06,0.05];

energy_nodes = logspace(logemin,logemax,NumNodes);
dsigmady = h5read(xsecfname, dxs_fname)';
RHSMatrix_NC = get_matrices(energy_nodes,energy_nodes,dsigmady);
RHSMatrix = cell(1,9);
for i=1:3
    sigma_array = h5read(xsecfname,sigma_fname{i})';
    RHSMatrix{i} = -diag(sigma_array) + RHSMatrix_NC;
    if flavor < 0 && i == 1
        %factor of 2 since there is 1 electron per two nucleons
        RHSMatrix{i} = RHSMatrix{i} - diag(get_glashow_total(energy_nodes)/2)...
            + get_glashow_partial(energy_nodes)/2;
    end
end

%tau energy loss in propogation
energy_tau = logspace(logemin,logemax,NumTau);
dsigma_crust=zeros(NumTau);
dsigma_mantle=zeros(NumTau);
dsigma_core=zeros(NumTau);
for i=1:length(elements)
    dsigma_crust=dsigma_crust+h5read(elossfname,['/epp/',elements{i}])'*frac_crust(i)/A(i);
    dsigma_crust=dsigma_crust+h5read(elossfname,['/pn/',elements{i}])'*frac_crust(i)/A(i);
    dsigma_mantle=dsigma_mantle+h5read(elossfname,['/epp/',elements{i}])'*frac_mantle(i)/A(i);
    dsigma_mantle=dsigma_mantle+h5read(elossfname,['/pn/',elements{i}])'*frac_mantle(i)/A(i);
    dsigma_core=dsigma_core+h5read(elossfname,['/epp/',elements{i}])'*frac_core(i)/A(i);
    dsigma_core=dsigma_core+h5read(elossfname,['/pn/',elements{i}])'*frac_core(i)/A(i);
end
%remove the diagonal entries which don't actually count
dsigma_crust=dsigma_crust-diag(diag(dsigma_crust));
dsigma_mantle=dsigma_mantle-diag(diag(dsigma_mantle));
dsigma_core=dsigma_core-diag(diag(dsigma_core));
RHSMatrix_crust=get_matrices_tau(Ei,energy_tau,dsigma_crust);
RHSMatrix_mantle=get_matrices_tau(Ei,energy_tau,dsigma_mantle);
RHSMatrix_core=get_matrices_tau(Ei,energy_tau,dsigma_core);
%total loss cross section
diffE=diff(Ei);
xsdiag_crust=sum(bsxfun(@times,diffE,dsigma_crust),2);
RHSMatrix_crust=-diag(xsdiag_crust)+RHSMatrix_crust;
xsdiag_mantle=sum(bsxfun(@times,diffE,dsigma_mantle),2);
RHSMatrix_mantle=-diag(xsdiag_mantle)+RHSMatrix_mantle;
xsdiag_core=sum(bsxfun(@times,diffE,dsigma_core),2);
RHSMatrix_core=-diag(xsdiag_core)+RHSMatrix_core;
RHSMatrix{4}=RHSMatrix_crust;
RHSMatrix{5}=RHSMatrix_mantle;
RHSMatrix{6}=RHSMatrix_core;

%tau production from nutau CC interactions
if flavor>0
    dtauCC=load('../../resources/CT14/differential_cross_sections/dxstau.dat');
else
    dtauCC=load('../../resources/CT14/differential_cross_sections/dxstaubar.dat');
end
RHSMatrix_43 = get_matrices(energy_nodes,energy_tau, dtauCC);
RHSMatrix{7}=RHSMatrix_43;

%%compute tau decay branching ratios
dndz_tau = zeros(NumTau,NumNodes);
dndz_emu = zeros(NumTau,NumNodes);
for i=1:NumTau
    Etau=energy_tau(i);
    for j=1:NumNodes
        Enutau=energy_nodes(j);
        z=Enutau/Etau;
        if z<1
           if flavor>0
            dndz_tau(i,j)=dntaudE(z,0,-1)+dntaudE(z,1,-1);
           else
            dndz_tau(i,j)=dntaudE(z,0,1)+dntaudE(z,1,1);
           end
           dndz_emu(i,j)=dnnottaudE(z);
        end
        dndz_tau(i,j)=dndz_tau(i,j)/Etau^2;
        dndz_emu(i,j)=dndz_emu(i,j)/Etau^2;
    end
end
RHSMatrix_34 = get_matrices(energy_tau,energy_nodes, dndz_tau);
RHSMatrix_14 = get_matrices(energy_tau,energy_nodes, dndz_emu);
%RHSMatrix_24 = RHSMatrix_14; %numu from tau decay
RHSMatrix{8}=RHSMatrix_34;%NA*c*ttau*rhobar/mtau is to be divided by later
RHSMatrix{9}=RHSMatrix_14;%NA*c*ttau*rhobar/mtau is to be divided by later

end

function  RHSMatrix = get_matrices(energy_in,energy_out,dsigmady)
[NumIN,NumOUT]=size(dsigmady);
DeltaE = diff(log(energy_in));
RHSMatrix = zeros(NumOUT,NumIN);
for i = 1:NumOUT
    for j = 2:NumIN
        RHSMatrix(i,j) = DeltaE(j-1)*dsigmady(j,i)*energy_in(j).^-1*energy_out(i).^2;
    end
end
end

function  RHSMatrix = get_matrices_tau(energy_edges,energy_nodes,dsigmady)
NumNodes = length(dsigmady);
DeltaE = diff(log(energy_edges));
RHSMatrix = zeros(NumNodes);
for i = 1:NumNodes
    for j = i+1:NumNodes
        RHSMatrix(i,j) = DeltaE(j)*dsigmady(j,i)*energy_nodes(j).^-1*energy_nodes(i).^2;
    end
end
end



function sig = get_glashow_total(energy_nodes)
Enu = energy_nodes;
GF = 1.16d-5;
hbarc=1.97d-14;
GW = 2.085;
MW = 80.385d0;
mmu=.106d0;
me=511.d-6;
selectron= 2.d0*me*Enu;
sig = 1.d0/3.d0*GF.^2.*selectron./pi.*(1.d0-(mmu.^2-me.^2)./selectron).^2 ...
    ./((1.d0-selectron/MW.^2).^2+GW.^2/MW.^2)*.676/.1057 * hbarc^2;
end

function dsig = get_glashow_partial(energy_nodes)
[Enuin, Enu] = meshgrid(energy_nodes,energy_nodes');
y = 1-Enu./Enuin;
GF = 1.16d-5;
hbarc=1.97d-14;
GW = 2.085;
MW = 80.385d0;
MZ = 91.18;
me=511.d-6;
s2t = 0.23;
gL =  s2t-0.5;
gR = s2t;
selectron= 2.d0*me*Enuin;
den = (1-selectron/MW^2).^2 + GW^2/MW^2;
t1 = gR^2./(1.+y.*selectron/MZ^2).^2;
t2 = gL./(1.+y.*selectron./MZ^2) + (1-selectron/MW^2)./den;
t3 = GW/MW./den;
dsig = GF^2*selectron/pi.*(t1 + (t2.^2+t3.^2).*(1-y).^2)*hbarc^2.*heavi(y);
dsig = dsig./Enuin;
end

function dndz = dntaudE(z,channel,pol)
%%tau- (from nutau) --> pol = -1; nutaubar --> pol = +1
if channel == 0  %electron channel, to be multiplied by 2
   g0 = 5./3. - 3.*z^2 + 4./3.*z^3;
   g1 = 1./3. - 3.*z^2 + 8./3.*z^3;
   dndz =  0.18*(g0+pol*g1);
elseif channel == 1 %hadron channel
    dndz = 0.;
    %1)pions
    r = 0.1395^2/1.776^2;
    if 1. > r + z
       g0 = 1./(1.-r); %*heaviside(1.d0-r-z)
       g1 = -(2.*z-1.+r)/(1.-r)^2; %*heaviside(1.d0-r-z)
       dndz = dndz + 0.12*(g0+pol*g1);
       %2) rho
       r = 0.77^2/1.776^2;
       if 1. > r+z
           g0 = 1./(1.-r); %*heaviside(1.d0-r-z)
           g1 = -(2.*z-1.+r)/(1.-r)*(1.-2*r)/(1+2*r); %*heaviside(1.d0-r-z)
           dndz = dndz  + 0.26*(g0+pol*g1);
           %3) a1 (yeah this is a thing!Why?)
           r = 1.26^2/1.776^2;
           if 1. > r+z
               g0 = 1./(1.-r); %*heaviside(1.d0-r-z)
               g1 = -(2.*z-1.+r)/(1.-r)*(1.-2*r)/(1+2*r); %*heaviside(1.d0-r-z)
               dndz = dndz + 0.13*(g0+pol*g1);
           end
       end
    end
    %4) X (everything else hadronic)
    if z > 0.3
    g0 = 1./0.30; %*heaviside(0.3-z)
    dndz = dndz + 0.13*g0;
    end
end
end

function dndz = dnnottaudE(z)
%nue or numu from tau decay
dndz = 0.18*(4.0-12.*z+12.*z^2-4.*z^3);
end


function he = heavi(y)
he = zeros(size(y));
[m,n]=size(y);
for i=1:m
    for j=1:n
        if y(i,j)>=0
            he(i,j)=1;
        else
            he(i,j)=0;
        end
    end
end
end
