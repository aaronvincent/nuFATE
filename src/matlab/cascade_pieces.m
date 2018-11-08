function [phi_1,energy_nodes] = cascade_pieces(varargin) %this is get_eigs
%%regeneration bits are included here and three nus+tau are solved at the 
%%same time. For treatment of on-the-spot decay of tau or without
%%regenerations, refer to the old code
%%input: flavor (only specify neutrino or antineutrino), length of path in
%%a sublayer-x (in cm), average density the
%%sublayer-rhobar (in g*cm^-3), input flux-phi_0, cross section file (optional)
%%output: output flux-phi_1, energy nodes-energy_nodes

flavor = varargin{1};
if flavor==0
    error(['You must specify a nonzero number for the flavor, a positive'...
    ' number for neutrinos or a negtive number for antineutrinos']);
end
x = varargin{2};
rhobar = varargin{3};
% defaults:
logemin = 3; %Min log energy (GeV) (do not touch unless you have recomputed cross sections)
logemax = 10; % Max log E
NumNodes = 200;

phi_0 = varargin{4}; %flux computed from the previous step
[m, n] = size(phi_0);
if m<n
    phi_0 = phi_0'; %turn into a column vector
end

if nargin >= 5
    xsecfname = varargin{5};
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


energy_nodes = logspace(logemin,logemax,NumNodes);
dsigmady = h5read(xsecfname, dxs_fname)';
RHSMatrix_NC = get_matrices(energy_nodes,dsigmady);
RHSMatrix = cell(1,4);
for i=1:3
    sigma_array = h5read(xsecfname,sigma_fname{i})';
    RHSMatrix{i} = -diag(sigma_array) + RHSMatrix_NC;
    if flavor < 0 && i == 1
        %factor of 2 since there is 1 electron per two nucleons
        RHSMatrix{i} = RHSMatrix{i} - diag(get_glashow_total(energy_nodes)/2)...
            + get_glashow_partial(energy_nodes)/2;
    end
end

%%compute tau decay branching ratios, can precompute and add to h5 file
%%to speed up
NA = 6.02e+23; %Avogadro's constant
c = 3.e+10; %speed of light in cm/s
mtau = 1.777; %tau mass in GeV
ttau = 2.906e-13; %tau lifetime in rest frame in s
RHSMatrix{4}=-diag(mtau./(NA*c*ttau*rhobar*energy_nodes)); %tau decay
dndz_tau = zeros(NumNodes, NumNodes);
dndz_emu = zeros(NumNodes, NumNodes);
for i=1:NumNodes
    Etau=energy_nodes(i);
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
dndz_tau=dndz_tau*mtau/(NA*c*ttau*rhobar);
dndz_emu=dndz_emu*mtau/(NA*c*ttau*rhobar);
RHSMatrix_34 = get_matrices(energy_nodes, dndz_tau); %nutau from tau decay
RHSMatrix_14 = get_matrices(energy_nodes, dndz_emu); %nue from tau decay
RHSMatrix_24 = RHSMatrix_14; %numu from tau decay
%tau production from nutau CC interactions
if flavor>0
    dtauCC=load('../../resources/CT14/differential_cross_sections/dxstau.dat');
else
    dtauCC=load('../../resources/CT14/differential_cross_sections/dxstaubar.dat');
end
RHSMatrix_43 = get_matrices(energy_nodes, dtauCC); %tau from nutau CC

%%construct the huge RHS Matrix
z0=zeros(NumNodes);
M=[RHSMatrix{1} z0 z0 RHSMatrix_14;
   z0 RHSMatrix{2} z0 RHSMatrix_24;
   z0 z0 RHSMatrix{3} RHSMatrix_34;
   z0 z0 RHSMatrix_43 RHSMatrix{4}];

%phi_0 = energy_nodes.^(2-g)';
[v,w] = eig(M);
ci = (v^-1)*phi_0;
w = diag(w);
t = rhobar*x*NA; %optical depth
phi_1 = v*(ci.*exp(w*t));

end

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
dsig = GF^2*selectron/pi.*(t1 + (t2.^2+t3.^2).*(1-y).^2)*hbarc^2.*heaviside(y);
dsig = dsig.*(1-y)./Enuin; %dy --> d\tilde E
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
