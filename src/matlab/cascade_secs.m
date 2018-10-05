function [w,v,ci,energy_nodes] = cascade_secs(varargin) %this is get_eigs

flavor = varargin{1};
% defaults:
g = 2;
logemin = 3;
logemax = 10;
NumNodes = 200;

if nargin >= 2
    g = varargin{2};
end
% if nargin >= 3 %you have changed the cross sections
%     logemin = varargin{3};
%     logemax = varargin{4};
%     NumNodes = varargin{5};
% end

if nargin >= 3
    xsecfname = varargin{3};
else
    xsecfname = '../../resources/NuFATECrossSections.h5';
end

if flavor==-1
    sigma_fname = '/total_cross_sections/nuebarxs';
elseif flavor == -2
    sigma_fname = '/total_cross_sections/numubarxs';
    
elseif flavor == 1
    sigma_fname = '/total_cross_sections/nuexs';
elseif flavor == 2
    sigma_fname = '/total_cross_sections/numuxs';
else
    error('flavor for secs must be +/- 1 or 2 ')
end
if flavor > 0
    dxs_fname = '/differential_cross_sections/dxsnu';
    sig3fname = '/total_cross_sections/nutauxs';
    secname = '/tau_decay_spectrum/secfull';
    regenname = '/tau_decay_spectrum/tfull';
else
    dxs_fname = '/differential_cross_sections/dxsnubar';
    sig3fname = '/total_cross_sections/nutaubarxs';
    secname = '/tau_decay_spectrum/secbarfull';
    regenname = '/tau_decay_spectrum/tbarfull';
end

energy_nodes = logspace(logemin,logemax,NumNodes);
[RHSMatrix] = get_RHS_matrices(energy_nodes,sigma_fname,sig3fname,dxs_fname,secname,regenname,xsecfname);




if flavor == -1 %add glashow piece
    Z = zeros(NumNodes);
    glashow_piece = (-diag(get_glashow_total(energy_nodes))+ get_glashow_partial(energy_nodes))/2.;
    RHSMatrix = RHSMatrix + [glashow_piece Z; Z Z];
end
phi_0 = energy_nodes.^(2-g)';
[v,w] = eig(RHSMatrix);
ci = (v^-1)*[phi_0; phi_0];
w = diag(w);

end

function  [RHSMatrix, sigma_array] = get_RHS_matrices(energy_nodes,sigma_fname,sig3fname,dxs_fname,secname,regenname,xsecfname)
NumNodes = length(energy_nodes);
h5flag = 1;
if h5flag
    sigma_array1 = h5read(xsecfname,sigma_fname)';
    sigma_array2 = h5read(xsecfname,sig3fname)';
    dsigmady = h5read(xsecfname,dxs_fname)';
    emuregen = h5read(xsecfname,secname)';
    tauregen = h5read(xsecfname,regenname)';
else
    sigma_array1 = load(['../data/' sigma_fname]);
    sigma_array2 = load(['../data/' sig3fname]);
    dsigmady = load(['../data/' dxs_fname]);
    emuregen = load(['../data/' secname]);
    tauregen = load(['../data/' regenname]);
end
DeltaE = diff(log(energy_nodes));
RHSMatrix1 = zeros(NumNodes);
RHSMatrix2 = zeros(NumNodes);
RHSMatrix3 = zeros(NumNodes);
RHSMatrix4 = zeros(NumNodes);


%     1: Regular NC bit
for i = 1:NumNodes
    for j = i+1:NumNodes
        RHSMatrix1(i,j) = DeltaE(j-1)*dsigmady(j,i)*energy_nodes(j).^-1*energy_nodes(i).^2;
    end
end

RHSMatrix1 = -diag(sigma_array1) + RHSMatrix1;
%    2: nue/numu secondary production
for i = 1:NumNodes
    for j = i+1:NumNodes
        RHSMatrix2(i,j) = DeltaE(j-1)*emuregen(j,i)*energy_nodes(j).^-1*energy_nodes(i).^2;
    end
end
%     3 is zero; 4 is tau regen
for i = 1:NumNodes
    for j = i+1:NumNodes
        RHSMatrix4(i,j) = DeltaE(j-1)*(dsigmady(j,i)+tauregen(j,i))*energy_nodes(j).^-1*energy_nodes(i).^2;
    end
end
RHSMatrix4 = -diag(sigma_array2) + RHSMatrix4;
RHSMatrix = [RHSMatrix1 RHSMatrix2; RHSMatrix3 RHSMatrix4];
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
s2t = 0.231;
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
