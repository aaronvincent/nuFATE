function [w,v,ci,energy_nodes] = cascade(varargin) %this is get_eigs
flavor = varargin{1};
% defaults:
g = 2;
logemin = 3; %Min log energy (GeV) (do not touch unless you have recomputed cross sections)
logemax = 10; % Max log E
NumNodes = 200;

%flavor=3;
%nargin=1;

if nargin >= 2
    g = varargin{2};
end
% if nargin >= 3 %you shouldn't be specifying only one of these:
%     logemin = varargin{3};
%     logemax = varargin{4};
%     NumNodes = varargin{5};
% end

if nargin >= 3
    xsecfname = varargin{3};
else
    xsecfname = '../../../resources/NuFATECrossSections.h5';
end
%get cross section locations
if flavor==-1
    sigma_fname = '/total_cross_sections/nuebarxs';
elseif flavor == -2
    sigma_fname = '/total_cross_sections/numubarxs';
elseif flavor == -3
    sigma_fname = '/total_cross_sections/nutaubarxs';
elseif flavor == 1
    sigma_fname = '/total_cross_sections/nuexs';
elseif flavor == 2
    sigma_fname = '/total_cross_sections/numuxs';
elseif flavor == 3
    sigma_fname = '/total_cross_sections/nutauxs';
end
if flavor > 0
    dxs_fname = '/differential_cross_sections/dxsnu';
else
    dxs_fname = '/differential_cross_sections/dxsnubar';
end

energy_nodes = logspace(logemin,logemax,NumNodes);
[RHSMatrix, sigma_array] = get_RHS_matrices(energy_nodes,sigma_fname,dxs_fname,xsecfname);
if flavor ==-3
    [RHregen,~] = get_RHS_matrices(energy_nodes,sigma_fname,'/tau_decay_spectrum/tbarfull',xsecfname);
    RHSMatrix = RHSMatrix + RHregen;
elseif flavor == 3
    [RHregen,~] = get_RHS_matrices(energy_nodes,sigma_fname,'/tau_decay_spectrum/tfull',xsecfname);
    RHSMatrix = RHSMatrix + RHregen;
elseif flavor == -1 %add the glashow piece
    sigma_array = sigma_array + get_glashow_total(energy_nodes)/2; %factor of 2 since there is 1 electron per two nucleons
    RHSMatrix = RHSMatrix + get_glashow_partial(energy_nodes)/2;
end



phi_0 = energy_nodes.^(2-g)';
[v,w] = eig(-diag(sigma_array) + RHSMatrix);
ci = (v^-1)*phi_0;
w = diag(w);

end

function  [RHSMatrix, sigma_array] = get_RHS_matrices(energy_nodes,sigma_fname,dxs_fname,xsecfname)
h5flag = 1; %this flag tells the code to read the HDF5 table, rather than a plain text file. It can be turned off here for legacy/comparison purposes
NumNodes = length(energy_nodes);
if h5flag
    %note the transpose ('): this is due to the way the h5 file is indexed
    sigma_array = h5read(xsecfname,sigma_fname)';
    dsigmady = h5read(xsecfname, dxs_fname)';    
    
else
    sigma_array = load(['data' sigma_fname '.dat']);
    dsigmady = load(['data' dxs_fname '.dat']);
end
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
