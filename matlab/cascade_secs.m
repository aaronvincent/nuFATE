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
if nargin >= 3 %you shouldn't be specifying only one of these:
    logemin = varargin{3};
    logemax = varargin{4};
    NumNodes = varargin{5};
end
if flavor==-1
    sigma_fname = '../data/nuebarxs.dat';
elseif flavor == -2
    sigma_fname = '../data/numubarxs.dat';

elseif flavor == 1
    sigma_fname = '../data/nuexs.dat';
elseif flavor == 2
    sigma_fname = '../data/numuxs.dat';
else
    error('flavor for secs must be +/- 1 or 2 ')
end
if flavor > 0
        dxs_fname = '../data/dxsnu.dat';
        sig3fname = '../data/nutauxs.dat';
        secname = '../data/secfull.dat';
        regenname = '../data/tfull.dat';
else
        dxs_fname = '../data/dxsnubar.dat';
        sig3fname = '../data/nutaubarxs.dat';
        secname = '../data/secbarfull.dat';
        regenname = '../data/tbarfull.dat';
end

energy_nodes = logspace(logemin,logemax,NumNodes);
[RHSMatrix] = get_RHS_matrices(energy_nodes,sigma_fname,sig3fname,dxs_fname,secname,regenname);


phi_0 = energy_nodes.^(2-g)';
[v,w] = eig(RHSMatrix);
ci = (v^-1)*[phi_0; phi_0];
w = diag(w); 

end

function  [RHSMatrix, sigma_array] = get_RHS_matrices(energy_nodes,sigma_fname,sig3fname,dxs_fname,secname,regenname)
    NumNodes = length(energy_nodes);
    sigma_array1 = load(sigma_fname);
    sigma_array2 = load(sig3fname);
    dsigmady = load(dxs_fname);
    emuregen = load(secname);
    tauregen = load(regenname);
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
