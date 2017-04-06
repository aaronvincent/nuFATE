function [w,v,ci,energy_nodes] = cascade(varargin) %this is get_eigs
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
elseif flavor == -3
    sigma_fname = '../data/nutaubarxs.dat';
elseif flavor == 1
    sigma_fname = '../data/nuexs.dat';
elseif flavor == 2
    sigma_fname = '../data/numuxs.dat';
elseif flavor == 3
    sigma_fname = '../data/nutauxs.dat';
end
if flavor > 0
    dxs_fname = '../data/dxsnu.dat';
else
    dxs_fname = '../data/dxsnubar.dat';
end

energy_nodes = logspace(logemin,logemax,NumNodes);
[RHSMatrix, sigma_array] = get_RHS_matrices(energy_nodes,sigma_fname,dxs_fname);
if flavor ==-3
        [RHregen,~] = get_RHS_matrices(energy_nodes,sigma_fname,'../data/tbarfull.dat');
        RHSMatrix = RHSMatrix + RHregen;
elseif flavor == 3
        [RHregen,~] = get_RHS_matrices(energy_nodes,sigma_fname,'../data/tfull.dat');
        RHSMatrix = RHSMatrix + RHregen;
end

phi_0 = energy_nodes.^(2-g)';
[v,w] = eig(-diag(sigma_array) + RHSMatrix);
ci = (v^-1)*phi_0;
w = diag(w); 

end

function  [RHSMatrix, sigma_array] = get_RHS_matrices(energy_nodes,sigma_fname,dxs_fname)
    NumNodes = length(energy_nodes);
    sigma_array = load(sigma_fname);
    dsigmady = load(dxs_fname);
    DeltaE = diff(log(energy_nodes));
    RHSMatrix = zeros(NumNodes);
    for i = 1:NumNodes
        for j = i+1:NumNodes
            RHSMatrix(i,j) = DeltaE(j-1)*dsigmady(j,i)*energy_nodes(j).^-1*energy_nodes(i).^2;
        end
    end
end
