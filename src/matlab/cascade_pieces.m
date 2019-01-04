function [phi_1,energy_nodes] = cascade_pieces(x,rhobar,radius,phi_0,RHSMatrix,energynodes) 
%%this is get_eigs
%%regeneration bits are included here and three nus+tau are solved at the 
%%same time. For treatment of on-the-spot decay of tau or without
%%regenerations, refer to the old code
%%input: length of path in a sublayer-x (in cm), average density the 
%%sublayer-rhobar (in g*cm^-3), input flux-phi_0, radius of layer-radius, 
%%RHS matrix elements-RHSMatrix, energy nodes-energy_nodes
%%output: output flux-phi_1, energy nodes-energy_nodes


NA = 6.02e+23; %Avogadro's constant
c = 3.e+10; %speed of light in cm/s
mtau = 1.777; %tau mass in GeV
ttau = 2.906e-13; %tau lifetime in rest frame in s
%ttau = 2.906e-14;
%define crust, mantle and core, it depends on the modeling of the Earth
REarth=6371.;
nodes_earth=[0., 1221., 3480., 5721., 5961., 6347., 6356.,6368., REarth];
Rcrust=nodes_earth(6);
Rmantle=nodes_earth(3);

% defaults:
NumNodes = 200;

[m, n] = size(phi_0);
if m<n
    phi_0 = phi_0'; %turn into a column vector
end

energy_nodes=energynodes;
%restore matrix elements
RHSMatrix_11=RHSMatrix{1,1};
RHSMatrix_22=RHSMatrix{1,2};
RHSMatrix_33=RHSMatrix{1,3};
if radius>=Rcrust%crust
    RHSMatrix_44=RHSMatrix{1,4};
elseif radius>=Rmantle%mantle
    RHSMatrix_44=RHSMatrix{1,5};
else%core
    RHSMatrix_44=RHSMatrix{1,6};
end
%RHSMatrix_44=0;
RHSMatrix_44=RHSMatrix_44-diag(mtau./(NA*c*ttau*rhobar*energy_nodes)); %tau decay
RHSMatrix_43=RHSMatrix{1,7};
RHSMatrix_34=RHSMatrix{1,8}*mtau/(NA*c*ttau*rhobar);
RHSMatrix_14=RHSMatrix{1,9}*mtau/(NA*c*ttau*rhobar);
RHSMatrix_24=RHSMatrix_14;


%%construct the huge RHS Matrix
z0=zeros(NumNodes);
M=[RHSMatrix_11 z0 z0 RHSMatrix_14;
   z0 RHSMatrix_22 z0 RHSMatrix_24;
   z0 z0 RHSMatrix_33 RHSMatrix_34;
   z0 z0 RHSMatrix_43 RHSMatrix_44];
% M=[RHSMatrix_11 z0 z0 z0;
%    z0 RHSMatrix_22 z0 z0;
%    z0 z0 RHSMatrix_33 RHSMatrix_34;
%    z0 z0 RHSMatrix_43 RHSMatrix_44];

%phi_0 = energy_nodes.^(2-g)';
[v,w] = eig(M);
ci = (v^-1)*phi_0;
w = diag(w);
t = rhobar*x*NA; %optical depth
phi_1 = v*(ci.*exp(w*t));

end
