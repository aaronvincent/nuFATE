function phi_1 = cascade_ode(theta,phi_0,RHSMatrix,energynodes,energytau,varargin) 
%%this is get_eigs
%%regeneration bits are included here and three nus+tau are solved at the 
%%same time. For treatment of on-the-spot decay of tau or without
%%regenerations, refer to the old code
%%input: length of path in a sublayer-x (in cm), average density the 
%%sublayer-rhobar (in g*cm^-3), input flux-phi_0, radius of layer-radius, 
%%RHS matrix elements-RHSMatrix, energy nodes-energy_nodes
%%output: output flux-phi_1, energy nodes-energy_nodes
if nargin > 5
    relerr = varargin{1};
else %default
    relerr = 1e-3;
end
REarth=6371.;%radius of the Earth
[m, n] = size(phi_0);
if m<n
    phi_0 = phi_0'; %turn into a column vector
end
energy_nodes = energynodes;
energy_tau = energytau;
xmax=2*REarth*cos(theta)*1e5;%in cm
odefun=@(x,phi) getdphidx(theta,x,phi,RHSMatrix,energy_nodes,energy_tau)*phi;
options=odeset('RelTol',relerr,'NonNegative',1,...
    'Jacobian',@(x,phi) getdphidx(theta,x,phi,RHSMatrix,energy_nodes,energy_tau));
[~,phisol]=ode15s(odefun,[0,xmax],phi_0,options);
phi_1=phisol(end,:)';
end

function dphidx = getdphidx(theta,x,phi,RHSMatrix,energy_nodes,energy_tau)
NA = 6.02e+23; %Avogadro's constant
c = 3.e+10; %speed of light in cm/s
mtau = 1.777; %tau mass in GeV
ttau = 2.906e-13; %tau lifetime in rest frame in s
%define crust, mantle and core, it depends on the modeling of the Earth
REarth=6371.*1e5;%radius of the Earth in cm
nodes_earth=[0., 1221., 3480., 5721., 5961., 6347., 6356.,6368., REarth]*1e5;
Rcrust=nodes_earth(6);
Rmantle=nodes_earth(3);
cth = cos(theta);
xtoR=x/REarth;
r=sqrt(1 + xtoR.^2 - 2.*xtoR*cth);
radius=REarth*r;
rho=rhoearth(r)*1e-3;%in g/cm^3

NumNodes = length(energy_nodes);
NumTau = length(energy_tau);
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
RHSMatrix_44=RHSMatrix_44-diag(mtau./(c*ttau*rho*NA*energy_tau)); %tau decay
RHSMatrix_43=RHSMatrix{1,7};
RHSMatrix_34=RHSMatrix{1,8}*mtau/(c*ttau*rho*NA);
RHSMatrix_14=RHSMatrix{1,9}*mtau/(c*ttau*rho*NA);
RHSMatrix_24=RHSMatrix_14;
%%construct the huge RHS Matrix
z0=zeros(NumNodes);
z0tau=zeros(NumTau,NumNodes);
M=[RHSMatrix_11 z0 z0 RHSMatrix_14;
   z0 RHSMatrix_22 z0 RHSMatrix_24;
   z0 z0 RHSMatrix_33 RHSMatrix_34;
   z0tau z0tau RHSMatrix_43 RHSMatrix_44];
dphidx=M*rho*NA;
end