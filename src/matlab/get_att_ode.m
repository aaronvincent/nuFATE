function [r_e, r_mu, r_tau, tau] = get_att_ode(flavor,g,zenith,E,varargin)
theta = (180-zenith)/180*pi;
logE = log10(E);
r_e = ones(size(E));
r_mu = r_e;
r_tau = r_e;
tau = zeros(size(E));
if zenith <= 90
    return
elseif zenith > 180
    error('Unphysical zenith angle, must be in between 90 and 180!');
end
[RHSMatrix,energy_nodes,energy_tau] = init_ode(flavor);
NumNodes = length(energy_nodes);
NumTau = length(energy_tau);
phi_nu = energy_nodes.^(2-g)';
%%need a prefactor to avoid negative solutions, default value works for
%%g<=4, increase it if you want sofer spectrum
if nargin > 4
    prefactor = varargin{1};
else %default
    prefactor = 1e70;
end


if nargin > 5
    relerr = varargin{2};
else %default
    relerr = 1e-3;
end

phi_0 = [phi_nu; phi_nu; phi_nu; zeros(NumTau,1)]*prefactor;
phi_out = cascade_ode(theta,phi_0,RHSMatrix,energy_nodes,energy_tau,relerr);
att_e=phi_out(1:NumNodes)./phi_0(1:NumNodes);
att_mu=phi_out(NumNodes+1:2*NumNodes)./phi_0(NumNodes+1:2*NumNodes);
att_tau=phi_out(2*NumNodes+1:3*NumNodes)./phi_0(2*NumNodes+1:3*NumNodes);
tau=phi_out(3*NumNodes+1:end)./energy_tau'.^2;
r_e = interp1(log10(energy_nodes),att_e,logE);
r_mu = interp1(log10(energy_nodes),att_mu,logE);
r_tau = interp1(log10(energy_nodes),att_tau,logE);
tau = interp1(log10(energy_tau), tau, logE)/prefactor;
end