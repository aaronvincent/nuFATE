function [r_e, r_mu, r_tau, tau] = get_att_pieces(flavor,g,zenith,E,num)
[RHSMatrix,energy_nodes] = init_pieces(flavor);
NumNodes = length(energy_nodes);
phi_nu = energy_nodes.^(2-g)';
phi_0 = [phi_nu; phi_nu; phi_nu; zeros(NumNodes,1)];
[x,rho,radii]=earth_pieces(zenith,num);
phi_in=phi_0;
for i=1:length(x)
%for i=1:1    
    phi_out=cascade_pieces(x(i),rho(i),radii(i),phi_in,RHSMatrix,energy_nodes);
    phi_in=phi_out;
end
att_e=phi_out(1:NumNodes)./phi_0(1:NumNodes);
att_mu=phi_out(NumNodes+1:2*NumNodes)./phi_0(NumNodes+1:2*NumNodes);
att_tau=phi_out(2*NumNodes+1:3*NumNodes)./phi_0(2*NumNodes+1:3*NumNodes);
tau=phi_out(3*NumNodes+1:4*NumNodes)./energy_nodes'.^2;
logE = log10(E);
r_e = interp1(log10(energy_nodes),att_e,logE);
r_mu = interp1(log10(energy_nodes),att_mu,logE);
r_tau = interp1(log10(energy_nodes),att_tau,logE);
tau = interp1(log10(energy_nodes), tau, logE);
end