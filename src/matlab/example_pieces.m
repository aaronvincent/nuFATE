%NuFate examples for matlab

%% get the spectrum for a single direction
g = 2; %spectral index
flav = -1; %negative sign means antineutrinos, only sign matters
zenith = 130; %choose a zenith.
num = 10; %number of subdivisions
E = logspace(3,10,500);

logemin = 3; %Min log energy (GeV) (do not touch unless you have recomputed cross sections)
logemax = 10; % Max log E
NumNodes = 200;
energy_nodes = logspace(logemin,logemax,NumNodes);
[att_e,att_mu,att_tau,~]=get_att_pieces(flav,g,energy_nodes,zenith,E,num);


figure
loglog(E,att_e,'--r',E,att_mu,'--g',E,att_tau,'--b','linewidth',2);
hold on;
%compare with the original nuFate result
f = [-3 -2 -1]; %array of flavors. 1,2,3 = e, mu, tau; negative sign means antineutrinos
% E = 1e5
for i = 1:length(f)
    flavor = f(i);
    if abs(flavor) ~= 3
        [w,v,ci,energy_nodes] = cascade_secs(flavor,g);
        %[w,v,ci,energy_nodes] = cascade(flavor,g);
    else
        [w,v,ci,energy_nodes] = cascade(flavor,g);
    end
    A = get_att_value(g,w,v,ci,energy_nodes,zenith,E');
    p(i) = loglog(E,A,'linewidth',2);
    hold on    
end
%L = legend([p(3) p(2) p(1)],'$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$');
L = legend('$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$',...
    '$\bar \nu_e$ nuFate','$\bar \nu_\mu$ nuFate','$\bar \nu_\tau$ nuFate');
set(L,'fontsize',20,'interpreter','latex')
legend boxoff
xlabel('$E_\nu$ (GeV)','interpreter','latex','fontsize',20)
ylabel('Attenuation','interpreter','latex','fontsize',20)
title(['zenith = ' num2str(zenith) '$^\circ$'],'fontsize',20,'interpreter','latex')

%compare average fluxes
ctvec = linspace(0,-1,10); %cos-theta-vec
ave_e = zeros(1,length(E));
ave_mu = zeros(1,length(E));
ave_tau = zeros(1,length(E));
for i = 1:length(ctvec)
    zenith=acos(ctvec(i))/pi*180;
    [att_e,att_mu,att_tau,~]=get_att_pieces(flav,g,energy_nodes,zenith,E,num);
    ave_e = ave_e + att_e;
    ave_mu = ave_mu + att_mu;
    ave_tau = ave_tau + att_tau;
end
ave_e = ave_e/length(ctvec);
ave_mu = ave_mu/length(ctvec);
ave_tau = ave_tau/length(ctvec);

figure
semilogx(E,ave_e,'--r',E,ave_mu,'--g',E,ave_tau,'--b','linewidth',2);



