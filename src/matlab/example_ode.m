%ODE examples for matlab
clear all
close all

%% get the spectrum for a single direction
g = 2; %spectral index
flav = -1; %plus(minus) sign for (anti)neutrinos, only sign matters
zenith = 130; %choose a zenith.

E = logspace(3,10,500);

%prefactor to avoid negative solutions, absent to use default value 1e60
prefactor = 1e70;
%relative error of the solutions, absent to use default value 1e-4
relerr = 1e-3;
[nuesol,numusol,nutausol,tausol]=get_att_ode(flav,g,zenith,E',prefactor,relerr);


figure
%plot attenuation
loglog(E,nuesol,'-r',E,numusol,'-g',E,nutausol,'-b','linewidth',2);
L = legend('$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$');
set(L,'fontsize',20,'interpreter','latex')
legend boxoff
xlabel('$E_\nu$ (GeV)','interpreter','latex','fontsize',20)
ylabel('Attenuation','interpreter','latex','fontsize',20)
title(['zenith = ' num2str(zenith) '$^\circ$'],'fontsize',20,'interpreter','latex')

L = legend('$\bar \nu_e$ ODE','$\bar \nu_\mu$ ODE','$\bar \nu_\tau$ ODE');
set(L,'fontsize',20,'interpreter','latex','location','southwest');


figure
%plot tau flux
loglog(E,tausol','-k','linewidth',2);
title(['zenith = ' num2str(zenith) '$^\circ$'],'fontsize',20,'interpreter','latex');
xlabel('$E_\tau$ (GeV)','interpreter','latex','fontsize',20)
ylabel('$\tau$ flux','interpreter','latex','fontsize',20);