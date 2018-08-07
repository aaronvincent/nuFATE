%NuFate examples for matlab

%% get the spectrum for a single direction
g = 2; %spectral index
f = [-3 -2 -1]; %array of flavors. 1,2,3 = e, mu, tau; negative sign means antineutrinos
zenith = 130; %choose a zenith.


E = logspace(4,9,500);
% E = 1e5
for i = 1:length(f)
    flavor = f(i);
    if abs(flavor) ~= 3
        [w,v,ci,energy_nodes] = cascade_secs(flavor,g);
    else
        [w,v,ci,energy_nodes] = cascade(flavor,g);
    end
    A = get_att_value(g,w,v,ci,energy_nodes,zenith,E');
    p(i) = loglog(E,A,'linewidth',2);
    hold on    
end
L = legend([p(3) p(2) p(1)],'$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$');
set(L,'fontsize',20,'interpreter','latex')
legend boxoff
xlabel('$E_\nu$ (GeV)','interpreter','latex','fontsize',20)
ylabel('Attenuation','interpreter','latex','fontsize',20)
title(['zenith = ' num2str(zenith) '$^\circ$'],'fontsize',20,'interpreter','latex')

%% Angle-averaged spectrum 
includesecondaries = 1; %0 to ignore secondary mu and e neutrinos from tau decay
figure
for i = 1:length(f)
    flavor = f(i)
    [energy_nodes, phisol] = plot_average_attenuation(flavor,g,includesecondaries);
    p(i) = semilogx(energy_nodes,phisol,'linewidth',2);
    hold on
end
L = legend([p(3) p(2) p(1)],'$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$');
set(L,'fontsize',20,'interpreter','latex')
legend boxoff
xlabel('$E_\nu$ (GeV)','interpreter','latex','fontsize',20)
ylabel('angle-averaged attenuation','interpreter','latex','fontsize',20)
title('Upgoing Antineutrinos','fontsize',20,'interpreter','latex')