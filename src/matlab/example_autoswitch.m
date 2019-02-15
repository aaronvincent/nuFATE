%example to switch between calssical ODE and nuFATE

%% get the spectrum for a single direction
g = 2; %spectral index
flav = -1; %plus(minus) sign for (anti)neutrinos, only sign matters for ODE
zenith = 120; %choose a zenith.

E = logspace(3,10,500);
%specify an angle to switch between classical nuFATE and ODE solver, 
%below that angle use ODE solver while above that angle use nuFATE instead to speed up
switchangle = 110;

if zenith <= switchangle %use ODE
    %prefactor to avoid negative solutions, absent to use default value 1e70
    prefactor = 1e70;
    %relative error of the solutions, absent to use default value 1e-3
    relerr = 1e-3;
    [nuesol,numusol,nutausol,tausol]=get_att_ode(flav,g,zenith,E',prefactor,relerr);
else %use classical nuFATE
    f = [-1 -2 -3]; %array of flavors. 1,2,3 = e, mu, tau; negative sign means antineutrinos
    for i = 1:length(f)
        flavor = f(i);
        if abs(flavor) ~= 3
            [w,v,ci,energy_nodes] = cascade_secs(flavor,g);
        else
            [w,v,ci,energy_nodes] = cascade(flavor,g);
        end
        if i == 1
            nuesol = get_att_value(g,w,v,ci,energy_nodes,zenith,E');
        elseif i == 2
            numusol = get_att_value(g,w,v,ci,energy_nodes,zenith,E');
        else
            nutausol = get_att_value(g,w,v,ci,energy_nodes,zenith,E');
        end
    end  
end


figure
%plot attenuation
loglog(E,nuesol,'-r',E,numusol,'-g',E,nutausol,'-b','linewidth',2);
L = legend('$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$');
set(L,'fontsize',20,'interpreter','latex')
legend boxoff
xlabel('$E_\nu$ (GeV)','interpreter','latex','fontsize',20)
ylabel('Attenuation','interpreter','latex','fontsize',20)
title(['zenith = ' num2str(zenith) '$^\circ$'],'fontsize',20,'interpreter','latex')

L = legend('$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$');
set(L,'fontsize',20,'interpreter','latex','location','southwest');
