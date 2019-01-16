
clear all
close all


%NuFate examples for matlab

%% get the spectrum for a single direction
g = 2; %spectral index
flav = -1; %negative sign means antineutrinos, only sign matters
zenith = 91; %choose a zenith.
num = 10; %number of subdivisions
E = logspace(3,10,500);
att=zeros(3,length(E));

[att(1,:),att(2,:),att(3,:),tau]=get_att_pieces(flav,g,zenith,E,num);

figure
loglog(E,att(1,:),'--r',E,att(2,:),'--g',E,att(3,:),'--b','linewidth',2);
xlabel('$E_\nu$ (GeV)','interpreter','latex','fontsize',20)
ylabel('Attenuation','interpreter','latex','fontsize',20)
title(['zenith = ' num2str(zenith) '$^\circ$'],'fontsize',20,'interpreter','latex')

L = legend('$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$');
set(L,'fontsize',20,'interpreter','latex','location','southwest');

%tau flux (E^2 factor removed)
figure
tau=tau';
loglog(E,tau','--k');
title('\tau flux');
xlabel('$E_\tau$ (GeV)','interpreter','latex','fontsize',20)
ylabel('\phi_\tau');

%compare with nuFate
figure
loglog(E,att(1,:),'--r',E,att(2,:),'--g',E,att(3,:),'--b','linewidth',2);
%semilogx(E,att_e,'--r',E,att_mu,'--g',E,att_tau,'--b','linewidth',2);
hold on;
f = [-1 -2 -3]; %array of flavors. 1,2,3 = e, mu, tau; negative sign means antineutrinos
color={'r','g','b'};
Areg=zeros(length(E),3);
for i = 1:length(f)
    flavor = f(i);
    if abs(flavor) ~= 3
       [w,v,ci,energy_nodes] = cascade_secs(flavor,g);
        %[w,v,ci,energy_nodes] = cascade(flavor,g);
    else
        [w,v,ci,energy_nodes] = cascade(flavor,g);
    end
    Areg(:,i) = get_att_value(g,w,v,ci,energy_nodes,zenith,E');
    p(i) = loglog(E,Areg(:,i),color{i},'linewidth',2);
    hold on    
end
L = legend('$\bar \nu_e$','$\bar \nu_\mu$','$\bar \nu_\tau$',...
    '$\bar \nu_e$ nuFate','$\bar \nu_\mu$ nuFate','$\bar \nu_\tau$ nuFate','location','southwest');
set(L,'fontsize',20,'interpreter','latex')
legend boxoff
xlabel('$E_\nu$ (GeV)','interpreter','latex','fontsize',20)
ylabel('Attenuation','interpreter','latex','fontsize',20)
title(['zenith = ' num2str(zenith) '$^\circ$'],'fontsize',20,'interpreter','latex')
