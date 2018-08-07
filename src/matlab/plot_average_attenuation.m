function [energy_nodes, phisol] = plot_average_attenuation(flavor,g,varargin)
%note: this is for UPGOING neutrinos. you should change ctvec limits for
%alls-ky. 
if nargin > 2 && abs(flavor) ~= 3
    secflag = varargin{1};
else
    secflag = 0;
end
if secflag
[w,v,ci,energy_nodes] = cascade_secs(flavor,g);
else
[w,v,ci,energy_nodes] = cascade(flavor,g);
end
Na = 6.022e23;
phisol = 0;
ctvec = linspace(0,1); %cos-theta-vec
for i = 1:length(ctvec)
tvec(i) = get_optical_depth(acos(ctvec(i)))*Na*1000/100^2;
bigphi = v*(ci.*exp(w*tvec(i)));
phi = bigphi(1:length(energy_nodes));
phisol = phisol + phi.*energy_nodes.^(g-2)';
end    
phisol = phisol/length(tvec);
if nargin > 3
%     this means  plot
semilogx(energy_nodes,phisol);
end
end
function t = get_optical_depth(theta)
    cth = cos(theta);
    xmax = 2*cth;
    Re = 6371;
    if xmax < 0
        t = 0;
    else
    rho = @(x)rhoearth(sqrt(1 + x.^2 - 2.*x*cth));
    t = integral(rho,0,xmax)*1000*Re;
    end
end