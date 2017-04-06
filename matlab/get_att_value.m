function A = get_att_value(g,w,v,ci,energy_nodes,zenith,E)

    Na = 6.022e23;
    theta = (180-zenith)/180*pi;
    logE = log10(E);
    t = get_optical_depth(theta)*Na*1000/100^2;
    phisol = v*(ci.*exp(w*t));
    phisol = phisol(1:length(energy_nodes)).*energy_nodes.^(g-2)'; %this step is in case secondaries were included
    A = interp1(log10(energy_nodes),phisol,logE);
end

function t = get_optical_depth(theta);
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