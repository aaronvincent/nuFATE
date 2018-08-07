function rho = rhoearth(x)

Re = 6371;
rho = x*0;
for i = 1:length(x)
    r = x(i)*Re;
if r < 1221. 
        p1 =  -0.0002177;
        p2 =  -4.265e-06;
        p3 =   1.309e+04;
elseif r < 3480 
        p1 =  -0.0002409;
        p2 =      0.1416;
        p3 =   1.234e+04;
elseif r < 5721 
        p1 =  -3.764e-05;
        p2 =     -0.1876;
        p3 =        6664;
elseif r < 5961 
        p1 = 0.;
        p2 =      -1.269;
        p3 =   1.131e+04;
elseif r < 6347 
        p1 = 0.;
        p2 = -.725;
        p3 = 7887.;
elseif r < 6356 
        p1 = 0;
        p2 = 0;
        p3 = 2900;
elseif r < 6368 
        p1 = 0;
        p2 = 0;
        p3 = 2600;
else 
        p1 = 0;
        p2 = 0;
        p3 = 1020;
end

rho(i) = p1*r.^2 + p2*r + p3;
end