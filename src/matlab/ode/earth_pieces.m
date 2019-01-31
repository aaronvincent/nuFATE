function [x,rhobar,radii] = earth_pieces(zenith,varargin)%radii to be tested
%%input: zenith angle-theta in degrees, number of subdivisions-num 
%%(optional, default to 10)
%%output: list of lengths of the path in the sublayers-x (in cm), 
%%the corresponding average densities in the sublayers-rhobar(in g/cm^3)
%%the smaller radius of the sublayer-radii (in km)
REarth=6371.;

%for test
% zenith=92;
% num=10;
% nargin=1;

theta=zenith/180.*pi;
if nargin>=2
    num=varargin{1};
else
    num=10;
end
%for downgoing neutrinos, path in the earth is 0
if zenith<=90.
    x=0;
    rhobar=2.7;
    radii=0;
else    

rmin=REarth*sin(theta);
%STW105 model, can change according to other earth models
nodes_earth=[0., 1221., 3480., 5721., 5961., 6347., 6356.,6368., REarth];
%%determine how many layers the neutrinos go through
layer=length(nodes_earth);
for i=length(nodes_earth):-1:1
    if rmin>=nodes_earth(i)
        layer=i;
        break;
    end
end
nodes=nodes_earth(layer:end);
%%assign weight to each layer and divide them according to the user
%%specified number
diff=zeros(length(nodes)-1);
w=zeros(length(diff),1);
for i=1:length(nodes)-1
    diff(i)=getrho(nodes(i))/getrho(nodes(i+1)-0.1)-1;
    w(i)=diff(i)*(nodes(i+1)-nodes(i));
end
if max(w)>0 %not constant density
    wt=w/sum(w);
    [m,ind]=max(wt);
    div=floor(num*wt);
    div(ind)=num-sum(div)+div(ind);
    nodes_div=zeros(length(nodes)+num,1);
    j=1;
    for i=1:length(nodes)-1
        nodes_div(j)=nodes(i);
        if div(i)==0
            j=j+1;
        else
            l=j;
            j=j+1;
            delta=(nodes(i+1)-nodes(i))/double(div(i)+1);
            for k=1:div(i)
                nodes_div(l+k)=nodes_div(l)+delta*double(k);
                j=j+1;
            end
        end
    end
    nodes_div(j)=nodes(end);
else %constant density, don't divide!
    nodes_div=nodes';
end
%%computer the average density in each sublayer
rhobar=zeros(length(nodes_div)-1,1);
for i=1:length(nodes_div)-1
    p=getps(nodes_div(i));
    rhobar(i)=(3./5.*p(1)*(nodes_div(i+1)^5-nodes_div(i)^5)+...
        3./4.*p(2)*(nodes_div(i+1)^4-nodes_div(i)^4))/...
        (nodes_div(i+1)^3-nodes_div(i)^3)+p(3);
end
rhobar=rhobar*1e-3; %from kg*m^-3 to g*cm^-3
%%compute the array of path length x in different sublayers
layer=length(nodes_div);
for i=length(nodes_div):-1:1
    if rmin>=nodes_div(i)
        layer=i+1;
        break;
    end
end
x=zeros(length(nodes_div)-layer+1,1);
if layer==length(nodes_div)
    x=2*REarth*cos(pi-theta);
else
    x(1)=sqrt(nodes_div(layer)^2-REarth^2*sin(theta)^2);
    j=2;
    for i=layer:length(nodes_div)-1
        x(j)=sqrt(nodes_div(i+1)^2-REarth^2*sin(theta)^2)-...
            sqrt(nodes_div(i)^2-REarth^2*sin(theta)^2);
        j=j+1;
    end
end
rhobar=rhobar(layer-1:end);
if length(x)<2
    radii=nodes_div(layer-1);
else
    radii=nodes_div(layer-1:end-1);
    radii=[flipud(radii);radii(2:end)];
    x=[flipud(x(2:end));2*x(1);x(2:end)]; %the order matters!
    rhobar=[flipud(rhobar);rhobar(2:end)];
end
x=x*1e5; %from km to cm
%done!

end
end


function rho = getrho(r)
    if r < 1221.
        p1 = -0.0002177;
        p2 = -4.265e-06;
        p3 = 1.309e+04;
    elseif r < 3480
        p1 = -0.0002409;
        p2 = 0.1416;
        p3 = 1.234e+04;
    elseif r < 5721
        p1 = -3.764e-05;
        p2 = -0.1876;
        p3 = 6664;
    elseif r < 5961
        p1 = 0.;
        p2 = -1.269;
        p3 = 1.131e+04;
    elseif r < 6347
        p1 = 0.;
        p2 = -0.725;
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

    rho = p1 * r^2 + p2 * r + p3;
end


function p = getps(r)
    if r < 1221.
        p1 = -0.0002177;
        p2 = -4.265e-06;
        p3 = 1.309e+04;
    elseif r < 3480
        p1 = -0.0002409;
        p2 = 0.1416;
        p3 = 1.234e+04;
    elseif r < 5721
        p1 = -3.764e-05;
        p2 = -0.1876;
        p3 = 6664;
    elseif r < 5961
        p1 = 0.;
        p2 = -1.269;
        p3 = 1.131e+04;
    elseif r < 6347
        p1 = 0.;
        p2 = -0.725;
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

    p=[p1, p2, p3];
end