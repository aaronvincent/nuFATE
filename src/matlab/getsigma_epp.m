%%define functions
function dsigma=getsigma_epp(rho,y,E,Z)
global mtau me alpha2 re2
ksi=getksi(y,rho);
beta=getbeta(y);
Ytau=getYtau(beta,rho,ksi);
Ye=getYe(beta,rho,ksi);
Ltau=getLtau(Z,rho,ksi,Ytau,E,y);
Le=getLe(Z,rho,ksi,Ye,E,y);
phitau=getphitau(beta,rho,ksi,Ltau);
phie=getphie(beta,rho,ksi,Le);
zeta=getzeta(E,Z);
dsigma=2/3/pi*Z*(Z+zeta)*alpha2*re2*(1-y)/y*(phie+me^2/mtau^2*phitau);

end

function B=getB(Z)
ZB=[1 202.4
2 151.9
3 159.9
4 172.3
5 177.9
6 178.3
7 176.6
8 173.4
9 170.0
10 165.8
11 165.8
12 167.1
13 169.1
14 170.8 
15 172.2
16 173.4
17 174.3
18 174.8
19 175.1
20 175.6
21 176.2
22 176.8
26 175.8
29 173.1
32 173.0
35 173.5
42 175.9
50 177.4
53 178.6
74 177.6
82 178.0
92 179.8];
if ~ismember(Z,ZB(:,1))
    B=182.7;
else
    B=ZB(find(ZB(:,1)==Z),2);
end

end

function ksi=getksi(y,rho)
global mtau me
ksi=(mtau*y/(2*me))^2*(1-rho.^2)/(1-y);
end
function beta=getbeta(y)
beta=y^2/2/(1-y);
end
function Ytau=getYtau(beta,rho,ksi)
Ytau=(4+rho.^2+3*beta*(1+rho.^2))./...
    ((1+rho.^2)*(3/2+2*beta).*log(3+ksi)+1-3/2*rho.^2);
end
function Ye=getYe(beta,rho,ksi)
Ye=(5-rho.^2+4*beta*(1+rho.^2))./...
    (2*(1+3*beta)*log(3+1./ksi)-rho.^2-2*beta*(2-rho.^2));
end
function Ltau=getLtau(Z,rho,ksi,Ytau,E,y)
global mtau me
B=getB(Z);
Ltau=log(2/3*mtau/me*B*Z^(-2/3)./...
    (1+2*me*sqrt(exp(1))*B*Z^(-1/3)*(1+ksi).*(1+Ytau)./(E*y*(1-rho.^2))));
end
function Le=getLe(Z,rho,ksi,Ye,E,y)
global mtau me
B=getB(Z);
Le=log(B*Z^(-1/3)*sqrt((1+ksi).*(1+Ye))./...
    (1+2*me*sqrt(exp(1))*B*Z^(-1/3)*(1+ksi).*(1+Ye)./(E*y*(1-rho.^2))))-...
    1/2*log(1+(3*me*Z^(1/3)/2/mtau)^2*(1+ksi).*(1+Ye));
end
function phitau=getphitau(beta,rho,ksi,Ltau)
phitau=((1+rho.^2)*(1+3/2*beta)-1./ksi*(1+2*beta).*(1-rho.^2)).*log(1+ksi)+...
    ksi.*(1-rho.^2-beta)./(1+ksi)+(1+2*beta)*(1-rho.^2);
phitau=phitau.*Ltau;
end
function phie=getphie(beta,rho,ksi,Le)
phie=((2+rho.^2)*(1+beta)+ksi.*(3+rho.^2)).*log(1+1./ksi)+...
    (1-rho.^2-beta)./(1+ksi)-(3+rho.^2);
phie=phie.*Le;
end
function zeta=getzeta(E,Z)
global gamma1 gamma2 mtau
zeta=(0.073*log(E/mtau/(1+gamma1*Z^(2/3)*E/mtau))-0.026)./...
    (0.058*log(E/mtau/(1+gamma2*Z^(1/3)*E/mtau))-0.014);
end

