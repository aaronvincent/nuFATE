%%define functions
function dsigma=getsigma_pn(Q2,y,E,Z,A)
global R cp1 cp2 cp3 ap1 ap2 ap3 bp1 bp2 bp3 cr1 cr2 cr3 ar1 ar2 ar3
global br1 br2 br3 m02 mp2 mr2 alpha2 M mtau
x=Q2/(2*M*E*y);
a=geta(A,x);
P=getP(x);
t=gett(Q2);
ar=getf(t,ar1,ar2,ar3);
br=getf(t,br1,br2,br3);
cr=getf(t,cr1,cr2,cr3);
bp=getf(t,bp1,bp2,bp3);
cp=getg(t,cp1,cp2,cp3);
ap=getg(t,ap1,ap2,ap3);
W2=getW2(Q2,E,y);
xr=getxi(Q2,W2,mr2);
xp=getxi(Q2,W2,mp2);
f2R=getF2i(x,cr,xr,ar,br);
f2P=getF2i(x,cp,xp,ap,bp);
f2p=Q2./(Q2+m02).*(f2P+f2R);
F2=getF2(a,A,Z,P,f2p);
dsigma=4*pi*alpha2./Q2.^2.*F2/y.*(1-y-M*x*y/(2*E)+...
    (1-2*mtau^2./Q2)*y^2.*(1+4*M^2*x.^2./Q2)/(2*(1+R)));


end

function t=gett(Q2)
global Q02 L2
t=log(log((Q2+Q02)/L2)/log(Q02/L2));
end
function f=getf(t,f1,f2,f3) 
f=f1+f2*t.^f3;
end
function g=getg(t,g1,g2,g3) 
g=g1+(g1-g2)*(1./(1+t.^g3)-1);
end
function xi=getxi(Q2,W2,m2)
global M
xi=(Q2+m2)./(Q2+m2+W2-M^2);
end
function W2=getW2(Q2,E,y)
global M
W2=M^2+2*M*E*y-Q2;
end
function F2i=getF2i(x,ci,xi,ai,bi) 
F2i=ci.*xi.^ai.*(1-x).^bi;
end
function P=getP(x) 
P=1-1.85*x+2.45*x.^2-2.35*x.^3+x.^4;
end
function F2=getF2(a,A,Z,P,F2p)
F2=a.*(Z+(A-Z)*P).*F2p;
end
function a=geta(A,x)
if x<0.0014
    a=A^-0.1;
elseif x<0.04
    a=A.^(0.069*log10(x)+0.097);
else
    a=1;
end
end
