function [phi,F]=quadLimbDark(p,ap,P,i,gamma1,gamma2,n,percentOfOrbit)

% Based on Mandel & Agol's Model for Quadratic Limb Darkening
% Mandel K. & Agol E. 2002, ApJ, 580, L171; please cite this
% paper if you make use of this in your research.  Also, a
% thanks to Gil Nachmani who wrote this routine would be appreciated.

% [phi,F] are the observed phase and relative flux
% p is the planet's radius in units of the star's radius (Rs)
% ap is the planet's orbital radius in units of Rs, assuming zero eccentricity
% P is the planet's period in days
% i is the inclination of the orbit in degrees
% gamma1, gamma2 are the Quadratic Limb Darkening coefficients:
% I(r) = 1-gamma1*(1-mu)-gamma2*(1-mu)^2, where: mu=cos(th)=sqrt(1-r^2);
% E.g. gamma1=0.296, gamma2=0.34 for HD209458
% n is the number of phase points in the resulting lightcurve
% percentOfOrbit is the percentage of orbital phase to be used for the flux
% calculations, e.g. for full orbit use 100(%).

% Gil Nachmani, April 2011

if nargin<8; percentOfOrbit=150*(p+1)/ap/2/pi; end
t=linspace(-P*percentOfOrbit/100,P*percentOfOrbit/100,n); phi=t/P;
Z=ap*(sin(2*pi/P*t).^2+(cos(pi/180*i)*cos(2*pi/P*t)).^2).^(1/2);

c1=0; c2=gamma1+2*gamma2; c3=0; c4=-gamma2; % I(r) = 1-sum(cn*(1-mu^(n/2)))
c0=1-c1-c2-c3-c4; ohmega = c0/(0+4)+c1/(1+4)+c2/(2+4)+c3/(3+4)+c4/(4+4);

j=0; 
for z=Z; j=j+1; 
    a=(z-p)^2; b=(z+p)^2; k=sqrt((1-a)/4/z/p); q=p^2-z^2; 
    k1=acos((1-p^2+z^2)/2/z); k0=acos((p^2+z^2-1)/2/p/z);
       
    %Evaluating lam_e
    if 1+p<z || abs(phi(j))>(p+1)/ap/2/pi;
        lam_e = 0;
    elseif abs(1-p)<z && z<=1+p;
        lam_e = 1/pi*(p^2*k0+k1-1/2*sqrt(4*z^2-(1+z^2-p^2)^2));
    elseif z<=1-p && z>p-1;
        lam_e = p^2;
    elseif z<=p-1;
        lam_e = 1;
    end
    
    %Evaluating lam_d and eta_d
    if z>=1+p || p==0 || abs(phi(j))>(p+1)/ap/2/pi; %1
        lam_d = 0; eta_d=0; 
    elseif z>=1/2+abs(p-1/2) && z<1+p; %2
        lam_d = lam1(p,z,a,b,k,q); eta_d = eta1(p,z,a,b,k1,k0);
    elseif p<0.5 && z>p && z<1-p; %3
        lam_d = lam2(p,z,a,b,k,q); eta_d = eta2(p,z);
    elseif p<0.5 && z==1-p; %4
        lam_d = lam5(p); eta_d = eta2(p,z);
    elseif p<0.5 && z==p; %5
        lam_d = lam4(p); eta_d = eta2(p,z);
    elseif p==0.5 && z==0.5; %6
        lam_d = 1/3-4/pi/9; eta_d = 3/32; 
    elseif p>0.5 && z==p; %7
        lam_d = lam3(p); eta_d = eta1(p,z,a,b,k1,k0); 
    elseif p>0.5 && z>=abs(1-p) && z<p; %8
        lam_d = lam1(p,z,a,b,k,q); eta_d = eta1(p,z,a,b,k1,k0); 
    elseif p<1 && z>0 && z<=1/2-abs(p-1/2); %9
        lam_d = lam2(p,z,a,b,k,q); eta_d = eta2(p,z); 
    elseif p<1 && z==0; %10
        lam_d = lam6(p); eta_d=eta2(p,z); 
    elseif p>1 && z<=p-1; %11
        lam_d = 0; eta_d=1/2; 
    end
F(j) = 1 - 1/(4*ohmega)*( (1-c2)*lam_e + c2*(lam_d+2/3*heaviside(p-z)) - c4*eta_d );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lam=lam1(p,z,a,b,k,q)
lam = 1/9/pi/sqrt(p*z) * ( ((1-b)*(2*b+a-3)-3*q*(b-2))*K(k) + 4*p*z*(z^2+7*p^2-4)*E(k)-3*q/a*PI((a-1)/a,k) );
function lam=lam2(p,z,a,b,k,q)
lam = 2/9/pi/sqrt(1-a) * ( (1-5*z^2+p^2+q^2)*K(1/k) + (1-a)*(z^2+7*p^2-4)*E(1/k)-3*q/a*PI((a-b)/a,1/k) );
function lam=lam3(p)
lam = 1/3 + 16*p/9/pi*(2*p^2-1)*E(1/2/p) - (1-4*p^2)*(3-8*p^2)/9/pi/p*K(1/2/p);
function lam=lam4(p)
lam = 1/3 + 2/9/pi*(4*(2*p^2-1)*E(2*p)+(1-4*p^2)*K(2*p));
function lam=lam5(p)
lam = 2/3/pi*acos(1-2*p) - 4/9/pi*(3+2*p-8*p^2)*sqrt(p*(1-p))-2/3*heaviside(p-1/2);
function lam=lam6(p)
lam = -2/3*(1-p^2)^(3/2);
function eta=eta1(p,z,a,b,k1,k0)
eta = 1/2/pi*(k1+2*eta2(p,z)*k0-1/4*(1+5*p^2+z^2)*sqrt((1-a)*(b-1)));
function eta=eta2(p,z)
eta = p^2/2*(p^2+2*z^2);
function f=E(k)
f = ellipticE(k.^2);
%f = mfun('EllipticE',k);
function f=PI(n,k) 
%f = mfun('EllipticPi',n,k);
f = ellipticPi(n,k.^2);
function f=K(k)
f = ellipticK(k.^2);
%f = mfun('EllipticK',k);


    
