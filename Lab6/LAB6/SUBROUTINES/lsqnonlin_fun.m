function out = lsqnonlin_fun(coeff)
load lsqnonlin.mat x y 
p = coeff(1);
ap = coeff(2);
gamma1 = coeff(3);
gamma2 = coeff(4);
P = 3.522;
i = coeff(5);
[phi,F] = quadLimbDark(p,ap,P,i,gamma1,gamma2,100,100);

y2 = interp1(phi,F,x);

out = y-y2;

end
