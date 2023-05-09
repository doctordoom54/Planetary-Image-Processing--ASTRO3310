clc
clear
data = importdata('plot.tbl');
p = 3.522499;
x_pf = mod(data.data(:,2),p);
y_f = data.data(:,3);
[minval,minind] = min(y_f);
x_pf = x_pf - x_pf(minind);
x_pf = x_pf./p;
y_f = y_f +1;
plot(x_pf,y_f,'c.')
r_td = sqrt(1-min(y_f));
coeff_innit = [0.1,7,0.3,0.3,0.83];
coeff_low = [0.01,5,0.2,0.07,80];
coeff_high = [0.15,10,0.6,0.36,90]; 