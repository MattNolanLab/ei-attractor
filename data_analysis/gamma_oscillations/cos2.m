t = 0:0.001:20;

om1 = 2*pi;
om2 = om1 + 0.1;
phi1 = pi/10;
phi2 = 0;

x1 = cos(om1.*t + phi1);
x2 = cos(om2.*t + phi2);

plot(t, x1+x2);