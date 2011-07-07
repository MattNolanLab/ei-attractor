t = 0:0.01:10;

om = 2*pi;
phi1 = 0;
phi2 = phi1 + 2*pi/10;

x1 = cos(om.*t + phi1);
x2 = cos(om.*t + phi2);

plot(t, x1+x2);