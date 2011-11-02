Ji = -10;
Je = 1;
tau_i = 0.01;
tau_e = 0.01;
D = 0.01;

om = 0:0.01:1000;

%me_iw = (1 + exp(-i*2*om*D)*Je*Ji ./ (1 + i*om*tau_i)) ./ (1 + i*om*tau_e);
me_iw = 1 ./ (1 + i*om*tau_e - Ji*exp(-i*om*D));

plot(om, abs(me_iw));
