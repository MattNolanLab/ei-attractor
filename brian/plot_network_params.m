lambda = 13;
beta = 3/lambda^2;
gamma = 1.05*beta;
a = 1;

abs_x = -20:0.1:20;
plot(abs_x, a*exp(-gamma*abs_x.^2) - exp(-beta*abs_x.^2));