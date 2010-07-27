lambda = 20;
beta = 3/lambda^2;
gamma = 1.05*beta;
a = 1;
connMult = 4;

abs_x = -40:0.1:40;
fun = a*exp(-gamma*abs_x.^2) - exp(-beta*abs_x.^2);
plot(abs_x, connMult * fun);