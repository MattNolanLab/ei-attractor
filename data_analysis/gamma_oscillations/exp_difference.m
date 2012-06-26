t = 0:0.01:40;
t_start = 0;

tau_dec = 5;
tau_rise = 1;
tau1 = tau_dec;
tau2 = tau_rise*tau_dec / (tau_rise + tau_dec);
B = 1/((tau2/tau1)^(tau_rise/tau1) - (tau2/tau1)^(tau_rise/tau2));
N = 1.358;

plot(t, B*(exp(-(t-t_start)/tau1) - exp(-(t-t_start)/tau2)));