% Simulate integrate and fire neuron
clear all;
close all;

Ne = 800
Ni = 200
N = Ne + Ni

% All variables are in basic units, i.e. s, volt, etc.
taum_e = 10e-3;
taum_i = 2e-3;

taue = 1e-3;
taui = 5e-3;
Vt = -59e-3;
Vr = -62e-3;
El = -60e-3;
we = 100e-3 / N
wi = 50e-3 / N

e_sparseness = 0.75
i_sparseness = 0.75

Ie = 1.4e-3
Ii = 0;

% Noise normalized per time unit (ms)
noise_sigma = 0.5e-3 / 1e-3;

% Euler settings
dt = 0.1e-3  % 0.1 ms


% Build excitatory neurons state
Ve = Vr + (Vt - Vr) * rand(Ne, 1);
Vi = Vr + (Vt - Vr) * rand(Ni, 1);

ge = zeros(Ni, 1);
gi = zeros(Ne, 1);

spikeMon_e = {};
spikeMon_i = {};

Emon_i = 200;
Imon_i = 100;
Vmon_e = [];
Vmon_i = [];
Vmon_t = [];

% Setup connections Mij: j --> i
% Assuming constant and uniform excitatory/inhibitory weights
Me = rand(Ni, Ne);
Mi = rand(Ne, Ni);
Me = double(Me <= e_sparseness);
Mi = double(Mi <= i_sparseness);

% Simulation
T = 2.5;
times = 0:dt:T;

spikeRecord_e = zeros(Ne, size(times, 2));

display 'Simulation running...'
t = 0;

fired_e = zeros(Ne, 1);
fired_i = zeros(Ni, 1);
f_i = 1;
t_spike = [];

t_i = 1;
for t = times
    Vmon_e = [Vmon_e Ve(Emon_i)];
    Vmon_i = [Vmon_i Vi(Imon_i)];
    Vmon_t = [Vmon_t t];
    

    % Check if neurons fired and add to syn. conductances
    
    gi = gi + Mi*fired_i * wi;
    ge = ge + Me*fired_e * we;
    
    dVe = dt * 1/taum_e * (El - Ve - gi + Ie);
    dVi = dt * 1/taum_i * (El - Vi + ge + Ii);
    
    dge = dt * -1/taue * ge;
    dgi = dt * -1/taui * gi;
    
    Ve = Ve + dVe + dt*noise_sigma*randn(Ne, 1);
    Vi = Vi + dVi + dt*noise_sigma*randn(Ni, 1);
    ge = ge + dge;
    gi = gi + dgi;
    
    fired_e = Ve > Vt;
    fired_i = Vi > Vt;
    
    Ve(fired_e) = Vr;
    Vi(fired_i) = Vr;
    
    if (nnz(fired_e) || nnz(fired_i))
        spikeMon_e{f_i} = find(fired_e);
        spikeMon_i{f_i} = find(fired_i);
        t_spike(f_i) = t;
        f_i = f_i + 1;
    end
    
    spikeRecord_e(:, t_i) = double(fired_e);
   
    %t = t + dt;
    
    t_i = t_i + 1;
end

x_lim = [2 2.5];

% Plot the results
figure('Position', [800 528 1200 800]);
subplot(5, 1, 1);
for it = 1:size(t_spike, 2)
    plot(spikeMon_e{it}*0 + t_spike(it), spikeMon_e{it}, '.');
    hold on;
end
title('Pyramidal neurons');
ylabel('Neuron number');
xlim(x_lim);
ylim([1 Ne]);

subplot(5, 1, 2);
plot(times, getFiringRate(spikeRecord_e, dt, 0.005));
ylabel('Firing rate (Hz)');
xlim(x_lim);


subplot(5, 1, 3);
plot(Vmon_t, Vmon_e*1000);
%title('Pyramidal neuron');
ylabel('Vm (mV)');
xlim(x_lim);


subplot(5, 1, 4);
for it = 1:size(t_spike, 2)
    plot(spikeMon_i{it}*0 + t_spike(it), spikeMon_i{it}, '.');
    hold on;
end

title('Interneurons');
ylabel('Neuron number');
xlim(x_lim);
ylim([1 Ni]);


subplot(5, 1, 5);
plot(Vmon_t, Vmon_i*1000);
%title('Interneuron');
ylabel('Vm (mV)');
xlabel('Time (s)');
xlim(x_lim);

hold off;

%figure();


