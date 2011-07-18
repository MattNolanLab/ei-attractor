% Simulate integrate and fire neuron
clear all;
close all;

% All variables are in basic units, i.e. s, volt, etc.
opt.Ne = 800;
opt.Ni = 200;
N = opt.Ne + opt.Ni;


% Excitatory cells
opt.taum_e = 9.3e-3;
opt.taue = 2e-3;
opt.El_e = -68.5e-3;
opt.Vt_e = -50.0e-3;
opt.Vr_e = -60.0e-3;
opt.e_sparseness = 0.75;
opt.Ie = 30.0e-3;
%we = 1/(taue*1000) * 700e-3 / N; % normalize the weight by time constant to inject constant charge
opt.we = 120e-3 / N;


% Inhibitory cell
opt.taum_i = 10e-3;
opt.taui = 5e-3;
opt.El_i = -60e-3;
opt.Vt_i = -50e-3;
opt.Vr_i = -58e-3;
opt.i_sparseness = 0.75;
opt.Ii = 5e-3;
%wi = 1/(taui*1000) * 600e-3 / N;
opt.wi = 0e-3 / N;



% Noise normalized per time unit (ms)
opt.noise_sigma = 0.05e-3 / 1e-3;



% Euler settings
opt.dt = 0.1e-3  % 0.1 ms
dt = opt.dt;


% Firing rate sliding window length
opt.rateWindowLen = 0.005; %ms
rateWindowLen = opt.rateWindowLen;

% Vm monitor, neuron index
opt.Emon_i = 100;
opt.Imon_i = 100;

% simulation time
opt.T = 0.5;


[spikeRecord_e, spikeRecord_i, spikeTimes, Vmon, times] = simulateEI(opt);
spikeMon_e = spikeTimes.e;
spikeMon_i = spikeTimes.i;
t_spike = spikeTimes.t_spike;


firingRate = getFiringRate(spikeRecord_e, dt, rateWindowLen);


x_lim = [0 opt.T];

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
ylim([1 opt.Ne]);

subplot(5, 1, 2);
plot(times, firingRate);
ylabel('Firing rate (Hz)');
xlim(x_lim);


subplot(5, 1, 3);
plot(Vmon.t, Vmon.e*1000);
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
ylim([1 opt.Ni]);


subplot(5, 1, 5);
plot(Vmon.t, Vmon.i*1000);
%title('Interneuron');
ylabel('Vm (mV)');
xlabel('Time (s)');
xlim(x_lim);

hold off;


%
% Plot firing rate fft
%
figure();
[Y f NFFT] = fourierTrans(firingRate, dt);
Y_abs = 2*abs(Y(1:NFFT/2+1));
plot(f,Y_abs);
xlim([0 200]);

