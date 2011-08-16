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
opt.Ie = 18.56e-3;
%we = 1/(taue*1000) * 700e-3 / N; % normalize the weight by time constant to inject constant charge
opt.we = 140e-3 / N;


% Inhibitory cell
opt.taum_i = 10e-3;
opt.taui = 5e-3;
opt.El_i = -60e-3;
opt.Vt_i = -50e-3;
opt.Vr_i = -58e-3;
opt.i_sparseness = 0.75;
opt.Ii = 6e-3;
%wi = 1/(taui*1000) * 600e-3 / N;
opt.wi = 15e-3 / N;


% Noise normalized per time unit (ms)
opt.noise_sigma = 0.1e-3 / 1e-3;



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
opt.T = 2;

tic;

[spikeRecord_e, spikeRecord_i, Vmon, times] = simulateEI(opt);
firingRate_e = getFiringRate(spikeRecord_e, dt, rateWindowLen);
spikeCell_e = spikeRecordToSpikeCell(spikeRecord_e, times);
spikeCell_i = spikeRecordToSpikeCell(spikeRecord_i, times);

elapsed_time = toc;
display(['simulation time:' num2str(elapsed_time)]);


x_lim = [0 opt.T];

% Plot the results
figure('Position', [800 528 1200 800]);
subplot(5, 1, 1);
spikeCellRasterPlot(spikeCell_e, '.');
title('Pyramidal neurons');
ylabel('Neuron number');
xlim(x_lim);
ylim([1 opt.Ne]);

subplot(5, 1, 2);
plot(times, firingRate_e);
ylabel('Firing rate (Hz)');
xlim(x_lim);


subplot(5, 1, 3);
plot(Vmon.t, Vmon.e*1000);
%title('Pyramidal neuron');
ylabel('Vm (mV)');
xlim(x_lim);


subplot(5, 1, 4);
spikeCellRasterPlot(spikeCell_i, '.');
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
[Y f NFFT] = fourierTrans(firingRate_e, dt);
Y_abs = 2*abs(Y(1:NFFT/2+1));
plot(f,Y_abs);
xlim([0 200]);

