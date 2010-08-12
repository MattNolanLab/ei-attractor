% Plot membrane potential and conductance of a stably firing neuron

conductance_freq_fileName = 'conductance_freq_firing.eps';

%load results/job23502_2010-08-07T19-37-40_output.mat

opt = parseOptions(options);
sheet_size = opt.sheet_size;


close all;
%hold off;

dt = 0.0001; %s
dt_rat = 0.02; % s - population response
delta_t = 0.25; % s - size of window for firing rate estimation
firingPop_t = 5;

SNList_nID = 11;
neuronID = SNList(SNList_nID);

mVolt = 1000;
nS = 1e9;

outputDir = '../../thesis/src/fig/';
    
figure;
plot(SNMonitor_times, SNMonitor_values(SNList_nID, :)*mVolt, 'LineWidth', 2);


% ------------------------------------------------
% Plot the spiking neuron
% ------------------------------------------------
%figure;


