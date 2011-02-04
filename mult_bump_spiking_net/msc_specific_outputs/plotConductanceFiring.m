% Plot membrane potential and conductance of a stably firing neuron

conductance_fileName = 'conductance_firing.eps';

%load results/job23601_2010-08-10T16-50-40_output.mat

opt = parseOptions(options);
sheet_size = opt.sheet_size;


close all;
%hold off;

dt = 0.0001; %s
dt_rat = 0.02; % s - population response
delta_t = 0.25; % s - size of window for firing rate estimation
firingPop_t = 5;

SNList_nID = 64;
neuronID = SNList(SNList_nID);

mVolt = 1000;
nS = 1e9;

outputDir = '../../thesis/src/fig/';
    
figure;
plot(SNMonitor_times, SNMonitor_values(SNList_nID, :)*mVolt, 'LineWidth', 2);


% ------------------------------------------------
% Plot the spiking neuron and its conductance
% ------------------------------------------------
fontSize = 16;
figure('Position', [900, 300, 900, 500]);

x_limits = [7 8];
y_limits = [-61 -53];

start_i = x_limits(1)/dt + 1;
end_i = x_limits(2)/dt + 1;
signal = SNMonitor_values(SNList_nID, start_i:end_i);

subplot(2, 3, [1 2], 'FontSize', fontSize, 'XMinorTick', 'on');
hold on;
plot(SNMonitor_times(start_i:end_i) - SNMonitor_times(start_i), signal*mVolt);
ylabel('V_{membrane} (mV)');
%xlim(x_limits);
%ylim(y_limits);
box on;
title 'A';


subplot(2, 3,  [4 5], 'FontSize', fontSize, 'XMinorTick', 'on');
hold on;
plot(SNgMonitor_times(start_i:end_i) - SNMonitor_times(start_i), SNgMonitor_values(SNList_nID, start_i:end_i)*nS);

xlabel('Time (s)');
ylabel('Syn. conductance (nS)');
%xlim(x_limits);
%ylim(y_limits);
box on;



subplot(2, 3, [3 6], 'FontSize', fontSize);
populationResponseFigure
title 'B';

annotation(gcf(),'arrow',[0.3178 0.3167],[0.727 0.356], 'Color', [0 0.8 0]);


set(gcf,'PaperPositionMode','auto');
print('-depsc2', [outputDir conductance_fileName]);
