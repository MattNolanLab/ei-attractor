% Illustrate conductance dynamics of the network
conductance_freq_fileName = 'conductance_freq_static.eps';

%load results/job23502_2010-08-07T19-37-40_output.mat

opt = parseOptions(options);
sheet_size = opt.sheet_size;


close all;
%hold off;

dt = 0.0001; %s
delta_t = 0.25; % s - size of window for firing rate estimation
firingPop_t = 5; % sec
dt_rat = 0.02; %sec

SNList_nID = 51
neuronID = SNList(SNList_nID)

mVolt = 1000;
nS = 1e9;

outputDir = '../../thesis/src/fig/';
    
figure;
plot(SNMonitor_times, SNMonitor_values(SNList_nID, :)*mVolt);

% ------------------------------------------------
% Plot membrane potential and conductance
% ------------------------------------------------
fontSize = 16;
figure('Position', [900, 300, 1100, 500]);

x_limits = [7 10];
y_limits = [-61 -53];

start_i = x_limits(1)/dt + 1;
end_i = x_limits(2)/dt + 1;
signal = SNMonitor_values(SNList_nID, start_i:end_i);

subplot(7, 3, [1 2 4 5 7 8], 'FontSize', fontSize, 'XMinorTick', 'on', 'YMinorTick', 'on');
hold on;
plot(SNMonitor_times(start_i:end_i), signal*mVolt);
plot([SNMonitor_times(start_i) SNMonitor_times(end_i)], [mean(signal)*mVolt mean(signal)*mVolt]);
ylabel('V_{membrane} (mV)');
xlim(x_limits);
%ylim(y_limits);
box on;
title 'A';

subplot(7, 3, [13 14 16 17 19 20], 'FontSize', fontSize, 'XMinorTick', 'on');
hold on;
plot(SNgMonitor_times, SNgMonitor_values(SNList_nID, :)*nS);

xlabel('Time (s)');
ylabel('Syn. conductance (nS)');
xlim(x_limits);
%ylim(y_limits);
box on;
title 'C'

%annotation(gcf(),'line',[0.6613 0.6613],[0.721 0.03474]);

%print('-depsc2', [outputDir 'subthreshold_conductance.eps']);

% ------------------------------------------------
% Plot the Fourier transform of membrane potential
% ------------------------------------------------
[Y, f, NFFT] = fourierTrans(signal, dt);

% Plot single-sided amplitude spectrum.
%figure;
subplot(7, 3, [15 18 21], 'FontSize', fontSize, 'XMinorTick', 'on');
hold on;
plot(f,2*abs(Y(1:NFFT/2+1))*1000);
axis tight;
xlim([0 40]);
xlabel('Frequency (Hz)');
ylabel('Power (\muW)');
box on;
title 'D';

% ------------------------------------------------
% Plot the recording position
% ------------------------------------------------
subplot(7, 3, [3 6 9], 'FontSize', fontSize);
populationResponseFigure
title 'B';

%set(gcf,'PaperPositionMode','auto');
%print('-depsc2', [outputDir conductance_freq_fileName]);


% ------------------------------------------------
% Plot the autocorrelation of the signal used
% ------------------------------------------------

% figure;
% normSignal = signal - mean(signal);
% corr = xcorr(normSignal);
% corrTimes = 0:numel(signal)*2 - 2;
% plot((corrTimes - numel(signal)+1)*dt, corr);
% xlim([0 numel(signal)*dt]);


% ------------------------------------------------
% Plot the spiking neuron
% ------------------------------------------------
%figure;


