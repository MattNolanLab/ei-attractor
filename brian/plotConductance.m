% Illustrate conductance dynamics of the network
conductance_freq_fileName = 'conductance_freq_static.eps';


close all;
%hold off;

dt = 0.0001; %s

SNList_nID = 19;
neuronID = SNList(SNList_nID);

mVolt = 1000;
nS = 1e9;

outputDir = '../../thesis/src/fig/';
    
figure;
plot(SNMonitor_times, SNMonitor_values(SNList_nID, :)*mVolt, 'LineWidth', 2);

% ------------------------------------------------
% Plot membrane potential and conductance
% ------------------------------------------------
fontSize = 18;
figure('Position', [900, 300, 1100, 500]);

x_limits = [8 9];
y_limits = [-61 -53];

subplot(3, 3, [1 2 4 5], 'FontSize', fontSize, 'XMinorTick', 'on');
hold on;
plot(SNMonitor_times, SNMonitor_values(SNList_nID, :)*mVolt, 'LineWidth', 2);
ylabel('V_{membrane} (mV)');
xlim(x_limits);
%ylim(y_limits);

subplot(3, 3, [7 8], 'FontSize', fontSize, 'XMinorTick', 'on');
hold on;
plot(SNgMonitor_times, SNgMonitor_values(SNList_nID, :)*nS, 'LineWidth', 2);

xlabel('Time (s)');
ylabel('Syn. conductance (nS)');
xlim(x_limits);
%ylim(y_limits);

%annotation(gcf(),'line',[0.6613 0.6613],[0.721 0.03474]);

%print('-depsc2', [outputDir 'subthreshold_conductance.eps']);

% ------------------------------------------------
% Plot the Fourier transform of membrane potential
% ------------------------------------------------
Fs = 1/dt;
start_i = x_limits(1)/dt + 1;
end_i = x_limits(2)/dt + 1;
signal = SNMonitor_values(SNList_nID, start_i:end_i);
signal = signal - mean(signal);
sL = numel(signal);

NFFT = 2^nextpow2(sL); % Next power of 2 from length of y
Y = fft(signal,NFFT)/sL;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
%figure;
subplot(3, 3, [3  6 9], 'FontSize', fontSize, 'XMinorTick', 'on');
hold on;
plot(f,2*abs(Y(1:NFFT/2+1))*1000, 'LineWidth', 2);
xlim([0 40]);
xlabel('Frequency (Hz)');
ylabel('Power (mW)');
%title(['Frequency response of membrane potential recording (' num2str(x_limits(1)) ' - ' num2str(x_limits(2)) ' s)']);

set(gcf,'PaperPositionMode','auto');
print('-depsc2', [outputDir conductance_freq_fileName]);


% ------------------------------------------------
% Plot the autocorrelation of the signal used
% ------------------------------------------------

figure;
normSignal = signal - mean(signal);
corr = xcorr(normSignal);
corrTimes = 0:numel(signal)*2 - 2;
plot((corrTimes - numel(signal)+1)*dt, corr);
xlim([0 numel(signal)*dt]);



