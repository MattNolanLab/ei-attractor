% Print current clamp data
close all;
clearvars -except c001* c002* c003* c004*

outputDir = '.';

pA = 1e12;
spike_th = 0;

fontSize = 18;
x_lim = [0 20];

times = c001_Time;
dt = times(2) - times(1);
sig_st = c003_Current_2;
sig_in = c002_Current_1;


figure('Position', [800 0 800 600]);
subplot(2,1,1, 'FontSize', fontSize);
plot(times, sig_st*pA, 'r');
%hold on;
ylabel('Current (pA)');
xlim(x_lim);
box off;

% [spike_t spike_id] = spikeDetect(sig_st, times, spike_th);
% plot(times(spike_id), sig_st(spike_id)*mV, 'ok');
% hold off;

subplot(2, 1, 2, 'FontSize', fontSize);
plot(times, sig_in*pA, 'b');
ylabel('Current (pA)');
xlim(x_lim);
box off;
xlabel('Time (s)');



set(gcf,'PaperPositionMode','auto');
print('-depsc2', sprintf('%s/stellate_interneuron_firing_rate.eps', ...
        outputDir));    
    
    
% Estimate voltage slope of stellate cells and interneurons
st_t2 = 4;
st_t1 = 3;

st_dV = sig_st(fix(st_t2/dt)+1) - sig_st(fix(st_t1/dt)+1);
st_dVdt = st_dV/(st_t2 - st_t1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Interneuron currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot just interneuron currents and low pass filtered version
figure('Position', [800 0 800 600]);
subplot(1,1,1, 'FontSize', fontSize);
%plot(times, sig_in*pA, 'b');
leg{1} = 'Original';
hold all;
leg_it = 2;
for cutoff = [4000]
    Fs = 1/dt;
    filt_N = 4;
    %cutoff = 25;
    [b, a] = butter(filt_N, cutoff/(Fs/2), 'low');
    sig_in_filt = filtfilt(b, a, sig_in);
    plot(times, sig_in_filt*pA, 'LineWidth', 2);
    hold all;
    xlabel('Time (s)');
    ylabel('Current (pA)');
    leg{leg_it} = sprintf('cutoff: %d Hz', cutoff);
    leg_it = leg_it + 1;
    
end
%legend(leg);
xlim([11.1 11.15]);
%grid minor;
%title('Interneuron');


hold on;

% Establish noise variance
noise_t_start = 0;
noise_t_end = 2;
noise_std = std(sig_in(fix(noise_t_start/dt+1):fix(noise_t_end/dt+1)));

% Find positions of peaks of EPSCs;
[epsc_id epsc_size epsc_start_id]= findEPSCs(sig_in_filt, noise_std);


plot(epsc_id*dt, sig_in(epsc_id)*pA, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(epsc_start_id*dt, sig_in(epsc_start_id)*pA, 'ko', 'MarkerSize', 10);

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'fig_epsc_detection.eps');


figure('Position', [800 0 600 800]);

%N_log = N;
%N_log(N_log == 0) = 1e-10;
sample_size = 10000;

logn_fit = lognfit(epsc_size*pA);
exp_fit = expfit(epsc_size*pA);
subplot(3, 1, 1, 'FontSize', fontSize);
[N_orig X_orig] = hist(epsc_size*pA, 100);
bar(X_orig, N_orig/numel(epsc_size));
xlim([0 1200]);
title('Data');
subplot(3, 1, 2, 'FontSize', fontSize);
[N_logn X_logn] = hist(lognrnd(logn_fit(1), logn_fit(2), 1, sample_size), 100);
bar(X_logn, N_logn/sample_size);
xlim([0 1200]);
title('Log normal');
subplot(3, 1, 3, 'FontSize', fontSize);
[N_exp, X_exp] = hist(exprnd(exp_fit, 1, sample_size), 100);
bar(X_exp, N_exp/sample_size);
xlabel('EPSC size (pA)');
xlim([0 1200]);
title('Exponential');

MU = logn_fit(1);
SIGMA = logn_fit(2);
M = exp(MU + SIGMA^2/2)
V = exp(2*MU + SIGMA^2) * (exp(SIGMA^2) - 1)

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'fig_epsc_size_distribution.eps');


%plot(X, N , 'o', X, 500*exp(-1/120*X));
%title('Interneuron');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Stellate currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot just interneuron currents and low pass filtered version
figure();
plot(times, sig_st*pA, 'r');
leg{1} = 'Original';
hold all;
leg_it = 2;
for cutoff = [4000]
    Fs = 1/dt;
    filt_N = 4;
    %cutoff = 25;
    [b, a] = butter(filt_N, cutoff/(Fs/2), 'low');
    sig_st_filt = filtfilt(b, a, sig_st);
    plot(times, sig_st_filt*pA, 'LineWidth', 1);
    hold all;
    xlabel('Time (s)');
    ylabel('Current (pA)');
    leg{leg_it} = sprintf('cutoff: %d Hz', cutoff);
    leg_it = leg_it + 1;
    
end
legend(leg);
%xlim([11.1 11.15]);
grid minor;
title('Stellate cell');


hold on;

% Establish noise variance
noise_t_start = 0;
noise_t_end = 2;
noise_std = std(sig_st(fix(noise_t_start/dt+1):fix(noise_t_end/dt+1)));

% Find positions of peaks of EPSCs;
[ipsc_id ipsc_size ipsc_start_id]= findIPSCs(sig_st_filt, noise_std);


plot(ipsc_id*dt, sig_st(ipsc_id)*pA, 'ro');
plot(ipsc_start_id*dt, sig_st(ipsc_start_id)*pA, 'ko');

% figure();
% [N X] = hist(epsc_size*pA, 100);
% 
% N_log = N;
% N_log(N_log == 0) = 1;
% P = polyfit(X, log(N_log), 1);
% exp_fit = exp(P(2) + P(1)*X);
% plot(X, N , 'o', X, 500*exp(-1/120*X));
