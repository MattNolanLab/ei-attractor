% Print current clamp data
close all;
clearvars -except c001* c002* c003* c004*

outputDir = '.';

mV = 1e3;
spike_th = 0;

fontSize = 18;
x_lim = [0 25];

times = c001_Time;
dt = times(2) - times(1);
sig_st = c003_Membrane_Voltage_2;
sig_in = c002_Membrane_Voltage_1;


figure('Position', [800 0 1000 600]);
subplot(3,1,1, 'FontSize', fontSize);
plot(times, sig_st*mV, 'r');
%hold on;
%xlabel('Time (s)');
ylabel('V_m (mV)');
xlim(x_lim);

% [spike_t spike_id] = spikeDetect(sig_st, times, spike_th);
% plot(times(spike_id), sig_st(spike_id)*mV, 'ok');
% hold off;

subplot(3, 1, 2, 'FontSize', fontSize);
plot(times, sig_in*mV, 'b');
ylabel('V_m (mV)');
xlim(x_lim);


subplot(3, 1, 3, 'FontSize', fontSize);
spikeRecord_st = sparse(1, numel(times));
spikeRecord_st(spikeDetect(sig_st, times, spike_th)) = 1;

spikeRecord_in = sparse(1, numel(times));
spikeRecord_in(spikeDetect(sig_in, times, spike_th)) = 1;

win_len = fix(1/dt);
noverlap = win_len/2;
[R POS] = slidingFiringRate([spikeRecord_st; spikeRecord_in], win_len, noverlap);
R = R / dt;

plot(POS*dt, R(:, 1), 'r', POS*dt, R(:, 2), 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
xlim(x_lim);

set(gcf,'PaperPositionMode','auto');
print('-depsc2', sprintf('%s/stellate_interneuron_firing_rate.eps', ...
        outputDir));    
    
    
% Estimate voltage slope of stellate cells and interneurons
st_t2 = 4;
st_t1 = 3;

st_dV = sig_st(fix(st_t2/dt)+1) - sig_st(fix(st_t1/dt)+1);
st_dVdt = st_dV/(st_t2 - st_t1)

