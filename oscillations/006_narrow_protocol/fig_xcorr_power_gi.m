% Plot cross correlation + gamma power/spectrograms

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '010';

nParam  = size(results, 1);
nTrials = size(results, 2);

pA = 1e12;
fontSize = 14;

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);

Nwe = size(results(1,1).opt.we_vec, 2);
Nsp = size(results(1,1).opt.sparseness_vec, 2);

sp_vec = results(1,1).opt.sparseness_vec;
we_vec = results(1,1).opt.we_vec;


F = 10:1:100;
sampling_rate = 1e4;
spec_win_len = 2000;
noverlap = spec_win_len /2;
plot_dB = true;

dt = results(1,1).opt.dt;
nPeriods = 1;
t_start = 1;
t_start_i = fix(t_start/dt) + 1;

xcorr_win_len_t = 0.5;
xcorr_win_len = fix(xcorr_win_len_t / dt);

% Simulation trial selection
sp = 0.05;
we = 1.5e-11;
par_it = getItFromSpWeVecs(sp, we, results, 1e-9);
trial_it = 1;

Ne_it = [1 2];

res = results(par_it, trial_it);

gi_sig1 = res.Vmon.gi(Ne_it(1), :)*pA;
gi_sig2 = res.Vmon.gi(Ne_it(2), :)*pA;

xcorr = slidingCorrelation(gi_sig1, gi_sig2, xcorr_win_len);


figure('Position', [800 800 1000 700]);
subplot(4, 1, [1 2], 'FontSize', fontSize);
plot(res.times, [gi_sig1; gi_sig2]);
ylabel('Syn. current (pA)');

subplot(4, 1, 3, 'FontSize', fontSize);
plot(res.times,  xcorr);
ylabel('Correlation coeff.');

subplot(4, 1, 4, 'FontSize', fontSize);
plotSpectrogramIntoAx(gi_sig1, spec_win_len, noverlap, F, sampling_rate, plot_dB);
xlim([0 15]);


