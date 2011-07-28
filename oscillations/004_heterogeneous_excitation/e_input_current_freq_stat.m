% Process input current simulation experiment
% Heterogeneous excitation to inhibitory and excitatory neurons
close all;
clearvars -except results;

path('..', path);

%load e_input_current_output_19-Jul-2011;

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;
nPar = 1;

dc_ratio = 1/20; % asynchronous mode detection
win_len = 0.002;


t_start = 1.5;
t_end   = 2.5;

for par_it = 1:nPar
    for trial_it = 1:nTrials
        trial_it;
        res = results(par_it, trial_it);
        
        t_start_i = t_start/res.opt.dt + 1;
        t_end_i   = t_end/res.opt.dt + 1;
    
        %firingRate_e = getFiringRate(res.spikeRecord_e(:, t_start_i:t_end_i), res.opt.dt, win_len);
        firingRate_e = res.firingRate_e(t_start_i:t_end_i);
        opt = res.opt;

        
        % Population frequency
%         [Y f NFFT] = fourierTrans(firingRate_e, opt.dt);
%         Y_abs = 2*abs(Y(1:NFFT/2+1));
%         Y_abs = Y_abs.^2;
% 
%         
%         [maxF maxFI] = max(Y_abs);
        fmax(trial_it, par_it) = getPopOscillationFreq(firingRate_e, dc_ratio, opt.dt);

        
        spikeRecord_e = res.spikeRecord_e(:, t_start_i:t_end_i);
        spikeRecord_i = res.spikeRecord_i(:, t_start_i:t_end_i);        
        times = res.times(t_start_i:t_end_i);

        spikeCnt_i = sum(spikeRecord_i');
        
        % mean firing rates of neurons in this trial
        mfr_T = t_end - t_start;

        
        e_mfr_all(:, trial_it, par_it) = full(sum(spikeRecord_e')/mfr_T);
        i_mfr_all(:, trial_it, par_it) = full(sum(spikeRecord_i')/mfr_T);
        e_mfr(trial_it, par_it) = mean(e_mfr_all(:, trial_it, par_it));
        i_mfr(trial_it, par_it) = mean(i_mfr_all(:, trial_it, par_it));
    end
end

nan_fmax = mean(fmax);
nan_fmax(find(isnan(mean(fmax)))) = 0;
nan_fmax(find(nan_fmax > 0)) = nan;

% Print the population and excitatory cells frequency depending on input
% parameter
for par_it = 1:nPar
    Ie_max(par_it) = results(par_it, 1).opt.Ie_max;
end

figure('Position', [840 800 800 500]);
subplot(1, 1, 1, 'FontSize', fontSize);
hold on;
plot_h = errorbar([Ie_max*1000; Ie_max*1000; Ie_max*1000]', ...
    [mean(fmax); mean(e_mfr); mean(i_mfr)]', ...
    [std(fmax); std(e_mfr); std(i_mfr)]', ...
    '-o', 'LineWidth', 1);
hold on;
plot(Ie_max*1000, nan_fmax, '-o', 'LineWidth', 1);

%errorbar_tick(plot_h, 80);
%errorbar(Ie*1000, mean(e_mfr), std(e_mfr));
xlabel('Max. input drive (mV)');
ylabel('Frequency (Hz)');
legend('Oscillation', 'E firing rate', 'I firing rate', 'Location', 'SouthEast');
axis tight;


% % Create histograms of maximum oscillation frequency for each parameter
% % value
% sp_cols = 5;
% sp_rows = ceil(nPar/sp_cols);
% figure('Position', [840 800 1100 1000]);
% fontSize = 14;
% for par_it = 1:nPar
%     subplot(sp_rows, sp_cols, par_it, 'FontSize', fontSize);
%     hist(fmax(:, par_it));
%     title(sprintf('Ie = %.2f mV', results(par_it, 1).opt.Ie*1000));
% end
% 
% 
% % Histograms of average firing frequency
% sp_cols = 5;
% sp_rows = ceil(nTrials/sp_cols);
% figure('Position', [840 800 1100 1000]);
% fontSize = 14;
% par_it = 5;
% for trial_it = 1:nTrials
%     subplot(sp_rows, sp_cols, trial_it, 'FontSize', fontSize);
%     hist(i_mfr_all(:, trial_it, par_it));
%     title(sprintf('f = %.2f Hz', fmax(trial_it, par_it)));
% end
% set(gcf,'PaperPositionMode','auto');
% %print('-depsc2', sprintf('output/2011-07-19/e_input_current_interneuron_firing_rates_trials_Ie_%.3f.eps', results(par_it, 1).opt.Ie*1000));
% 
