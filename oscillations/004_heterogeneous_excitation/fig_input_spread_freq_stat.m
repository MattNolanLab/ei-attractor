% Process input spread simulation experiment
% Heterogeneous excitation to inhibitory and excitatory neurons
close all;
clearvars -except results;

path('../include/', path);

%load e_input_current_output_19-Jul-2011;

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;
nPar = nParam;

dc_ratio = 1/15; % asynchronous mode detection
win_len = 0.002;

autoCorrNPeaks = 10;


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
        [fmax(trial_it, par_it) coherence(trial_it, par_it)] = getPopOscFreqAutoCorr(firingRate_e, opt.dt, autoCorrNPeaks);

        
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

% nan_fmax = mean(fmax, 1);
% nan_fmax(find(isnan(mean(fmax, 1)))) = 0;
% nan_fmax(find(nan_fmax > 0)) = nan;

% Print the population and excitatory cells frequency depending on input
% parameter
D = results(1, 1).opt.D;  % radius of the microcircuit
is_vec = results(1, 1).opt.input_spread_vec / D;

figure('Position', [840 800 800 900]);
subplot(5, 1, [1 2 3 4], 'FontSize', fontSize);
grid on;
plot_h = errorbar([is_vec; is_vec; is_vec]', ...
    [mean(fmax, 1); mean(e_mfr, 1); mean(i_mfr, 1)]', ...
    [std(fmax, 0, 1); std(e_mfr, 0, 1); std(i_mfr, 0, 1)]', ...
    '-o', 'LineWidth', 1);

xlabel('Normalized input spread');
ylabel('Frequency (Hz)');
legend('Oscillation', 'E firing rate', 'I firing rate', 'Location', 'SouthEast');
axis tight;
xlim([opt.input_spread_vec(1)/D opt.input_spread_vec(end)/D]);


subplot(5,1,[5], 'FontSize', fontSize);
plot(is_vec, mean(coherence), '-o');
ylabel('Coherence');
xlabel('Input drive (mV)');
legend('E population', 'Location', 'SouthEast');
grid on;
