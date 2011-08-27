% Process input current simulation experiment
close all;
clearvars -except results;

path('../include', path);

%load e_input_current_output_19-Jul-2011;

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;
nPar = 31;

freqEstnPeaks = 10;


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

        
        [fmax(trial_it, par_it) coherence(trial_it, par_it)] = getPopOscFreqAutoCorr(firingRate_e, opt.dt, freqEstnPeaks);

        
        spikeRecord_e = res.spikeRecord_e(:, t_start_i:t_end_i);
        spikeRecord_i = res.spikeRecord_i(:, t_start_i:t_end_i);        

        spikeCnt_i = sum(spikeRecord_i');
        
        % mean firing rates of neurons in this trial
        mfr_T = t_end - t_start;
        
        e_mfr_all(:, trial_it, par_it) = full(sum(spikeRecord_e')/mfr_T);
        i_mfr_all(:, trial_it, par_it) = full(sum(spikeRecord_i')/mfr_T);
        e_mfr(trial_it, par_it) = mean(e_mfr_all(:, trial_it, par_it));
        i_mfr(trial_it, par_it) = mean(i_mfr_all(:, trial_it, par_it));
    end
end


% Print the population and excitatory cells frequency depending on input
% parameter
for par_it = 1:nPar
    Ie(par_it) = results(par_it, 1).opt.Ie;
end

figure('Position', [840 800 800 900]);
subplot(5, 1, [1 2 3 4], 'FontSize', fontSize);
hold on;
plot_h = errorbar([Ie*1000; Ie*1000; Ie*1000]', ...
    [mean(fmax); mean(e_mfr); mean(i_mfr)]', ...
    [std(fmax); std(e_mfr); std(i_mfr)]', ...
    '-o', 'LineWidth', 1);

ylabel('Frequency (Hz)');
legend('E pop oscillation', 'E firing rate', 'I firing rate', 'Location', 'NorthWest');
xlim([Ie(1) Ie(end)]*1000);
grid on;
box on;

subplot(5,1,[5], 'FontSize', fontSize);
plot(Ie*1000, mean(coherence), '-o');
ylabel('Coherence');
xlabel('Input drive (mV)');
legend('E population', 'Location', 'SouthWest');
xlim([Ie(1) Ie(end)]*1000);
grid on;
