close all;
clearvars -except results;

outputDir = 'output_local/';
outputFileRef = '008';
fontSize = 16;

par_it = 21;

nParam  = size(results, 1);
nTrials = size(results, 2);

ddC_eps = -0.005;
min_peak_ratio = 1/10;


parfor trial_it = 1:nTrials
    figure('Position', [800 0 1400 1000], 'Visible', 'off');
    subplot(1,1,1, 'FontSize', fontSize);
    
    res = results(par_it,trial_it);

    signal = res.firingRate_e;
    s_len = numel(signal);
    C = xcorr(signal, 'coeff');

    plot_C = C(s_len:end);
    plot_times = [0:numel(plot_C)-3] * res.opt.dt * 1000;
    dC = diff(plot_C);
    ddC = diff(dC);
    plot(plot_times, [plot_C(1:end-2)]); %; dC(1:end-1); diff(plot_C, 2)]);
    hold on;
    
    dC_zero = sign(dC(2:end) .* dC(1:end-1)) == -1;
    peak_i = find(dC_zero & ddC < ddC_eps);
    C_peak = plot_C(peak_i);
    
    % remove peaks which are insignificant compared to maximal peak
    max_peak = max(C_peak);
    peak_i(find(C_peak < min_peak_ratio)) = [];
    C_peak = plot_C(peak_i);
    
    plot(plot_times(peak_i), C_peak, 'o');
    
    xlim([0 200]);
    ylim([0 1]);
    xlabel('Autocorrelation time shift (ms)');
    ylabel('Normalized correlation');
    
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('%s/%s_e_input_current_autocorr_Ie_%.3fmV_trial_%.3d.eps', outputDir, outputFileRef, res.opt.Ie*1000, trial_it));

end