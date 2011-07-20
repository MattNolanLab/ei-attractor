% Create frequency plot of excitatory population responses at specified Ie
% value
close all;
clearvars -except results;

fontSize = 14;


%load e_input_current_output_19-Jul-2011;

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;

Ie_all = [18.6 18.8 19.2 19.6 19.8] .* 1e-3; 
find_eps = 1e-9;

for Ie = Ie_all
    par_it = 1;
    while par_it <= nParam
        if abs(results(par_it, 1).opt.Ie - Ie) < find_eps
            break;
        end;
        par_it = par_it + 1;
    end
    if par_it > nParam
        warning(sprintf('Could not find Ie = %f mV', Ie*1000));
        continue;
    end;

    t_start = 1.5;
    t_end   = 2.5;
    f_lim = [0 200];

    sp_cols = 5;
    sp_rows = ceil(nTrials/sp_cols);

    % Plot big figure with all the trials
    figure('Position', [600 0 1500 800], 'Visible', 'off');
    for trial_it = 1:nTrials
        res = results(par_it, trial_it);

        t_start_i = t_start/res.opt.dt + 1;
        t_end_i   = t_end/res.opt.dt + 1;

        firingRate_e = res.firingRate_e(t_start_i:t_end_i);
        opt = res.opt;


        % Population frequency
        [Y f NFFT] = fourierTrans(firingRate_e, opt.dt);
        Y_abs = 2*abs(Y(1:NFFT/2+1));
        Y_abs = Y_abs.^2;

        subplot(sp_rows, sp_cols, trial_it, 'FontSize', fontSize);
        plot(f,Y_abs);
        xlim(f_lim);
        %xlabel('Frequency (Hz)');
        %ylabel('Power (Hz^2)');
        %title(['Excitatory population power spectrum. Ie = ' num2str(opt.Ie) 'mV, trial: ' num2str(trial_it)]);
        %title(['Ie = ' num2str(opt.Ie) 'mV, trial: ' num2str(trial_it)]);

        F_arr(trial_it, :) = Y_abs;
    end

    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('output/2011-07-19/e_input_current_power_spectrum_detail_Ie_%.3fmV.eps', opt.Ie*1000));


    figure('Visible', 'off');
    subplot(1, 1, 1, 'FontSize', fontSize);
    plot(f, mean(F_arr));
    xlim(f_lim);
    xlabel('Frequency (Hz)');
    ylabel('Power (Hz^2)');

    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('output/2011-07-19/e_input_current_power_spectrum_average_Ie_%.3fmV.eps', opt.Ie*1000));
end