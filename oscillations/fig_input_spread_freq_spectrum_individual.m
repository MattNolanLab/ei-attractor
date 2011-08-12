% Create frequency plot of excitatory population responses at specified
% input spread value
close all;
clearvars -except results;

fontSize = 14;


%load e_input_current_output_19-Jul-2011;
outputDir = 'output';

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;

spread_all = 0:0.25:5;
find_eps = 1e-9;
D = results(1, 1).opt.D;

for spread = spread_all
    par_it = 1;
    while par_it <= nParam
        if abs(results(par_it, 1).opt.input_spread/D - spread) < find_eps
            break;
        end;
        par_it = par_it + 1;
    end
    if par_it > nParam
        warning(sprintf('Could not find spread = %f', spread));
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
    print('-depsc2', sprintf('%s/e_input_spread_power_spectrum_detail_spread_%.3fmV.eps', outputDir, opt.input_spread/D));


    figure('Visible', 'off');
    subplot(1, 1, 1, 'FontSize', fontSize);
    plot(f, mean(F_arr));
    xlim(f_lim);
    xlabel('Frequency (Hz)');
    ylabel('Power (Hz^2)');

    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('%s/e_input_spread_power_spectrum_average_spread_%.3fmV.eps', outputDir, opt.input_spread/D));
    
end