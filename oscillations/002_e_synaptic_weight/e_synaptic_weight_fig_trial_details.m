% Synaptic weight simulations
% Plot detailed responses for each trial

close all;
clearvars -except results;

path('..', path);

fontSize = 14;


%load e_input_current_output_19-Jul-2011;

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;

we_all = results(1, 1).opt.we_vec;
find_eps = 1e-9;

for we = we_all
    par_it = 1;
    while par_it <= nParam
        if abs(results(par_it, 1).opt.we - we) < find_eps
            break;
        end;
        par_it = par_it + 1;
    end
    if par_it > nParam
        warning(sprintf('Could not find we = %f mV', we*1000));
        continue;
    end;

    t_start = 1.5;
    t_end   = 2.5;
    f_lim = [0 200];

    
    % Plot detailed responses of a randomly selected trial
    res = results(par_it, 1);
    opt = res.opt;
    
    plot_options.x_lim = [1.5 2.5];
    plot_options.fontSize = fontSize;
    
    figure('Position', [800 0 1400 1000]);
    plotTrialDetails(results, par_it, 1, plot_options);
end