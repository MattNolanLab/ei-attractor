% Create population detailed plot for all trials, specified by Ie
% value
close all;
clearvars -except results;

path('..', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

spread_all = 0:0.25:5;
find_eps = 1e-9;
D = results(1, 1).opt.D;
dt = results(1,1).opt.dt;

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

    T_i.start = t_start/dt + 1;
    T_i.end   = t_end/dt + 1;
    
    parfor trial_it = 1:nTrials
        res = results(par_it, trial_it);

        
        % Plot detailed responses of a randomly selected trial
        x_lim = [1.5 2.5];
        opt = res.opt;

        figure('Position', [800 0 1400 1000], 'Visible', 'off');
        fig_populationDetail(res, x_lim, T_i, fontSize);
        

        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('%s/002_e_input_spread_population_detail_spread_%3.3f_trial_%.3d.eps', outputDir, opt.input_spread/D, trial_it));
    end
end