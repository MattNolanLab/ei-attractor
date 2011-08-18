% Create population detailed plot for all trials, specified by Ie
% value
close all;
clearvars -except results;

path('..', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '003';

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

spread_all = [0.5 2 4 5 10];
find_eps = 1e-9;
D = results(1, 1).opt.D;
dt = results(1,1).opt.dt;

t_start = 58;
t_end   = 60;
f_lim = [0 200];


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

    T_i.start = t_start/dt + 1;
    T_i.end   = t_end/dt + 1;
    
    parfor trial_it = 1:nTrials
        res = results(par_it, trial_it);

        
        % Plot detailed responses of a randomly selected trial
        x_lim = [t_start t_end];
        opt = res.opt;

        figure('Position', [800 0 1400 1000], 'Visible', 'off');
        fig_populationDetail(res, x_lim, T_i, fontSize);
        

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'Painters');
        print('-depsc2', sprintf('%s/%s_e_input_spread_population_detail_spread_%3.3f_trial_%.3d.eps', ...
            outputDir, outputNum, opt.input_spread/D, trial_it));
        
        
        % Cv of excitatory and inhibitory neurons
        min_stat = 50;
        nBins = 20;
        figure('Position', [900 400 1000 400], 'Visible', 'off');
        fig_Cv_stat(res, min_stat, nBins, fontSize);
        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('%s/%s_e_input_spread_Cv_spread_%3.3f_trial_%.3d.eps', ...
            outputDir, outputNum, opt.input_spread/D, trial_it));    

    end
end