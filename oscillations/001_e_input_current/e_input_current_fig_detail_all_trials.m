% Create population detailed plot for all trials, specified by Ie
% value
close all;
clearvars -except results;


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local/';
outputFileRef = '009';

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

NCells = 25;

Ie_all = [19.45] * 1e-3;
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

    x_lim = [2 2.2];

    
    parfor trial_it = 1:nTrials
        res = results(par_it, trial_it);
        opt = res.opt;

        % Plot detailed responses of a randomly selected trial
        figure('Position', [800 0 400 900], 'Visible', 'off');
        fig_populationDetail(res, x_lim, nan, fontSize, NCells);
        
        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('%s/%s_e_input_current_population_detail_Ie_%.3fmV_trial_%.3d.eps', outputDir, outputFileRef, opt.Ie*1000, trial_it));
    end
end