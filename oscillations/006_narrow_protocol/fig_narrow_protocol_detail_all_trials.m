% Create population detailed plot for all trials, specified by Ie
% value
close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local/';
outputNum = '012';

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

dt = results(1,1).opt.dt;

N_spikes = 25;

t_start = 0;
t_end   = 10;
f_lim = [0 200];


parfor par_it = 1:nParam
    T_i = struct();
    T_i.start = t_start/dt + 1;
    T_i.end   = t_end/dt + 1;
    
    for trial_it = 1:nTrials
        res = results(par_it, trial_it);

        x_lim = [t_start t_end];
        opt = res.opt;

        figure('Position', [0 0 1680 1050], 'Visible', 'off');
        fig_populationDetail(res, x_lim, T_i, fontSize, N_spikes);
        

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'Painters');
        print('-depsc2', sprintf('%s/%s_population_detail_par_it_%.3d_trial_%.3d.eps', ...
            outputDir, outputNum, par_it, trial_it));
        
        
%         % Cv of excitatory and inhibitory neurons
%         min_stat = 50;
%         nBins = 20;
%         figure('Position', [900 400 1000 400], 'Visible', 'off');
%         fig_Cv_stat(res, min_stat, nBins, fontSize);
%         set(gcf,'PaperPositionMode','auto');
%         print('-depsc2', sprintf('%s/%s_e_input_spread_Cv_spread_%3.3f_trial_%.3d.eps', ...
%             outputDir, outputNum, opt.input_spread/D, trial_it));    

    end
end