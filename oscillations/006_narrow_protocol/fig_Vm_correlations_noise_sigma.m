close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '003';

nParam  = size(results, 1);
nTrials = size(results, 2);

spread_all = [0.5 2 5 10];
find_eps = 1e-9;
D = results(1, 1).opt.D;
dt = results(1,1).opt.dt;

N_spikes = 25;

t_start = 5;
t_end   = 7;
f_lim = [0 200];

%par_it = 1;
trial_it = 1;


for par_it = 1:numel(results(1,1).opt.noise_sigma_vec)
    res = results(par_it, trial_it);
    opt = res.opt;

    x_lim = [4 10];

    ti_start = x_lim(1)/opt.dt + 1;
    ti_end = x_lim(2)/opt.dt + 1;
    win_len = 0.5 / dt; %s

    plot_opt.fontSize = 16;
    plot_opt.x_lim = x_lim;

    
    parfor it1 = 1:size(res.Vmon.e, 1)
        for it2 = it1+1:size(res.Vmon.e, 1)
            figure('Position', [0 0 1600 800], 'Visible', 'off');

            Vm_correlations_EE(res, it1, it2, ti_start, ti_end, win_len, plot_opt);

            set(gcf,'PaperPositionMode','auto');    
            print('-depsc2', sprintf('%s/%s_narrow_ramp_2Vm_exc_noise_%3.3f_trial_%.3d_nid_%.3d_%.3d.eps', ...
                    outputDir, outputNum, opt.noise_sigma*1000, trial_it, res.opt.Emon_i(it1), res.opt.Emon_i(it2)));    
        end
    end


    parfor it1 = 1:size(res.Vmon.i, 1)
        for it2 = it1+1:size(res.Vmon.i, 1)
            figure('Position', [0 0 1600 800], 'Visible', 'off');

            Vm_correlations_II(res, it1, it2, ti_start, ti_end, win_len, plot_opt);

            set(gcf,'PaperPositionMode','auto');    
            print('-depsc2', sprintf('%s/%s_narrow_ramp_2Vm_inh_noise_%3.3f_trial_%.3d_nid_%.3d_%.3d.eps', ...
                    outputDir, outputNum, opt.noise_sigma*1000, trial_it, res.opt.Emon_i(it1), res.opt.Emon_i(it2)));    
        end
    end
end