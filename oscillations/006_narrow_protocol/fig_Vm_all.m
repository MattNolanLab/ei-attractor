close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '001';

fontSize = 16;

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

par_it = 1;
trial_it = 1;

x_lim = [5 7];

res = results(par_it, trial_it);
opt = res.opt;

parfor it = 1:size(res.Vmon.e, 1);
    figure('Position', [800 0 1000 350], 'Visible', 'off');

    subplot(1, 1, 1, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.e(it, :)*1000);
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    box on;
    axis tight;
    xlim(x_lim);

    set(gcf,'PaperPositionMode','auto');    
    print('-depsc2', sprintf('%s/%s_narrow_spiking_exc_spread_%3.3f_trial_%.3d_nid_%.3d.eps', ...
            outputDir, outputNum, opt.input_spread/D, trial_it, res.opt.Emon_i(it)));    
end


parfor it = 1:size(res.Vmon.i, 1);
    figure('Position', [800 0 1000 350], 'Visible', 'off');

    subplot(1, 1, 1, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.i(it, :)*1000);
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    box on;
    axis tight;
    xlim(x_lim);

    set(gcf,'PaperPositionMode','auto');    
    print('-depsc2', sprintf('%s/%s_narrow_spiking_inh_spread_%3.3f_trial_%.3d_nid_%.3d.eps', ...
            outputDir, outputNum, opt.input_spread/D, trial_it, res.opt.Imon_i(it)));    
end