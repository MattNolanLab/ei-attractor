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


res = results(par_it, trial_it);
opt = res.opt;

x_lim = [4 8];
ti_start = x_lim(1)/opt.dt + 1;
ti_end = x_lim(2)/opt.dt + 1;
win_len = 0.5 / dt; %s


parfor it1 = 1:size(res.Vmon.e, 1)
    for it2 = it1+1:size(res.Vmon.e, 1)
        figure('Position', [0 0 1600 800], 'Visible', 'off');
        
        s1 = res.Vmon.e(it1, ti_start:ti_end);
        s2 = res.Vmon.e(it2, ti_start:ti_end);

        subplot(3, 1, 1, 'FontSize', fontSize);
        plot(res.Vmon.t(ti_start:ti_end), [s1; s2]*1000);
        ylabel('Vm (mV)');
        box on;
%        axis tight;
        xlim(x_lim);
        title(sprintf('Excitatory neurons %d:%d', res.opt.Emon_i(it1), res.opt.Emon_i(it2)));
        
        subplot(3, 1, 2, 'FontSize', fontSize);
        plot(res.Vmon.t(ti_start:ti_end-win_len), slidingCorrelation(s1, s2, win_len));
        %xlabel('Time (s)');
        ylabel('Corr. coefficient');
        box on;
        xlim(x_lim);
        
        subplot(3, 1, 3, 'FontSize', fontSize);
        plot(res.times(ti_start:ti_end), res.firingRate_i(ti_start:ti_end));
        ylabel('Firing rate (Hz)');
        xlim(x_lim);
        box on;
        xlabel('Time (s)');
        title('Interneurons');



        set(gcf,'PaperPositionMode','auto');    
        print('-depsc2', sprintf('%s/%s_narrow_2Vm_exc_%3.3f_trial_%.3d_nid_%.3d_%.3d.eps', ...
                outputDir, outputNum, opt.input_spread/D, trial_it, res.opt.Emon_i(it1), res.opt.Emon_i(it2)));    
    end
end


parfor it1 = 1:size(res.Vmon.i, 1)
    for it2 = it1+1:size(res.Vmon.i, 1)
        figure('Position', [0 0 1600 800], 'Visible', 'off');
        
        s1 = res.Vmon.i(it1, ti_start:ti_end);
        s2 = res.Vmon.i(it2, ti_start:ti_end);

        subplot(3, 1, 1, 'FontSize', fontSize);
        plot(res.Vmon.t(ti_start:ti_end), [s1; s2]*1000);
        ylabel('Vm (mV)');
        box on;
%        axis tight;
        xlim(x_lim);
        title(sprintf('Interneurons neurons %d:%d', res.opt.Imon_i(it1), res.opt.Imon_i(it2)));
        
        subplot(3, 1, 2, 'FontSize', fontSize);
        plot(res.Vmon.t(ti_start:ti_end-win_len), slidingCorrelation(s1, s2, win_len));
        %xlabel('Time (s)');
        ylabel('Corr. coefficient');
        box on;
        xlim(x_lim);
        
        subplot(3, 1, 3, 'FontSize', fontSize);
        plot(res.times(ti_start:ti_end), res.firingRate_i(ti_start:ti_end));
        ylabel('Firing rate (Hz)');
        xlim(x_lim);
        box on;
        xlabel('Time (s)');

        
        set(gcf,'PaperPositionMode','auto');    
        print('-depsc2', sprintf('%s/%s_narrow_2Vm_inh_%3.3f_trial_%.3d_nid_%.3d_%.3d.eps', ...
                outputDir, outputNum, opt.input_spread/D, trial_it, res.opt.Emon_i(it1), res.opt.Emon_i(it2)));    
    end
end