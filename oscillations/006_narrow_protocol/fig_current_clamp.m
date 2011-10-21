% Plot lag of onset of spiking between interneurons and stellate cells, as
% a function of excitatory synaptic strength x density (sparsity)

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '009';

nParam  = size(results, 1);
nTrials = size(results, 2);

mA = 1000;
fontSize = 14;

spread_all = [0.5 2 5 10];

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);

Nwe = size(results(1,1).opt.we_vec, 2);
Nsp = size(results(1,1).opt.sparseness_vec, 2);

par_it = 20;
trial_it = 1;

Ni_it = 1;
Ne_it = 1;

res = results(par_it, trial_it);

ge_sig = res.Vmon.ge(Ni_it, :)*mA;
gi_sig = -res.Vmon.gi(Ne_it, :)*mA;


% Time 1
figure('Position', [800 800 1000 600]);
x_lim = [7.8 8];
subplot(2,1,1, 'FontSize', fontSize);
plot(res.times, gi_sig);
xlabel('Time (s)');
ylabel('Syn. current (mA)');
title('Stellate cell');
xlim(x_lim);

subplot(2, 1, 2, 'FontSize', fontSize);
plot(res.times, ge_sig);
xlabel('Time (s)');
ylabel('Syn. current (mA)');
title('Interneuron');
xlim(x_lim);

set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
print('-depsc2', sprintf('%s/%s_gi_ge_currents_time1.eps', ...
        outputDir, outputNum));    


% Time 2
figure('Position', [800 800 1000 600]);
x_lim = [14 14.2];
subplot(2,1,1, 'FontSize', fontSize);
plot(res.times, gi_sig);
xlabel('Time (s)');
ylabel('Syn. current (mA)');
title('Stellate cell');
xlim(x_lim);

subplot(2, 1, 2, 'FontSize', fontSize);
plot(res.times, ge_sig);
xlabel('Time (s)');
ylabel('Syn. current (mA)');
title('Interneuron');
xlim(x_lim);

set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
print('-depsc2', sprintf('%s/%s_gi_ge_currents_time2.eps', ...
        outputDir, outputNum));    



for Ne_it = 1%:numel(res.opt.Emon_i)
    figure('Position', [800 800 1000 600]);

    gi_sig = - res.Vmon.gi(Ne_it, :)*mA;

    plot_opt.fontSize = 14;
    plot_opt.F = 10:1:100;
    plot_opt.sampling_rate = 1e4;
    plot_opt.win_len = 3000;
    
    which_cell = 'stellate';
    
    
    plot_opt.plot_dB = false;
    plot_opt.x_lim = [0 10];
    
    [Y F T P P_plot] = plot_g_freq(gi_sig, res.times, which_cell, plot_opt);

    set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
    print('-depsc2', sprintf('%s/%s_gi_spectrogram_N%.3d.eps', ...
            outputDir, outputNum, Ne_it));    

        
    figure('Position', [800 800 1000 600]);
    %gap = 0.3;
    f_ax = subplot(3,1,[1 2], 'FontSize', fontSize);
    
    [fmax fmax_i] = max(P_plot);
    plot(T, F(fmax_i))
    
    P_ax = subplot(3, 1, 3, 'FontSize', fontSize);
    plot(T, fmax);
    
    f_pos = get(f_ax, 'Position');
    P_pos = get(P_ax, 'Position');

    %f_pos(2) = P_pos(2) + P_pos(4) + gap;
    %set(f_ax, 'Position', f_pos, 'XTick', []);
    xlabel('Time (s)');    

end


figure('Position', [800 800 1000 600]);
which_cell = 'interneuron';
plot_g_freq(res.Vmon.ge(Ni_it, :)*mA, res.times, which_cell, plot_opt);



