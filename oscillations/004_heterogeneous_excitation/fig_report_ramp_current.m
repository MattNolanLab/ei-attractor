% Create population detailed plot for all trials, specified by Ie
% value
close all;
clearvars -except results;

path('../include', path);

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

N = 25;


x_lim = [4.1 4.7];


%
% Plot detailed responses of selected trial
%

figure('Position', [800 0 1000 800]);

par_it = 1;
trial_it = 2;  

    res = results(par_it, trial_it);
    opt = res.opt;


    subplot(6, 1, [1 2], 'FontSize', fontSize);
    spikeCellRasterPlot(res.spikeCell_e, '.', N);
    title('Principal cells');
    %ylabel('Neuron no.');
    xlim(x_lim);
    ylim([1 N]);
    box on;
    set(gca, 'Xtick', []);
    set(gca, 'Ytick', []);
%        set(gca, 'Position', [0.13 0.754 0.85 0.154]);


    subplot(6, 1, 3, 'FontSize', fontSize);
    plot(res.times, res.firingRate_e);
    ylabel('Firing rate (Hz)');
    box on;
    set(gca, 'Xtick', []);
    axis tight;
    xlim(x_lim);
%        set(gca, 'Position', [0.13 0.505 0.85 0.097336845583677]);


    subplot(6, 1, [4 5], 'FontSize', fontSize);
    spikeCellRasterPlot(res.spikeCell_i, '.', N);
    title('Interneurons');
    %ylabel('Neuron no.');
    xlim(x_lim);
    ylim([1 N]);
    set(gca, 'Xtick', []);
    set(gca, 'Ytick', []);
    box on;


    subplot(6, 1, 6, 'FontSize', fontSize);
    plot(res.times, getFiringRate(res.spikeRecord_i, res.opt.dt, res.opt.rateWindowLen));
    ylabel('Firing rate (Hz)');
    xlim(x_lim);
    box on;
    xlabel('Time (s)');


    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', 'output/final_ramp_current.eps');

    
 
    % Two excitatory neurons and an inhibitory one
    e_nid = [100 200];
    i_nid = 100;
    
    figure('Position', [800 0 1000 700]);

    subplot(2, 1, 1, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.e(1:2, :)*1000);
    ylabel('Vm (mV)');
    box on;
    axis tight;
    xlim(x_lim);
    title('Excitatory cells');

    subplot(2, 1, 2, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.i(1:2, :)*1000);
    xlabel('Time (s)');
    ylabel('Vm (mV)');
    box on;
    axis tight;
    xlim(x_lim);
    title('Inhibitory cells');

    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', 'output/final_ramp_current_Vm_record.eps');
    