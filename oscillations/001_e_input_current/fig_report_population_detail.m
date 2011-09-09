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

N = 25;

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


    % Plot detailed responses of a randomly selected trial
    figure('Position', [800 0 850 1050]);

    trial_it = 1;  

        res = results(par_it, trial_it);
        opt = res.opt;



        subplot(5, 2, 1, 'FontSize', fontSize);
        spikeCellRasterPlot(res.spikeCell_e, '.', N);
        title('Principal cells');
        %ylabel('Neuron no.');
        xlim(x_lim);
        ylim([1 N]);
        box on;
        set(gca, 'Xtick', []);
        set(gca, 'Ytick', []);
        set(gca, 'Position', [0.13 0.754 0.335 0.154]);


        subplot(5, 2, 3, 'FontSize', fontSize);
        plot(res.times, res.firingRate_e);
        ylabel('Firing rate (Hz)');
        xlim(x_lim);
        box on;
        set(gca, 'Xtick', []);
        set(gca, 'Position', [0.13 0.634748272458045 0.334659090909091 0.097336845583677]);

        subplot(5, 2, 5, 'FontSize', fontSize);
        plot(res.Vmon.t, res.Vmon.e*1000);
        ylabel('Vm (mV)');
        box on;
        axis tight;
        xlim(x_lim);
        set(gca, 'Xtick', []);
        set(gca, 'Position',[0.13 0.50542941757157 0.334659090909091 0.101678183613031]);

        subplot(5, 2, 7, 'FontSize', fontSize);
        spikeCellRasterPlot(res.spikeCell_i, '.', N);
        title('Interneurons');
        %ylabel('Neuron no.');
        xlim(x_lim);
        ylim([1 N]);
        set(gca, 'Xtick', []);
        set(gca, 'Ytick', []);

        box on;
        set(gca, 'Position',[0.13 0.281 0.334659090909091 0.143139190523197]);

        subplot(5, 2, 9, 'FontSize', fontSize);
        plot(res.Vmon.t, res.Vmon.i*1000);
        ylabel('Vm (mV)');
        xlabel('Time (s)');
        box on;
        axis tight;
        xlim(x_lim);
        set(gca, 'Position',[0.13 0.156303499422758 0.334659090909091 0.102]);
    

    trial_it = 17;  

        res = results(par_it, trial_it);
        opt = res.opt;


        subplot(5, 2, 2, 'FontSize', fontSize);
        spikeCellRasterPlot(res.spikeCell_e, '.', N);
        title('Principal cells');
        %ylabel('Neuron no.');
        xlim(x_lim);
        ylim([1 N]);
        box on;
        set(gca, 'Xtick', []);
        set(gca, 'Ytick', []);
        set(gca, 'Position', [0.55 0.754 0.335 0.154]);


        subplot(5, 2, 4, 'FontSize', fontSize);
        plot(res.times, res.firingRate_e);
        %ylabel('Firing rate (Hz)');
        xlim(x_lim);
        box on;
        set(gca, 'Xtick', []);
        set(gca, 'Position', [0.55 0.634748272458045 0.334659090909091 0.097336845583677]);

        subplot(5, 2, 6, 'FontSize', fontSize);
        plot(res.Vmon.t, res.Vmon.e*1000);
        %ylabel('Vm (mV)');
        box on;
        axis tight;
        xlim(x_lim);
        set(gca, 'Xtick', []);
        set(gca, 'Position',[0.55 0.50542941757157 0.334659090909091 0.101678183613031]);

        subplot(5, 2, 8, 'FontSize', fontSize);
        spikeCellRasterPlot(res.spikeCell_i, '.', N);
        title('Interneurons');
        %ylabel('Neuron no.');
        xlim(x_lim);
        ylim([1 N]);
        set(gca, 'Xtick', []);
        set(gca, 'Ytick', []);

        box on;
        set(gca, 'Position',[0.55 0.281342546890424 0.334659090909091 0.143139190523197]);

        subplot(5, 2, 10, 'FontSize', fontSize);
        plot(res.Vmon.t, res.Vmon.i*1000);
        %ylabel('Vm (mV)');
        xlabel('Time (s)');
        box on;
        axis tight;
        xlim(x_lim);
        set(gca, 'Position',[0.55 0.156303499422758 0.334659090909091 0.102]);
        
        
        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('output/report_spiking_details.eps', outputDir, outputFileRef, opt.Ie*1000, trial_it));
end