% Histograms of input current to populations and their firing rates, at specified
% input spread value
close all;
clearvars -except results;

fontSize = 14;


%load e_input_current_output_19-Jul-2011;
outputDir = 'output';

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;

spread_all = 0:0.25:5;
find_eps = 1e-9;
D = results(1, 1).opt.D;
nBins = 50;

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
    
    res = results(par_it, 2);

    % Excitatory input histogram
    figure('Position', [800 500 800 400], 'Visible', 'off');
    subplot(1, 2, 1, 'FontSize', fontSize);
    hist(res.opt.Ie*1000, nBins);
    axis tight;
    xlabel('Excitatory input drive (mV)');
    ylabel('No. of neurons');
    xlim([0 res.opt.Ie_max*1000]);
    
    subplot(1, 2, 2, 'FontSize', fontSize);
    hist(sum(res.spikeRecord_e')/res.opt.T, nBins);
    axis tight;
    xlabel('Excitatory firing rate (Hz)');

    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('%s/input_spread_Ie_hist_spread_%.3f.eps', outputDir, res.opt.input_spread/res.opt.D));
    

    
    % Inhibitory input histogram
    figure('Position', [800 500 800 400], 'Visible', 'off');
    subplot(1, 2, 1, 'FontSize', fontSize);
    hist(res.opt.Ii*1000, nBins);
    axis tight;
    xlabel('Inhibitory input drive (mV)');
    ylabel('No. of neurons');
    xlim([0 res.opt.Ii_max*1000]);

    subplot(1, 2, 2, 'FontSize', fontSize);
    hist(sum(res.spikeRecord_i')/res.opt.T, nBins);
    axis tight;
    xlabel('Inhibitory firing rate (Hz)');
    
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('%s/input_spread_Ii_hist_spread_%.3f.eps', outputDir, res.opt.input_spread/res.opt.D));
    

    
end