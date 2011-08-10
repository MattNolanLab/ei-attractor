function plotTrialDetails(results, par_it, trial_it, plot_options)
    % Raster plot and Vm of population and recorded neurons
   
    res = results(par_it, trial_it);
    opt = res.opt;
    
    x_lim = plot_options.x_lim;
    fontSize = plot_options.fontSize;
    
    subplot(5, 1, 1, 'FontSize', fontSize);
    spikeCellRasterPlot(res.spikeCell_e, '.');
    title('Pyramidal neurons');
    ylabel('Neuron number');
    xlim(x_lim);
    ylim([1 opt.Ne]);

    subplot(5, 1, 2, 'FontSize', fontSize);
    plot(res.times, res.firingRate_e);
    ylabel('Firing rate (Hz)');
    xlim(x_lim);


    subplot(5, 1, 3, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.e*1000);
    %title('Pyramidal neuron');
    ylabel('Vm (mV)');
    xlim(x_lim);


    subplot(5, 1, 4, 'FontSize', fontSize);
    spikeCellRasterPlot(res.spikeCell_i, '.');
    title('Interneurons');
    ylabel('Neuron number');
    xlim(x_lim);
    ylim([1 opt.Ni]);


    subplot(5, 1, 5, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.i*1000);
    %title('Interneuron');
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    xlim(x_lim);

    hold off;
end