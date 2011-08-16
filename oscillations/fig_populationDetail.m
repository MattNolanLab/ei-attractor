function fig_populationDetail(res, x_lim, T_i, fontSize)
    %firingRate_e = res.firingRate_e(T_i.start:T_i.end);
    opt = res.opt;
    
    hold off;

    subplot(5, 1, 1, 'FontSize', fontSize);
    spikeCellRasterPlot(res.spikeCell_e, '.');
    title('Pyramidal neurons');
    ylabel('Neuron number');
    xlim(x_lim);
    ylim([1 opt.Ne]);
    box on;

    subplot(5, 1, 2, 'FontSize', fontSize);
    plot(res.times, res.firingRate_e);
    ylabel('Firing rate (Hz)');
    xlim(x_lim);
    box on;

    subplot(5, 1, 3, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.e*1000);
    %title('Pyramidal neuron');
    ylabel('Vm (mV)');
    xlim(x_lim);
    box on;

    subplot(5, 1, 4, 'FontSize', fontSize);
    spikeCellRasterPlot(res.spikeCell_i, '.');
    title('Interneurons');
    ylabel('Neuron number');
    xlim(x_lim);
    ylim([1 opt.Ni]);
    box on;

    subplot(5, 1, 5, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.i*1000);
    %title('Interneuron');
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    xlim(x_lim);
    box on;
    

end