function fig_populationDetail(res, x_lim, T_i, fontSize, N)
    opt = res.opt;
    
    hold off;

%     subplot(6, 1, 1, 'FontSize', fontSize);
%     spikeCellRasterPlot(res.spikeCell_e, '.', N);
%     title('Pyramidal neurons');
%     %ylabel('Neuron no.');
%     xlim(x_lim);
%     ylim([1 N]);
%     box on;

    subplot(4, 1, 1, 'FontSize', fontSize);
    plot(res.times, res.firingRate_e);
    ylabel('Firing rate (Hz)');
    xlim(x_lim);
    box on;

    subplot(4, 1, 2, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.e(1, :)*1000);
    ylabel('Vm (mV)');
    box on;
    axis tight;
    xlim(x_lim);

%     subplot(6, 1, 4, 'FontSize', fontSize);
%     spikeCellRasterPlot(res.spikeCell_i, '.', N);
%     title('Interneurons');
%     %ylabel('Neuron no.');
%     xlim(x_lim);
%     ylim([1 N]);
%     box on;
%     
    subplot(4, 1, 3, 'FontSize', fontSize);
    plot(res.times, res.firingRate_i);
    ylabel('Firing rate (Hz)');
    xlim(x_lim);
    box on;


    subplot(4, 1, 4, 'FontSize', fontSize);
    plot(res.Vmon.t, res.Vmon.i(1, :)*1000);
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    box on;
    axis tight;
    xlim(x_lim);
    

end