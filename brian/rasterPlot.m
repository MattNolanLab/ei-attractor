% Print population firing raster plot

opt = parseOptions(options);
sheet_size = opt.sheet_size;

close all;
figure;
hold on;

for x_i = 0:(sheet_size-1)
    for y_i = 0:(sheet_size-1)
        neuronID = y_i*sheet_size + x_i;
        neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
        nID = neuronID + zeros(size(neuronSpikes));
        plot(neuronSpikes, nID, '.');
    end
end

hold off;
