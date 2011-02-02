% Print population firing raster plot in 3d

opt = parseOptions(options);
sheet_size = opt.sheet_size;
timeLimit = 2; % sec

%close all;
figure;
hold on;

for x_i = 0:(sheet_size-1)
    for y_i = 0:(sheet_size-1)
        neuronID = y_i*sheet_size + x_i;
        neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
        spikeData = neuronSpikes(find(neuronSpikes <= timeLimit));
        X = 50 + zeros(size(spikeData));
        Y = y_i + zeros(size(spikeData));
        plot3(spikeData, X, Y, '.', 'MarkerSize', 6);
    end
end

hold off;
