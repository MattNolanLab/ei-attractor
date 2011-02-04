% Print population firing raster plot

opt = parseOptions(options);
sheet_size = opt.sheet_size;

%close all;
figure;
hold on;

for x_i = 0:(sheet_size-1)
    for y_i = 0:(sheet_size-1)
        neuronID = y_i*sheet_size + x_i;
        neuronSpikes = spikeCell{neuronID+1};
        nID = neuronID + zeros(size(neuronSpikes));
        plot(neuronSpikes, nID, '.', 'MarkerSize', 6);
    end
end

hold off;
