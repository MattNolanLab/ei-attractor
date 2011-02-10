function rasterPlot(spikeCell, neuronIDs, flat)
% RASTERPLOR Produce a figure containing raster plot of specified neurons
% (or all neurons is neuronIDs == [])
%
% Use the flat parameter if you want to display Y axis as 1:N instead of
% actual neuron numbers

%opt = parseOptions(options);
%sheet_size = opt.sheet_size;

%close all;
figure;
hold on;

if neuronIDs == -1
    neuronIDs = 1:numel(spikeCell);
end

if flat
    for it = 1:numel(neuronIDs)
        neuronSpikes = spikeCell{neuronIDs(it)};
        nID = it + zeros(size(neuronSpikes));
        plot(neuronSpikes, nID, '.', 'MarkerSize', 6);
    end
else
    for it = neuronIDs
        neuronSpikes = spikeCell{it};
        nID = it + zeros(size(neuronSpikes));
        plot(neuronSpikes, nID, '.', 'MarkerSize', 6);
    end
end

hold off;

xlabel('Time (s)');
ylabel('Neuron no.');

end