function [spikeHist] = createSpikeHistCell(neuronIDs, spikeTimes, dt, startTime, endTime)
% ------------------------------------------------------------------------
% Preprocesses spiking data (as cell array) into spikeHist value

N = numel(neuronIDs); % Number of neurons
edges = linspace(startTime, endTime, (endTime-startTime)/dt + 1);
spikeHist = zeros(N, numel(edges));

for it = 1:N
        neuronSpikes = spikeTimes{neuronIDs(it)};
        e_size = size(neuronSpikes, 1);
        if (e_size ~= 0)
            spikeHist(it, :) = histc(neuronSpikes, edges);
        end
    end
end
