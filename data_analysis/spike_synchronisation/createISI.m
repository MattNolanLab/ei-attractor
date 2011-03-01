function [isi] = createISI(spikeCell)
% CREATEISI preprocess spike data into ISI values

for neuronID = 1:numel(spikeCell)
    neuronSpikes = spikeCell{neuronID};
    e_size = size(neuronSpikes, 1);
    if (e_size == 0)
        isi{neuronID} = []; 
    else
        isi{neuronID} = neuronSpikes(2:end) - neuronSpikes(1:end-1);
    end
end

end