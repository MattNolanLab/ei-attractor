% ------------------------------------------------------------------------
% Preprocesses spikeMonitor_times_nXXX into isi - array of isi for all
% neurons
% 
% BEWARE - this is a script - side effects
% ------------------------------------------------------------------------

% sheet_size must be defined

sheet_size = double(sheet_size);

for x_i = 0:(sheet_size-1)
    for y_i = 0:(sheet_size-1)
        neuronID = y_i*sheet_size + x_i;            
        neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
        e_size = size(neuronSpikes, 1);
        if (e_size == 0)
            isi{neuronID+1} = []; 
        else
            isi{neuronID+1} = neuronSpikes(2:end) - neuronSpikes(1:end-1);
        end
    end
end