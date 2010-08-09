% ------------------------------------------------------------------------
% Preprocesses spikeMonitor_times_nXXX into spikeHist value
% 
% BEWARE - this is a script - side effects
% ------------------------------------------------------------------------

% sheet_size and dt_track must be defined

% Preprocess spiking data: firing time --> histogram of firing, for
% each neuron
edges = linspace(0, endTime, endTime/dt_track);
spikeHist = zeros(sheet_size^2, numel(edges));
for x_i = 0:(sheet_size-1)
    for y_i = 0:(sheet_size-1)
        neuronID = y_i*sheet_size + x_i;            
        neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
        e_size = size(neuronSpikes, 1);
        if (e_size == 0)
            spikeHist(neuronID+1, :) = 0;
        else
            spikeHist(neuronID+1, :) = histc(neuronSpikes, edges);
        end
    end
end
