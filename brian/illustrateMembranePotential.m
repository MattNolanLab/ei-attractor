% Illustrate how membrane potential changes with the rat moving

%load results/membrane_potential/job23101_2010-08-05T09-54-38_output.mat

dt_rat = 0.02;  % sec
dt_track = 0.1; % sec
delta_t = 0.5;  % sec

opt = parseOptions(options);
sheet_size = opt.sheet_size;

endTime = 1200; %sec

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
    
