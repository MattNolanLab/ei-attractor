function [firingRate] = getNeuronFiringRate(neuronSpikeHist, startTime, endTime, dt_track, delta_t)
% Compute firing rate of a single neuron, assuming neuronSpikeHist starts
% at 0

end_i = size(neuronSpikeHist, 2);

it = 1;
for t = startTime:dt_track:endTime
    hist_i = fix(t/dt_track + 1);
    
    nbins_pre = fix(min(t/dt_track, delta_t/dt_track/2));
    nbins_post = fix(min((end_i-hist_i), delta_t/dt_track/2));

    % take the window at position specified by t and sum up
    % number of spikes
    bins_i = hist_i-nbins_pre:hist_i+nbins_post;
    bins = neuronSpikeHist(bins_i);

    firingRate(it) = sum(bins)/(dt_track*numel(bins));
    it = it + 1;
end
    

end