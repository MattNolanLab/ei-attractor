function [spike_id] = findSpikes(Vm, threshold)
% FINDSPIKES simple spike detection algorithm. Finds all maxima of the
% signal sig which have Vm > threshold
    
    [xmax imax xmin imin] = extrema(Vm);
    
    false_spikes = find(xmax <= threshold);
    imax(false_spikes) = [];
    spike_id = imax;
end