function [firingRate] = computeFiringRate(spikeTimes, startTime, endTime, dt, delta_t)
    % spikeTimes - an array of spike times to estimate the rate for
    % endTime - end time of the data to compute firing rate for
    % dt - time resolution
    % delta_t - sliding window resolution

    nel = numel(spikeTimes);
    
    times = startTime:dt:endTime;
    firingRate = zeros(1, numel(times));
    
    if nel == 0
        return;
    end
    
    maxTime = spikeTimes(nel);
    
    
    t_i = 1;
    for t = times
        numSpikes = numel(find((spikeTimes > t - delta_t/2) & (spikeTimes < t + delta_t/2)));
        firingRate(t_i) = numSpikes/delta_t;
        t_i = t_i + 1;
    end
end