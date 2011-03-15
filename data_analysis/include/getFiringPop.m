function [firingPop] = getFiringPop(spikeHist, t, dt_track, delta_t)
    % Compute population firing rate at time t, do not assume any
    % topographical organisation of neurons.
    
    %sheet_size = sqrt(size(spikeHist, 1));
    
    end_i = size(spikeHist, 2);
    hist_i = fix(t/dt_track + 1);    
    nbins_pre = fix(min(t/dt_track, delta_t/dt_track/2));
    nbins_post = fix(min((end_i-hist_i), delta_t/dt_track/2));
    
    N = size(spikeHist, 1);
    for nID = 1:N
        % take the window at position specified by t and sum up
        % number of spikes
        bins = spikeHist(nID, hist_i-nbins_pre:hist_i+nbins_post);
        firingPop(nID) = sum(bins)/delta_t;
    end
end