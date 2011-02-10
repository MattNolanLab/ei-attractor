function [firingPop] = getFiringPop(spikeHist, t, dt_track, delta_t)
    % Compute population firing rate at time t
    sheet_size = sqrt(size(spikeHist, 1));
    
    end_i = size(spikeHist, 2);
    hist_i = fix(t/dt_track + 1);    
    nbins_pre = fix(min(t/dt_track, delta_t/dt_track/2));
    nbins_post = fix(min((end_i-hist_i), delta_t/dt_track/2));
    
    for x_i = 0:(sheet_size-1)
        for y_i = 0:(sheet_size-1)
            neuronID = y_i*sheet_size + x_i;

            % take the window at position specified by t and sum up
            % number of spikes
            bins = spikeHist(neuronID+1, hist_i-nbins_pre:hist_i+nbins_post);

            firingPop(x_i+1, y_i+1) = sum(bins)/delta_t;
        end
    end

    firingPop = firingPop';    
end