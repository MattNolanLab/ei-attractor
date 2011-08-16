function [firingRate] = getFiringRate(spikeHist, dt, win_len)
    % Compute population firing rate at each time interval of the
    % population vector

    N = size(spikeHist, 1);  % Number of neurons in population
    s_hist = sum(spikeHist);
    nbins_half = fix(win_len/dt/2);
    
    firingRate = zeros(1, size(spikeHist, 2));
    

    end_i = size(spikeHist, 2);
  
    for t_i = 1:end_i
        bins_start = max(t_i - nbins_half, 1);
        bins_end   = min(t_i + nbins_half, end_i);

        firingRate(t_i) = sum(s_hist(bins_start:bins_end))/win_len/N;
    end
end