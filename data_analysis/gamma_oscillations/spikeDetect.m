function [spike_ids] = spikeDetect(sig, times, threshold)
    % SPIKEDETECT simple spike detection algorithm, which takes signal
    % 'sig', with sampling rate 'dt', finds local maxima, which are then
    % compared to threshold. The extrema which exceed threshold are
    % considered spikes

    d2sig_eps = 0;
    
    dsig = diff(sig);
    d2sig = diff(dsig);
    
    dsig_zero = sign(dsig(2:end) .* dsig(1:end-1)) == -1;
    peak_i = find(dsig_zero & d2sig < d2sig_eps) + 1;
   
    % Remove local extrema which are less than threshold
    spike_ids = peak_i(find(sig(peak_i) > threshold));
    spike_times = times(spike_ids);
end