function [R POS] = slidingFiringRate(spikeRecord, win_len, noverlap)
    % Computes number of spikes per window length, sliding noverlap points
    % each iteration

    rec_end = size(spikeRecord, 2);
    
    pos_N = 1:noverlap:rec_end;
    R = zeros(numel(pos_N), size(spikeRecord, 1));
    POS = zeros(1, numel(pos_N));
    
    it = 1;
    for pos = pos_N;
      
        R(it, :) = full(sum(spikeRecord(:, pos:min(pos+win_len-1, rec_end))'));
        POS(it) = pos;
        
        it = it + 1;
    end
    
    R = R / win_len;
end