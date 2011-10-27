function [N, E POS] = slidingFiringRate(spikeRecord, win_len, noverlap, nBins)

    rec_end = size(spikeRecord, 2);
    
    pos_N = 1:noverlap:rec_end;
    N = zeros(numel(pos_N), nBins+1);
    POS = zeros(1, numel(pos_N));
    
    edges = linspace(0, 100, nBins+1);
    
    it = 1;
    for pos = pos_N;
      
        N(it, :) = histc(full(sum(spikeRecord(:, pos:min(pos+win_len-1, rec_end))')), edges);
        POS(it) = pos;
        
        it = it + 1;
    end
    
    E = edges;
end