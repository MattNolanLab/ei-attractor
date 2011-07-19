function [spikeCell] = spikeRecordToSpikeCell(spikeRecord, times)
    % Process spikeRecord entity (binary array of spikes for each
    N = size(spikeRecord, 1);
    
    for it = 1:N
        spikeCell{it} = times(spikeRecord(it, :) == 1);
    end
end