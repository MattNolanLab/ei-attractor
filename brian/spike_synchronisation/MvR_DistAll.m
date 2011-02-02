function [D] = MvR_DistAll(neuronIDs, spikeCell, tc, dt, startTime, endTime, spikeNumThreshold)
% Compute MvR distance of all spike train pairs in the cell array set
% with time constant t_c
% For more informations cf. MCW van Rossum, A Novel spike distance, Neural
% computation, 2001

% Exponential
expo_dur = 10*tc; %sec
e_t = 0:dt:expo_dur;
expo = exp(-(e_t)/tc);


N = numel(neuronIDs);
D = nan * ones(N);

for it1 = 1:N
    it1
    if (numel(spikeCell{neuronIDs(it1)}) < spikeNumThreshold)
        continue;
    end
    response1 = createSpikeHistCell(neuronIDs(it1), spikeCell, dt, startTime, endTime);
    
    for it2 = it1:N
        if (numel(spikeCell{neuronIDs(it2)}) < spikeNumThreshold)
            continue;
        end

        response2 = createSpikeHistCell(neuronIDs(it2), spikeCell, dt, startTime, endTime);
        D(it1, it2) = 1/tc * trapz((response1 - response2).^2)*dt;
    end
end