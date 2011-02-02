function [D] = MvR_distAll(spikeCell, tc, dt)
% Compute MvR distance of all spike train pairs in the cell array set
% with time constant t_c
% For more informations cf. MCW van Rossum, A Novel spike distance, Neural
% computation, 2001

% Timing properties for response functions
%dt = 0.1;
%startTime = 0;
%endTime = max(spikeCell{nid1}(end), spikeCell{nid2}(end));

% Exponential properties
expo_dur = 2*tc; %sec

% create spike response functions for both neurons
%responseFunc = createSpikeHistCell([nid1 nid2], spikeCell, dt, startTime, endTime);

% Create exponential
e_t = 0:dt:expo_dur;
expo = exp(-(e_t)/tc);

res1 = conv(train1, expo);
res2 = conv(train2, expo);

D = 1/tc * trapz((res(1, :) - res(2, :)).^2)*dt

end