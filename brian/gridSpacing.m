% Track blobs in the population firing patterns with different lambda_net
% values and make an estimation.

clear all;

firingTime = 4; %sec
delta_t = 1; % Time window (sec)

folder = 'results/static_wave/lambda_net/';
files = dir([folder '/job21000*.mat'])

nFiles = numel(files);

%firingPop = [];
%opts = [];

% Firstly process firing rates from all specified files
for it = 1:nFiles
    [firingPop(:, :, it) opts(it)] = getFiringPopFromFile([folder files(it).name], firingTime, delta_t);
end

firingPop