% Track blobs in the population firing patterns with different lambda_net
% values and make an estimation.

%clear all;
close all;
hold off;

firingTime = 4; %sec
delta_t = 1; % Time window (sec)

folder = 'results/static_wave/lambda_net/';
files = dir([folder '/job213*.mat'])

nFiles = numel(files);
fileStep = 1;

edgeThreshold = 5;

%firingPop = [];
%opts = [];

fileInd = 1:fileStep:nFiles; 

% Firstly process firing rates from all specified files
%  for it = fileInd
%      disp(['Processing file: ' files(it).name]);
%      [firingPop(:, :, (it-1)/fileStep+1) opts((it-1)/fileStep+1)] = getFiringPopFromFile([folder files(it).name], firingTime, delta_t);
%  end

sheet_size = opts(1).sheet_size;

estDist = [];

% Compute distances of all pairs of blobs from all the files
for it = 1:nFiles; 
    [r, c] = trackBlobs(firingPop(:, :, it));
    
    % Find blobs which are very near the edge, and remove them
    edge_i = find(r < edgeThreshold | r > sheet_size - edgeThreshold);
    r(edge_i) = [];
    c(edge_i) = [];
    edge_i = find(c < edgeThreshold | c > sheet_size - edgeThreshold);
    r(edge_i) = [];
    c(edge_i) = [];

    rep_r = repmat(r', 1, size(r, 2));
    rep_c = repmat(c', 1, size(c, 2));

    r = repmat(r, size(r, 2), 1);
    c = repmat(c, size(c, 2), 1);

    dist = sqrt((r-rep_r).^2 + (c-rep_c).^2);

    % Now sort the distances and take the second minimal value (the first
    % one is always zero) and average these distances
    sortDist = sort(dist);
    estDist(it) = mean(sortDist(2, :))
end

createFigureSpacing([opts.lambda_net], estDist);

printDir = '../../thesis/src/fig/';
printFile = 'blobSpacingPlot.eps';
print('-depsc', [printDir printFile]);