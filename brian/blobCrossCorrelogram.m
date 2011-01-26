% Create cross correlograms of pairs of neurons within one blob
% Assuming that blobs are stable


% createSpikeHist parameters
startTime = 10;
endTime = 200;
dt_track = 0.1;
delta_t = 0.25; % Should be this value.

preprocess = true;

if preprocess
    disp 'Creating spike histogram';
    createSpikeHist;  % script - computes spikeHist for current file        
    createSpikeCell;  % script - transforms spike times to a cell array
end


% Get positions of blobs and neuron numbers around the centers
blobEstTime = 20; % sec
firingPop = getFiringPop(spikeHist, blobEstTime, dt_track, delta_t);
[r, c] = trackBlobs(firingPop);
r = fix(r)
c = fix(c)
%blobCenterNID = r*sheet_size + c;

for blobID = 1:numel(r); 
    blobCent_r = r(blobID);
    blobCent_c = c(blobID);
    blobCenterNID = blobCent_r*sheet_size + blobCent_c;
    rad = 2;
    
    [neigh_r, neigh_c] = getBlobNeighborhood(blobCent_r, blobCent_c, rad);
    neigh_NID = neigh_r*sheet_size + neigh_c
    
    % Draw cross correlation histograms of spike times (1st neuron with others)
    figure('Position', [0 0 1680 1050], 'Visible', 'off');
    nplots = sqrt(numel(neigh_NID));
    for it = 1:numel(neigh_NID)
        subplot(nplots, nplots, it);
        neuronNum = neigh_NID(it)
        if (neuronNum < 1 || neuronNum > sheet_size^2)
            continue;
        end
        crossCorrelation(spikeCell{blobCenterNID}, spikeCell{neuronNum}, ...
            0.0075, 0.25, 200);
        title(['n' int2str(blobCenterNID) ', n' ...
            int2str(neuronNum)]);
    end
    
    fileName = ['cross_corr_n' int2str(blobCenterNID(1)) '_neighborhood'];
    
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', ['results/fig/spike_statistics/within_blob_corr/unstable_net/' fileName]);
end