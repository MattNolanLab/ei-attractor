% Create cross correlograms of pairs of neurons within one blob
% Assuming that blobs are stable


% createSpikeHist parameters
startTime = 10;
endTime = 200;
dt_track = 0.1;
delta_t = 0.25; % Should be this value.

preprocess = false;

if preprocess
    disp 'Creating spike histogram';
    createSpikeHist;  % script - computes spikeHist for current file        
    createSpikeCell;  % script - transforms spike times to a cell array
end


% Get positions of blobs and neuron numbers around the centers
blobEstTime = 10; % sec
firingPop = getFiringPop(spikeHist, blobEstTime, dt_track, delta_t);
[r, c] = trackBlobs(firingPop);
blobCenterNID = fix(r)*sheet_size + fix(c);


% Draw cross correlation histograms of spike times (1st neuron with others)
for it = 1:numel(blobCenterNID)
    figure;
    neuronNum = blobCenterNID(it)
    crossCorrelation(spikeCell{blobCenterNID(3)}, spikeCell{neuronNum}, ...
        0.0075, 0.25, 200);
    
    fileName = ['cross_corr_n' int2str(blobCenterNID(1)) '_n' int2str(neuronNum)];
    
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', ['results/fig/spike_statistics/' fileName]);
end