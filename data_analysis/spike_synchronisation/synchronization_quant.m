path('..', path);

neuronIDs = 1000:2000;
tc = 0.050;
dt = 0.001;
startT = 0;
endT = 200;
spikeCntThreshold = 10;

% startTime = 0;
% endTime = 200;
% dt_track = 0.1;
% delta_t = 0.25; % Should be this value.
% 
% preprocess = true;
% 
% if preprocess
%     disp 'Creating spike histogram';
%     spikeHist = createSpikeHistCell(1:sheet_size^2, spikeCell, dt_track, startTime, endTime);  % script - computes spikeHist for current file        
% end
% 
% % Get positions of blobs and neuron numbers around the centers
% blobEstTime = startT; % sec
% firingPop = getFiringPop(spikeHist, blobEstTime, dt_track, delta_t);
% [r, c] = trackBlobs(firingPop);
% r = fix(r)
% c = fix(c)
% %blobCenterNID = r*sheet_size + c;
% 

centerID = 2400;
r = fix(centerID/sheet_size);
c = mod(centerID, sheet_size);

for blobID = 1; 
    blobCent_r = r(blobID);
    blobCent_c = c(blobID);
    blobCenterNID = blobCent_r*sheet_size + blobCent_c;
    rad = 10;
    
    [neigh_r, neigh_c] = getBlobNeighborhood(blobCent_r, blobCent_c, rad);
    neigh_NID = neigh_r*sheet_size + neigh_c
end


D = MvR_DistAll(neigh_NID, spikeCell, tc, dt, startT, endT, spikeCntThreshold);
D(find(D == 0)) = nan;

pcolor(D'); shading flat; colorbar

D_nnan = D(find(~isnan(D)));
avg_D = mean(D_nnan)
sum_D = sum(D_nnan)