% Specify a neuron ID and its neighborhood radius, to create a figure of a
% synchronization index (defined by a MvR distance)
%
% Just now, we assume that population response doesn't move

path('../include/', path);
close all;

opt = parseOptions(options);

% MvR distance parameters
tc = 0.025;
expo_dur = 10*tc; %sec
dt = 0.001;
startT = 0;
endT = 200;
ISI = 0.1; % average isi
rad = 26; % Neighborhood radius - in neural units

syncWinT = 2; % sec
syncWin = syncWinT/dt;

spikeCntThreshold = syncWinT/ISI/2;


% The value of venterListID has been determined from freq. spectra of Vm recordings of the
% neurons see script figureMembraneVFreq.m

nID = 4728; % job40000 - estimated from raster plot
%nID = 7170; % job40003
center_r = fix(nID/sheet_size);
center_c = mod(nID, sheet_size);

[neigh_r, neigh_c] = getBlobNeighborhood(center_r, center_c, rad, sheet_size);
neigh_NID = neigh_r*sheet_size + neigh_c + 1;


% Remove neurons which have spike count < threshold
it_del = [];
for it = 1:numel(neigh_NID)
    spikes = spikeCell{neigh_NID(it)};
    if numel(find(spikes >= startT & spikes <= endT)) < (endT-startT)/ISI/2;
        it_del = [it_del it];
    end
end
neigh_NID(it_del) = [];

% Check for the neighborhood neurons
figure(1);
plot(mod(neigh_NID-1, sheet_size), fix((neigh_NID-1)./sheet_size), '.');
display('Displaying the neighborhood plot. Press any key to continue.');
pause;


% Create the spike response function for the neighborhood neurons and
% convolve it with the exponential
spikeHist = createSpikeHistCell(neigh_NID, spikeCell, dt, startT, endT);
spikeHist = spikeHist;
sz = size(spikeHist);
N = sz(1)

e_t = 0:dt:expo_dur;
expo = exp(-(e_t)/tc);

response = [];
for it = 1:N
    response(it, :) = conv(spikeHist(it, :), expo);
end

pairTimedSync = [];
it = 1;
timeSteps = 1:syncWin:sz(2);

for it1 = 1:N
    it1
    for it2 = it1:N
        fmg = (response(it1, :) - response(it2, :)).^2;
        real_t = startT;
        for t = 1:numel(timeSteps)
            n1_spikes = spikeCell{neigh_NID(it1)};
            n2_spikes = spikeCell{neigh_NID(it2)};
            n1_scnt = numel(find(n1_spikes >= real_t & n1_spikes <= real_t+syncWinT));
            n2_scnt = numel(find(n2_spikes >= real_t & n2_spikes <= real_t+syncWinT));
            if (n1_scnt < spikeCntThreshold || n2_scnt < spikeCntThreshold)
                pairTimedSync(it, t) = nan;
            else            
                spikeCnt_coeff = (n1_scnt + n2_scnt)/2;

                pairTimedSync(it, t) = trapz(fmg(t:t+syncWin)) / spikeCnt_coeff;

                
            end
            real_t = real_t + syncWinT;
        end        
        it = it + 1;

    end
end

pairTimedSync = pairTimedSync/tc*dt;
meanTimedSync = [];
for t = 1:numel(timeSteps)
    validIDs = find(~isnan(pairTimedSync(:, t)) & pairTimedSync(:, t) ~= 0);
    meanTimedSync(t) = mean(pairTimedSync(validIDs, t));
    stdTimedSync(t) = std(pairTimedSync(validIDs, t));
end

figure(2);
plot((timeSteps-1)/syncWin*syncWinT, meanTimedSync);
ylim([0 1]);
xlabel('Time (s)');
ylabel('Mean D^2_{MvR} of pairs (norm. by (N+M)/2');
title(['Time dependent sync. measure. Blob and around. t_c = ' num2str(tc) ' s.']);


% -------------------------------------------------------------------------
% Raster plot of neighborhood neurons
% -------------------------------------------------------------------------
figure(3);
rasterPlot(spikeCell, neigh_NID, true);
title(['Raster plot of neighborhood neurons of neuron (r: ' ...
    int2str(center_r) ', c: ' int2str(center_c) ')']);
xlim([startT endT]);



