% Specify a neuron ID and its neighborhood radius, to create a figure of
%  1) Cross correlation histograms of relevant neurons (which fire
%     sufficiently strong)
%  2) MvR overall distance between these neurons
%
% To detect synchronisations between neurons projecting onto the specified
% neuron. Assuming that we have a multiple bump solution here.
%
% In this specific case, job 23500 has been used, which contains membrane
% potential recording of a full row of neurons. For setups without membrane
% potential recordings, it is necessary to specify neuron ID directly

path('../include/', path);
%load ../../mult_bump_spiking_net/results/job23500_2010-08-07T19-37-39_output.mat;

close all;

opt = parseOptions(options);

% MvR distance parameters
tc = 0.025;
dt = 0.001;
startT = 0;
endT = 1;
spikeCntThreshold = 5;

% Correlogram start and end time setting
corrStartT = 0; % sec
corrEndT = corrStartT + 10; % sec

% Correlogram parameters
corr_dt = 0.02; %sec
corr_range = 0.5;
corr_T = 10;

% Neighborhood radius - in neural units
rad = 13;


% The value of venterListID has been determined from freq. spectra of Vm recordings of the
% neurons see script figureMembraneVFreq.m
%centerListID = 51;
%nID = double(SNList(centerListID));
nID = 5330; % job40000 - estimated from raster plot
%nID = 7170; % job40003
center_r = fix(nID/sheet_size);
center_c = mod(nID, sheet_size);

[neigh_r, neigh_c] = getBlobNeighborhood(center_r, center_c, rad, sheet_size);
neigh_NID = neigh_r*sheet_size + neigh_c + 1;


% Remove neurons which have spike count < threshold
it_del = [];
for it = 1:numel(neigh_NID)
    spikes = spikeCell{neigh_NID(it)};
    if numel(find(spikes >= startT & spikes <= endT)) < spikeCntThreshold
        it_del = [it_del it];
    end
end
neigh_NID(it_del) = [];

% -------------------------------------------------------------------------
% Raster plot of neighborhood neurons
% -------------------------------------------------------------------------
rasterPlot(spikeCell, neigh_NID, true);
title(['Raster plot of neighborhood neurons of neuron (r: ' ...
    int2str(center_r) ', c: ' int2str(center_c) ')']);
xlim([startT endT]);


% -------------------------------------------------------------------------
% Draw cross correlation histograms of spike times (1st neuron with others)
% -------------------------------------------------------------------------
% plotPx = 50;
% fontSize = 8;
% nplots = numel(neigh_NID)
% figure('Position', [0 0 nplots*plotPx*1.5 nplots*plotPx], 'Visible', 'off');
% 
% cnt = 1;
% for nid1 = neigh_NID
%     nid1
%     for nid2 = neigh_NID
%         subplot(nplots, nplots, cnt, 'FontSize', fontSize);
%         %neuronNum = neigh_NID(it)
% 
%         % Setup spikes from the specified interval
%         spikes1 = spikeCell{nid1};
%         spikes1 = spikes1(find(spikes1 >= corrStartT & spikes1 <= corrEndT));
%         spikes2 = spikeCell{nid2};
%         spikes2 = spikes2(find(spikes2 >= corrStartT & spikes2 <= corrEndT));
% 
%         crossCorrelation(spikes1, spikes2, corr_dt, corr_range, corr_T);
%         title(['n' int2str(nid2) ', n' int2str(nid2)]);
%         
%         cnt = cnt + 1;
%     end
% end
% 
% set(gcf,'PaperPositionMode','auto');
% fileName = 'silentNeuronCrossCorrHist.eps';
% print('-depsc2', ['fig/' fileName]);


% -------------------------------------------------------------------------
% Draw cross correlation histograms of spike times (1st neuron with others)
% -------------------------------------------------------------------------
D = MvR_DistAll(neigh_NID, spikeCell, tc, dt, startT, endT, spikeCntThreshold);
D(find(D == 0)) = nan;

figure;
pcolor(D'); shading flat; colorbar
xlabel('Neuron no.');
ylabel('Neuron no.');
title('MvR distance, normalized by (N+M)/2');

D_nnan = D(find(~isnan(D)));
avg_D = mean(D_nnan)
sum_D = sum(D_nnan)

hist_precision = 40;
figure;
hist(D_nnan, hist_precision);
xlabel('D_{MvR}');
ylabel('Count (pairs)');
title('Histogram of distances of pairs of neurons');
