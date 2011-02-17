% Specify a neuron ID and its neighborhood radius, to create a figure of a
% synchronization index (defined by a MvR distance)
%
% Just now, we assume that population response doesn't move

path('../include/', path);
close all;

opt = parseOptions(options);

% MvR distance parameters
tc = 0.025;
dt = 0.001;
startT = 0;
endT = 200;
ISI = 0.1; % average isi
rad = 20; % Neighborhood radius - in neural units

syncWinT = 2; % sec


% The value of venterListID has been determined from freq. spectra of Vm recordings of the
% neurons see script figureMembraneVFreq.m

nID = 4728; % job40000 - estimated from raster plot
%nID = 7170; % job40003


[pairTimedSync, meanTimedSync, stdTimedSync, neigh_NID] = ...
    timeDependentSync(spikeCell, tc, dt, startT, endT, ISI, rad, nID, syncWinT);

sz = size(pairTimedSync);


% Check for the neighborhood neurons
figure(1);
plot(mod(neigh_NID-1, sheet_size), fix((neigh_NID-1)./sheet_size), '.');
%display('Displaying the neighborhood plot. Press any key to continue.');
%pause;


figure(2);
set(gcf, 'Position', [883 528 600 500]);
% -------------------------------------------------------------------------
% Plot time dependent histogram map of sync. coefficient
% -------------------------------------------------------------------------
subplot(2, 1, 1);
[timedHist, binCenters] = hist(pairTimedSync, 40);
pcolor(startT:syncWinT:endT, binCenters, timedHist);
xlim([startT endT]);
shading flat;
%colorbar;
xlabel('Time (s)');
ylabel('D^2_{MvR}');
title('Time dependent proportion of pairs of spikes given a sync. coeff');


% -------------------------------------------------------------------------
% Plot time dependent average sync. coefficient
% -------------------------------------------------------------------------

subplot(2, 1, 2);
plot(startT:syncWinT:endT, meanTimedSync);
xlim([startT endT]);
ylim([0 1]);
xlabel('Time (s)');
ylabel('Mean D^2_{MvR}');
title(['Time dependent sync. measure. Blob and around. t_c = ' num2str(tc) ' s.']);


% -------------------------------------------------------------------------
% Raster plot of neighborhood neurons
% -------------------------------------------------------------------------
figure(3);
rasterPlot(spikeCell, neigh_NID, true);
title(['Raster plot of neighborhood neurons of neuron (r: ' ...
    int2str(center_r) ', c: ' int2str(center_c) ')']);
xlim([startT endT]);



