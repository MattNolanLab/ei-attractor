% Specify a neuron ID and its neighborhood radius, to create a figure of a
% synchronization index (defined by a MvR distance)
%
% Just now, we assume that population response doesn't move

path('../include/', path);
close all;
clear all;

jobNums = 40000:40040;
dataFolder = '../../../central_data_store/simulation_data/spike_synchronization/';
outputFolder = '../../../central_data_store/data_analysis/spike_synchronisation/time_dependent_sync/SyncAndDrift/'
load([dataFolder 'driftsResults_40000-40099.mat'], 'blobTracks_c', 'blobTracks_r');


% MvR distance parameters
tc = 0.025;
dt = 0.001;
startT = 10;
endT = 200;
ISI = 0.1; % average isi
rad = 20; % Neighborhood radius - in neural units

syncWinT = 2; % sec



% The value of venterListID has been determined from freq. spectra of Vm recordings of the
% neurons see script figureMembraneVFreq.m

nID = 4728; % job40000 - estimated from raster plot
%nID = 7170; % job40003

spikes_runs = {};
for jobNum = jobNums
    d = dir([dataFolder 'job' num2str(jobNum) '*.mat'])
    % Load file which contains the tracking data and extract the spike
    % histogram and blob positions
    d(end).name
    clear spikeCell;
    load([dataFolder d(end).name], 'spikeCell');      
    %opt = parseOptions(options);
    spikes_runs{jobNum-jobNums(1)+1} = spikeCell;
end

parfor jobNum = jobNums
    cell_id = jobNum - jobNums(1) + 1

    [pairTimedSync, meanTimedSync, stdTimedSync, neigh_NID] = ...
        timeDependentSync(spikes_runs{cell_id}, tc, dt, startT, endT, ISI, rad, nID, syncWinT);

    sz = size(pairTimedSync);

    figure('Position', [883 528 800 900]); %, 'Visible', 'off');
    % -------------------------------------------------------------------------
    % Plot time dependent histogram map of sync. coefficient
    % -------------------------------------------------------------------------
    cb_ax = subplot(12, 1, [1 2 3 4]);
    [timedHist, binCenters] = hist(pairTimedSync, 40);
    pcolor(startT:syncWinT:endT, binCenters, timedHist);
    xlim([startT endT]);
    shading flat;
    colorbar('NorthOutside');
    xlabel('Time (s)');
    ylabel('D^2_{MvR}');
    title('Time dependent proportion of pairs of spikes given a sync. coeff');


    % -------------------------------------------------------------------------
    % Plot time dependent average sync. coefficient
    % -------------------------------------------------------------------------

    subplot(12, 1, [6 7 8]);
    plot(startT:syncWinT:endT, meanTimedSync);
    xlim([startT endT]);
    ylim([0 1]);
    xlabel('Time (s)');
    ylabel('Mean D^2_{MvR}');
    title(['Time dependent sync. measure. Blob and around. t_c = ' num2str(tc) ' s.']);

%     % Check for the neighborhood neurons
%     figure(2);
%     plot(mod(neigh_NID-1, sheet_size), fix((neigh_NID-1)./sheet_size), '.');



    % % -------------------------------------------------------------------------
    % % Raster plot of neighborhood neurons
    % % -------------------------------------------------------------------------
    % figure(3);
    % rasterPlot(spikeCell, neigh_NID, true);
    % title(['Raster plot of neighborhood neurons of neuron (r: ' ...
    %     int2str(center_r) ', c: ' int2str(center_c) ')']);
    % xlim([startT endT]);
    % 


    % -------------------------------------------------------------------------
    % Now produce plot of drift as a third subplot to compare with sync
    % Requires loading of preprocess tracked blobs
    % -------------------------------------------------------------------------
    trackStartTime = 10;
    trackEndTime = 200;

    dt_track = 0.1;
    delta_t = 0.25; % Should be this value.

    trackTimeAxis = trackStartTime:dt_track:trackEndTime - dt_track;
    
    
    job_id = jobNum - jobNums(1) + 1;
    subplot(12, 1, [10 11 12]);
    plot(trackTimeAxis, sqrt(blobTracks_r(:, job_id).^2 + blobTracks_c(:, job_id).^2));
    xlim([startT endT]);
    
    xlabel('Time (s)');
    ylabel('Drift [neurons]');
    title('Time dependent drift (Euc. dist. from zero)');
    
    
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', [outputFolder 'SyncAndDrift_job' num2str(jobNum) '.eps']);
end




