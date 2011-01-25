% Plot of average maximum drift in both directions, over 200 s of
% simulation.
%
% All the input files should contain whole population spike responses of
% the rat without any velocity input, i.e. the rat not moving
%
% Algorithm:
%  1. jobNums specifies range of job numbers to use
%  2. if preprocess is true, the files are read from disk from directory
%     specified by 'folder' and relative drifts extracted from each file
%  3. Drifts are stored in 'blobTracks_r/c' array
%  4. Then a 2D plot is created by plotting maximum drift in x and y
%     direction, plus a detailed plot of how the blobs in population
%     response drift, created from run which is most erroneous

close all;

startTime = 10;
endTime = 200;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.

preprocess = false;
    
jobNums = 40000:40099;

% Preprocess tracking data if necessary
if (preprocess == true)
    %folder = 'results/static_wave/zero_velocity/';
    folder = 'results/';

    nFiles = numel(jobNums);
    
    spikeHist_arr = [];
    blobTracks_r = [];
    blobTracks_c = [];

    for f_it = 1:nFiles
        d = dir([folder 'job' num2str(jobNums(f_it)) '*.mat'])
        % Load file which contains the tracking data and extract the spike
        % histogram and blob positions
        d(end).name
        load([folder d(end).name]);       
        opts = parseOptions(options);
        sheet_size = opts.sheet_size;
        disp 'Creating spike histogram';
        createSpikeHist;  % script - computes spikeHist for current file        
        spikeHist_arr(:, :, f_it) = spikeHist;
        
        disp 'Tracking blobs'
        [blobTracks_r(:, f_it) blobTracks_c(:, f_it)] = ...
            trackPopulationDrift(startTime, endTime, ...
            spikeHist_arr(:, :, f_it), dt_track, delta_t);        
    end
    
    
    for n_it = 0:sheet_size^2-1
        clear(['spikeMonitor_times_n' num2str(n_it)]);
    end
    clear spikeHist

    save('-v7.3', ['results/driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end)) '.mat']);
    
else  
    % Assuming results/driftResults.mat has been loaded
end



% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------
%plotTrack_r = abs(blobTracks_r);
%plotTrack_c = abs(blobTracks_c);

max_r = max(blobTracks_r);
max_c = max(blobTracks_c);
min_r = min(blobTracks_r);
min_c = min(blobTracks_c);

max_min_r = [max_r; min_r];
max_min_c = [max_c; min_c];

[absMax_r absMax_ri] = max(abs(max_min_r));
[absMax_c absMax_ci] = max(abs(max_min_c));

[absMax_rr absMax_rri] = max(absMax_r);


fontSize = 14;
figure('Position', [883 528 800 420]);
times = 0:dt_track:endTime-startTime-delta_t/2;
subplot(2, 5, [1 2 3], 'FontSize', fontSize);
plot(times, blobTracks_c(1:numel(times), absMax_rri), 'b');
ylabel('X drift (neurons)');
subplot(2, 5, [6 7 8], 'FontSize', fontSize);
plot(times, blobTracks_r(1:numel(times), absMax_rri), 'b');
%ylim([-5 25]);

xlabel('Time (s)');
ylabel('Y drift (neurons)');
%axis normal;

subplot(1, 10, [8 9 10], 'FontSize', fontSize);
plot(max_min_c(absMax_ci, :), max_min_r(absMax_ri, :), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
xlabel('X drift (neurons)');
ylabel('Y drift (neurons)');
hold on;
plot(max_min_c(absMax_ci(absMax_rri), absMax_rri), max_min_r(absMax_ri(absMax_rri), absMax_rri), 'sb', 'MarkerSize', 30);
%plot(max_min_c(absMax_ci(absMax_rri), absMax_rri), max_min_r(absMax_ri(absMax_rri), absMax_rri), 'ob', 'MarkerFaceColor', 'b');
%quiver(max_c(max_ri), max_r(max_ri), -5, -5);
hold off;
axis equal;
%ylim([0 24]);
%xlim([0 14]);


set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'results/maximumDrifts_Burak-Fiete-prefDirs.eps');
