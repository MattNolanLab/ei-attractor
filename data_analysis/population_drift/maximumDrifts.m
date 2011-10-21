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
path('../include/', path);
close all;
clear all;

startTime = 10;
endTime = 200;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.


folder = 'data/000_002_NoPrefDirs_RandomInit/';
jobNums = 21900:21999;


folder = params.folder; 
jobNums = params.jobNums;

    parfor f_it = 1:nFiles
        jobNums(f_it)
        d = dir([folder 'job' num2str(jobNums(f_it)) '*.mat'])
        if (numel(d) == 0)
            warning(['Could not open job no.' num2str(jobNums(f_it))]);
            continue;
        end
        % Load file which contains the tracking data and extract the spike
        % histogram and blob positions
        d(end).name
        %clear spikeCell
        dataLoad = load([folder d(end).name], 'spikeCell');       
        spikeCell = dataLoad.spikeCell;
        %opts = parseOptions(dataLoad.options);
        disp 'Creating spike histogram';
        spikeHist = createSpikeHistCell(1:numel(spikeCell), spikeCell, ...
            dt_track, 0, endTime);
        
        disp 'Tracking blobs'
        [blobTracks_r(:, f_it) blobTracks_c(:, f_it)] = ...
            trackPopulationDrift(startTime, endTime, spikeHist, ...
            dt_track, delta_t);        
    end
    
    
    clear spikeHist spikeCell

% Preprocess tracking data if necessary
if (params.preprocess == true)
    [blobTracks_r, blobTracks_c] = preprocessPopulationDrift(params, startTime, endTime, dt_track, delta_t);
else  
    load([folder '/driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end))]);
end



% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------

fontSize = 16;
figure('Position', [883 528 800 600]);
subplot(1, 1, 1, 'FontSize', fontSize);
for it = 1:size(blobTracks_r, 2)
    hold on;
    plot(max_min_c(absMax_ci(it), it), max_min_r(absMax_ri(it), it), '.', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
end
xlabel('X drift (neurons)');
ylabel('Y drift (neurons)');
%title(['Maximum drifts in X and Y directions, ' ...
%    'bump NO initialisation, lambda\_net = 18']);
title('Maximum drifts in X and Y directions, Burak&Fiete preferred directions');
axis equal;

set(gcf,'PaperPositionMode','auto');
print('-depsc2', ...
   'output/000_002_NoPrefDirs_RandomInit/noPrefDirs_RandomInit-maximumDrifts.eps');


% ------------------------------------------------------------------------
% Plot high drifts in both directions
% ------------------------------------------------------------------------
highDriftThreshold_c = 5;
highDriftThreshold_r = 10;

high_id_c = find(absMax_c > highDriftThreshold_c);
high_id_r = find(absMax_r > highDriftThreshold_r);

figure('Position', [0 0 800, 500]);
subplot(2, 1, 1);
if (numel(high_id_c) >0)
    plot(startTime:dt_track:endTime-dt_track, blobTracks_c(:, high_id_c));
end
ylabel('X Drift (neurons)');
title(['Detailed drifts of populations that reach absmax > ' num2str(highDriftThreshold_c)]);
for it = 1:numel(high_id_c)
    leg_c{it} = ['job ' num2str(high_id_c(it) - 1)];
end
legend(leg_c, 'Location', 'SouthEast');

subplot(2, 1, 2);
if (numel(high_id_r) > 0)
    plot(startTime:dt_track:endTime-dt_track, blobTracks_r(:, high_id_r));
end
xlabel('Time (s)');
ylabel('Y Drift (neurons)');
title(['Detailed drifts of populations that reach absmax > ' num2str(highDriftThreshold_r)]);
for it = 1:numel(high_id_c)
    leg_c{it} = ['job ' num2str(high_id_c(it) - 1)];
end
legend(leg_c, 'Location', 'SouthEast');
