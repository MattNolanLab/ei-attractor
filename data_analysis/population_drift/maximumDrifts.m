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


driftsParams;  % include parameters
params = SpacingDrifts_009.lambda_net_18;


folder = params.folder; 
jobNums = params.jobNums;


% Preprocess tracking data if necessary
if (params.preprocess == true)
    [blobTracks_r, blobTracks_c] = preprocessPopulationDrift(params, startTime, endTime, dt_track, delta_t);
else  
    load([folder '/driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end))]);
end



% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------

fontSize = 14;
figure('Position', [883 528 800 420]);
subplot(1, 1, 1, 'FontSize', fontSize);
    
[absMax_r, absMax_c] = plotMaximumDrifts(gca, params, blobTracks_r, blobTracks_c);

set(gcf,'PaperPositionMode','auto');
print('-depsc2', params.output);


% ------------------------------------------------------------------------
% Plot high drifts in both directions
% ------------------------------------------------------------------------
highDriftThreshold_c = 5;
highDriftThreshold_r = 5;

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
for it = 1:numel(high_id_r)
    leg_c{it} = ['job ' num2str(high_id_r(it) - 1)];
end
legend(leg_c, 'Location', 'SouthEast');
