% Plot of maximum drift in both directions, over 200 s of
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
%  4. NS directions (job21800) along with Burak and Fiete pref. dirs into one
%     plot

path('../include/', path);
close all;
clear all;

startTime = 10;
endTime = 200;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.


driftsParams;  % include parameters

par{1} = prefDirs_NS_000_005;
par{2} = Burak_Fiete_PrefDirs_000_001;

par{1}.MarkerFaceColor = [0 0 1];
par{2}.MarkerFaceColor = [1 0 0];


fontSize = 16;
figure('Position', [883 528 800 600]);

for par_it = 1:numel(par)
    params = par{par_it};

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
    subplot(1, 1, 1, 'FontSize', fontSize);
    box on;
    plotMaximumDrifts(gca, params, blobTracks_r, blobTracks_c);
end

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'output/report_2011-03-25-mult_bump_spiking_net/fig/compare_NS_vs_Burak.eps');


