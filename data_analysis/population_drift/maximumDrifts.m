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
endTime = 30;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.


% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
folders.noisy_membrane_potential_001 = 'data/001_noisy_membrane_potential/sigma_0_01/';

Burak_Fiete_PrefDirs_000_001.folder = 'data/000_001_Burak_Fiete_PrefDirs/';
Burak_Fiete_PrefDirs_000_001.jobNums = 40000:40099;
Burak_Fiete_PrefDirs_000_001.title = '';
Burak_Fiete_PrefDirs_000_001.preprocess = false;
Burak_Fiete_PrefDirs_000_001.xlim = [-20 20];
Burak_Fiete_PrefDirs_000_001.ylim = [-15 15];
Burak_Fiete_PrefDirs_000_001.output = 'output/000_001_Burak_Fiete_PrefDirs/maximumDrifts-40000-40099.eps';


bump_initialized_002.folder = 'data/002_bump_initialized/';
bump_initialized_002.jobNums = 2100:2199;
bump_initialized_002.title = '';
bump_initialized_002.preprocess = false;
bump_initialized_002.xlim = [-20 20];
bump_initialized_002.ylim = [-15 15];
bump_initialized_002.output = 'output/002_bump_initialized/maximumDrifts-2100-2199.eps';

timeNoise_noInitNoise_003.sigma_0_1.folder = 'data/003_timeNoise_noInitNoise/sigma_0_1/';
timeNoise_noInitNoise_003.sigma_0_1.jobNums = 20200:20299;
timeNoise_noInitNoise_003.sigma_0_1.title = '';
timeNoise_noInitNoise_003.sigma_0_1.preprocess = false;
timeNoise_noInitNoise_003.sigma_0_1.xlim = [-10 10];
timeNoise_noInitNoise_003.sigma_0_1.ylim = [-10 10];
timeNoise_noInitNoise_003.sigma_0_1.output = 'output/003_timeNoise_noInitNoise/noise_sigma_0_1-maximumDrifts.eps';

StartFromEL_006.lambda_net_20.folder = 'data/006_StartFromEL/lambda_net_20';
StartFromEL_006.lambda_net_20.jobNums = 30100:30199;
StartFromEL_006.lambda_net_20.title = '';
StartFromEL_006.lambda_net_20.preprocess = false;
StartFromEL_006.lambda_net_20.xlim = [-10 10];
StartFromEL_006.lambda_net_20.ylim = [-10 10];
StartFromEL_006.lambda_net_20.output = 'output/006_StartFromEL/lambda_net_20-maximumDrifts.eps';

SpacingDrifts_009.lambda_net_13.folder = 'data/009_SpacingDrifts/lambda_net_13/';
SpacingDrifts_009.lambda_net_13.jobNums = 90000:90099;
SpacingDrifts_009.lambda_net_13.title = '';
SpacingDrifts_009.lambda_net_13.preprocess = false;
SpacingDrifts_009.lambda_net_13.xlim = [-10 10];
SpacingDrifts_009.lambda_net_13.ylim = [-10 10];
SpacingDrifts_009.lambda_net_13.output = 'output/009_SpacingDrifts/lambda_net_13-maximumDrifts.eps';

SpacingDrifts_009.lambda_net_14.folder = 'data/009_SpacingDrifts/lambda_net_14/';
SpacingDrifts_009.lambda_net_14.jobNums = 90100:90199;
SpacingDrifts_009.lambda_net_14.title = '';
SpacingDrifts_009.lambda_net_14.preprocess = false;
SpacingDrifts_009.lambda_net_14.xlim = [-10 10];
SpacingDrifts_009.lambda_net_14.ylim = [-10 10];
SpacingDrifts_009.lambda_net_14.output = 'output/009_SpacingDrifts/lambda_net_14-maximumDrifts.eps';

SpacingDrifts_009.lambda_net_16.folder = 'data/009_SpacingDrifts/lambda_net_16/';
SpacingDrifts_009.lambda_net_16.jobNums = 90200:90299;
SpacingDrifts_009.lambda_net_16.title = '';
SpacingDrifts_009.lambda_net_16.preprocess = false;
SpacingDrifts_009.lambda_net_16.xlim = [-10 10];
SpacingDrifts_009.lambda_net_16.ylim = [-10 10];
SpacingDrifts_009.lambda_net_16.output = 'output/009_SpacingDrifts/lambda_net_16-maximumDrifts.eps';

SpacingDrifts_009.lambda_net_18.folder = 'data/009_SpacingDrifts/lambda_net_18/';
SpacingDrifts_009.lambda_net_18.jobNums = 90300:90399;
SpacingDrifts_009.lambda_net_18.title = '';
SpacingDrifts_009.lambda_net_18.preprocess = false;
SpacingDrifts_009.lambda_net_18.xlim = [-10 10];
SpacingDrifts_009.lambda_net_18.ylim = [-10 10];
SpacingDrifts_009.lambda_net_18.output = 'output/009_SpacingDrifts/lambda_net_18-maximumDrifts.eps';

% -------------------------------------------------------------------------
% End parameters
% -------------------------------------------------------------------------
params = SpacingDrifts_009.lambda_net_18;


folder = params.folder; %'data/001_noisy_membrane_potential/sigma_0_01/';
jobNums = params.jobNums;

% Preprocess tracking data if necessary
if (params.preprocess == true)

    nFiles = numel(jobNums);
    
    spikeHist_arr = [];
    blobTracks_r = zeros(numel(startTime:dt_track:endTime)-1, numel(jobNums));
    blobTracks_c = blobTracks_r;

    for f_it = 1:nFiles
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

    save('-v7.3', [folder 'driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end)) '.mat']);
    
else  
    load([folder '/driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end))]);
end



% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------
max_r = max(blobTracks_r);
max_c = max(blobTracks_c);
min_r = min(blobTracks_r);
min_c = min(blobTracks_c);

max_min_r = [max_r; min_r];
max_min_c = [max_c; min_c];

[absMax_r absMax_ri] = max(abs(max_min_r));
[absMax_c absMax_ci] = max(abs(max_min_c));

[absMax_rr absMax_rri] = max(absMax_c);


fontSize = 14;
figure('Position', [883 528 800 420]);


subplot(1, 1, 1, 'FontSize', fontSize);
for it = 1:size(blobTracks_r, 2)
    hold on;
    plot(max_min_c(absMax_ci(it), it), max_min_r(absMax_ri(it), it), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
end
xlabel('X drift (neurons)');
ylabel('Y drift (neurons)');

title(params.title);
% title(['Maximum drifts in X and Y directions, ' ...
%     'bump NO initialisation, noise\_sigma = 0.01 mV']);
axis equal;
xlim(params.xlim);
ylim(params.ylim);


set(gcf,'PaperPositionMode','auto');
print('-depsc2', params.output);


% ------------------------------------------------------------------------
% Plot high drifts in both directions
% ------------------------------------------------------------------------
highDriftThreshold_c = 6;
highDriftThreshold_r = 6;

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
