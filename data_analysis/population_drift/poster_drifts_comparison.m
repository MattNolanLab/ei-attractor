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

fontSize = 22;
abs_shift = 0.0;


driftsParams;  % include parameters
% ------------------------------------------------------------------------
% 
% High Drift
%
% ------------------------------------------------------------------------
params = Burak_Fiete_PrefDirs_000_001;
params.MarkerFaceColor = 'b';

folder = params.folder; 
jobNums = params.jobNums;

load([folder '/driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end))]);



% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------

figure('Position', [600 0 1050 900]);
subplot(3, 3, 1, 'FontSize', fontSize);
box on;
    
[absMax_r_high, absMax_c_high] = plotMaximumDrifts(gca, params, blobTracks_r, blobTracks_c);
xlabel('');
pos = get(gca, 'Position');
pos(1) = pos(1) - abs_shift*pos(3);
set(gca, 'Position', pos);

% ------------------------------------------------------------------------
% Plot high drift magnitude (sqrt(c^2 + r^2))
% ------------------------------------------------------------------------
highDriftMagnitudeTh = 11;

driftMag = sqrt(blobTracks_c.^2 + blobTracks_r.^2);
maxMag = max(driftMag);

high_id = find(maxMag > highDriftMagnitudeTh);

subplot(3, 3, 2:3, 'FontSize', fontSize);
if (numel(high_id) >0)
    plot(startTime:dt_track:endTime-dt_track, driftMag(:, high_id));
end
%ylabel('Drift magnitude (nrns)');
%title(['Detailed drifts of populations that reach absmax > ' num2str(highDriftMagnitudeTh)]);
ylim([0 20]);



% ------------------------------------------------------------------------
% 
% Low Drift
%
% ------------------------------------------------------------------------
params = bump_initialized_002;
params.MarkerFaceColor = [0 0.5 0];

folder = params.folder; 
jobNums = params.jobNums;

load([folder '/driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end))]);

fontSize = 22;


% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------
subplot(3, 3, 4, 'FontSize', fontSize);    
set(gca, 'FontSize', fontSize);
[absMax_r_low, absMax_c_low] = plotMaximumDrifts(gca, params, blobTracks_r, blobTracks_c);
box on;
pos = get(gca, 'Position');
pos(1) = pos(1) - abs_shift*pos(3);
set(gca, 'Position', pos);
xlabel('');




% ------------------------------------------------------------------------
% Plot high drift magnitude (sqrt(c^2 + r^2))
% ------------------------------------------------------------------------
highDriftMagnitudeTh = 2;

driftMag = sqrt(blobTracks_c.^2 + blobTracks_r.^2);
maxMag = max(driftMag);

high_id = find(maxMag > highDriftMagnitudeTh);

subplot(3, 3, [5 6], 'FontSize', fontSize);
if (numel(high_id) >0)
    plot(startTime:dt_track:endTime-dt_track, driftMag(:, high_id));
end
%ylabel('Drift magnitude (nrns)');
%title(['Detailed drifts of populations that reach absmax > ' num2str(highDriftMagnitudeTh)]);
ylim([0 20]);
%xlabel('Time (s)');




% ------------------------------------------------------------------------
% 
% Asymmetric drift
%
% ------------------------------------------------------------------------
params = NoisyNetwork_NoPrefDirs_004.sigma_0_1;
params.MarkerFaceColor = [0.5 0 0];
params.xlim = [-20 20];
params.ylim = [-20 20];

folder = params.folder; 
jobNums = params.jobNums;

load([folder '/driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end))]);

fontSize = 22;


% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------
subplot(3, 3, 7, 'FontSize', fontSize);    
set(gca, 'FontSize', fontSize);
[absMax_r_asym, absMax_c_asym] = plotMaximumDrifts(gca, params, blobTracks_r, blobTracks_c);
box on;
pos = get(gca, 'Position');
pos(1) = pos(1) - abs_shift*pos(3);
set(gca, 'Position', pos);




% ------------------------------------------------------------------------
% Plot high drift magnitude (sqrt(c^2 + r^2))
% ------------------------------------------------------------------------
highDriftMagnitudeTh = 9;

driftMag = sqrt(blobTracks_c.^2 + blobTracks_r.^2);
maxMag = max(driftMag);

high_id = find(maxMag > highDriftMagnitudeTh);

subplot(3, 3, [8 9], 'FontSize', fontSize);
if (numel(high_id) >0)
    plot(startTime:dt_track:endTime-dt_track, driftMag(:, high_id));
end
%ylabel('Drift magnitude (nrns)');
%title(['Detailed drifts of populations that reach absmax > ' num2str(highDriftMagnitudeTh)]);
ylim([0 20]);
xlabel('Time (s)');


set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'poster_drifts_high_low.eps');


% ------------------------------------------------------------------------
% 
% Bar plot
%
% ------------------------------------------------------------------------
figure();
subplot(1,1,1, 'FontSize', fontSize);
absmax = [mean(absMax_c_high) mean(absMax_c_low) mean(absMax_c_asym); ...
    mean(absMax_r_high) mean(absMax_r_low) mean(absMax_r_asym)];
bar(absmax);
map(1, :) = [0     0 0.5];
map(2, :) = [0   0.5   0];
map(3, :) = [0.5   0   0];
colormap(map);
set(gca, 'XTickLabel', ['X'; 'Y']);
%ylabel('Maximal drift (neurons)');

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'poster_drifts_high_low_mean_bar.eps');



