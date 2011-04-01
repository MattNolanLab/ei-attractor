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
%  4. Plot for different lambda_net parameters (13, 14, 16 18, 20)

path('../include/', path);
close all;
clear all;

startTime = 10;
endTime = 200;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.


driftsParams;  % include parameters

par{1} = SpacingDrifts_009.lambda_net_13;
par{2} = SpacingDrifts_009.lambda_net_14;
par{3} = SpacingDrifts_009.lambda_net_16;
par{4} = SpacingDrifts_009.lambda_net_18;
par{5} = StartFromEL_006.lambda_net_20;

fontSize = 14;
figure('Position', [883 528 1000 600]);

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
    subplot(2, 3, par_it, 'FontSize', fontSize);
    plotMaximumDrifts(gca, params, blobTracks_r, blobTracks_c);
    title(['\lambda_{net} = ' num2str(params.lambda_net)]);

    if (par_it ~= 1 && par_it ~= 4)
        ylabel('');
    end
    if (par_it < 3)
        xlabel('');
    end
end

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'output/009_SpacingDrifts/lambda_net_drifts.eps');


