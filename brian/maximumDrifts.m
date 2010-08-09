% Plot of average maximum drift in both directions, over 200 s of
% simulation withou any velocity input


startTime = 10;
endTime = 200;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.

preprocess = true;
    
jobNums = 21600:21679;

% Preprocess tracking data if necessary
if (preprocess == true)
    folder = 'results/static_wave/zero_velocity/';

    nFiles = numel(jobNums);
    
    spikeHist_arr = [];
    blobTracks_r = [];
    blobTracks_c = [];

    for f_it = 1:nFiles
        d = dir([folder 'job' num2str(jobNums(f_it)) '*.mat'])
        % Load file which contains the tracking data and extract the spike
        % histogram and blob positions
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

    save('-v7.3', 'results/driftsResults.mat');
    
else
end



% ------------------------------------------------------------------------
% Plot  both drifts in x and y directions
% ------------------------------------------------------------------------