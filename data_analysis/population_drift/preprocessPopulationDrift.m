function [blobTracks_r, blobTracks_c] = preprocessPopulationDrift(params, startTime, endTime, dt_track, delta_t)
    % PREPROCESSPOPULATIONDRIFT take a set of parameters (folder, job numbers,
    % etc.), load files one by one (or in parallel) and output the tracking data
    % for each of the files. The tracking data are stored column-wise per each
    % processed file

    folder = params.folder; 
    jobNums = params.jobNums;

    nFiles = numel(jobNums);
    
    spikeHist_arr = [];
    %blobTracks_r = zeros(numel(startTime:dt_track:endTime)-1, numel(jobNums));
    %blobTracks_c = blobTracks_r;

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

    save('-v7.3', [folder 'driftsResults_' num2str(jobNums(1)) '-' num2str(jobNums(end)) '.mat']);
end
