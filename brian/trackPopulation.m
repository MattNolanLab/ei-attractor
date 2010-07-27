% Track the shift of population response throughout the simulation

%close all;

    opt = parseOptions(options);
    sheet_size = opt.sheet_size;

    %sheet_size = double(sheet_size);
    dt_rat = 0.02; % sec
    dt_track = 0.1; % sec; dt for the tracking algorithm
    delta_t = 0.5; % sec
    startTime = 0;
    endTime = 148; % sec
    
    firingPop = zeros(sheet_size, sheet_size);
    
    % Set the initial blob position to the center of the response
    blobPos_x = sheet_size/2;
    blobPos_y = sheet_size/2;
    
    % Preprocess spiking data: firing time --> histogram of firing, for
    % each neuron
    edges = linspace(0, endTime, endTime/dt_track);
    spikeHist = zeros(sheet_size^2, numel(edges));
    for x_i = 0:(sheet_size-1)
        for y_i = 0:(sheet_size-1)
            neuronID = y_i*sheet_size + x_i;            
            neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
            e_size = size(neuronSpikes, 1);
            if (e_size == 0)
                spikeHist(neuronID+1, :) = e_size;
            else
                spikeHist(neuronID+1, :) = histc(neuronSpikes, edges);
            end
        end
    end
    
    
    blobPos_x = [];
    blobPos_y = [];
    
    currBlobPos_x = sheet_size/2;
    currBlobPos_y = sheet_size/2;
    for t = startTime+delta_t/dt_track/2:dt_track:endTime-delta_t/dt_track/2
        for x_i = 0:(sheet_size-1)
            for y_i = 0:(sheet_size-1)
                neuronID = y_i*sheet_size + x_i;
                
                % take the window at position specified by t and sum up
                % number of spikes
                hist_i = fix(t/dt_track + 1);
                nbins = fix(delta_t/dt_track/2);
                bins = spikeHist(neuronID+1, hist_i-nbins:hist_i+nbins);

                firingPop(x_i+1, y_i+1) = sum(bins)/delta_t;
            end
        end

        firingPop = firingPop';

        % Simply threshold the population response to segment the image
        % This should easily work, since the blobs are coherent
        firingThr = 0.35;
        thrFiringPop = double(im2bw(firingPop/max(max(firingPop)), firingThr));

        [r, c] = trackBlobs(thrFiringPop);
        
        % Find the nearest blob to the last position
        [minDist min_i] = min((r-currBlobPos_y).^2 + (c-currBlobPos_x).^2);
        currBlobPos_x = c(min_i);
        currBlobPos_y = r(min_i);
        blobPos_x = [blobPos_x currBlobPos_x];
        blobPos_y = [blobPos_y currBlobPos_y];
        t
    end
    
