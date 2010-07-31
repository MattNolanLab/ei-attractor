function trackPopulation(fileName)
    % Track the shift of population response throughout the simulation
    fileName
    
    saveDir = 'results/fig/';
    [start_i end_i] = regexp(fileName, 'job\d+_\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_', 'Start', 'End');
    fileBase = fileName(start_i:end_i-1);
    
    load(fileName);

    opt = parseOptions(options);
    optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];   

    sheet_size = opt.sheet_size;

    %sheet_size = double(sheet_size);
    dt_rat = 0.02; % sec
    dt_track = 0.1; % sec; dt for the tracking algorithm
    delta_t = 0.5; % sec
    startTime = 0;
    endTime = 150; % sec
    
    saveFig = true;
    
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
                spikeHist(neuronID+1, :) = 0;
            else
                spikeHist(neuronID+1, :) = histc(neuronSpikes, edges);
            end
        end
    end
    
    
    blobPos_x = [0];
    blobPos_y = [0];
    
    t_start = startTime+delta_t/2;
    t_end = endTime-delta_t/2;
    
    firingPop = getFiringPop(spikeHist, t_start, dt_track, delta_t);
    [currBlobPos_r currBlobPos_c] = trackBlobs(firingPop);
    
    for t = t_start:dt_track:t_end
        firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
        [r, c] = trackBlobs(firingPop);
        
        % Find the nearest blob to the last position
        curr_r = repmat(currBlobPos_r', 1, size(r, 2));
        curr_c = repmat(currBlobPos_c', 1, size(c, 2));
        
        r = repmat(r, size(currBlobPos_r, 2), 1);
        c = repmat(c, size(currBlobPos_c, 2), 1);
        
        % Don't change order of r - blobDist_r neither for c. It is
        % important for minimum in the next step
        dist = (r-curr_r).^2 + (c-curr_c).^2;
        [minDist min_i] = min(dist);
        [sortDist sort_i] = sort(minDist);
        
        nBlobs = size(sort_i, 2);
        quasiMedian_i = round(nBlobs/2);
        
        new_r = r(1, sort_i(quasiMedian_i-2:quasiMedian_i+2));
        new_c = c(1, sort_i(quasiMedian_i-2:quasiMedian_i+2));
        
        old_r = currBlobPos_r(min_i(sort_i(quasiMedian_i-2:quasiMedian_i+2)));
        old_c = currBlobPos_c(min_i(sort_i(quasiMedian_i-2:quasiMedian_i+2)));

        r_shift = mean(new_r - old_r);
        c_shift = mean(new_c - old_c);
        
        blobPos_x = [blobPos_x (blobPos_x(end) + c_shift)];
        blobPos_y = [blobPos_y (blobPos_y(end) + r_shift)];
        
        currBlobPos_r = r(1, :);
        currBlobPos_c = c(1, :);
        t
    end
    
    figure('Visible', 'off');
    plot(blobPos_x, blobPos_y);
    %xlim([0 sheet_size]); ylim([0 sheet_size]);
    axis square;
    
    xlabel('Neuron no.');
    ylabel('Neuron no.');
    title(['Blob position: threshold:' num2str(opt.threshold) ', connMult: ' num2str(opt.connMult) ', alpha:' num2str(opt.alpha) ',lambda_{net}: ' num2str(opt.lambda_net)]);
    
    if saveFig
        popPlotFile =  [saveDir fileBase '_tracking.eps'];
        print('-depsc', popPlotFile);
    end

    
end
