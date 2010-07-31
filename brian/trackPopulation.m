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
    
    
    blobPos_r = [0];
    blobPos_c = [0];
    
    t_start = startTime+delta_t/2;
    t_end = endTime-delta_t/2;
    
    firingPop = getFiringPop(spikeHist, t_start, dt_track, delta_t);
    [currBlobPos_r currBlobPos_c] = trackBlobs(firingPop);
    
    % Find a blob nearest to the center of sheet
    c_r = sheet_size/2;
    c_c = sheet_size/2;
    [cMin cMin_i] = min((currBlobPos_r - c_r).^2 + (currBlobPos_c - c_c).^2);
    
    currBlobPos_r = currBlobPos_r(cMin_i);
    currBlobPos_c = currBlobPos_c(cMin_i);
    
    shift_r = currBlobPos_r - blobPos_r(end);
    shift_c = currBlobPos_c - blobPos_c(end);
    
    for t = t_start:dt_track:t_end
        firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
        [r, c] = trackBlobs(firingPop);
        
        % Find the blob nearest to the last position
        %curr_r = repmat(currBlobPos_r', 1, size(r, 2));
        %curr_c = repmat(currBlobPos_c', 1, size(c, 2));
        
        %r = repmat(r, size(currBlobPos_r, 2), 1);
        %c = repmat(c, size(currBlobPos_c, 2), 1);
        
        % Don't change order of r - blobDist_r neither for c. It is
        % important for minimum in the next step
        
        dist = (r-currBlobPos_r).^2 + (c-currBlobPos_c).^2;
        [minDist min_i] = min(dist);
        
        %[sortDist sort_i] = sort(minDist);        
        %nBlobs = size(sort_i, 2);
        %quasiMedian_i = round(nBlobs/2);
        
        new_r = r(min_i);
        new_c = c(min_i);
        
        %old_r = currBlobPos_r(min_i(sort_i(quasiMedian_i-2:quasiMedian_i+2)));
        %old_c = currBlobPos_c(min_i(sort_i(quasiMedian_i-2:quasiMedian_i+2)));

        %r_shift = mean(new_r - old_r);
        %c_shift = mean(new_c - old_c);
        
        
        blobPos_r = [blobPos_r (new_r - shift_r)];
        blobPos_c = [blobPos_c (new_c - shift_c)];

        % Check if new position of the blob is too near the edge, if so
        % choose new blob near the center and reset shift
        minEdgeDist =  10;
        if (sheet_size - new_r < minEdgeDist || new_r < minEdgeDist || ...
            sheet_size - new_c < minEdgeDist || new_c < minEdgeDist)
            disp('change of blob');
            [cMin cMin_i] = min((r - c_r).^2 + (c - c_c).^2);
            currBlobPos_r = r(cMin_i);
            currBlobPos_c = c(cMin_i);
            
            shift_r = currBlobPos_r - blobPos_r(end);
            shift_c = currBlobPos_c - blobPos_c(end);
        else
            currBlobPos_r = new_r;
            currBlobPos_c = new_c;
        end
        
        %t
    end
    
    figure('Visible', 'off');
    plot(blobPos_c, blobPos_r);
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
