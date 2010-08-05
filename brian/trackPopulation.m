function trackPopulation(fileName, oldFormat, startTime, endTime)
    % Track the shift of population response throughout the simulation
    fileName
    
    saveDir = 'results/fig/';
    [start_i end_i] = regexp(fileName, 'job\d+_\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_', 'Start', 'End');
    fileBase = fileName(start_i:end_i-1);
    
    load(fileName);
    
    if (oldFormat)
        pos_x = ratData_pos_x;
        pos_y = ratData_pos_y;
    end
    


    opt = parseOptions(options);
    optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];   

    sheet_size = opt.sheet_size;

    %sheet_size = double(sheet_size);
    dt_rat = 0.02; % sec
    dt_track = 0.1; % sec; dt for the tracking algorithm
    delta_t = 0.5; % sec
    %startTime = 0;
    %endTime = 150; % sec
    
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
    
    t_start = startTime;
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
        
        % Find blob nearest to the last position and record its position
        dist = (r-currBlobPos_r).^2 + (c-currBlobPos_c).^2;
        [minDist min_i] = min(dist);
        
        new_r = r(min_i);
        new_c = c(min_i);
                
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
    
    % Integrate the velocity of the rat and velocity of the neural pattern
    % to obtain scaling factor
    startInt = 10; %sec
    endInt = 12; %sec
    
    startBlobPos_i = (startInt-startTime)/dt_track + 1; % compensate for the time shift of start of tracking
    endBlobPos_i   = (endInt-startTime)/dt_track + 1;
    
    startRatPos_i = startInt/dt_rat + 1; % rat timestamps begin at 0
    endRatPos_i = endInt/dt_rat + 1;
    
    blobDist = sqrt((blobPos_r(endBlobPos_i) - blobPos_r(startBlobPos_i))^2 + ...
        (blobPos_c(endBlobPos_i) - blobPos_c(startBlobPos_i))^2);
    
    ratDist = sqrt((pos_x(startRatPos_i) - pos_x(endRatPos_i))^2 + ...
        (pos_y(startRatPos_i) - pos_y(endRatPos_i))^2);
    
    % Scale factor is the fraction of how much the rat traveled in the
    % estimation period to how many neurons the population firing pattern
    % traveled
    scaleFactor = ratDist/blobDist
    
    % It should suffice to multiply the neuron position by the scale factor
    % to obtain the trajectory in cm
    blobPos_cm_r = pos_y(1) + blobPos_r * scaleFactor;
    blobPos_cm_c = pos_x(1) + blobPos_c * scaleFactor;
    
    startRatDrift_i = ceil(t_start/dt_rat + 1);
    endRatDrift_i = ceil(t_end/dt_rat + 1);
    nDrift = startRatDrift_i:dt_track/dt_rat:endRatDrift_i;
    drift = sqrt((blobPos_cm_r(1:numel(nDrift))' - pos_y(nDrift)).^2 + ...
        (blobPos_cm_c(1:numel(nDrift))' - pos_x(nDrift)).^2);
    
    fontSize = 14;
    
    figure('Visible', 'off');
    subplot(1, 1, 1, 'FontSize', fontSize);
    times = t_start:dt_track:t_end;
    plot(times(1:numel(nDrift)), drift, 'k');
    
    xlabel('Time (s)');
    ylabel('Drift (cm)');
    %title(['Estimated drift: threshold: ' num2str(opt.threshold) ', connMult: ' num2str(opt.connMult) ', alpha: ' num2str(opt.alpha) ',lambda_{net}: ' num2str(opt.lambda_net)]);
    
    
%     figure(2);
%     plot(blobPos_cm_c, blobPos_cm_r);
%     xlim([0 sheet_size]); ylim([0 sheet_size]);
%     axis equal;    
    
    if saveFig
        popPlotFile =  [saveDir fileBase '_tracking.eps'];
        print('-depsc', popPlotFile);
    end
end
