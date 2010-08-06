%function trackPopulation(fileName, oldFormat, startTime, endTime)
% Track the shift of population response throughout the simulation

close all;
    

    % Load file and process? If not - assuming the tracking has already
    % been done and results saved to tracking*.mat file
    loadFlag = false;
    
    jobId = 22204;
    d = dir(['results/job' num2str(jobId) '*.mat']);
    fileName = ['results/' d(end).name]
    
    oldFormat = false;
    startTime = 4;
    endTime = 1200;
    
    saveDir = 'results/fig/';
    [start_i end_i] = regexp(fileName, 'job\d+_\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_', 'Start', 'End');
    fileBase = fileName(start_i:end_i-1);
    
    if loadFlag == true
    
        load(fileName);

        if (oldFormat)
            pos_x = ratData_pos_x;
            pos_y = ratData_pos_y;
        end



        opt = parseOptions(options);
        optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];   

        sheet_size = opt.sheet_size;

        dt_rat = 0.02; % sec
        dt_track = 0.02; % sec; dt for the tracking algorithm
        delta_t = 0.25; % sec

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

        updateId = 0;
        nUpdate = 100;
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

            if (mod(updateId, nUpdate) == 0)
                t
            end
            updateId = updateId + 1;
        end
    end
        
    % Integrate the velocity of the rat and velocity of the neural pattern
    % to obtain scaling factor
    startInt = startTime + 10; %sec
    endInt = startInt+2; %sec
    
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
    blobPos_cm_r = pos_y(t_start/dt_rat + 1) + blobPos_r * scaleFactor;
    blobPos_cm_c = pos_x(t_start/dt_rat + 1) + blobPos_c * scaleFactor;
    
    startRatDrift_i = ceil(t_start/dt_rat) + 1;
    endRatDrift_i = ceil(t_end/dt_rat) + 1;
    nDrift = startRatDrift_i:dt_track/dt_rat:endRatDrift_i;
    drift = sqrt((blobPos_cm_r(1:numel(nDrift))' - pos_y(nDrift)).^2 + ...
        (blobPos_cm_c(1:numel(nDrift))' - pos_x(nDrift)).^2);
    
    
    fontSize = 15;
    figure('Position', [1000, 300, 500, 500]);
    subplot(2, 2, [3 4], 'FontSize', fontSize);
    times = (1:numel(nDrift)) * dt_track + t_start;
    plot(times, drift, 'k');
    
    xlabel('Time (s)');
    ylabel('Drift (cm)');
    
    %---------------------------------------------
    % Plot single neuron rate and spike responses
    %---------------------------------------------
    neuronNum = sheet_size^2 / 2;
    neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);
    h = 5.0;  % cm
    arenaDiam = 180;   % cm

    subplot(2, 2, 1, 'FontSize', fontSize);
    plotSpikes_xy(neuronSpikes, pos_x, pos_y, dt_rat, neuronNum);
    %xlabel('Rat position [cm]');
    %ylabel('Rat position [cm]');
    xlim([-arenaDiam/2 arenaDiam/2]);
    ylim([-arenaDiam/2 arenaDiam/2]);
    set(gca(), 'XTick', [-arenaDiam/2 arenaDiam/2]);
    set(gca(), 'YTick', [-arenaDiam/2 arenaDiam/2]);
    %set(gca(), 'XTick', []);
    %set(gca(), 'YTick', []);

    
    subplot(2, 2, 2, 'FontSize', fontSize);
    plotSNResponse(neuronSpikes, pos_x, pos_y, arenaDiam, h, dt_rat, neuronNum);
    xlim([-arenaDiam/2 arenaDiam/2]);
    ylim([-arenaDiam/2 arenaDiam/2]);
    set(gca(), 'XTick', [-arenaDiam/2 arenaDiam/2]);
    set(gca(), 'YTick', [-arenaDiam/2 arenaDiam/2]);
    %set(gca(), 'XTick', []);
    %set(gca(), 'YTick', []);

     
    if saveFig
        popPlotFile =  [saveDir fileBase '_tracking.eps'];
        print('-depsc', popPlotFile);
        
        for n_it = [0:neuronNum - 1 neuronNum+1:sheet_size^2-1]
            clear(['spikeMonitor_times_n' num2str(n_it)]);
        end
        clear spikeHist
        
        save('-v7.3', ['results/tracking_' d(end).name]);
    end
%end
