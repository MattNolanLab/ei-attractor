%function trackPopulation(fileName, oldFormat, startTime, endTime)
% Track the shift of population response throughout the simulation
path('../include', path);

close all;
clear all;
    

    % Load file and process? If not - assuming the tracking has already
    % been done and results saved to tracking*.mat file
    loadFlag = true;
    
    jobId = 40106;
    folder = 'simulation_data/000_003_Burak_Fiete_Path_Integ/';
    d = dir([folder 'job' num2str(jobId) '*.mat']);
    fileName = [folder d(end).name]
    
    oldFormat = false;
    
    dt_rat = 0.02; % sec
    dt_track = 0.1; % sec; dt for the tracking algorithm
    delta_t = 0.25; % sec

    startTime = 10; t_start = startTime;
    endTime = 1200; t_end = endTime - delta_t/2;
    
    
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

        saveFig = true;

        %firingPop = zeros(sheet_size, sheet_size);
        %createSpikeHist  % a script - side effects
        %[blobPos_r blobPos_c] = trackPopulationDrift(startTime, endTime, spikeHist, dt_track, delta_t);
    end
        
    % Integrate the velocity of the rat and velocity of the neural pattern
    % to obtain scaling factor
%     startInt = startTime + 10; %sec
%     endInt = startInt+2; %sec
%     
%     startBlobPos_i = (startInt-startTime)/dt_track + 1; % compensate for the time shift of start of tracking
%     endBlobPos_i   = (endInt-startTime)/dt_track + 1;
%     
%     startRatPos_i = startInt/dt_rat + 1; % rat timestamps begin at 0
%     endRatPos_i = endInt/dt_rat + 1;
%     
%     blobDist = sqrt((blobPos_r(endBlobPos_i) - blobPos_r(startBlobPos_i))^2 + ...
%         (blobPos_c(endBlobPos_i) - blobPos_c(startBlobPos_i))^2);
%     
%     ratDist = sqrt((pos_x(startRatPos_i) - pos_x(endRatPos_i))^2 + ...
%         (pos_y(startRatPos_i) - pos_y(endRatPos_i))^2);
%     
%     % Scale factor is the fraction of how much the rat traveled in the
%     % estimation period to how many neurons the population firing pattern
%     % traveled
%     scaleFactor = ratDist/blobDist
%     
%     % It should suffice to multiply the neuron position by the scale factor
%     % to obtain the trajectory in cm
%     blobPos_cm_r = pos_y(t_start/dt_rat + 1) + blobPos_r * scaleFactor;
%     blobPos_cm_c = pos_x(t_start/dt_rat + 1) + blobPos_c * scaleFactor;
%     
%     startRatDrift_i = ceil(t_start/dt_rat) + 1;
%     endRatDrift_i = ceil(t_end/dt_rat) + 1;
%     nDrift = startRatDrift_i:dt_track/dt_rat:endRatDrift_i;
%     drift = sqrt((blobPos_cm_r(1:numel(nDrift))' - pos_y(nDrift)).^2 + ...
%         (blobPos_cm_c(1:numel(nDrift))' - pos_x(nDrift)).^2);
%     
%     
      fontSize = 15;
%     figure('Position', [1000, 300, 600, 500]);
%     subplot(2, 2, [3 4], 'FontSize', fontSize);
%     times = (1:numel(nDrift)) * dt_track + t_start;
%     plot(times, drift, 'k');
%     
%     xlabel('Time (s)');
%     ylabel('Drift (cm)');
%     title 'C';
%     
    %---------------------------------------------
    % Plot single neuron rate and spike responses
    %---------------------------------------------
    neuronNum = sheet_size^2 / 2 + sheet_size/2;
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
    title 'A';

    
    subplot(2, 2, 2, 'FontSize', fontSize);
    plotSNResponse(neuronSpikes, pos_x, pos_y, arenaDiam, h, dt_rat, neuronNum);
    xlim([-arenaDiam/2 arenaDiam/2]);
    ylim([-arenaDiam/2 arenaDiam/2]);
    set(gca(), 'XTick', [-arenaDiam/2 arenaDiam/2]);
    set(gca(), 'YTick', [-arenaDiam/2 arenaDiam/2]);
    %set(gca(), 'XTick', []);
    %set(gca(), 'YTick', []);
    title 'B';

     
    if saveFig
        popPlotFile =  [folder 'data_analysis/' fileBase '_tracking.eps'];
        set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
        print('-depsc', popPlotFile);
        
        for n_it = [0:neuronNum - 1 neuronNum+1:sheet_size^2-1]
            clear(['spikeMonitor_times_n' num2str(n_it)]);
        end
        clear spikeHist
        
        save('-v7.3', [folder 'data_analysis/tracking_' d(end).name]);
    end
%end
