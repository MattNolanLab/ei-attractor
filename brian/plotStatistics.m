function [spikeFig rateFig] = plotStatistics(fileName, neuronIDs)

    fileName
    neuronIDs

    errCode = 0;

    % Print statistics about the data

    close all;
    %clear all;

    % dir = 'results/fig/';
    % params = 's64_alpha0.01';
    % spikePlotFile = [dir 'spikeMap_' params '.eps'];
    % ratePlotFile  = [dir 'rateMap_'  params '.eps'];
    %loadDir = 'results/';
    saveDir = 'results/fig/';
    [start_i end_i] = regexp(fileName, '\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_job\d+_', 'Start', 'End');
    fileBase = fileName(start_i:end_i-1);
    
    load(fileName);

    [opt.sheet_size opt.alpha] = parseOptions(options);
    optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];
    
    sheet_size = opt.sheet_size;


    saveFig = true;
    %saveFig = false;


    dt_rat = 0.02; % sec
    %endTime = 1200 - dt_rat; % sec
    

    h = 5.0;  % cm
    arenaDiam = 180;   % cm

    for neuronNum = neuronIDs
        disp(['Processing data for neuron ' int2str(neuronNum) '...']);
        
        spikePlotFile = [saveDir fileBase optStr '_n' int2str(neuronNum) '_spikeMap.eps'];
        ratePlotFile =  [saveDir fileBase optStr '_n' int2str(neuronNum) '_rateMap.eps'];

        neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);
        %neuronSpikes = rat11015_t5c3_timeStamps;

        spikeFig = figure('Visible', 'off');
        plotSpikes_xy(neuronSpikes, ratData_pos_x, ratData_pos_y, dt_rat, neuronNum);
        if saveFig
            print('-depsc', spikePlotFile);
        end

        rateFig = figure('Visible', 'off');
        plotSNResponse(neuronSpikes, ratData_pos_x, ratData_pos_y, arenaDiam, h, dt_rat, neuronNum);
        if saveFig
            print('-depsc', ratePlotFile);
        end
    end

    
    %
    % Print population firing rate at specified time
    %
    startTime = 3;
    endTime = startTime; % sec
    delta_t = 1; % sec
    
    firingPop = zeros(sheet_size, sheet_size);
    
    for x_i = 0:(sheet_size-1)
        for y_i = 0:(sheet_size-1)
            neuronID = y_i*sheet_size + x_i;
            neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
            firingRate = computeFiringRate(neuronSpikes, startTime, endTime, dt_rat, delta_t);
            
            firingPop(x_i+1, y_i+1) = firingRate(numel(firingRate));
        end
    end

    popFig = figure('Visible', 'off');
    pcolor(0:sheet_size-1,0:sheet_size-1,firingPop);

    colorbar;
    axis square tight;
    shading interp;    
    
    if saveFig
        popPlotFile =  [saveDir fileBase optStr '_popRateMap.eps'];
        print('-depsc', popPlotFile);
    end
end

