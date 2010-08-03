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
    [start_i end_i] = regexp(fileName, 'job\d+_\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_', 'Start', 'End');
    fileBase = fileName(start_i:end_i-1);
    
    load(fileName);

    opt = parseOptions(options);
    optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];
    
    sheet_size = opt.sheet_size;


    printSN = true;
    printPop = false;
    saveFig = true;
    %saveFig = false;


    dt_rat = 0.02; % sec
    %endTime = 1200 - dt_rat; % sec
    

    h = 5.0;  % cm
    arenaDiam = 180;   % cm

    if printSN
        for neuronNum = neuronIDs
            disp(['Processing data for neuron ' int2str(neuronNum) '...']);

            spikePlotFile = [saveDir fileBase optStr '_n' int2str(neuronNum) '_spikeMap.eps'];
            ratePlotFile =  [saveDir fileBase optStr '_n' int2str(neuronNum) '_rateMap.eps'];

            neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);
            %neuronSpikes = rat11015_t5c3_timeStamps;

            spikeFig = figure('Visible', 'off');
            plotSpikes_xy(neuronSpikes, pos_x, pos_y, dt_rat, neuronNum);
            title(['Spike plot - SN - l: ' num2str(opt.l) ', threshold:' num2str(opt.threshold) ', connMult: ' num2str(opt.connMult) ', alpha:' num2str(opt.alpha) ',lambda_{net}: ' num2str(opt.lambda_net)]);
            %title(['Spike plot - SN - taum: ' num2str(opt.taum) ', taui: ' num2str(opt.taui)]);
            %title(['Spike plot - SN - lambda_net: ' num2str(opt.lambda_net)]);
            if saveFig
                print('-depsc', spikePlotFile);
            end

            rateFig = figure('Visible', 'off');
            plotSNResponse(neuronSpikes, pos_x, pos_y, arenaDiam, h, dt_rat, neuronNum);
            title(['Rate plot - SN - l: ' num2str(opt.l) ', threshold:' num2str(opt.threshold) ', connMult: ' num2str(opt.connMult) ', alpha:' num2str(opt.alpha) ',lambda_{net}: ' num2str(opt.lambda_net)]);
            %title(['Rate plot - SN - taum: ' num2str(opt.taum) ', taui: ' num2str(opt.taui)]);
            %title(['Rate plot - SN - lambda_net: ' num2str(opt.lambda_net)]);
            if saveFig
                print('-depsc', ratePlotFile);
            end
        end
    end

    
    if printPop
        disp 'Processing population firing rate'
        %
        % Print population firing rate at specified time
        %
        startTime = 4;
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
        %popFig = figure();
        pcolor(0:sheet_size-1,0:sheet_size-1,firingPop);
        title(['Population rate plot - l: ' num2str(opt.l) ', threshold:' num2str(opt.threshold) ', connMult: ' num2str(opt.connMult) ', alpha:' num2str(opt.alpha) ',lambda_{net}: ' num2str(opt.lambda_net)]);
        %title(['Population rate plot - taum: ' num2str(opt.taum) ', taui: ' num2str(opt.taui)]);
        %title(['Population rate plot - SN - lambda_net: ' num2str(opt.lambda_net)]);

        colorbar;
        axis square tight;
        shading interp;    

        if saveFig
            popPlotFile =  [saveDir fileBase optStr '_popRateMap.eps'];
            print('-depsc', popPlotFile);
        end
    end
end

