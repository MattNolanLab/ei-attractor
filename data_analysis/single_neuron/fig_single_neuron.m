%function trackPopulation(fileName, oldFormat, startTime, endTime)
% Track the shift of population response throughout the simulation
path('../include', path);

close all;
clear all;
    

    % Load file and process? If not - assuming the tracking has already
    % been done and results saved to tracking*.mat file
    %loadFlag = true;
    
    jobId = 40100;
    folder = 'simulation_data/000_003_Burak_Fiete_Path_Integ/';
    d = dir([folder 'job' num2str(jobId) '*.mat']);
    fileName = [folder d(end).name]
    
    saveFig = true;
    
    oldFormat = false;
    
    dt_rat = 0.02; % sec
    dt_track = 0.1; % sec; dt for the tracking algorithm
    delta_t = 0.25; % sec

    startTime = 10; t_start = startTime;
    endTime = 1200; t_end = endTime - delta_t/2;
    
    
    saveDir = 'results/fig/';
    [start_i end_i] = regexp(fileName, 'job\d+_\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_', 'Start', 'End');
    fileBase = fileName(start_i:end_i-1);
    
    load(fileName);
    opt = parseOptions(options);
    optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];   
    sheet_size = opt.sheet_size;

    if (oldFormat)
        pos_x = ratData_pos_x;
        pos_y = ratData_pos_y;
    end


    %---------------------------------------------
    % Plot single neuron rate and spike responses
    %---------------------------------------------
    neuronNum = sheet_size^2 / 2 + sheet_size/2;
    neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);
    h = 5.0;  % cm
    arenaDiam = 180;   % cm
    
    fontSize = 18;

    figure('Position', [300 200 1300 600]);
    subplot(1, 2, 1, 'FontSize', fontSize);
    plotSpikes_xy(neuronSpikes, pos_x, pos_y, dt_rat, neuronNum);
    %xlabel('Rat position [cm]');
    %ylabel('Rat position [cm]');
    xlim([-arenaDiam/2 arenaDiam/2]);
    ylim([-arenaDiam/2 arenaDiam/2]);
%    set(gca(), 'XTick', [-arenaDiam/2 arenaDiam/2]);
%    set(gca(), 'YTick', [-arenaDiam/2 arenaDiam/2]);
    %set(gca(), 'XTick', []);
    %set(gca(), 'YTick', []);
    title 'A';

    
    subplot(1, 2, 2, 'FontSize', fontSize);
    plotSNResponse(neuronSpikes, pos_x, pos_y, arenaDiam, h, dt_rat, neuronNum);
    xlim([-arenaDiam/2 arenaDiam/2]);
    ylim([-arenaDiam/2 arenaDiam/2]);
%    set(gca(), 'XTick', [-arenaDiam/2 arenaDiam/2]);
%    set(gca(), 'YTick', [-arenaDiam/2 arenaDiam/2]);
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
