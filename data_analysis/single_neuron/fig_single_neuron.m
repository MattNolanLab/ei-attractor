%function trackPopulation(fileName, oldFormat, startTime, endTime)
% Track the shift of population response throughout the simulation
path('../include', path);

close all;
clear all;
    

    % Load file and process? If not - assuming the tracking has already
    % been done and results saved to tracking*.mat file
    %loadFlag = true;
    
    jobId = 2385;
    chpntId = 1;
    folder = '../../../central_data_store/simulation_data/007_mult_bump_path_integ/';
    d = dir([folder 'job' num2str(jobId) '*.mat']);
    fileName = [folder d(chpntId).name]
    
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
    neuronNum = 1;
    neuronSpikes = spikeCell{neuronNum};
    h = 5.0;  % cm
    arenaDiam = 180;   % cm
    
    fontSize = 25;

    figure('Position', [100 200 1000 300]);
    subplot(1, 3, 1, 'FontSize', fontSize);
    plotSpikes_xy(neuronSpikes, pos_x, pos_y, dt_rat, neuronNum);
    xlim([-arenaDiam/2 arenaDiam/2]);
    ylim([-arenaDiam/2 arenaDiam/2]);
    %title 'A';
    axis equal tight;

    
    subplot(1, 3, 2, 'FontSize', fontSize);
    [firingHist xedges yedges] = plotSNResponse(neuronSpikes, pos_x, pos_y, arenaDiam, h, dt_rat, neuronNum);
    xlim([-arenaDiam/2 arenaDiam/2]);
    ylim([-arenaDiam/2 arenaDiam/2]);
    axis off;
    axis equal tight;
     

    %figure('Position', [100 200 1500 600]);    
    subplot(1, 3, 3, 'FontSize', fontSize);
    edges = linspace(-arenaDiam, arenaDiam, 2*arenaDiam/h + 1);
    
    firingRate = firingHist;
    firingRate(isnan(firingRate)) = 0;
    autocorr = xcorr2(firingRate - mean(reshape(firingRate, 1, numel(firingRate))));
    [X Y] = meshgrid(edges);
    autocorr(sqrt(X.^2 + Y.^2) > arenaDiam) = nan;
    pcolor(edges, edges, autocorr);
    shading interp;
    hold on;
    plot([30 80], [-200 -200], 'LineWidth', 7, 'Color', 'k');
    hold off;
    axis off;
    axis equal tight;
    
    
    
%    if saveFig
        set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
        print('-depsc', sprintf('%s/job%.4d_single_neuron_response.eps', ...
            folder, jobId));
        print('-dpng', sprintf('%s/job%.4d_single_neuron_response.png', ...
            folder, jobId));
        
%    end

    
    figure('Position', [800 200 500 500]);
    subplot(10, 1, 1:9, 'FontSize', fontSize);
    [G crossCorr angles] = cellGridness(firingHist, xedges);
    plot(angles, crossCorr, 'LineWidth', 2);
    xlabel('Rotation angle (deg)');
    ylabel('Correlation');
    axis square tight;
    pos = get(gca, 'Position');
    pos(1) = pos(1) + 0.05*pos(3);
    set(gca, 'Position', pos);


    set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
    print('-depsc', sprintf('%s/job%.4d_single_neuron_autocorr_corr.eps', ...
        folder, jobId));
