% Print coloured 2d histogram of the population response at different
% time snapshots

function timePopulationPlot(fileName)

close all;


saveDir = 'results/fig/';
[start_i end_i] = regexp(fileName, 'job\d+_\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_', 'Start', 'End');
fileBase = fileName(start_i:end_i-1);

load(fileName);

opt = parseOptions(options);
sheet_size = opt.sheet_size;

for t = 4:6:150;
    disp(['Plotting population response for time t = ' num2str(t)]);
    %sheet_size = double(sheet_size);
    dt_rat = 0.02; % sec
    delta_t = 0.8; % sec
    startTime = t;
    endTime = startTime; % sec

    firingPop = zeros(sheet_size, sheet_size);

    for x_i = 0:(sheet_size-1)
        for y_i = 0:(sheet_size-1)
            neuronID = y_i*sheet_size + x_i;
            neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
            firingRate = computeFiringRate(neuronSpikes, startTime, endTime, dt_rat, delta_t);

            firingPop(x_i+1, y_i+1) = firingRate(numel(firingRate));
        end
    end

    %histmat = hist2(id_x, id_y, xedges, yedges);
    f = figure('Visible', 'off');
    pcolor(0:sheet_size-1,0:sheet_size-1,firingPop');
    xlabel('Neuron number');
    ylabel('Neuron number');
    title(['Population response: t: ' num2str(t)]);

    colorbar;
    axis square tight;
    shading interp;
    
    popPlotFile =  [saveDir fileBase '_timedPopRateMap_t' num2str(t, '%.3d') '.eps'];
    print('-depsc', popPlotFile);
end

end