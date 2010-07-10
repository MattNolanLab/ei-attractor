% Print coloured 2d histogram of the population response
%function printPopulationResponse(sheet_size, startTime, endTime)

%close all;

    sheet_size = parseOptions(options);

    %sheet_size = double(sheet_size);
    dt_rat = 0.02; % sec
    delta_t = 0.8; % sec
    startTime = 60;
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
    %figure;
    pcolor(0:sheet_size-1,0:sheet_size-1,firingPop);

    colorbar;
    axis square tight;
    shading interp;
%end