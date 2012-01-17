% Print coloured 2d histogram of the population response
% function printPopulationResponse(sheet_size, startTime, endTime)

    close all;
    path('../include', path);

    path('../include', path);

    opt = parseOptions(options);
    sheet_size = opt.sheet_size;

    %sheet_size = double(sheet_size);
    dt_rat = 0.02; % sec
    delta_t = 0.25; % sec
    startTime = 4;
    endTime = startTime; % sec
    
    firingPop = zeros(sheet_size, sheet_size);
    
    for x_i = 0:(sheet_size-1)
        for y_i = 0:(sheet_size-1)
            neuronID = y_i*sheet_size + x_i;
            neuronSpikes = spikeCell{neuronID + 1}; %eval(['spikeMonitor_times_n' num2str(neuronID)]);
            firingRate = computeFiringRate(neuronSpikes, startTime, endTime, dt_rat, delta_t);
            
            firingPop(x_i+1, y_i+1) = firingRate(numel(firingRate));
        end
    end
    
    firingPop = firingPop';

    %histmat = hist2(id_x, id_y, xedges, yedges);
    figure(3);
    hold off;
    pcolor(0:sheet_size-1,0:sheet_size-1,firingPop);
    axis square tight;
    colorbar;
    shading flat;
    colormap(jet);
    
    
    % Print the population response vector 3D
    figure
    surf(firingPop);
    shading interp
    view(-52, 80);
    axis off;
