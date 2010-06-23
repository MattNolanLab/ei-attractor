% Print coloured 2d histogram of the population response
%function printPopulationResponse(sheet_size, startTime, endTime)

    sheet_size = 40;
    startTime = 1100;
    endTime = startTime+1;
    
    xedges = 0:sheet_size-1;
    yedges = 0:sheet_size-1;
    
    id_x = [];
    id_y = [];

    for neuronID = 0:sheet_size^2 - 1
        neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
        ids = find(neuronSpikes >= startTime & neuronSpikes <= endTime);
        
        id_x = [id_x (ids*0 + mod(neuronID, sheet_size))'];
        id_y = [id_y (ids*0 + floor(neuronID/sheet_size))'];
    end

    histmat = hist2(id_x, id_y, xedges, yedges);
    pcolor(xedges,yedges,histmat);

    colorbar;
    axis square tight;
    shading interp;
%end