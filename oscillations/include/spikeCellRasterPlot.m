function [] = spikeCellRasterPlot(spikeCell, markerString)
    % plot raster plot of spike times into the current figure
    for it = 1:size(spikeCell, 2);
        plot(spikeCell{it}, ones(1, size(spikeCell{it}, 2))*it, markerString);
        hold on;
    end
end
