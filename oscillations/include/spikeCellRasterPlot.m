function [] = spikeCellRasterPlot(spikeCell, markerString, N)
    if nargin < 3
        N = inf;
    end
    
    % plot raster plot of spike times into the current figure
    for it = 1:min(N, size(spikeCell, 2));
        plot(spikeCell{it}, ones(1, size(spikeCell{it}, 2))*it, markerString);
        hold on;
    end
end
