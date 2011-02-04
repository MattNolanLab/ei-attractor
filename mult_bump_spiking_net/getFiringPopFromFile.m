% Load a file and compute population firing rate at specified time t
function [firingPop, opt] = getFiringPopFromFile(fileName, t, delta_t)

    [start_i end_i] = regexp(fileName, 'job\d+_\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\d_', 'Start', 'End');
    fileBase = fileName(start_i:end_i-1);
    
    load(fileName);    

    opt = parseOptions(options);
    sheet_size = opt.sheet_size;
    
    firingPop = zeros(sheet_size, sheet_size);
    
    for x_i = 0:(sheet_size-1)
        for y_i = 0:(sheet_size-1)
            neuronID = y_i*sheet_size + x_i;
            neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);          
            firingPop(x_i+1, y_i+1) = numel(find((neuronSpikes > t - delta_t/2) & (neuronSpikes < t + delta_t/2)));
        end
    end
    
    firingPop = firingPop';
end
