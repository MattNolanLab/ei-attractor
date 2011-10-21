% Convert all files specified by jobNums into spikeCell representation
jobNums = 21900:21999;
dataFolder = '../../../central_data_store/simulation_data/000_002_NoPrefDirs_RandomInit/';

for jobNum = jobNums
    d = dir([dataFolder 'job' num2str(jobNum) '*.mat'])
    % Load file which contains the tracking data and extract the spike
    % histogram and blob positions
    d(end).name
    load([dataFolder d(end).name]);      
    opt = parseOptions(options);
    
    createSpikeCell;
    clear spikeMonitor_times*;
    
    save([dataFolder d(end).name]);
end
