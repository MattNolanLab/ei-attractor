% Print statistics about the data

close all;
clear all;

% dir = 'results/fig/';
% params = 's64_alpha0.01';
% spikePlotFile = [dir 'spikeMap_' params '.eps'];
% ratePlotFile  = [dir 'rateMap_'  params '.eps'];
loadDir = 'results/';
saveDir = 'results/fig/';
fileBase = '2010-07-02T22-04-45';

loadFile = [loadDir fileBase '_output.mat'];
load(loadFile);

[opt.sheet_size opt.alpha] = parseOptions(options);
optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];
spikePlotFile = [saveDir fileBase optStr '_spikeMap.eps'];
ratePlotFile =  [saveDir fileBase optStr '_rateMap.eps'];



saveFig = true;
%saveFig = false;


dt_rat = 0.02; % sec
endTime = 1200 - dt_rat; % sec
delta_t = 1; % sec

h = 5.0;  % cm
arenaDiam = 180;   % cm

for i = 801
    neuronNum = i-1;

    neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);
    %neuronSpikes = rat11015_t5c3_timeStamps;

    figure(1);
    plotSpikes_xy(neuronSpikes, ratData_pos_x, ratData_pos_y, dt_rat, neuronNum);

    figure(2);
    plotSNResponse(neuronSpikes, ratData_pos_x, ratData_pos_y, arenaDiam, h, dt_rat, neuronNum);
    
    if saveFig
        figure(1);
        print('-depsc', spikePlotFile);
        figure(2);
        print('-depsc', ratePlotFile);
    end
end
