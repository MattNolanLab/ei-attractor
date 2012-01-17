% Compute gridness score of all files specified by jobNums
path('../include', path);

close all;
clear all;
    

params.jobNums = 2300:2399;
params.chkpntId = 1;
params.folder = '../../../central_data_store/simulation_data/007_mult_bump_path_integ/';

params.preprocess = false;


neuronNum = 1;
h = 5.0;  % cm
arenaDiam = 180;   % cm
dt_rat = 0.02; % sec


if params.preprocess
    [G crossCorr angles] = cellGridnessScoreFromFiles(params, ...
        arenaDiam, h, dt_rat, neuronNum);
    
    save('-v7.3', sprintf('%s/gridnessResults_%d-%d.mat', ...
        params.folder, params.jobNums(1), params.jobNums(end)));
else
    load(sprintf('%s/gridnessResults_%d-%d.mat', ...
        params.folder, params.jobNums(1), params.jobNums(end)));
end

fontSize = 20;

figure('Position', [800 200 500 450]);
subplot(10,1,1:9, 'FontSize', fontSize);

hist(G, 20);
xlabel('Gridness score');
ylabel('Frequency');
colormap([0 0.5 0]);

set(gcf(), 'PaperPositionMode', 'auto');
print('-depsc', sprintf('%s/gridnessHistogram_%d-%d.eps', ...
    params.folder, params.jobNums(1), params.jobNums(end)));


