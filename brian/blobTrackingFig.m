% Illustrate the blob tracking algorithm

close all;


preprocess = false;
%clear all;
%load results/path_integration/sheet_size/job22205_2010-08-05T19-05-54_output.mat

dt_rat = 0.02; % sec
dt_track = 0.1; % sec; dt for the tracking algorithm
delta_t = 0.25; % sec

startTime = 10; t_start = startTime;
endTime = 450; t_end = endTime - delta_t/2;


opt = parseOptions(options);
optStr = ['_s' num2str(opt.sheet_size) '_alpha' num2str(opt.alpha)];   

sheet_size = opt.sheet_size;

if preprocess == true
    createSpikeHist  % a script - side effects
end

t = 100;
firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
[r, c, thrFiringPop, segFiringPop] = trackBlobs(firingPop);

fontSize = 14;
figure('Position', [680 660 1000 280]);

subplot(1, 4, 1, 'FontSize', fontSize);
pcolor(firingPop);
colormap(gca(), 'Jet');
axis square;
shading flat;
set(gca(), 'XTick', []);
set(gca(), 'YTick', []);
title 'A'

subplot(1, 4, 2, 'FontSize', fontSize);
pcolor(thrFiringPop);
axis square;
shading flat;
set(gca(), 'XTick', []);
set(gca(), 'YTick', []);
title 'B'

subplot(1, 4, 3, 'FontSize', fontSize);
pcolor(thrFiringPop);
hold on;
plot(c, r, 'og', 'LineWidth', 3, 'MarkerSize', 2)
axis square;
shading flat;
set(gca(), 'XTick', []);
set(gca(), 'YTick', []);
title 'C'


t = 419;
firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
[r, c, thrFiringPop, segFiringPop] = trackBlobs(firingPop);

subplot(1, 4, 4, 'FontSize', fontSize);
pcolor(firingPop);
hold on;
plot(c, r, 'ok', 'LineWidth', 3, 'MarkerSize', 2)
axis square;
shading flat;
set(gca(), 'XTick', []);
set(gca(), 'YTick', []);
title 'D'

set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
print -depsc2 '../../thesis/src/fig/blobTracking.eps';
