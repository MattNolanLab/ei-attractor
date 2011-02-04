% Print the spike and SN rate plot produced by Hafting et al. data

close all;
load '../../data/hafting_et_al_2005/Hafting_Fig2c_Trial1_SN_thesis.mat';

dt_rat = 0.02; % sec
h = 5.0;  % cm
arenaDiam = 180;   % cm
neuronSpikes = rat11015_t5c1_timeStamps;
neuronNum = -10;

fontSize = 14;
figure('Position', [1012 500 640 420]);

subplot(1, 2, 1, 'FontSize', fontSize);
plotSpikes_xy(neuronSpikes, pos_x, pos_y, dt_rat, neuronNum);
xlabel('X (cm)');
ylabel('Y (cm)');
xlim([-arenaDiam/2 arenaDiam/2]);
ylim([-arenaDiam/2 arenaDiam/2]);
set(gca(), 'XTick', [-arenaDiam/2 arenaDiam/2]);
set(gca(), 'YTick', [-arenaDiam/2 arenaDiam/2]);
title 'A';

subplot(1, 2, 2, 'FontSize', fontSize);
plotSNResponse(neuronSpikes, pos_x, pos_y, arenaDiam, h, dt_rat, neuronNum);
xlabel('X (cm)');
ylabel('Y (cm)');
xlim([-arenaDiam/2 arenaDiam/2]);
ylim([-arenaDiam/2 arenaDiam/2]);
set(gca(), 'XTick', [-arenaDiam/2 arenaDiam/2]);
set(gca(), 'YTick', [-arenaDiam/2 arenaDiam/2]);
title 'B';

colorbar('peer',gca(), [0.9399 0.2786 0.03402 0.4786], 'FontSize', fontSize);

set(gcf,'PaperPositionMode','auto');
set(gcf(), 'Renderer', 'painters');
print('-depsc2', '../../thesis/src/fig/haftingSN.eps');