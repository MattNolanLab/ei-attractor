% Print rat trajectory and instantaneous velocity of the rat

close all;

dt_rat = 0.02; % sec

% ------------------------------------------------------------------------
% Plot the original trajectory of the rat
load '../../data/hafting_et_al_2005/original/Hafting_Fig2c_Trial1.mat'

figure('Position', [600 500 10150 420]);

fontSize = 19;
subplot(1, 2, 1, 'FontSize', fontSize);
plot(pos_x, pos_y, 'k');
axis square;
xlabel('Horizontal rat position (cm)');
ylabel('Vertical rat position (cm)');
title('A');


subplot(2, 2, 2, 'FontSize', fontSize);
pos_x_diff = pos_x(1:end-1) - pos_x(2:end);
pos_y_diff = pos_y(1:end-1) - pos_y(2:end);
pos_timeStamps_diff = pos_timeStamps(2:end) - pos_timeStamps(1:end-1);

speed = sqrt(pos_x_diff.^2 + pos_y_diff.^2) ./ pos_timeStamps_diff;

plot(pos_timeStamps(1:end-1), speed, 'k');
xlim([0 1200]);

ylabel('Rat speed (cm/s)');
title('B');


clear pos_x pos_y pos_timeStamps
load '../../data/hafting_et_al_2005/Hafting_Fig2c_Trial1_preprocessed.mat'

subplot(2, 2, 4, 'FontSize', fontSize);
pos_x_diff = pos_x(1:end-1) - pos_x(2:end);
pos_y_diff = pos_y(1:end-1) - pos_y(2:end);
pos_timeStamps_diff = pos_timeStamps(2:end) - pos_timeStamps(1:end-1);

speed = sqrt(pos_x_diff.^2 + pos_y_diff.^2) ./ pos_timeStamps_diff;

plot(pos_timeStamps(1:1200/dt_rat), speed(1:1200/dt_rat), 'k');
xlim([0 1200]);

xlabel('Time (s)');
ylabel('Rat speed (cm/s)');
title('');

set(gcf,'PaperPositionMode','auto');
print('-depsc2', '../../thesis/src/fig/ratTrajectory.eps');