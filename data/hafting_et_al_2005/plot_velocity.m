% Print velocity in X and Y direction
close all
clear all
load Hafting_Fig2c_Trial1_preprocessed
%load rat_data_original


folder = 'output';

dt = 0.02;
fs = 1/dt;
stopf = 0.5; % Hz
filtord = 4;

remove_1st_4s = true;

if remove_1st_4s
    pos_x(1:4/dt) = [];
    pos_y(1:4/dt) = [];
    pos_timeStamps(end-4/dt+1:end) = [];
end

[B, A] = butter(filtord, stopf / (fs/2));

fontSize = 14;
histBins = 30;

totalTime = 1200; %s
pos_x = pos_x(1:totalTime/dt);
pos_y = pos_y(1:totalTime/dt);
pos_timeStamps = pos_timeStamps(1:totalTime/dt);

vel_x = diff(pos_x)/dt;
vel_y = diff(pos_y)/dt;


% Original hafting velocities
figure
subplot(2, 1, 1, 'FontSize', fontSize);
plot(pos_timeStamps(1:end-1), vel_x)
xlabel('Time (s)');
ylabel('Rat velocity in X (cm/s)');
legend(['Mean: ' num2str(mean(abs(vel_x)))]);
%xlim([0 300]);

subplot(2,1,2, 'FontSize', fontSize);
plot(pos_timeStamps(1:end-1), vel_y)
xlabel('Time (s)');
ylabel('Rat velocity in Y (cm/s)');
legend(['Mean: ' num2str(mean(abs(vel_y)))]);
%xlim([0 300]);

max(abs(vel_x))
max(abs(vel_y))

set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
print('-depsc', sprintf('%s/velocities_original.eps', ...
    folder));


figure
subplot(1, 1, 1, 'FontSize', fontSize);
plot(pos_x, pos_y);
xlabel('X (cm)')
ylabel('Y (cm)')
title('Rat trajectory: original');

set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
print('-depsc', sprintf('%s/trajectory_original.eps', ...
    folder));


figure
title('Velocity distribution: original');
subplot(1, 2, 1, 'FontSize', fontSize);
hist(vel_x, histBins);
xlabel('X Velocity (cm/s)');
ylabel('Count');
%xlim([-150 150]);
subplot(1, 2, 2, 'FontSize', fontSize);
hist(vel_y, histBins);
xlabel('Y velocity (cm/s)');
ylabel('Count');
%xlim([-150 150]);

set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
print('-depsc', sprintf('%s/velocities_hist_original.eps', ...
    folder));



% Low pass filtered
vel_x_new = filtfilt(B, A, vel_x);
vel_y_new = filtfilt(B, A, vel_y);
pos_x_new = pos_x(1) + cumsum(vel_x_new)*dt;
pos_y_new = pos_y(1) + cumsum(vel_y_new)*dt;

figure
title('Velocity clipped');
subplot(2, 1, 1, 'FontSize', fontSize);
plot(pos_timeStamps(1:end-1), vel_x_new)
%xlim([0 300]);
xlabel('Time (s)');
ylabel('Rat velocity in X (cm/s)');
legend(['Mean: ' num2str(mean(abs(vel_x_new)))]);


subplot(2,1,2, 'FontSize', fontSize);
plot(pos_timeStamps(1:end-1), vel_y_new)
xlabel('Time (s)');
ylabel('Rat velocity in Y (cm/s)');
legend(['Mean: ' num2str(mean(abs(vel_y_new)))]);

%xlim([0 300]);

set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
print('-depsc', sprintf('%s/velocities_lowpass.eps', ...
    folder));


figure
subplot(1, 1, 1, 'FontSize', fontSize);
plot(pos_x_new, pos_y_new);
xlabel('X (cm)')
ylabel('Y (cm)')
title('Rat trajectory: low pass filtered');

set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
print('-depsc', sprintf('%s/trajectory_lowpass.eps', ...
    folder));


figure
title('Velocity distribution: low pass filtered');
subplot(1, 2, 1, 'FontSize', fontSize);
hist(vel_x_new, histBins);
xlabel('X Velocity (cm/s)');
ylabel('Count');
xlim([-150 150]);
subplot(1, 2, 2, 'FontSize', fontSize);
hist(vel_y_new, histBins);
xlabel('Y velocity (cm/s)');
ylabel('Count');
xlim([-150 150]);

set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
print('-depsc', sprintf('%s/velocities_hist_lowpass.eps', ...
    folder));

% % Save original data
% save('rat_trajectory_original.mat', 'pos_timeStamps', 'pos_x', 'pos_y', 'dt');
% 
% % Save low pass filtered data
% pos_timeStamps = pos_timeStamps(1:end-1);
% pos_x = pos_x_new;
% pos_y = pos_y_new;
% save('rat_trajectory_lowpass.mat', 'pos_timeStamps', 'pos_x', 'pos_y', 'dt')
