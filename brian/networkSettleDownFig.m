% Print a figure illustrating how the network settles down into an
% attractor state
close all;

%clear all;
load results/static_wave/zero_velocity/job23602_2010-08-10T16-50-34_output.mat

opt = parseOptions(options);
sheet_size = opt.sheet_size;

%preprocess = false;


dt_rat = 0.02; % sec
dt_track = 0.1; % sec; dt for the tracking algorithm
delta_t = 0.25; % sec

startTime = 10; t_start = startTime;
endTime = 450; t_end = endTime - delta_t/2;


% ------------------------------------------------------------------------
% Print network which settles down successfully
% ------------------------------------------------------------------------

%if preprocess == true
    createSpikeHist  % a script - side effects
%end


fontSize = 16;
figure('Position', [680 660 1000 400]);

%t = [3 4 5 6];
t = [3.75 4 4.25 4.5 4.75];

for t_i = 1:numel(t);
    firingPop = getFiringPop(spikeHist, t(t_i), dt_track, delta_t);
    %[r, c, thrFiringPop, segFiringPop] = trackBlobs(firingPop);

    subplot(2, numel(t), t_i, 'FontSize', fontSize);
    pcolor(firingPop);
    %colormap(gca(), 'Jet');
    axis square;
    shading flat;
    set(gca(), 'XTick', []);
    set(gca(), 'YTick', []);
    title(['t = ' num2str(t(t_i)) ' s']);
    
    if (t_i == 1)
        ylabel('A    ' , 'Rotation', 0);
    end
end


% ------------------------------------------------------------------------
% Print network which doesn't settle down
% ------------------------------------------------------------------------
load results/static_wave/zero_velocity/job21600_2010-08-09T18-52-45_output.mat
opt = parseOptions(options);
sheet_size = opt.sheet_size;

createSpikeHist; % a script - side effects

old_t = numel(t);
t = [10 50 75 100 125];
if (numel(t) ~= old_t)
    error 'numbers of columns don''t match'
end

for t_i = 1:numel(t);
    firingPop = getFiringPop(spikeHist, t(t_i), dt_track, delta_t);
    %[r, c, thrFiringPop, segFiringPop] = trackBlobs(firingPop);

    subplot(2, numel(t), t_i + numel(t), 'FontSize', fontSize);
    pcolor(firingPop);
    %colormap(gca(), 'Jet');
    axis square;
    shading flat;
    set(gca(), 'XTick', []);
    set(gca(), 'YTick', []);
    title(['t = ' num2str(t(t_i)) ' s']);
    
    if (t_i == 1)
        ylabel('B    ', 'Rotation', 0);
    end

end

set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
print -depsc2 '../../thesis/src/fig/networkSettleDown.eps';