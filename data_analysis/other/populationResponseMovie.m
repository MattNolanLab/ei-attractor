% Create population response movie
close all;

opt = parseOptions(options);

startTime = 0;
endTime = 70;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.


disp 'Creating spike histogram';
spikeHist = createSpikeHistCell(1:numel(spikeCell), spikeCell, ...
    dt_track, 0, endTime);

figure(1);
set(gcf, 'Visible', 'off');

disp 'Printing frames';
it = 1;
for t = startTime:dt_track:endTime
    t
    firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
    image(firingPop, 'Parent', gca);
    axis equal tight;
    colorbar;
    
    M(it) = getframe(gcf);
    
    it = it+1;
end

clear spikeHist;
