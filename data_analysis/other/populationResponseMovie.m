% Create population response movie
path('../include', path);

startTime = 0;
endTime = 10;

dt_track = 0.1;
delta_t = 0.25; % Should be this value.

folder = 'simulation_data/008_BlobSpacing/';
jobNums = 80000:80019;


nFiles = numel(jobNums);    
for f_it = 1:nFiles
    
    close all;
    clear M;

    jobNums(f_it)
    d = dir([folder 'job' num2str(jobNums(f_it)) '*.mat'])
    if (numel(d) == 0)
        warning(['Could not open job no.' num2str(jobNums(f_it))]);
        continue;
    end
    % Load file which contains the tracking data and extract the spike
    % histogram and blob positions
    d(end).name
    %clear spikeCell
    dataLoad = load([folder d(end).name], 'spikeCell');       
    spikeCell = dataLoad.spikeCell;
    %opts = parseOptions(dataLoad.options);

    disp 'Creating spike histogram';
    spikeHist = createSpikeHistCell(1:numel(spikeCell), spikeCell, ...
        dt_track, 0, endTime);

    figure(2);
    %set(gcf, 'Visible', 'off');

    disp 'Printing frames';
    it = 1;
    for t = startTime:dt_track:endTime
        t
        firingPop = getFiringPop(spikeHist, t, dt_track, delta_t);
        firingPop = reshape(firingPop, sheet_size, sheet_size)';
        image(firingPop, 'Parent', gca);
        axis equal tight;
        colorbar;

        M(it) = getframe(gcf);

        it = it+1;
    end

    outFile = ['output/008_BlobSpacing/job' num2str(jobNums(f_it)) '_movie.avi'];
    movie2avi(M, outFile, 'FPS', 10);
    
    clear spikeHist;
end