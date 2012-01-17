function [G crossCorr angles] = cellGridnessScoreFromFiles(params, ...
    arenaDiam, h, dt_rat, neuronID)
    % CELLGRIDNESSSCOREFROMFILES Load files and compute gridness score,
    % rotational Correlation per each file for cell specified by neuronID
    % Each row of G crossCorr and angles contains data of one instance


folder = params.folder; 
jobNums = params.jobNums;
fileId = params.chkpntId;

nFiles = numel(jobNums);

G = zeros(nFiles, 1)*nan;

parfor f_it = 1:nFiles
    jobNums(f_it)
    d = dir([folder 'job' num2str(jobNums(f_it)) '*.mat'])
    if (numel(d) < fileId)
        warning(['Could not open job no.' num2str(jobNums(f_it))]);
        continue;
    end
    % Load file which contains the tracking data and extract the spike
    % histogram and blob positions
    d(fileId).name
    %clear spikeCell
    dataLoad = load([folder d(fileId).name], 'spikeCell', 'pos_x', ...
        'pos_y');       
    spikeCell = dataLoad.spikeCell;
    pos_x = dataLoad.pos_x;
    pos_y = dataLoad.pos_y;
    
    [firingHist xedges yedges] = SNSpatialRate(spikeCell{neuronID}, ...
        pos_x, pos_y, arenaDiam, h, dt_rat);

    [score C ang] = cellGridness(firingHist, xedges);
    G(f_it) = score;
    crossCorr(f_it, :) = C;
    angles(f_it, :) = ang;
    

    fontSize = 18;
    figure('Position', [800 200 700 600], 'Visible', 'off');
    subplot(10, 1, 1:9, 'FontSize', fontSize);
    plot(ang, C, 'LineWidth', 2);
    xlabel('Rotation angle (deg)');
    ylabel('Correlation');
    axis square tight;

    set(gcf(), 'PaperPositionMode', 'auto', 'Renderer', 'painters');
    print('-depsc', sprintf('%s/job%.4d_single_neuron_autocorr_corr.eps', ...
        folder, jobNums(f_it)));

end