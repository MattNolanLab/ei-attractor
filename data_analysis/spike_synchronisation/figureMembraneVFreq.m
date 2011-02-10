% Plot a huge figure with frequency spectra of membrane potentials of
% neurons - all of the recorded ones in the batch

close all;

opt = parseOptions(options);
sheet_size = opt.sheet_size;

N = numel(SNList); % Number of recorded neurons
%N = 9;

% Monitor dt - can be different from sim_dt
dt = 0.0001; %s


mVolt = 1000;
nS = 1e9;
    
x_limits = [7 10];
y_limits = [-61 -53];
start_i = x_limits(1)/dt + 1;
end_i = x_limits(2)/dt + 1;

    
fontSize = 10;
figure('Position', [0, 0, 1680, 3000], 'Visible', 'off');
plotEdgeC = 5;
plotEdgeR = ceil(N/plotEdgeC);

for list_id = 1:N
    list_id
    signal = SNMonitor_values(list_id, start_i:end_i);
    
    [Y, f, NFFT] = fourierTrans(signal, dt);

    % Plot single-sided amplitude spectrum.
    subplot(plotEdgeR, plotEdgeC, list_id, 'FontSize', fontSize, 'XMinorTick', 'on');
    plot(f,2*abs(Y(1:NFFT/2+1))*1000);
    axis tight;
    xlim([0 40]);
    %xlabel('');
    %ylabel('Power (\muW)');
    title(['N\_' int2str(SNList(list_id)) '(' int2str(list_id) '); X: [Hz], Y: [\muW]']);
    box on;
end

set(gcf,'PaperPositionMode','auto');
fileName = 'membraneVFreq.eps'
print('-depsc2', ['fig/' fileName]);