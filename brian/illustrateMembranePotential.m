% Illustrate how membrane potential changes with the rat moving

close all;

load results/membrane_potential/job23101_2010-08-05T09-54-38_output.mat
%load results/job23101_2010-08-05T09-54-38_output.mat

dt_rat = 0.02;  % sec
dt_track = 0.1; % sec
delta_t = 0.5;  % sec

opt = parseOptions(options);
sheet_size = opt.sheet_size;

endTime = 55; %sec

% Preprocess spiking data: firing time --> histogram of firing, for
% each neuron
edges = linspace(0, endTime, endTime/dt_track + 1);
spikeHist = zeros(sheet_size^2, numel(edges));
for x_i = 0:(sheet_size-1)
    for y_i = 0:(sheet_size-1)
        neuronID = y_i*sheet_size + x_i;            
        neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronID)]);
        e_size = size(neuronSpikes, 1);
        if (e_size == 0)
            spikeHist(neuronID+1, :) = 0;
        else
            spikeHist(neuronID+1, :) = histc(neuronSpikes, edges);
        end
    end
end
    
SNList_nID = 39;
neuronID = double(SNList(SNList_nID));


fontSize = 14;
miliSec = 1000; % ms multiplier
x_limits = [8 16]; %sec

f = figure('Position', [1, 1, 700, 700]);

% Plot membrane potential at the top and firing rate at the bottom
subplot(18, 1, 8:12, 'TickDir', 'out', 'FontSize', fontSize);
hold all;
plot(SNMonitor_times, SNMonitor_values(SNList_nID, :)*miliSec, 'LineWidth', 2);
axis tight;
xlim(x_limits);
box on;
ylabel('V_{membrane} (mV)');

nFiringRate = getNeuronFiringRate(spikeHist(neuronID, :), 0, endTime, dt_track, delta_t);
subplot(18, 1, 15:18, 'TickDir', 'out', 'FontSize', fontSize);
hold all;
plot(0:dt_track:endTime, nFiringRate, 'LineWidth', 2);
axis tight;
xlim(x_limits);
box on;
xlabel('Time (s)');
ylabel('Firing rate (Hz)');


firingPop(:, :, 1) = getFiringPop(spikeHist, 9, dt_track, delta_t);
firingPop(:, :, 2) = getFiringPop(spikeHist, 12, dt_track, delta_t);
firingPop(:, :, 3) = getFiringPop(spikeHist, 15, dt_track, delta_t);

for it = 1:3
    sPlot_i = (it-1)*4;
    subplot(18, 11, sPlot_i + [1 2 3   12 13 14   23 24 25   34 35 36    45 46 47   56 57 58], 'TickDir', 'out', 'FontSize', fontSize);
    pcolor(1:sheet_size,1:sheet_size,firingPop(:, :, it));
    axis square tight;
    shading flat;

    set(gca, 'XTick', [1 sheet_size]);
    set(gca, 'YTick', [1 sheet_size]);
    
    recN_x = mod(neuronID, sheet_size) + 1;
    recN_y = floor(neuronID/sheet_size) + 1;
    hold all;
    plot(recN_x+1/5*sheet_size, recN_y+1/5*sheet_size, '.', 'Color', [0 1 0], 'MarkerSize', 20);
    plot([recN_x recN_x+1/5*sheet_size], [recN_y recN_y + 1/5*sheet_size], 'Color', [0 1 0], 'LineWidth', 2); 
end

set(gcf,'PaperPositionMode','auto');
fileName = '../../thesis/src/fig/membraneFiring.eps';
print('-depsc2', fileName);



