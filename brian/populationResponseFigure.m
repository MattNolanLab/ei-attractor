% Script, not function!, to plot population response

% Print coloured 2d histogram of the population response
%opt = parseOptions(options);
%sheet_size = opt.sheet_size;

%sheet_size = double(sheet_size);
%dt_rat = 0.02; % sec
%delta_t = 0.25; % sec

firingPop = zeros(sheet_size, sheet_size);

for x_i = 0:(sheet_size-1)
    for y_i = 0:(sheet_size-1)
        neuronID = y_i*sheet_size + x_i;
        neuronSpikes = eval(['spikeMonitor_times_n' num2str(neuronID)]);
        firingRate = computeFiringRate(neuronSpikes, firingPop_t, firingPop_t, dt_rat, delta_t);

        firingPop(x_i+1, y_i+1) = firingRate(numel(firingRate));
    end
end

firingPop = firingPop';

%histmat = hist2(id_x, id_y, xedges, yedges);
%figure(3);
hold off;
pcolor(1:sheet_size,1:sheet_size,firingPop);
axis square tight;
shading flat;
set(gca(), 'XTick', [1 sheet_size], 'YTick', [1 sheet_size]);

%SNList_nID = 38;
drawPin(double(SNList(SNList_nID)), sheet_size, [1 1 0]);
