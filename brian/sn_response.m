% Plot histogram of firing as a function of 2d arena space

dt_rat = 0.02; % sec
endTime = 1200 - dt_rat; % sec
delta_t = 1; % sec

arenaDiam = 180; %cm
precision = 180; % cm

%sheet_size = 40;
neuronNum = 800;

neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);

xedges = linspace(-arenaDiam/2, arenaDiam/2, precision);
yedges = linspace(-arenaDiam/2, arenaDiam/2, precision);

% Transform an array of spike times of a SN to an array of x and y
% positions
firingRate = computeFiringRate(neuronSpikes, endTime, dt_rat, delta_t);

firingHist = zeros(numel(xedges), numel(yedges));
% for each firing rate find indices in the histogram of x and y positions
for t_i = 1:numel(firingRate)
    x_ind = find(histc(ratData_pos_x(t_i), xedges) == 1);
    y_ind = find(histc(ratData_pos_y(t_i), yedges) == 1);
    
    if firingHist(x_ind, y_ind) < firingRate(t_i)
        firingHist(x_ind, y_ind) = firingRate(t_i);
    end
end

%histmat = hist2(spike_pos_x, spike_pos_y, xedges, yedges);
%plot(ratData_pos_x, ratData_pos_y);
%hold all;
pcolor(xedges,yedges,firingHist');

colorbar;
axis square tight;
shading interp;

xlabel('Rat position [cm]');
ylabel('Rat position [cm]');
