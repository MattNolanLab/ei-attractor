% Plot histogram of firing as a function of 2d arena space

arenaDiam = 180; %cm
precision = 180; % cm

sheet_size = 40;
neuronNum = 800;

neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);

xedges = linspace(-arenaDiam/2, arenaDiam/2, precision);
yedges = linspace(-arenaDiam/2, arenaDiam/2, precision);

% Transform an array of spike times of a SN to an array of x and y
% positions
spike_pos_x = [];
spike_pos_y = [];

spike_time_ind = floor(neuronSpikes / 0.02) + 1;
for ind = 1:numel(neuronSpikes)
    spike_pos_x = [spike_pos_x;ratData_pos_x(spike_time_ind(ind))];
    spike_pos_y = [spike_pos_y;ratData_pos_y(spike_time_ind(ind))];
end

histmat = hist2(spike_pos_x, spike_pos_y, xedges, yedges);
%plot(ratData_pos_x, ratData_pos_y);
%hold all;
pcolor(xedges,yedges,histmat);

colorbar;
axis square tight;
shading interp;
