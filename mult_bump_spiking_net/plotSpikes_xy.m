function plotSpikes_xy(neuronSpikes, ratData_pos_x, ratData_pos_y, dt_rat, neuronID)
    % Plot the spike positions of a single neuron
    
    %dt_rat = 0.02; % sec
    %endTime = 1200 - dt_rat; % sec

    %arenaDiam = 180;   % cm
    %precision = 180/h; % cm

    %neuronNum = 800;
    %neuronSpikes = eval(['spikeMonitor_times_n' int2str(neuronNum)]);


    neuronPos_i = floor(neuronSpikes/dt_rat) + 1;  % indices into the position vector
    neuronPos_x = ratData_pos_x(neuronPos_i);
    neuronPos_y = ratData_pos_y(neuronPos_i);

    %fig = figure;
    plot(ratData_pos_x, ratData_pos_y);
    hold all;
    plot(neuronPos_x, neuronPos_y, '.r', 'MarkerSize', 6);
    hold off;
    axis square tight;
end
