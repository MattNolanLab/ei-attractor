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
    m_i = max(neuronPos_i);

    %fig = figure;
    plot(ratData_pos_x(1:m_i), ratData_pos_y(1:m_i));
    hold on;
    plot(neuronPos_x, neuronPos_y, '.r', 'MarkerSize', 10);
    plot([30 80], [-100 -100], 'LineWidth', 7, 'Color', 'k');
    hold off;
    axis off;
    axis equal tight;
end
