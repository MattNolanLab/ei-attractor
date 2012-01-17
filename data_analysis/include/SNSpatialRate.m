function [firingHist xedges yedges] = SNSpatialRate(neuronSpikes, ratData_pos_x, ratData_pos_y, arenaDiam, h, dt_rat)
    % Plot histogram of firing as a function of 2d arena space
    % Use spatial smoothing algorithm from Hafting et al., 2005

    precision = arenaDiam/h;
    xedges = linspace(-arenaDiam/2, arenaDiam/2, precision + 1);
    yedges = linspace(-arenaDiam/2, arenaDiam/2, precision + 1);

    firingHist = zeros(numel(xedges), numel(yedges)) * nan;

    for x_i = 1:numel(xedges)
        for y_i = 1:numel(xedges)
            x = xedges(x_i);
            y = yedges(y_i);
            isNearTrack = numel(find(sqrt((ratData_pos_x-x).^2 + (ratData_pos_y-y).^2) <= h)) > 0;

            if isNearTrack
                normConst = dt_rat*trapz(gaussianFilter(sqrt((ratData_pos_x-xedges(x_i)).^2 + (ratData_pos_y-yedges(y_i)).^2), h));
                neuronPos_i = floor(neuronSpikes/dt_rat) + 1;  % indices into the position vector
                neuronPos_x = ratData_pos_x(neuronPos_i);
                neuronPos_y = ratData_pos_y(neuronPos_i);

                spikes = sum(gaussianFilter(sqrt((neuronPos_x-xedges(x_i)).^2 + (neuronPos_y - yedges(y_i)).^2), h));
                firingHist(x_i, y_i) = spikes/normConst;
            end
        end
    end

    firingHist = firingHist';
end