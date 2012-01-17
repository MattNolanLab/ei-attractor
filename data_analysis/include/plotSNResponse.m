function [firingHist xedges yedges] = plotSNResponse(neuronSpikes, ratData_pos_x, ratData_pos_y, arenaDiam, h, dt_rat, neuronID)
    % Plot histogram of firing as a function of 2d arena space
    % Use spatial smoothing algorithm from Hafting et al., 2005

    [firingHist xedges yedges] = SNSpatialRate(neuronSpikes, ...
        ratData_pos_x, ratData_pos_y, arenaDiam, h, dt_rat);
    pcolor(xedges,yedges,firingHist);
    hold on;
    plot([30 80], [-100 -100], 'LineWidth', 7, 'Color', 'k');
    hold off;
    axis equal;
    axis off;    %colorbar;
    shading interp;
    colormap(gca, 'jet');
end