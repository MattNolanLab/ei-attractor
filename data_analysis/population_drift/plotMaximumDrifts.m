function [absMax_r, absMax_c] = plotMaximumDrifts(ax, params, blobTracks_r, blobTracks_c)
    % PLOTMAXIMUMDRIFTS plots maximum drifts given by blobTracks variables into axis ax, parameterized by params input

    max_r = max(blobTracks_r);
    max_c = max(blobTracks_c);
    min_r = min(blobTracks_r);
    min_c = min(blobTracks_c);
    
    max_min_r = [max_r; min_r];
    max_min_c = [max_c; min_c];
    
    [absMax_r absMax_ri] = max(abs(max_min_r));
    [absMax_c absMax_ci] = max(abs(max_min_c));
    
    [absMax_rr absMax_rri] = max(absMax_c);
    
    
    for it = 1:size(blobTracks_r, 2)
        hold on;
        plot(ax, max_min_c(absMax_ci(it), it), max_min_r(absMax_ri(it), it), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 2);
    end
    xlabel('X drift (neurons)');
    ylabel('Y drift (neurons)');
    
    title(params.title);
    axis equal;
    xlim(params.xlim);
    ylim(params.ylim);
end
