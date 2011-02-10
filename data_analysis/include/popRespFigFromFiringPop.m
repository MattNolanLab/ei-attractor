function popRespFigFromFiringPop(firingPop, sheet_size, pinID)
    % Draw a simple pcolor graphs into current axes
    
    hold off;
    pcolor(1:sheet_size,1:sheet_size,firingPop);
    axis square tight;
    shading flat;
    set(gca(), 'XTick', [1 sheet_size], 'YTick', [1 sheet_size]);

    
    if (pinID ~= false)
        drawPin(pinID, sheet_size, [1 1 0]);
    end
end