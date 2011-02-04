function drawPin(neuronID, sheet_size , color)
% Draw a coloured pin with tip at x,y position

    recN_x = mod(neuronID, sheet_size) + 1;
    recN_y = floor(neuronID/sheet_size) + 1;
    hold on;
    
    %plot([recN_x recN_x+1/5*sheet_size], [recN_y recN_y + 1/5*sheet_size], 'Color', color, 'LineWidth', 2); 
    %plot(recN_x+1/5*sheet_size, recN_y+1/5*sheet_size, 'o', 'Color', color, 'MarkerSize', 6, 'LineWidth', 5);
    quiver(recN_x+1/5*sheet_size, recN_y + 1/5*sheet_size, -1/5*sheet_size, -1/5*sheet_size, 'MaxHeadSize',2, 'Color', color, 'LineWidth', 1);
    %hold off;
end