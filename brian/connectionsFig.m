% print connections + topography

close all;

sheet_size = sqrt(size(connections, 1));

%subplot(1, 3, 1);
%torus(15, 40, 50);
%axis square;
%title 'A';
%colormap(gca(), 'gray');

%grid off;
%set(gca
%set(gca(), 'Visible', 'off');

neuronNum = sheet_size^2/2 + sheet_size/2;
fontSize = 14;

% Plot outgoing connection weights of a neuron
subplot(2, 2, 1, 'FontSize', fontSize);
pcolor(reshape(connections(neuronNum, :), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'A';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);

% Plot input connection weights of a neuron
subplot(2, 2, 2, 'FontSize', fontSize);
pcolor(reshape(connections(:, neuronNum), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'B';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);


% A neuron at the edge of the sheet
neuronNum = 8*sheet_size + 8;

subplot(2, 2, 3, 'FontSize', fontSize);
pcolor(reshape(connections(neuronNum, :), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'C';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);

subplot(2, 2, 4, 'FontSize', fontSize);
pcolor(reshape(connections(:, neuronNum), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'D';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);


set(gcf,'PaperPositionMode','auto');
set(gcf(), 'Renderer', 'painters');
print('-depsc2', '../../thesis/src/fig/connections.eps');