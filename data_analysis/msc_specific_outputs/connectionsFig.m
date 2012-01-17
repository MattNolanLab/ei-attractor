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
fontSize = 16;

figure('Position', [800 500 900 420]);

% Plot outgoing connection weights of a neuron
subplot(2, 4, 3, 'FontSize', fontSize);
pcolor(reshape(connections(neuronNum, :), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'B';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);
colorbar;

% Plot input connection weights of a neuron
subplot(2, 4, 4, 'FontSize', fontSize);
pcolor(reshape(connections(:, neuronNum), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'C';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);


% A neuron at the edge of the sheet
neuronNum = 8*sheet_size + 8;

subplot(2, 4, 7, 'FontSize', fontSize);
pcolor(reshape(connections(neuronNum, :), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'D';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);

subplot(2, 4, 8, 'FontSize', fontSize);
pcolor(reshape(connections(:, neuronNum), sheet_size, sheet_size)');
colormap(gca(), 'default');
axis square;
shading flat;
title 'E';
set(gca(), 'XTick', [1 sheet_size]);
set(gca(), 'YTick', [1 sheet_size]);


% ------------------------------------------------------------------------
% print weighing function for different a parameters
% ------------------------------------------------------------------------
connMult = 20;
lambda = 60;
beta = 3/lambda^2; gamma = 1.05*beta;

abs_x = -40:0.01:40;

a = 1;
fun = a*exp(-gamma*abs_x.^2) - exp(-beta*abs_x.^2);
subplot(1, 2, 1, 'FontSize', fontSize);
hold all;
plot(abs_x, connMult * fun, 'LineWidth', 1);

a = 1.01;
fun = a*exp(-gamma*abs_x.^2) - exp(-beta*abs_x.^2);
plot(abs_x, connMult * fun, 'LineWidth', 1);
title 'A';

xlabel('Distance between neurons');
ylabel('Connection weight (nS)');

legend ('a = 1', 'a = 1.01');


set(gcf,'PaperPositionMode','auto');
set(gcf(), 'Renderer', 'painters');
print('-depsc2', '../../thesis/src/fig/connections.eps');