% print torus to figure

close all;
subplot(1, 1, 1);
torus(20, 50, 50);

set(gca(), 'Visible', 'off');
set(gcf,'PaperPositionMode','auto');
print('-deps2', '../../thesis/src/fig/torus.eps');