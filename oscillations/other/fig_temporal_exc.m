% print a simple figure of excitatory temporal input pattern
close all;

figure('Position', [800 200 500 250]);
subplot(1, 1, 1, 'FontSize', 16);
t = 0:0.01:10;
plot(t, [t/10; sin(2*pi*t) + 1; zeros(1, numel(t)) + 1], 'LineWidth', 1);
box on;
xlabel('Time (s)');
ylabel('Input current');
set(gca, 'YTick', []);
set(gca, 'XTick', []);

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'fig_temporal_exc.eps'); 