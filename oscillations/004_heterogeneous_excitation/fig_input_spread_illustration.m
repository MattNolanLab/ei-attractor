% Input spread illustration

close all;

Ie_max =  900e-12;
sigma_vec = [0.5 1 2 5 10];
thr = 420e-12;

pA = 1e12;

d = -1:0.001:1;

fontSize = 16;
subplot(10,1,1:9, 'FontSize', fontSize);

it = 1;
for sigma = sigma_vec
    plot(d, Ie_max*exp(-d.^2./sigma^2)*pA, 'LineWidth', 2);
    hold all;
    leg{it} = sprintf('\\sigma = %.1f', sigma);
    it = it+1;
end

plot([-1 1], [thr thr]*pA, '--k', 'LineWidth', 1);
xlabel('Relative distance from centre');
ylabel('Excitatory current (pA)');
%ylim([10 Ie_max+0.5]);
legend(leg, 'Location', 'south', 'FontSize', fontSize-2);


set(gcf,'PaperPositionMode','auto', 'Renderer', 'Painters');
print('-depsc2', 'fig_input_spread_illustration.eps');