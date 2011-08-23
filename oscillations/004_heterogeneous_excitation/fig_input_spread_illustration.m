% Input spread illustration

close all;

Ie_max =  24;
sigma_vec = [0.5 1 2 5 10];
thr = 18.5;

d = -1:0.001:1;

fontSize = 16;
subplot(1,1,1, 'FontSize', fontSize);

it = 1;
for sigma = sigma_vec
    plot(d, Ie_max*exp(-d.^2./sigma^2), 'LineWidth', 2);
    hold all;
    leg{it} = sprintf('\\sigma = %.1f', sigma);
    it = it+1;
end

plot([-1 1], [thr thr], '--k', 'LineWidth', 1);
xlabel('Relative distance from centre');
ylabel('Input drive (mV)');
ylim([10 Ie_max+0.5]);
legend(leg, 'Location', 'south', 'FontSize', fontSize-2);
