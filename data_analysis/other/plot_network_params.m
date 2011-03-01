a = 1;
connMult = 20;

abs_x = -60:0.01:60;


fontSize = 14;
%close all;
figure('Position', [878 452 600 350]);

subplot(1, 1, 1, 'FontSize', fontSize);

lam = [30];

for l_it = 1:numel(lam)
    lambda = lam(l_it);
    beta = 3/lambda^2; gamma = 1.05*beta;

    fun = a*exp(-gamma*abs_x.^2) - exp(-beta*abs_x.^2);
    hold all;
    plot(abs_x, connMult * fun, 'LineWidth', 2);
end

xlabel('Distance between neurons');
ylabel('Connection strength (nS)');

legend(...
    ['\lambda_{net} = ' num2str(lam(1))], ...
    ['\lambda_{net} = ' num2str(lam(2))]);

legend1 = legend(gca(),'show');
set(legend1,'Position',[0.6457 0.1672 0.2312 0.2038]);



%set(gcf,'PaperPositionMode','auto');
%print('-depsc2', '../../thesis/src/fig/connectivityFunction.eps');