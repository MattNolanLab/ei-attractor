% Analyze optogenetic manipulation data in the Entorhinal cortex
close all;
sample_data = c003_Current_2;
dt = 1e-4;
pA = 1e12;
nA = 1e9;

fontSize = 14;

figure(1);
set(gca, 'FontSize', fontSize);
plot(c001_Time, sample_data*pA);
xlabel('Time (s)');
ylabel('Syn. current (pA)');

subsample_time = [2 10];
subsample_i = subsample_time*1/dt + 1

subsample_data = sample_data(subsample_i(1):subsample_i(2));
subsample_data = subsample_data - mean(subsample_data);

figure(2);
set(gca, 'FontSize', fontSize);
plot(c001_Time(subsample_i(1):subsample_i(2)), subsample_data*pA);
xlim([4.8 4.9]);
xlabel('Time (s)');
ylabel('Syn. current (pA)');

figure(3);
set(gca, 'FontSize', fontSize);
F = 30:2:150;
[Y, F, T, P] = spectrogram(subsample_data*pA, 1024, 940, F, 1e4);
% The following code produces the same result as calling
% spectrogram with no outputs:

%filt = find(P < 1e-1);
%P(filt) = nan;
plot_dB = false;

if (plot_dB)
    P_plot = 10*log10(abs(P));
    power_label = 'Power (dB)';
else
    P_plot = abs(P);
    power_label = 'Power';
end

min_P = min(min(P_plot));
max_P = max(max(P_plot));

surf(T,F,P_plot,'EdgeColor','none');
axis xy; axis tight; %view(0,90);
%zlim([min_P max_P]);
xlabel('Time');
ylabel('Frequency (Hz)');
zlabel(power_label);
colormap jet;
view(-90, -90)
%colorbar;

set(gcf,'PaperPositionMode','auto');
print('-depsc2', 'wide_c003_Current_2_spectrogram_non_dB.eps');


% figure(4);
% for it = 1:numel(T)
%     plot(F, P_plot(:, it));
%     xlabel('Frequency (Hz)');
%     ylabel(power_label);
%     ylim([min_P max_P]);
%     
%     M(it) = getframe(gcf);
% end
% 
% outFile = ['freq_movie_wide_c003_Current_2_non_dB.avi'];
% movie2avi(M, outFile, 'FPS', 10);
% 
