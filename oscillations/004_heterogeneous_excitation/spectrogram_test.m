% Compute spectrogram of population firing rate of neurons

par_it = 17;
trial_it = 25;
res = results(par_it, trial_it);
signal = results(par_it, trial_it).firingRate_e;

fontSize = 16;
subplot(1,1,1, 'FontSize', fontSize);

win = 512;
noverlap = 500;
F = 20:2:150;
[Y, F, T, P] = spectrogram(signal, win, noverlap, F, fix(1/res.opt.dt), 'yaxis');
% The following code produces the same result as calling
% spectrogram with no outputs:

%filt = find(P < 1e-1);
%P(filt) = nan;
plot_dB = true;

if (plot_dB)
    P_plot = 10*log10(abs(P));
    power_label = 'Power (dB)';
else
    P_plot = abs(P);
    power_label = 'Power';
end

min_P = min(min(P_plot));
max_P = max(max(P_plot));

surf(T,F,P_plot, 'EdgeColor', 'none');
axis xy; axis tight; %view(0,90);
%zlim([min_P max_P]);
xlabel('Time');
ylabel('Frequency (Hz)');
%zlabel(power_label);
colormap jet;
view(0, 90)
colorbar;
