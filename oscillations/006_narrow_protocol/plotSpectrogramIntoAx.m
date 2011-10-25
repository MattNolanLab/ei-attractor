function [Y F T P P_Plot] = ...
    plotSpectrogramIntoAx(sig, win_len, noverlap, F, sampling_rate, plot_dB)

    % Stellate cell
    [Y, F, T, P] = spectrogram(sig, win_len, noverlap, F, sampling_rate);
    % The following code produces the same result as calling
    % spectrogram with no outputs:


    if (plot_dB)
        P_plot = 10*log10(abs(P));
        power_label = 'Power (dB)';
    else
        P_plot = abs(P);
        power_label = 'Power';
    end

    surf(T,F,P_plot,'EdgeColor','none');
    axis xy; axis tight; %view(0,90);
    %zlim([min_P max_P]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    zlabel(power_label);
    colormap jet;
    view(0, 90)
end