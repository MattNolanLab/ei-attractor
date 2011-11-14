function plot2g_xcorr(sig1, sig2, times, p_o)
    % PLOT2G_XCORR
    % Plot sliding x-correlation of two signals with their power spectra
    % and power in gamma range
    
    %xcorr = slidingCorrelation(sig1, sig2, p_o.xcorr_win_len);

%     subplot(5, 1, [1 2], 'FontSize', p_o.fontSize);
%     plot(times, [sig1; sig2]);
%     ylabel('Syn. current (pA)');
%     xlim(p_o.x_lim);
%     title(p_o.title);

%     subplot(5, 1, 4:5, 'FontSize', p_o.fontSize);
%     plotSpectrogramIntoAx(sig1, p_o.spec_win_len, p_o.noverlap, p_o.F, p_o.sampling_rate, p_o.plot_dB);
%     xlim(p_o.x_lim);
%     set(gca, 'XTick', []);
%     xlabel('');

    subplot(2, 1, 2, 'FontSize', p_o.fontSize);
    [Y F T P P_plot] = plotSpectrogramIntoAx(sig2, p_o.spec_win_len, p_o.noverlap, p_o.F, p_o.sampling_rate, p_o.plot_dB);
    xlim(p_o.x_lim);

    
    subplot(2, 1, 1, 'FontSize', p_o.fontSize);
    [max_P max_P_i] = max(P);
    
    [ax] = plot(T, max_P, 'r', 'LineWidth', 2);
    ylabel('Power (pA^2/Hz)');
    set(gca, 'XTick', []);
    xlim(p_o.x_lim);
    
%     set(get(ax(1), 'YLabel'), 'String', 'Corr. coeff.');
%     set(get(ax(2), 'YLabel'), 'String', 'Power (pA^2)');
%     set(ax(1), 'LineWidth', 1);
%     set(ax(2), 'LineWidth', 1);
%     set(ax(1), 'XTick', [], 'FontSize', p_o.fontSize);
%     set(get(ax(1), 'XLabel'), 'FontSIze', p_o.fontSize);
%     set(get(ax(1), 'YLabel'), 'FontSIze', p_o.fontSize);
%     set(ax(2), 'FontSize', p_o.fontSize);
%     set(get(ax(2), 'YLabel'), 'FontSIze', p_o.fontSize);
%     set(ax(2), 'XTick', []);
%     %axes(ax(2)); axis tight;
%     xlim(ax(1), p_o.x_lim);
%     xlim(ax(2), p_o.x_lim);

    
end