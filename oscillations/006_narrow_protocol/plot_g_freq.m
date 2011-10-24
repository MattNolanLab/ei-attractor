function [Y F T P P_plot] = plot_g_freq(g_sig, times, which_cell, plot_opt)
    if (strcmp(which_cell, 'stellate'))
        top_title = 'Stellate cell';
    elseif strcmp(which_cell, 'interneuron')
        top_title = 'Interneuron';
    else
        error 'Unknown cell description';
    end

    %gi_sig = - res.Vmon.gi(Ne_it, :)*mA;

    g_ax = subplot(2,1, 1, 'FontSize', plot_opt.fontSize);
    plot(times, g_sig);
    %xlabel('Time (s)');
    set(gca, 'XTick', []);
    ylabel('Syn. current (pA)');
    title(top_title);
    xlim(plot_opt.x_lim);


    %set(gca, 'FontSize', fontSize);
    F = plot_opt.F;
    sampling_rate = plot_opt.sampling_rate;
    win_len = plot_opt.win_len;
    noverlap = win_len /2;
    [Y, F, T, P] = spectrogram(g_sig, win_len, noverlap, F, sampling_rate);
    % The following code produces the same result as calling
    % spectrogram with no outputs:

    %filt = find(P < 1e-1);
    %P(filt) = nan;
    plot_dB = plot_opt.plot_dB;

    if (plot_dB)
        P_plot = 10*log10(abs(P));
        power_label = 'Power (dB)';
    else
        P_plot = abs(P);
        power_label = 'Power';
    end

    min_P = min(min(P_plot));
    max_P = max(max(P_plot));


    spec_ax = subplot(2, 1, 2, 'FontSize', plot_opt.fontSize);
    surf(T,F,P_plot,'EdgeColor','none');
    axis xy; axis tight; %view(0,90);
    %zlim([min_P max_P]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    zlabel(power_label);
    colormap jet;
    view(0, 90)
    colorbar('Location', 'EastOutside');
    xlim(plot_opt.x_lim);

    % [x y width height]
    pos_g = get(g_ax, 'Position');
    pos_spec = get(spec_ax, 'Position');

    gap = 0.03;
    pos_g(3) = pos_spec(3);
    pos_g(2) = pos_spec(2) + pos_spec(4) + gap;
    set(g_ax, 'Position', pos_g);
    set(spec_ax, 'Position', pos_spec);
    
    set(gcf, 'Renderer', 'painters');
end