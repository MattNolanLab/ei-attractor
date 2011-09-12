function Vm_correlations_II(res, it1, it2, ti_start, ti_end, win_len, p_o)

    s1 = res.Vmon.i(it1, ti_start:ti_end);
    s2 = res.Vmon.i(it2, ti_start:ti_end);

    subplot(3, 1, 1, 'FontSize', p_o.fontSize);
    plot(res.Vmon.t(ti_start:ti_end), [s1; s2]*1000);
    ylabel('Vm (mV)');
    box on;
    axis tight;
    xlim(p_o.x_lim);
    title(sprintf('Interneurons %d:%d', res.opt.Imon_i(it1), res.opt.Imon_i(it2)));

    subplot(3, 1, 2, 'FontSize', p_o.fontSize);
    plot(res.Vmon.t(ti_start:ti_end-win_len), slidingCorrelation(s1, s2, win_len));
    %xlabel('Time (s)');
    ylabel('Corr. coefficient');
    box on;
    xlim(p_o.x_lim);

    subplot(3, 1, 3, 'FontSize', p_o.fontSize);
    plot(res.times(ti_start:ti_end), res.firingRate_i(ti_start:ti_end));
    ylabel('Firing rate (Hz)');
    xlim(p_o.x_lim);
    box on;
    xlabel('Time (s)');

end