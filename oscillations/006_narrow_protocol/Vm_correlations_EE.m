function Vm_correlations_EE(res, it1, it2, ti_start, ti_end, win_len, plot_opt)
    % plot sliding correlation coefficient of two excitatory neurons

    s1 = res.Vmon.e(it1, ti_start:ti_end);
    s2 = res.Vmon.e(it2, ti_start:ti_end);

    subplot(3, 1, 1, 'FontSize', plot_opt.fontSize);
    plot(res.Vmon.t(ti_start:ti_end), [s1; s2]*1000);
    ylabel('Vm (mV)');
    box on;
%        axis tight;
    xlim(plot_opt.x_lim);
    title(sprintf('Excitatory neurons %d:%d', res.opt.Emon_i(it1), res.opt.Emon_i(it2)));

    subplot(3, 1, 2, 'FontSize', plot_opt.fontSize);
    plot(res.Vmon.t(ti_start:ti_end-win_len), slidingCorrelation(s1, s2, win_len));
    %xlabel('Time (s)');
    ylabel('Corr. coefficient');
    box on;
    xlim(plot_opt.x_lim);

    subplot(3, 1, 3, 'FontSize', plot_opt.fontSize);
    plot(res.times(ti_start:ti_end), res.firingRate_i(ti_start:ti_end));
    ylabel('Firing rate (Hz)');
    xlim(plot_opt.x_lim);
    box on;
    xlabel('Time (s)');
    title('Interneurons');

end