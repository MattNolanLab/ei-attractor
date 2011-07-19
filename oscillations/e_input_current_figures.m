% Process input current simulation experiment
close all;

%load e_input_current_output_19-Jul-2011;

fontSize = 16;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;
nPar = 21;

for par_it = 1:nPar
    res = results(par_it, trial_it);
    
    spikeRecord_e = res.spikeRecord_e;
    spikeRecord_i = res.spikeRecord_i;
    Vmon = res.Vmon;
    times = res.times;
    firingRate_e = res.firingRate_e;
    spikeCell_e = res.spikeCell_e;
    spikeCell_i = res.spikeCell_i;
    opt = res.opt;
    
    
    figure('Position', [800 528 1200 800]);
    
    x_lim = [2 2.5];

    % Plot the results

    subplot(5, 1, 1);
    spikeCellRasterPlot(spikeCell_e, '.');
    title('Pyramidal neurons');
    ylabel('Neuron number');
    xlim(x_lim);
    ylim([1 opt.Ne]);

    subplot(5, 1, 2);
    plot(times, firingRate_e);
    ylabel('Firing rate (Hz)');
    xlim(x_lim);


    subplot(5, 1, 3);
    plot(Vmon.t, Vmon.e*1000);
    %title('Pyramidal neuron');
    ylabel('Vm (mV)');
    xlim(x_lim);


    subplot(5, 1, 4);
    spikeCellRasterPlot(spikeCell_i, '.');
    title('Interneurons');
    ylabel('Neuron number');
    xlim(x_lim);
    ylim([1 opt.Ni]);


    subplot(5, 1, 5);
    plot(Vmon.t, Vmon.i*1000);
    %title('Interneuron');
    ylabel('Vm (mV)');
    xlabel('Time (s)');
    xlim(x_lim);

    hold off;
    
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', ['output/2011-07-19/e_input_current_' num2str(par_it) '.eps']);



    %
    % Plot firing rate fft
    %
%     figure();
%     [Y f NFFT] = fourierTrans(firingRate_e, dt);
%     Y_abs = 2*abs(Y(1:NFFT/2+1));
%     plot(f,Y_abs);
%     xlim([0 200]);
end

t_start = 1.5;
t_end   = 2.5;

for par_it = 1:nPar
    for trial_it = 1:nTrials
        res = results(par_it, trial_it);
        
        t_start_i = t_start/res.opt.dt + 1;
        t_end_i   = t_end/res.opt.dt + 1;
    
        firingRate_e = res.firingRate_e(t_start_i:t_end_i);
        %spikeCell_e = res.spikeCell_e;
        %spikeCell_i = res.spikeCell_i;
        opt = res.opt;

        
        % Population frequency
        [Y f NFFT] = fourierTrans(firingRate_e, dt);
        Y_abs = 2*abs(Y(1:NFFT/2+1));
        Y_abs = Y_abs.^2;
%         figure('Visible', 'off');
%         plot(f,Y_abs);
%         xlim([0 200]);
%         xlabel('Frequency (Hz)');
%         ylabel('Power (Hz^2)');
%         title(['Excitatory population power spectrum. Ie = ' num2str(opt.Ie) 'mV, trial: ' num2str(trial_it)]);
%         
%         set(gcf,'PaperPositionMode','auto');
%         print('-depsc2', sprintf('output/2011-07-19/e_input_current_power_spectrum_par_it_%.4d_trial_%.4d', par_it, trial_it));

        
        [maxF maxFI] = max(Y_abs);
        fmax(trial_it, par_it) = f(maxFI);

        
        spikeRecord_e = res.spikeRecord_e(:, t_start_i:t_end_i);
        spikeRecord_i = res.spikeRecord_i(:, t_start_i:t_end_i);        
        times = res.times(t_start_i:t_end_i);

        spikeCnt_i = sum(spikeRecord_i');
        
        % mean firing rates of neurons in this trial
        mfr_T = t_end - t_start;
        e_mfr(trial_it, par_it) = mean(sum(spikeRecord_e'))/mfr_T;
        i_mfr(trial_it, par_it) = mean(sum(spikeRecord_i'))/mfr_T;
    end
end


% Print the population and excitatory cells frequency depending on input
% parameter
for par_it = 1:nPar
    Ie(par_it) = results(par_it, 1).opt.Ie;
end

figure('Position', [840 800 800 500]);
subplot(1, 1, 1, 'FontSize', fontSize);
hold on;
errorbar([Ie*1000; Ie*1000; Ie*1000]', ...
    [mean(fmax); mean(e_mfr); mean(i_mfr)]', ...
    [std(fmax); std(e_mfr); std(i_mfr)]', ...
    'LineWidth', 2);
%errorbar(Ie*1000, mean(e_mfr), std(e_mfr));
xlabel('Input drive (mV)');
ylabel('Frequency (Hz)');
legend('Population', 'E cells', 'I cells', 'Location', 'SouthEast');
axis tight;

