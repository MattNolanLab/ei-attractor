% Process input current simulation experiment
close all;

%load e_input_current_output_19-Jul-2011;

nParam  = size(results, 1);
nTrials = size(results, 2);

trial_it = 1;
nPar = 12;

% for par_it = 1:nPar
%     res = results(par_it, trial_it);
%     
%     spikeRecord_e = res.spikeRecord_e;
%     spikeRecord_i = res.spikeRecord_i;
%     Vmon = res.Vmon;
%     times = res.times;
%     firingRate_e = res.firingRate_e;
%     spikeCell_e = res.spikeCell_e;
%     spikeCell_i = res.spikeCell_i;
%     opt = res.opt;
%     
%     
%     figure('Position', [800 528 1200 800]);
%     
%     x_lim = [2 2.5];
% 
%     % Plot the results
% 
%     subplot(5, 1, 1);
%     spikeCellRasterPlot(spikeCell_e, '.');
%     title('Pyramidal neurons');
%     ylabel('Neuron number');
%     xlim(x_lim);
%     ylim([1 opt.Ne]);
% 
%     subplot(5, 1, 2);
%     plot(times, firingRate_e);
%     ylabel('Firing rate (Hz)');
%     xlim(x_lim);
% 
% 
%     subplot(5, 1, 3);
%     plot(Vmon.t, Vmon.e*1000);
%     %title('Pyramidal neuron');
%     ylabel('Vm (mV)');
%     xlim(x_lim);
% 
% 
%     subplot(5, 1, 4);
%     spikeCellRasterPlot(spikeCell_i, '.');
%     title('Interneurons');
%     ylabel('Neuron number');
%     xlim(x_lim);
%     ylim([1 opt.Ni]);
% 
% 
%     subplot(5, 1, 5);
%     plot(Vmon.t, Vmon.i*1000);
%     %title('Interneuron');
%     ylabel('Vm (mV)');
%     xlabel('Time (s)');
%     xlim(x_lim);
% 
%     hold off;
%     
%     set(gcf,'PaperPositionMode','auto');
%     print('-depsc2', ['output/2011-07-19/e_input_current_' num2str(par_it) '.eps');
% 
% 
% 
%     %
%     % Plot firing rate fft
%     %
% %     figure();
% %     [Y f NFFT] = fourierTrans(firingRate_e, dt);
% %     Y_abs = 2*abs(Y(1:NFFT/2+1));
% %     plot(f,Y_abs);
% %     xlim([0 200]);
% end


for par_it = 1:nPar
    for trial_it = 1:nTrials
        res = results(par_it, trial_it);
    
        spikeRecord_e = res.spikeRecord_e;
        spikeRecord_i = res.spikeRecord_i;
        Vmon = res.Vmon;
        times = res.times;
        firingRate_e = res.firingRate_e;
        spikeCell_e = res.spikeCell_e;
        spikeCell_i = res.spikeCell_i;
        opt = res.opt;

        
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

    end
end


% Print the population and excitatory cells frequency depending on input
% parameter
for par_it = 1:nPar
    Ie(par_it) = results(par_it, 1).opt.Ie
end

figure();
errorbar(Ie*1000, mean(fmax), std(fmax));
xlabel('Input drive (mV)');
ylabel('Frequency (Hz)');
legend('Population');

