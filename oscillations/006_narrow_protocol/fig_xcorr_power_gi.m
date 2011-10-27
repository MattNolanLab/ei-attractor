% Plot cross correlation + gamma power/spectrograms

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '011';

nParam  = size(results, 1);
nTrials = size(results, 2);

pA = 1e12;
fontSize = 14;

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);

Nwe = size(results(1,1).opt.we_vec, 2);
Nsp = size(results(1,1).opt.sparseness_vec, 2);

sp_vec = results(1,1).opt.sparseness_vec;
we_vec = results(1,1).opt.we_vec;

sp = 0.01;
we = 90e-12;


figure('Position', [800 800 1000 1000], 'Visible', 'off', 'Renderer', 'painters');
for par_it = getItFromSpWeVecs(sp, we, results, 1e-9)%1:nParam
    for trial_it = 1%:nTrials
        p_o = struct();
        p_o.F = 10:1:150;
        p_o.sampling_rate = 1e4;
        p_o.spec_win_len = 2000;
        p_o.noverlap = p_o.spec_win_len /2;
        p_o.plot_dB = true;
        p_o.fontSize = fontSize;

        dt = results(1,1).opt.dt;
        nPeriods = 1;
        t_start = 1;
        t_start_i = fix(t_start/dt) + 1;

        p_o.x_lim = [0 15];

        xcorr_win_len_t = 0.5;
        p_o.xcorr_win_len = fix(xcorr_win_len_t / dt);

        % Simulation trial selection
        %sp = 0.01;
        %we = 4e-11;

        res = results(par_it, trial_it);
        we = res.opt.we;
        sp = res.opt.e_sparseness;
        p_o.title = sprintf('sparseness: %f, we: %f pA, trial: %.3d',sp, we*pA, trial_it);
        %par_it = getItFromSpWeVecs(sp, we, results, 1e-9)


        Ne_it = [1 2];
        Ni_it = [1 2];


%         % Current onto stellate cells
%         gi_sig1 = res.Vmon.gi(Ne_it(1), :)*pA;
%         gi_sig2 = res.Vmon.gi(Ne_it(2), :)*pA;
%         
%         plot2g_xcorr(gi_sig1, gi_sig2, res.times, p_o);
% 
% 
%         % % Current onto interneurons
%         % figure('Position', [800 800 1000 1000]);
%         % ge_sig1 = res.Vmon.ge(Ni_it(1), :)*pA;
%         % ge_sig2 = res.Vmon.ge(Ni_it(2), :)*pA;
%         % plot2g_xcorr(ge_sig1, ge_sig2, res.times, p_o);
% 
% 
%         % firing rate histograms
%         set(gcf,'PaperPositionMode','auto');
%         print('-depsc2', sprintf('%s/%s_xcorr_power_par_it%.3d_trial_%.3d.eps', ...
%             outputDir, outputNum, par_it, trial_it));    
%         
%         
%         % Compare two different samples when power is low and high
%         % Time 1 - low gamma
%         figure('Position', [800 800 1000 600]);
%         x_lim = [10 10.2];
%         subplot(2,1,1, 'FontSize', fontSize);
%         plot(res.times, gi_sig1);
%         xlabel('Time (s)');
%         ylabel('Syn. current (pA)');
%         title('Stellate cell 1');
%         xlim(x_lim);
% 
%         subplot(2, 1, 2, 'FontSize', fontSize);
%         plot(res.times, gi_sig2);
%         xlabel('Time (s)');
%         ylabel('Syn. current (pA)');
%         title('Stellate cell 2');
%         xlim(x_lim);
% 
%         set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
%         print('-depsc2', sprintf('%s/%s_gi_gi_currents_par_it%.3d_trial_%.3d_low_gamma.eps', ...
%                 outputDir, outputNum, par_it, trial_it));    
%             
%         % Time 2 - high gamma
%         figure('Position', [800 800 1000 600]);
%         x_lim = [14 14.2];
%         subplot(2,1,1, 'FontSize', fontSize);
%         plot(res.times, gi_sig1);
%         xlabel('Time (s)');
%         ylabel('Syn. current (pA)');
%         title('Stellate cell 1');
%         xlim(x_lim);
% 
%         subplot(2, 1, 2, 'FontSize', fontSize);
%         plot(res.times, gi_sig2);
%         xlabel('Time (s)');
%         ylabel('Syn. current (pA)');
%         title('Stellate cell 2');
%         xlim(x_lim);
% 
%         set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
%         print('-depsc2', sprintf('%s/%s_gi_gi_currents_par_it%.3d_trial_%.3d_high_gamma.eps', ...
%                 outputDir, outputNum, par_it, trial_it));    

            
        % Plot a histogram of average firing of excitatory cells
        figure('Position', [800 800 1000 400]);
        subplot(1,1,1, 'FontSize', fontSize');
        avg_firing_nbins = 40;
        win_len = 1e4;
        noverlap = win_len/4;
        [h_N, h_E h_POS] = slidingFiringRate(res.spikeRecord_e, win_len, noverlap, avg_firing_nbins);
        pcolor((h_POS-1)*dt, h_E, h_N');

        title('Sliding histogram of firing rates');
        ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        c_h = colorbar;
        set(get(c_h, 'Title'), 'String', 'Neuron count', 'FontSize', fontSize);
        shading flat;
        xlim(p_o.x_lim);
        
        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-depsc2', sprintf('%s/%s_e_firing_rate_par_it%.3d_trial_%.3d.eps', ...
                outputDir, outputNum, par_it, trial_it));    

    end
end
    