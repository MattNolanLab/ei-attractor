% Plot cross correlation + gamma power/spectrograms

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '012';

nParam  = size(results, 1);
nTrials = size(results, 2);

pA = 1e12;
pS = 1e12;
fontSize = 17;

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);

%Nwe = size(results(1,1).opt.we_vec, 2);
%Nsp = size(results(1,1).opt.sparseness_vec, 2);

%sp_vec = results(1,1).opt.sparseness_vec;
%we_vec = results(1,1).opt.we_vec;

dt = results(1,1).opt.dt;

sp = 0.1;
we = 600e-12;


for par_it = 1%getItFromSpWeVecs(sp, we, results, 1e-9)%1:nParam
    for trial_it = 1%:nTrials
        p_o = struct();
        p_o.F = 10:1:150;
        p_o.sampling_rate = 1/dt;
        p_o.spec_win_len = fix(200e-3/dt);
        p_o.noverlap = 0.75*p_o.spec_win_len;
        p_o.plot_dB = false;
        p_o.fontSize = fontSize;

        dt = results(1,1).opt.dt;
        nPeriods = 1;
        t_start = 1;
        t_start_i = fix(t_start/dt) + 1;

        p_o.x_lim = [0 10];

        xcorr_win_len_t = 0.5;
        p_o.xcorr_win_len = fix(xcorr_win_len_t / dt);

        % Simulation trial selection
        %sp = 0.01;
        %we = 4e-11;

        res = results(par_it, trial_it);
        we = res.opt.we;
        sp = res.opt.e_sparseness;
        p_o.title = sprintf('Syn. input to stellate cell; sparseness: %f, we: %f pS, trial: %.3d',sp, we*pS, trial_it);
        %par_it = getItFromSpWeVecs(sp, we, results, 1e-9)


        Ne_it = [1 2];
        Ni_it = [1 2];

        Isyn_e1 = -res.Vmon.Isyn_e(Ne_it(1), :)*pA;
        Isyn_e2 = -res.Vmon.Isyn_e(Ne_it(2), :)*pA;
        Isyn_i1 = -res.Vmon.Isyn_i(Ni_it(1), :)*pA;
        Isyn_i2 = -res.Vmon.Isyn_i(Ni_it(2), :)*pA;


        % Plot all the currents in time
        figure('Position', [600 800 1000 500]);
        subplot(10,1,1:9, 'FontSize', fontSize);
        plot(res.times, Isyn_e1, 'r', res.times, Isyn_i1, 'b');
        xlabel('Time (s)');
        ylabel('Current (pA)');
        box off;
        xlim(p_o.x_lim);
        
        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_currents_whole_ste_int.eps', ...
            outputDir, outputNum, par_it, trial_it));    
        

        % Current onto stellate cells
        figure('Position', [600 800 1000 800], 'Renderer', 'painters');
        
        plot2g_xcorr(Isyn_e1, Isyn_e2, res.times, p_o);

        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_xcorr_power_ste.eps', ...
            outputDir, outputNum, par_it, trial_it));    
        print('-dpng', sprintf('%s/%s_par_it%.3d_trial_%.3d_xcorr_power_ste.png', ...
            outputDir, outputNum, par_it, trial_it));    

        
%         % Current onto interneurons
%         figure('Position', [800 800 1000 1000], 'Renderer', 'painters');
%         p_o.title = sprintf('Syn. input to interneurons');
%         plot2g_xcorr(Isyn_i1, Isyn_i2, res.times, p_o);
% 
%         set(gcf,'PaperPositionMode','auto');
%         print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_xcorr_power_int.eps', ...
%             outputDir, outputNum, par_it, trial_it));    
% 

        
        
        % Compare two different samples when power is low and high
        % Time 1 - low gamma
        figure('Position', [800 800 800 400]);
        x_lim = [6.5 6.7];
        subplot(10,1,1:9, 'FontSize', fontSize);
        plot(res.times, Isyn_e1, 'r', res.times, Isyn_i1, 'b');
        xlabel('Time (s)');
        ylabel('Syn. current (pA)');
        %axis tight;
        xlim(x_lim);
        box off;

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_gi_ge_currents_low_gamma.eps', ...
                outputDir, outputNum, par_it, trial_it));  
            
            
        % Time 2 - high gamma
        figure('Position', [800 800 800 400]);
        x_lim = [9 9.2];
        subplot(10,1,1:9, 'FontSize', fontSize);
        plot(res.times, Isyn_e1, 'r', res.times, Isyn_i1, 'b');
        xlabel('Time (s)');
        ylabel('Syn. current (pA)');
        %axis tight;
        xlim(x_lim);
        box off;

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_gi_ge_currents_high_gamma.eps', ...
                outputDir, outputNum, par_it, trial_it));    

            
        % Print currents onto stellate cells
        figure('Position', [800 800 800 400]);
        x_lim = [9 9.2];
        subplot(10,1,1:9, 'FontSize', fontSize);
        plot(res.times, Isyn_e1, 'r', res.times, Isyn_e2, 'k');
        xlabel('Time (s)');
        ylabel('Syn. current (pA)');
        %axis tight;
        xlim(x_lim);
        box off;

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_gi_gi_currents_high_gamma.eps', ...
                outputDir, outputNum, par_it, trial_it));    

            
            
        % Plot a histogram of average firing of excitatory cells
        maxF_slidingRate = 80;

        figure('Position', [800 800 1000 400]);
        subplot(10,1,2:9, 'FontSize', fontSize');
        avg_firing_nbins = 40;
        win_len = 1/dt;
        noverlap = win_len/4;
        [h_N, h_E h_POS] = slidingFiringRate(res.spikeRecord_e, win_len, ...
            noverlap, avg_firing_nbins, maxF_slidingRate);
        pcolor((h_POS-1)*dt, h_E, h_N');

        title('Sliding histogram of firing rates: Excitatory neurons');
        ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        c_h = colorbar;
        set(get(c_h, 'Title'), 'String', 'Neuron count', 'FontSize', fontSize);
        shading flat;
        xlim(p_o.x_lim);
        
        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-dpng', sprintf('%s/%s_par_it%.3d_trial_%.3d_e_firing_rate_st.png', ...
                outputDir, outputNum, par_it, trial_it));    

            
        % Plot a histogram of average firing of inhibitory cells
        figure('Position', [800 800 1000 400]);
        subplot(10,1,2:9, 'FontSize', fontSize');
        avg_firing_nbins = 40;
        win_len = 1/dt;
        noverlap = win_len/4;
        [h_N, h_E h_POS] = slidingFiringRate(res.spikeRecord_i, win_len, ...
            noverlap, avg_firing_nbins, maxF_slidingRate);
        pcolor((h_POS-1)*dt, h_E, h_N');

        title('Sliding histogram of firing rates: Interneurons');
        ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        c_h = colorbar;
        set(get(c_h, 'Title'), 'String', 'Neuron count', 'FontSize', fontSize);
        shading flat;
        xlim(p_o.x_lim);
        
        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-dpng', sprintf('%s/%s_par_it%.3d_trial_%.3d_e_firing_rate_int.png', ...
                outputDir, outputNum, par_it, trial_it));    
            

            
        % Plot excitation profile
        figure('Position', [800 800 1000 400]);
        subplot(10,1,1:9, 'FontSize', fontSize);

        t = 0:0.01:p_o.x_lim(2);
        plot(t, [res.opt.Ie_max/res.opt.T.*t; res.opt.Ii_max/res.opt.T.*t]*pA, 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Current (pA)');
        legend('Excitatory neuron', 'Interneuron', 'Location', 'NorthWest');
        title('External current injected to cells');

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_ramp_current_fig.eps', ...
                outputDir, outputNum, par_it, trial_it));    

    end
end
    