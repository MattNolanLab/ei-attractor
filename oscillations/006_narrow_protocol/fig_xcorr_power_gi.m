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
fontSize = 14;

Ne = size(results(1,1).spikeCell_e, 2);
Ni = size(results(1,1).spikeCell_i, 2);

Nwe = size(results(1,1).opt.we_vec, 2);
Nsp = size(results(1,1).opt.sparseness_vec, 2);

sp_vec = results(1,1).opt.sparseness_vec;
we_vec = results(1,1).opt.we_vec;

sp = 0.1;
we = 4080e-12;


for par_it = getItFromSpWeVecs(sp, we, results, 1e-9)%1:nParam
    for trial_it = 1%:nTrials
        p_o = struct();
        p_o.F = 10:1:150;
        p_o.sampling_rate = 1e4;
        p_o.spec_win_len = 2000;
        p_o.noverlap = p_o.spec_win_len /2;
        p_o.plot_dB = false;
        p_o.fontSize = fontSize;

        dt = results(1,1).opt.dt;
        nPeriods = 1;
        t_start = 1;
        t_start_i = fix(t_start/dt) + 1;

        p_o.x_lim = [0 20];

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


        % Current onto stellate cells
        figure('Position', [800 800 1000 1000], 'Renderer', 'painters');
        Isyn_e1 = res.Vmon.Isyn_e(Ne_it(1), :)*pA;
        Isyn_e2 = res.Vmon.Isyn_e(Ne_it(2), :)*pA;
        
        plot2g_xcorr(Isyn_e1, Isyn_e2, res.times, p_o);



        % firing rate histograms
        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_xcorr_power.eps', ...
            outputDir, outputNum, par_it, trial_it));    
        
        % Current onto interneurons
        figure('Position', [800 800 1000 1000], 'Renderer', 'painters');
        Isyn_i1 = res.Vmon.Isyn_i(Ni_it(1), :)*pA;
        Isyn_i2 = res.Vmon.Isyn_i(Ni_it(2), :)*pA;
        p_o.title = sprintf('Syn. input to interneurons');
        plot2g_xcorr(Isyn_i1, Isyn_i2, res.times, p_o);


        
        
        % Compare two different samples when power is low and high
        % Time 1 - low gamma
        figure('Position', [800 800 1000 600]);
        x_lim = [10 10.2];
        subplot(2,1,1, 'FontSize', fontSize);
        plot(res.times, Isyn_e1);
        xlabel('Time (s)');
        ylabel('Syn. current (pA)');
        title('Stellate cell 1');
        xlim(x_lim);

        subplot(2, 1, 2, 'FontSize', fontSize);
        plot(res.times, Isyn_e2);
        xlabel('Time (s)');
        ylabel('Syn. current (pA)');
        title('Stellate cell 2');
        xlim(x_lim);

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_gi_gi_currents_low_gamma.eps', ...
                outputDir, outputNum, par_it, trial_it));    
            
        % Time 2 - high gamma
        figure('Position', [800 800 1000 600]);
        x_lim = [14 14.2];
        subplot(2,1,1, 'FontSize', fontSize);
        plot(res.times, Isyn_e1);
        xlabel('Time (s)');
        ylabel('Syn. current (pA)');
        title('Stellate cell 1');
        xlim(x_lim);

        subplot(2, 1, 2, 'FontSize', fontSize);
        plot(res.times, Isyn_e2);
        xlabel('Time (s)');
        ylabel('Syn. current (pA)');
        title('Stellate cell 2');
        xlim(x_lim);

        set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
        print('-depsc2', sprintf('%s/%s_par_it%.3d_trial_%.3d_gi_gi_currents_high_gamma.eps', ...
                outputDir, outputNum, par_it, trial_it));    

            
        % Plot a histogram of average firing of excitatory cells
        maxF_slidingRate = 150;

        figure('Position', [800 800 1000 400]);
        subplot(10,1,2:9, 'FontSize', fontSize');
        avg_firing_nbins = 40;
        win_len = 1e4;
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
        win_len = 1e4;
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
        figure('Position', [800 800 1000 600]);
        subplot(1,1,1, 'FontSize', fontSize);

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
    