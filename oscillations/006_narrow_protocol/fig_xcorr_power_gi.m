% Plot cross correlation + gamma power/spectrograms

close all;
clearvars -except results;

path('../include', path);


%load e_input_current_output_19-Jul-2011;
outputDir = 'output_local';
outputNum = '010';

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

trial_it = 1;

figure('Position', [800 800 1000 1000], 'Visible', 'off', 'Renderer', 'painters');
for par_it = 1:size(results,1)
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
        p_o.title = sprintf('sparseness: %f, we: %f pA',sp, we*pA);
        %par_it = getItFromSpWeVecs(sp, we, results, 1e-9)


        Ne_it = [1 2];
        Ni_it = [1 2];


        % Current onto stellate cells
        gi_sig1 = res.Vmon.gi(Ne_it(1), :)*pA;
        gi_sig2 = res.Vmon.gi(Ne_it(2), :)*pA;
        %times = res.times;
        
        %subs = 10;
        %gi_sig1 = gi_sig1(1:subs:end);
        %gi_sig2 = gi_sig2(1:subs:end);
        %times = times(1:subs:end);
        plot2g_xcorr(gi_sig1, gi_sig2, res.times, p_o);


        % % Current onto interneurons
        % figure('Position', [800 800 1000 1000]);
        % ge_sig1 = res.Vmon.ge(Ni_it(1), :)*pA;
        % ge_sig2 = res.Vmon.ge(Ni_it(2), :)*pA;
        % plot2g_xcorr(ge_sig1, ge_sig2, res.times, p_o);


        % firing rate histograms
        set(gcf,'PaperPositionMode','auto');
        print('-depsc2', sprintf('%s/%s_xcorr_power_par_it%.3d.eps', ...
            outputDir, outputNum, par_it));    
end
    