% Simulate how synchronisation dynamics of a network changes as a function
% of network size
%
% This variant scales conductances as a function of network size
%clear all;
close all;

path('../include/', path);

outputDir = 'output_local';
outputNum = 'tmp';

fontSize = 16;

sim_flag = true;


nTrials = 1;
Ni_it = 1;
Ne_it = 1;
pA = 1e12;



if sim_flag
    N_vec = 1000;%[100:50:1000];
    NN = numel(N_vec);


    alpha_E_max = 5.40e-8;
    alpha_I_max = 4.5e-8;

    al_coeff_vec = 1; %[0.5:0.05:1.5];
    Nalpha = numel(al_coeff_vec);

    % N...outer, alpha...inner cycle
    parfor it = 1:NN*Nalpha
        display(sprintf('parameter %d out of %d', it, NN*Nalpha));
            opt = getStdCellProperties();


            opt.N = N_vec(fix((it-1)/Nalpha) + 1);
            opt.alpha_coeff = al_coeff_vec(mod(it-1, Nalpha) + 1);

            opt.Ne = fix(0.9 * opt.N);
            opt.Ni = fix(0.1 * opt.N);
            opt.N_ext = 1000;

            % All variables are in basic units, i.e. s, volt, etc.
            opt.e_sparseness = 0.4;
            opt.i_sparseness = 0.8;
            
            % Densities of external connections onto E/I populations
            opt.e_ext_density = 0.25;
            opt.i_ext_density = 0.25;
            
            % Firing rates of external connections (per ext. neuron)
            opt.r_ext = 40;


            opt.we = 1/opt.Ne * alpha_E_max / opt.e_sparseness * opt.alpha_coeff;
            opt.we_std = 600e-12;
            opt.refrac_e = opt.refrac_e_mean + opt.refrac_e_std*randn(opt.Ne, 1);


            opt.wi = 1/opt.Ni * alpha_I_max / opt.i_sparseness * opt.alpha_coeff; % pS
            opt.refrac_i = opt.refrac_i_mean + opt.refrac_i_std*randn(opt.Ni, 1);

            opt.we_ext_e = 1400e-12;
            opt.we_ext_i = opt.we_ext_e/4.5;

            % Noise (mV)
            opt.noise_sigma = 0.02e-3;
            opt.sigma_init_cond = 10e-3; % 10mV


            % Euler settings
            opt.dt = 0.05e-3  % 0.05 ms
            dt = opt.dt;


            % Firing rate sliding window length
            opt.rateWindowLen = 0.005; %ms
            rateWindowLen = opt.rateWindowLen;

            % Vm monitor, neuron index
            opt.Emon_i = [10 20];
            opt.Imon_i = [1  2];

            % simulation time
            opt.T = 2.5;


            param_i = 1;

            net_data = struct();

            %
            % Start loop
            %        
            [net_data.Me net_data.Mi] = MeMi(opt);

            F = zeros(1, nTrials) * nan;
            coh = zeros(1, nTrials) * nan;
            tmpresults = [];
            for trialNum = 1:nTrials
                trialNum;
                [spikeRecord_e, spikeRecord_i, Vmon, times] = simulateEIHomo(opt, net_data);

                tmpresults(trialNum).spikeRecord_e = spikeRecord_e;
                tmpresults(trialNum).spikeRecord_i = spikeRecord_i;
                tmpresults(trialNum).Me = net_data.Me;
                tmpresults(trialNum).Mi = net_data.Mi;
                tmpresults(trialNum).Vmon = Vmon;
                tmpresults(trialNum).times = times;


                tmpresults(trialNum).firingRate_e = getFiringRate(spikeRecord_e, dt, rateWindowLen);
                tmpresults(trialNum).firingRate_i = getFiringRate(spikeRecord_i, dt, rateWindowLen);
                tmpresults(trialNum).spikeCell_e = spikeRecordToSpikeCell(spikeRecord_e, times);
                tmpresults(trialNum).spikeCell_i = spikeRecordToSpikeCell(spikeRecord_i, times);

                tmpresults(trialNum).opt = opt;
                
                % Data processing part
                nPeriods = 1;
                t_start = 1;
                t_start_i = fix(t_start/dt) + 1;

                Isyn_e = Vmon.Isyn_e(Ni_it, t_start_i:end)*pA;
                [F(trialNum) coh(trialNum)] = getSigFreqCoherence(Isyn_e, opt.dt, nPeriods);
            end

            results_F(it, :) = F;
            results_coh(it, :) = coh;
            results(it, :) = tmpresults;

    end
    
    save('-v7.3', sprintf('%s_homogeneous_input_%s.mat', outputNum, datestr(now, 'yyyy-mm-dd_HH-MM-SS')));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% coherence_th = 0.2;
% 
% trial_it = 1;
% F = reshape(results_F(:, trial_it), Nalpha, NN)';
% coh = reshape(results_coh(:, trial_it), Nalpha, NN)';
% 
% F(coh < coherence_th) = nan;
% coh(coh < coherence_th) = nan;
% 
% figure('Position', [800 1050 900 1050]);
% [N_grid al_grid] = meshgrid(al_coeff_vec, N_vec);
% subplot(2,1,1, 'FontSize', fontSize);
% surf(al_grid, N_grid, F);
% xlabel('Network size');
% ylabel('Synaptic coupling (normalised)');
% zlabel('Frequency (Hz)');
% view(-64, 46);
% 
% subplot(2,1,2, 'FontSize', fontSize);
% surf(al_grid, N_grid, coh);
% xlabel('Network size');
% ylabel('Synaptic coupling (normalised)');
% zlabel('Coherence');
% view(-64, 46);
% 
% 
% set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
% print('-depsc2', sprintf('%s/%s_freq_coh_net_size_syn_strength.eps', ...
%         outputDir, outputNum));    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print samples of synaptic conductances for each trial
p_o.x_lim = [1.5 2.5];
p_o.figVisible = 'off';
p_o.f_x_lim = [0 400];
p_o.C_x_lim = [0 200];

for it = 1:NN*Nalpha
    N = N_vec(fix((it-1)/Nalpha) + 1);
    alpha_coeff = al_coeff_vec(mod(it-1, Nalpha) + 1);
    
    res = results(it, trial_it);
    dt = res.opt.dt;
    
    t_start_i = fix(p_o.x_lim(1)/dt) + 1;
    t_end_i = fix(p_o.x_lim(2)/dt) + 1;

    
    Isyn_e1 = -res.Vmon.Isyn_e(Ne_it(1), t_start_i:t_end_i)*pA;
    Isyn_i1 = -res.Vmon.Isyn_i(Ni_it(1), t_start_i:t_end_i)*pA;
    times = res.times(t_start_i:t_end_i);


    % Plot all the currents in time
    figure('Position', [600 800 1000 800], 'Visible', p_o.figVisible);
    subplot(3,2,1:2, 'FontSize', fontSize);
    plot(times, Isyn_e1, 'r', times, Isyn_i1, 'b');
    xlabel('Time (s)');
    ylabel('Current (pA)');
    box off;
    xlim(p_o.x_lim);
    title(sprintf('%d neurons; total syn. coupling: %.3f; A', N, alpha_coeff));
    
    % Power spectrum of inhibitory currents
    subplot(3, 2, 3, 'FontSize', fontSize);
    [e_Y e_f e_NFFT] = fourierTrans(Isyn_e1 - mean(Isyn_e1), dt);
    e_Y_pow = (2*abs(e_Y(1:e_NFFT/2+1))).^2;
    plot(e_f, e_Y_pow, 'r');
    xlabel('Frequency (Hz)');
    ylabel('Power (pA^2)');
    xlim(p_o.f_x_lim);
    title('B');

    % Power spectrum of excitatory currents
    subplot(3, 2, 4, 'FontSize', fontSize);
    [i_Y i_f i_NFFT] = fourierTrans(Isyn_i1 - mean(Isyn_i1), dt);
    i_Y_pow = (2*abs(i_Y(1:i_NFFT/2+1))).^2;
    plot(i_f, i_Y_pow, 'b');
    xlabel('Frequency (Hz)');
    ylabel('Power (pA^2)');
    xlim(p_o.f_x_lim);
    title('C');
    
    % Autocorrelation of inhibitory currents
    subplot(3, 2, 5, 'FontSize', fontSize);
    e_C = autoCorr(Isyn_e1);
    plot((0:dt:times(end)-times(1))*1000, e_C, 'r');
    xlabel('Time (ms)');
    ylabel('Correlation');
    xlim(p_o.C_x_lim);
    title('D');
    
    % Autocorrelation of inhibitory currents
    subplot(3, 2, 6, 'FontSize', fontSize);
    i_C = autoCorr(Isyn_i1);
    plot((0:dt:times(end)-times(1))*1000, i_C, 'b');
    xlabel('Time (ms)');
    ylabel('Correlation');
    xlim(p_o.C_x_lim);
    title('E');
    

    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('%s/%s_N%.4d_al%.3f_trial_%.3d_currents_whole_ste_int.eps', ...
        outputDir, outputNum, N, alpha_coeff, trial_it));    
end
    
