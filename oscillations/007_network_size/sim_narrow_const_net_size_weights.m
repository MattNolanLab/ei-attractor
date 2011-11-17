% Simulate how synchronisation dynamics of a network changes as a function
% of network size
%
% This variant scales conductances as a function of network size
%clear all;
close all;

path('../include/', path);

outputDir = 'output_local';
outputNum = '002';

fontSize = 16;

sim_flag = false;


nTrials = 2;
Ni_it = 1;
Ne_it = 1;
pA = 1e12;



if sim_flag
    N_vec = [100:50:1000];
    NN = numel(N_vec);

    %global_opt.we_std_vec = [600:10:1300] * 1e-12;
    %Nwe_std = numel(global_opt.we_std_vec);

    alpha_E_max = 5.40e-8;
    alpha_I_max = 3.96e-8;

    al_coeff_vec = [0.5:0.05:1.5];
    Nalpha = numel(al_coeff_vec);

    % N...outer, alpha...inner cycle
    parfor it = 1:NN*Nalpha
        display(sprintf('parameter %d out of %d', it, NN*Nalpha));
            opt = getStdCellProperties();



            opt.N = N_vec(fix((it-1)/Nalpha) + 1);
            opt.alpha_coeff = al_coeff_vec(mod(it-1, Nalpha) + 1);

            opt.Ne = fix(0.9 * opt.N);
            opt.Ni = fix(0.1 * opt.N);

            % All variables are in basic units, i.e. s, volt, etc.
            opt.e_sparseness = 0.4;
            opt.i_sparseness = 0.8;


            opt.Ie_0 = 900e-12;  % pA
            opt.Ie_max = 900e-12; % pA
            opt.we = 1/opt.Ne * alpha_E_max / opt.e_sparseness * opt.alpha_coeff;
            opt.we_std = 600e-12;
            opt.refrac_e = opt.refrac_e_mean + opt.refrac_e_std*randn(opt.Ne, 1);


            opt.Ii_0 = 200e-12; % pA
            opt.Ii_max = 200e-12; % pA
            opt.wi = 1/opt.Ni * alpha_I_max / opt.i_sparseness * opt.alpha_coeff; % pS
            opt.refrac_i = opt.refrac_i_mean + opt.refrac_i_std*randn(opt.Ni, 1);

            % Current distribution settings
            % Diameter of the activated area
            opt.D = 100e-6; % micrometers
            opt.input_spread = 2*opt.D;


            % Noise (mV)
            opt.noise_sigma = 0.02e-3;
            opt.sigma_init_cond = 10e-3; % 2mV


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


            % 
            % Create simulation results
            %
            opt.dists_e = opt.D*rand(opt.Ne, 1);
            opt.dists_i = opt.D*rand(opt.Ni, 1);    


            % Now assume neurons are uniformly distributed in the specified area,
            % generate their distances and input current according to the specified
            % spread function
            opt.Ie = gaussianSpread(opt.dists_e, opt.input_spread, opt.Ie_0);
            opt.Ii = gaussianSpread(opt.dists_i, opt.input_spread, opt.Ii_0);

            Ie_max = gaussianSpread(opt.dists_e, opt.input_spread, opt.Ie_max);
            Ii_max = gaussianSpread(opt.dists_i, opt.input_spread, opt.Ii_max);

            opt.dIe = (Ie_max - opt.Ie) / opt.T;
            opt.dIi = (Ii_max - opt.Ii) / opt.T;

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
                [spikeRecord_e, spikeRecord_i, Vmon, times] = simulateEIRamp(opt, net_data);

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
    
    save('-v7.3', sprintf('001_narrow_const_net_size_weights_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS')));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coherence_th = 0.2;

trial_it = 1;
F = reshape(results_F(:, trial_it), Nalpha, NN)';
coh = reshape(results_coh(:, trial_it), Nalpha, NN)';

F(coh < coherence_th) = nan;
coh(coh < coherence_th) = nan;

figure('Position', [800 1050 900 1050]);
[N_grid al_grid] = meshgrid(al_coeff_vec, N_vec);
subplot(2,1,1, 'FontSize', fontSize);
surf(al_grid, N_grid, F);
xlabel('Network size');
ylabel('Synaptic coupling (normalised)');
zlabel('Frequency (Hz)');
view(-64, 46);

subplot(2,1,2, 'FontSize', fontSize);
surf(al_grid, N_grid, coh);
xlabel('Network size');
ylabel('Synaptic coupling (normalised)');
zlabel('Coherence');
view(-64, 46);


set(gcf,'PaperPositionMode','auto', 'Renderer', 'painters');
print('-depsc2', sprintf('%s/%s_freq_coh_net_size_syn_strength.eps', ...
        outputDir, outputNum));    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print samples of synaptic conductances for each trial
p_o.figVisible = false;
p_o.x_lim = [0 2.5];
p_o.figVisible = 'off';

parfor it = 1:NN*Nalpha
    N = N_vec(fix((it-1)/Nalpha) + 1);
    alpha_coeff = al_coeff_vec(mod(it-1, Nalpha) + 1);
    
    res = results(it, trial_it);
    
    Isyn_e1 = -res.Vmon.Isyn_e(Ne_it(1), :)*pA;
    Isyn_i1 = -res.Vmon.Isyn_i(Ni_it(1), :)*pA;


    % Plot all the currents in time
    figure('Position', [600 800 1000 500], 'Visible', p_o.figVisible);
    subplot(10,1,1:9, 'FontSize', fontSize);
    plot(res.times, Isyn_e1, 'r', res.times, Isyn_i1, 'b');
    xlabel('Time (s)');
    ylabel('Current (pA)');
    box off;
    xlim(p_o.x_lim);
    title(sprintf('%d neurons; total syn. coupling: %.3f', N, alpha_coeff));

    set(gcf,'PaperPositionMode','auto');
    print('-depsc2', sprintf('%s/%s_N%.4d_al%.3f_trial_%.3d_currents_whole_ste_int.eps', ...
        outputDir, outputNum, N, alpha_coeff, trial_it));    
end
    
