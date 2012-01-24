% Simulate how synchronisation dynamics of a network changes as a function
% of network size
%
% This variant scales conductances as a function of network size
%clear all;
close all;

path('../include/', path);

outputDir = 'output_local';
outputNum = '000';

fontSize = 16;

sim_flag = true;


nTrials = 4;
Ni_it = 1;
Ne_it = 1;
pA = 1e12;



if sim_flag
    N_vec = [100];
    NN = numel(N_vec);


    alpha_E_max = 5.40e-8;
    alpha_I_max = 3.96e-8;

    al_coeff_vec = [1];
    Nalpha = numel(al_coeff_vec);

    % N...outer, alpha...inner cycle
    parfor it = 1:NN*Nalpha
        display(sprintf('parameter %d out of %d', it, NN*Nalpha));
            opt = getStdCellProperties();


            opt.N = N_vec(fix((it-1)/Nalpha) + 1);
            opt.alpha_coeff = al_coeff_vec(mod(it-1, Nalpha) + 1);

            opt.Ne = fix(0.9 * opt.N);
            opt.Ni = fix(0.1 * opt.N);
            opt.Ne_ext = 10;
            opt.Ni_ext = 10;

            % All variables are in basic units, i.e. s, volt, etc.
            opt.e_sparseness = 0.4;
            opt.i_sparseness = 0.8;
            
            % Densities of external connections onto E/I populations
            opt.e_ext_density = 1;
            opt.i_ext_density = 1;
            
            % Firing rates of external connections (per ext. neuron)
            opt.r_e_ext = 40;
            opt.r_i_ext = 40;


            opt.we = 1/opt.Ne * alpha_E_max / opt.e_sparseness * opt.alpha_coeff;
            opt.we_std = 600e-12;
            opt.refrac_e = opt.refrac_e_mean + opt.refrac_e_std*randn(opt.Ne, 1);


            opt.wi = 1/opt.Ni * alpha_I_max / opt.i_sparseness * opt.alpha_coeff; % pS
            opt.refrac_i = opt.refrac_i_mean + opt.refrac_i_std*randn(opt.Ni, 1);


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
            end
    end
    
    %save('-v7.3', sprintf('%s_homogeneous_input_%s.mat', outputNum, datestr(now, 'yyyy-mm-dd_HH-MM-SS')));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

