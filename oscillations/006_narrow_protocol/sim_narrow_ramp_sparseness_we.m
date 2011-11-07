% Determine lag between average first spike of interneurons and stellate
% cell
% This simulation produces data as a function of sparseness and synaptic
% strength

clear all;
close all;

path('../include/', path);

%results = [];

global_opt.Ne = 800;
global_opt.Ni = 200;
N = global_opt.Ne + global_opt.Ni;

global_opt.sparseness_vec = [0.1];
%global_opt.we_vec = [600 800 1000 1200 1400 1600] * 1e-12;
global_opt.taue_vec = [1:0.2:2] * 1e-3;


%Nspar = numel(global_opt.sparseness_vec);
%Nwe = numel(global_opt.we_vec);
Ntaue = numel(global_opt.taue_vec);

parfor it = 1:Ntaue
    %sparseness = global_opt.sparseness_vec(fix((it-1)/Nwe) + 1)
    %we = global_opt.we_vec(mod(it-1, Nwe) + 1)
    taue = global_opt.taue_vec(it);
    
    %for we = opt.we_vec
        opt = global_opt;

        % All variables are in basic units, i.e. s, volt, etc.

        % Excitatory cells
        opt.taum_e = 9.3e-3;
        opt.taue = taue;
        opt.El_e = -68.5e-3;
        opt.Vt_e = -50.0e-3;
        opt.Vr_e = opt.El_e;
        opt.Rm_e = 60e6; % MOhm
        opt.Ie_0 = 800e-12;  % pA
        opt.Ie_max = 800e-12; % pA
        %opt.we_vec = global_opt.we_vec;
        opt.we = 1000e-12;
        opt.refrac_e_mean = 40e-3;
        opt.refrac_e_std = 5e-3;
        opt.refrac_e = opt.refrac_e_mean + opt.refrac_e_std*randn(global_opt.Ne, 1);
        opt.refrac_e_g_inc = 1/opt.Rm_e/2;


        % Inhibitory cell
        opt.taum_i = 10e-3;
        opt.taui = 7e-3;
        opt.El_i = -60e-3;
        opt.Vt_i = -50e-3;
        opt.Vr_i = opt.El_i;
        opt.Rm_i = 80e6; % MOhm
        opt.Ii_0 = 125e-12; % pA
        opt.Ii_max = 125e-12; % pA
        opt.wi = 350e-12; % pS
        opt.refrac_i_mean = 7.5e-3; %msec
        opt.refrac_i_std  = 0.5e-3;
        opt.refrac_i = opt.refrac_i_mean + opt.refrac_i_std*randn(global_opt.Ni, 1);
        opt.refrac_i_g_inc = 1/opt.Rm_i;

        opt.spikeVm = 40e-3;
        
        % Reversal potentials
        opt.V_rev_e = 0e-3;
        opt.V_rev_i = -75e-3;

        opt.sparseness_vec = global_opt.sparseness_vec;
        opt.e_sparseness = 0.1;
        opt.i_sparseness = 0.8;
        
        opt.Vclamp = -50e-3;


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
        opt.Emon_i = [100 120 140 160];
        opt.Imon_i = [100 120 140 160];

        % simulation time
        opt.T = 5;


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

        nTrials = 1;
        param_i = 1;

        net_data = struct();
    
        %
        % Start loop
        %        
        [net_data.Me net_data.Mi] = MeMi(opt);

        tmpresults = [];
        for trialNum = 1:nTrials
            trialNum
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
        end

        results(it, :) = tmpresults;
    %end
end

clearvars -except results;
save('-v7.3', sprintf('013_narrow_const_sparseness_we_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS')));
