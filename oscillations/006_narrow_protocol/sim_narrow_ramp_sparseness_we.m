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

global_opt.sparseness_vec = [0.01 0.05 0.1:0.1:0.8];
global_opt.we_vec = [50 100 150 200:200:4000] *1e-3 / N;


Nspar = numel(global_opt.sparseness_vec);
Nwe = numel(global_opt.we_vec);

parfor it = 1:Nspar*Nwe
    sparseness = global_opt.sparseness_vec(fix((it-1)/Nwe) + 1)
    we = global_opt.we_vec(mod(it-1, Nwe) + 1)
    
    %for we = opt.we_vec
        opt = global_opt;

        % All variables are in basic units, i.e. s, volt, etc.

        % Excitatory cells
        opt.taum_e = 9.3e-3;
        opt.taue = 2e-3;
        opt.El_e = -68.5e-3;
        opt.Vt_e = -50.0e-3;
        opt.Vr_e = -60.0e-3;
        opt.Ie_0 = 40e-3;
        opt.Ie_max = 40e-3;
        opt.we_vec = global_opt.we_vec;
        opt.we = we;


        % Inhibitory cell
        opt.taum_i = 10e-3;
        opt.taui = 10e-3;
        opt.El_i = -60e-3;
        opt.Vt_i = -50e-3;
        opt.Vr_i = -58e-3;
        opt.Ii_0 = 9e-3;
        opt.Ii_max = 9e-3;
        opt.wi = 800e-3 / N;

        opt.spikeVm = 0;

        opt.sparseness_vec = global_opt.sparseness_vec;
        opt.e_sparseness = sparseness;
        opt.i_sparseness = 0.8;


        % Current distribution settings
        % Diameter of the activated area
        opt.D = 100e-6; % micrometers
        opt.input_spread = 2*opt.D;


        % Noise (mV)
        opt.noise_sigma = 0.02e-3;
        opt.sigma_init_cond = 2e-3; % 2mV


        % Euler settings
        opt.dt = 0.1e-3  % 0.1 ms
        dt = opt.dt;


        % Firing rate sliding window length
        opt.rateWindowLen = 0.005; %ms
        rateWindowLen = opt.rateWindowLen;

        % Vm monitor, neuron index
        opt.Emon_i = [100 120];
        opt.Imon_i = [100 120 140];

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
save('-v7.3', sprintf('010_narrow_const_sparseness_we_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS')));
