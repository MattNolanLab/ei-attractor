% Narrow stimulation protocol

clear all;
close all;

path('../include/', path);

% All variables are in basic units, i.e. s, volt, etc.
opt.Ne = 800;
opt.Ni = 200;
N = opt.Ne + opt.Ni;


% Excitatory cells
opt.taum_e = 9.3e-3;
opt.taue = 2e-3;
opt.El_e = -68.5e-3;
opt.Vt_e = -50.0e-3;
opt.Vr_e = -60.0e-3;
opt.e_sparseness = 0.75;
opt.Ie_0 = 18e-3;
opt.Ie_max = 20e-3;
opt.we = 400e-3 / N;


% Inhibitory cell
opt.taum_i = 10e-3;
opt.taui = 5e-3;
opt.El_i = -60e-3;
opt.Vt_i = -50e-3;
opt.Vr_i = -58e-3;
opt.i_sparseness = 0.75;
opt.Ii_0 = 0;
opt.Ii_max = 0;
opt.wi = 20e-3 / N;

opt.spikeVm = 0;



% Current distribution settings
% Diameter of the activated area
opt.D = 100e-6; % micrometers
opt.input_spread = 10*opt.D;


% Noise (mV)
opt.noise_sigma = 0.02e-3;
opt.sigma_init_cond = 2e-3; % 2mV


% Euler settings
opt.dt = 0.5e-3  % 0.5 ms
dt = opt.dt;


% Firing rate sliding window length
opt.rateWindowLen = 0.005; %ms
rateWindowLen = opt.rateWindowLen;

% Vm monitor, neuron index
opt.Emon_i = [100 110 120 130 140 150 200 300 400 500 600 700 800];
opt.Imon_i = [100 110 120 130 140 150 160 170 180 190];

% simulation time
opt.T = 8;

opt.input_spread_vec = [2]*opt.D;

% 
% Create simulation results
%
nTrials = 8;
param_i = 1;

% The same network structure for all simulations
[net_data.Me net_data.Mi] = MeMi(opt);

opt.dists_e = opt.D*rand(opt.Ne, 1);
opt.dists_i = opt.D*rand(opt.Ni, 1);    


for spread = opt.input_spread_vec

    opt.input_spread = spread;
    
    % Now assume neurons are uniformly distributed in the specified area,
    % generate their distances and input current according to the specified
    % spread function
    opt.Ie = gaussianSpread(opt.dists_e, opt.input_spread, opt.Ie_0);
    opt.Ii = gaussianSpread(opt.dists_i, opt.input_spread, opt.Ii_0);
    
    Ie_max = gaussianSpread(opt.dists_e, opt.input_spread, opt.Ie_max);
    Ii_max = gaussianSpread(opt.dists_i, opt.input_spread, opt.Ii_max);
    
    opt.dIe = (Ie_max - opt.Ie) / opt.T;
    opt.dIi = (Ii_max - opt.Ii) / opt.T;

    
    parfor trialNum = 1:nTrials
        trialNum
        [spikeRecord_e, spikeRecord_i, Vmon, times] = simulateEIRamp(opt, net_data);
        
        tmpresults(trialNum).spikeRecord_e = spikeRecord_e;
        tmpresults(trialNum).spikeRecord_i = spikeRecord_i;
        tmpresults(trialNum).Vmon = Vmon;
        tmpresults(trialNum).times = times;

        
        tmpresults(trialNum).firingRate_e = getFiringRate(spikeRecord_e, dt, rateWindowLen);
        tmpresults(trialNum).firingRate_i = getFiringRate(spikeRecord_i, dt, rateWindowLen);
        tmpresults(trialNum).spikeCell_e = spikeRecordToSpikeCell(spikeRecord_e, times);
        tmpresults(trialNum).spikeCell_i = spikeRecordToSpikeCell(spikeRecord_i, times);

        tmpresults(trialNum).opt = opt;
    end
    
    results(param_i, :) = tmpresults;
    
    param_i = param_i + 1;
end

clear tmpresults;
save('-v7.3', sprintf('001_narrow_ramp_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS')));