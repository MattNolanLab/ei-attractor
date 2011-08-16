% Test whether the high frequency simulation states are stable over a
% longer period --> generate lots of simulations


% Simulate integrate and fire neuron
clear all;
close all;

path('..', path);


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
opt.Ie = 18.5e-3;
%we = 1/(taue*1000) * 700e-3 / N; % normalize the weight by time constant to inject constant charge
opt.we = 180e-3 / N;


% Inhibitory cell
opt.taum_i = 10e-3;
opt.taui = 5e-3;
opt.El_i = -60e-3;
opt.Vt_i = -50e-3;
opt.Vr_i = -58e-3;
opt.i_sparseness = 0.75;
opt.Ii = 6e-3;
%wi = 1/(taui*1000) * 600e-3 / N;
opt.wi = 15e-3 / N;


% Noise normalized per time unit (ms)
opt.noise_sigma = 0.02e-3;
opt.sigma_init_cond = 2e-3; % 2mV


% Euler settings
opt.dt = 0.5e-3  % 0.1 ms
dt = opt.dt;


% Firing rate sliding window length
opt.rateWindowLen = 0.002; %ms
rateWindowLen = opt.rateWindowLen;

% Vm monitor, neuron index
opt.Emon_i = 100;
opt.Imon_i = 100;

% simulation time
opt.T = 60;


% 
% Create simulation results
%
nTrials = 300;

% The same network structure for all simulations
[net_data.Me net_data.Mi] = MeMi(opt);

param_i = 1;
for Ie = 19.5e-3
    opt.Ie = Ie;
    Ie
    
    parfor trialNum = 1:nTrials
        trialNum
        [spikeRecord_e, spikeRecord_i, Vmon, times] = simulateEI(opt, net_data);
        
        tmpresults(trialNum).spikeRecord_e = spikeRecord_e;
        tmpresults(trialNum).spikeRecord_i = spikeRecord_i;
        tmpresults(trialNum).Vmon = Vmon;
        tmpresults(trialNum).times = times;

        
        tmpresults(trialNum).firingRate_e = getFiringRate(spikeRecord_e, dt, rateWindowLen);
        tmpresults(trialNum).spikeCell_e = spikeRecordToSpikeCell(spikeRecord_e, times);
        tmpresults(trialNum).spikeCell_i = spikeRecordToSpikeCell(spikeRecord_i, times);

        tmpresults(trialNum).opt = opt;
    end
    
    results(param_i, :) = tmpresults;
    
    param_i = param_i + 1;
end

clear tmpresults;
save('-v7.3', sprintf('001_high_freq_stability_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS')));
