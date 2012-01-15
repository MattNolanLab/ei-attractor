% Generate different OrnStein-Uhlenbeck processes
%
% There are 2 time constants: fast (3ms) and slow (10ms) and for each
% protocol, these two OU processess are summed up together to produce one
% waveform

close all;
clear all;

pA = 1e12;

tau_fast = 3*1e-3;
tau_slow = 10*1e-3;

sigma_fast = [0.18 0.36];
sigma_slow = [0.18 0.25];

DC_bias = [0 0.02 0.03 0.06];
currentFactor = [250 350 450 550 650 750] * 1e-12;

dt = 0.05e-3;
T = 40;
silence_period = 3;

X = [];
silence = zeros(1, numel(0:dt:silence_period-dt));

for curr_f = 1:numel(currentFactor)
    it = 1;
    for sig = 1:2
         for bias = 1:4
            X_fast = OrnsteinUhlenbeckProcess(tau_fast, 0, ...
                sigma_fast(sig), dt, T);
            X_slow = OrnsteinUhlenbeckProcess(tau_slow, 0, ...
                sigma_slow(sig), dt, T);

            X = X_fast + X_slow + DC_bias(bias);
            X = X * currentFactor(curr_f);
            X = [silence X silence];
            t = [0:numel(X)-1] * dt;

            
            save(sprintf('OrnsteinUhlenbeck_%dpA_%d.mat', ...
                fix(currentFactor(curr_f)*pA), it), 'X', 't');

            it = it+1;
        end
    end
end