% Generate different OrnStein-Uhlenbeck processes
%
% There are 2 time constants: fast (3ms) and slow (10ms) and for each
% protocol, these two OU processess are summed up together to produce one
% waveform
% All current units are in pA

close all;
clear all;

noise_range = 1000;
noise_T = 20;

tau_fast = 3*1e-3;
tau_slow = 10*1e-3;

sigma_fast = [0.18 0.36];
sigma_slow = [0.18 0.25];
sig = 2;

DC_bias = [0.06 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
currentFactor = [550 650 750];

dt = 0.05e-3;
T = 40;

% Estimated kernel size - in order to use fast algorithm for kernel
% estimation, last ker_sz points of uniform noise input current must be 0
ker_sz = 50e-3/dt;


X = [];

for curr_f = 1:numel(currentFactor)
    it = 1;
    for bias = 1:numel(DC_bias)
        X_fast = OrnsteinUhlenbeckProcess(tau_fast, 0, ...
            sigma_fast(sig), dt, T);
        X_slow = OrnsteinUhlenbeckProcess(tau_slow, 0, ...
            sigma_slow(sig), dt, T);
    
        X = X_fast + X_slow + DC_bias(bias);
        X = X * currentFactor(curr_f);
        
        noise_uni = -noise_range + 2*noise_range * ...
            rand(1, noise_T/dt);
        noise_uni(end-ker_sz+1:end) = 0;
        
        X = [noise_uni X];
        t = [0:numel(X)-1] * dt;
    
        
        bias_abs = DC_bias(bias)*currentFactor(curr_f);
        fileName = sprintf('OrnsteinUhlenbeck_uni_%dpA_%dpA_DC_%dpA', ...
                fix(noise_range), fix(currentFactor(curr_f)), ...
                fix(bias_abs))
        save([fileName '.mat'], 'X', 't');
        dlmwrite([fileName '.txt'], [t' X'])
    
        it = it+1;
    end
end
