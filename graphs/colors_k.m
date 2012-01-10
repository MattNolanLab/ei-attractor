% Create color distributions based on average node degree parameter

close all;
clear all;

N = 100;
beta = 1;
kMax = 90;
nTrials = 25;

ks = 2:2:kMax;
C_all = zeros(nTrials, numel(ks));

for k_it = 1:numel(ks)
    k = ks(k_it)
    for trial_it = 1:nTrials
        M = WattsStrogatzGraph(N, k, beta);
        [C V] = greedyColouring(M);
        if (checkGraphColoring(V, M, C) == 0)
            error('Coloring failed: improper coloring produced!');
        end
        
        C_all(trial_it, k_it) = C;
    end
end

plot(ks, mean(C_all), 'o');