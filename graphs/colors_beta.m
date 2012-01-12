% Chromatic number of a Watts-Strogatz graph as a function of the rewiring
% probability

% Create color distributions based on average node degree parameter

close all;
clear all;

fontSize = 14;

N = 100;
%beta = 1;
k = N / 10;
nTrials = 100;

betas = 1;%0:0.1:1

C_all = zeros(nTrials, numel(betas));

for beta_it = 1:numel(betas)
    beta = betas(beta_it)
    for trial_it = 1:nTrials
        M = WattsStrogatzGraph(N, k, beta);
        [C V] = greedyColouring(M);
        if (checkGraphColoring(V, M, C) == 0)
            error('Coloring failed: improper coloring produced!');
        end
        
        C_all(trial_it, beta_it) = C;
    end
end

subplot(1,1,1, 'FontSize', fontSize);
plot(betas, mean(C_all), 'o-');
xlabel('\beta');
ylabel('Chromatic number (greedy coloring)');
