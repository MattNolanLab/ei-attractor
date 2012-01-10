function M = WatzStrogatzGraph(N, k, beta)
% WATZSTROGATZGRAPH Create a Watz-Strogatz graph
% Return connection matrix of an undirected graph according to the
% Watz-Strogatz definition (see Watz and Strogatz, Nature, 1998)
% N Number of nodes
% k Node degree; k/2 will be the number of edges
% beta rewiring probability

if (N <= 0)
    error 'Number of vertices must be > 0'
end
if (k > N || k < 0)
    error(sprintf('k must be in the range [0, %d], not %d', ...
        N, k));
end
if (mod(k, 2) == 1)
    error('k must be even');
end
if (beta < 0 || beta > 1)
    error('beta is probability, must be in [0, 1]');
end

M = zeros(N);

for it = 0:N-1
    ind = mod([it-k/2:it-1 it+1:it+k/2], N) + 1;
    M(it+1, ind) = 1;
end

% Rewiring
for k_it = 1:k/2
    for N_it = 1:N
        if (rand() < beta)
            row = M(N_it, :);
            row(N_it) = 1;
            nzero = find(row == 0);
            rewire_idx = nzero(randi(numel(nzero)));
            M(N_it, rewire_idx) = 1;
            M(rewire_idx, N_it) = 1;
            M(N_it, mod(N_it-1+k_it, N) + 1) = 0;
            M(mod(N_it-1+k_it, N) + 1, N_it) = 0;
        end
    end
end

end