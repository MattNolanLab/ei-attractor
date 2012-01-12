function M = shuffleGraph(M)
% Take an incidence matrix of an undirected graph and just shuffle the
% nodes

N = size(M, 1);

for it = 1:100*N
    swap_it1 = randi(N);
    swap_it2 = randi(N);
    
    tmpRow = M(swap_it1, :);
    M(swap_it1, :) = M(swap_it2, :);
    M(swap_it2, :) = tmpRow;
    
    tmpCol = M(:, swap_it1);
    M(:, swap_it1) = M(:, swap_it2);
    M(:, swap_it2) = tmpCol;
end



end