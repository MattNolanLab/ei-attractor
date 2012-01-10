function ok = checkGraphColoring(V, M, C)
% Check if the graph M contains proper coloring in V, with C the number of
% colors

ok = 1;

for it = 1:numel(V)
    conn_colors = M(it, :) .* V;
    if (nnz(conn_colors == V(it)))
        ok = 0;
        break;
    end
end


end