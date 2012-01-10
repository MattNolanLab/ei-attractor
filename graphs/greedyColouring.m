function [C V] = greedyColouring(M)
% GREEDYCOLOURING Calculate chromatic number of a graph given by M, using
% greedy algorithm

sz = size(M);
if (sz(1) ~= sz(2) || numel(sz) ~= 2)
    error('Incidence matrix mus be square!');
end

N = sz(1);

% Order the vertices by degree
degs = sum(M);
[sortDegs sortDegs_I] = sort(degs, 'descend');

% Assign color
% Colors are 1...C, zero means unassigned color.
% Therefore the chromatic number will be C
node_colors = zeros(1, N);
C = 0;

for unsort_it = 1:N
    it = sortDegs_I(unsort_it);
    
    avail_color = find((1 - M(it, :)) .* node_colors ~= 0);
    if (numel(avail_color) == 0)
        % new color
        C = C + 1;
        node_colors(it) = C;
    else
        new_c_flag = 1;
        for av_c_it = 1:numel(avail_color)
            if (nnz(M(it, :) .* node_colors == node_colors(avail_color(av_c_it))) == 0) 
                node_colors(it) = node_colors(avail_color(av_c_it));
                new_c_flag = 0;
                break;
            end
        end
        
        if(new_c_flag == 1)
            C = C + 1;
            node_colors(it) = C;
        end
    end
end

V = node_colors;
end