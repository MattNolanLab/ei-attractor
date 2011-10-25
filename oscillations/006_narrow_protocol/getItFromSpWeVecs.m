function [it] = getItFromSpWeVecs(sp_val, we_val, res, eps)

    eps = min(min(res(1,1).opt.we_vec), min(res(1,1).opt.sparseness_vec)) * eps;
    found = false;
    
    for it = 1:size(res, 1)
        if (abs(res(it, 1).opt.we - we_val) < eps && ...
                abs(res(it, 1).opt.e_sparseness - sp_val) < eps)
            found = true;
            break;
        end
    end
    
    if (~found)
        it = -1;
    end
end