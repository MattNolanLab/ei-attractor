function [Me, Mi] = MeMi(opt)
    % Setup connections Mij: j --> i
    % Assuming constant and uniform excitatory/inhibitory weights
    Me = rand(opt.Ni, opt.Ne);
    Mi = rand(opt.Ne, opt.Ni);
    Me = double(Me <= opt.e_sparseness);
    Mi = double(Mi <= opt.i_sparseness);
end