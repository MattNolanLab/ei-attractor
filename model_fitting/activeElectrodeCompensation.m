function [f] = activeElectrodeCompensation(Vm_rec, Iin, t)


dt = t(2) - t(1);
M = 10e-3/dt;

Iin = Iin*1e9;
Vm_rec = Vm_rec*1e3;

% First 3 seconds should contain Iin = 0, i.e. Vm_rest
Vm_rest = mean(Vm_rec);

% See Badel, L., Lefort, S., Brette, R., Petersen, C. C. H., Gerstner, W.,
% & Richardson, M. J. E. (2008). Dynamic I-V curves are reliable predictors
% of naturalistic pyramidal-neuron voltage traces J Neurophysiol, 99(2), 656?666.


Iin_dot = [diff(Iin)/dt 0];
Vm_rec_dot = [diff(Vm_rec)/dt 0];


A = zeros(M);
sigma_Iin = sigma(Iin, Iin);

for j = 0:M-1
    j
    for k = 0:M-1
        A(j+1,k+1) = X(Iin_dot, Iin_dot, j, k) - ...
            S(Iin, Iin_dot, j) * S(Iin, Iin_dot, k) / sigma_Iin;
    end
end
b(M+1) = mean(Vm_rec);

f = inv(A) * b;


% fit the kernel tail with an exponential
end