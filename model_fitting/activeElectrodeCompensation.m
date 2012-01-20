function [f] = activeElectrodeCompensation(Vm_rec, Iin, t)

V_rest_range = 50e-3;
spike_th = 0;
after_spike_t = 200e-3; % 200ms

dt = t(2) - t(1);
M = 4e-3/dt;


% First 3 seconds should contain Iin = 0, i.e. Vm_rest
Vm_rest = mean(Vm_rec(1:2.5/dt+1));

% Pick analysis points which are only near Vm_rest and at least 200ms after
% last spike
anal_ids = find(abs(Vm_rec - Vm_rest) < V_rest_range);
spike_ids = findSpikes(Vm_rec, spike_th);

for it = 1:numel(spike_ids)
    spike_t = (spike_ids(it) - 1)*dt;
    
    anal_t = (anal_ids - 1)*dt;
    
    dist = anal_t - spike_t;
    
    rm_id = find(dist > 0 & dist < after_spike_t);
    anal_ids(rm_id) = [];
end


% anal_ids contains indices into Vm_rec than will be used for the combined
% electrode+membrane filter estimation
% See Badel, L., Lefort, S., Brette, R., Petersen, C. C. H., Gerstner, W.,
% & Richardson, M. J. E. (2008). Dynamic I-V curves are reliable predictors
% of naturalistic pyramidal-neuron voltage traces J Neurophysiol, 99(2), 656?666.


Iin_dot = [diff(Iin)/dt 0];
Vm_rec_dot = [diff(Vm_rec)/dt 0];

% Now select only the analysed data
Iin_dot = Iin_dot(anal_ids);
Vm_rec_dot = Vm_rec_dot(anal_ids);
Iin = Iin(anal_ids);
%Vm_rec = Vm_rec(anal_ids);

A = zeros(M);
sigma_Iin = sigma(Iin, Iin);

for j = 0:M-1
    j
    for k = 0:M-1
        A(j+1,k+1) = X(Iin_dot, Iin_dot, j, k) - ...
            S(Iin, Iin_dot, j) * S(Iin, Iin_dot, k) / sigma_Iin;
    end
end

% F ... uncorrected membrane current
F = Vm_rec_dot - sigma(Vm_rec_dot, Iin)/sigma_Iin * Iin;
b = zeros(M, 1);
for j = 0:M-1
    b(j+1) = S(F, Iin_dot, j);
end

f = inv(A) * b;


end