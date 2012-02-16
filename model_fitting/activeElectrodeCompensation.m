function [f] = activeElectrodeCompensation(Vm_rec, Iin, t)


dt = t(2) - t(1);
M = 10e-3/dt;

Iin = Iin*1e9;
Vm_rec = Vm_rec*1e3;

A = zeros(M+1);
b = zeros(M+1, 1);

for i = 0:M-1
    i
    for j = 0:M-1
        i1 = [zeros(1, j) Iin(1:end-j)];
        i2 = [zeros(1, i) Iin(1:end-i)];
        A(i+1, j+1) = mean(i1 .* i2);
    end
end

for i = 0:M-1
    tmp = mean([zeros(1, i) Iin(1:end-i)]);
    A(i+1, M+1) = tmp;
    A(M+1, i+1) = tmp;
end

A(M+1, M+1) = 1;

for i = 0:M-1
    Itmp = [zeros(1, i) Iin(1:end-i)];
    b(i+1) = mean(Vm_rec .* Itmp);
end
b(M+1) = mean(Vm_rec);

f = inv(A) * b;


% fit the kernel tail with an exponential
end