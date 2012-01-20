% Dynamic I-V curve of a neuron
close all;
clear leg;

mV = 1e3;
pA = 1e12;
nA = 1e9;

fontSize = 14;

dt = c001_Time(2) - c001_Time(1);

i_start = 0/dt+1;
i_end = 10/dt + 1;

Iin = c003_Current_2(i_start:i_end);
Vm = c002_Membrane_Voltage_2(i_start:i_end);
t = c001_Time(i_start:i_end);

C = 186*1e-12; %pF

% cutoff = 5000;
% Fs = 1/dt;
% filt_N = 4;
% [b, a] = butter(filt_N, cutoff/(Fs/2), 'low');
% Vm = filtfilt(b, a, Vm);


dVm = diff(Vm);
Vm_mean = mean(Vm);

Im = Iin(1:end-1) - C*dVm/dt;

x_lim = [0 4];

figure('Position', [800 0 1000 1000]);

subplot(3,1,1, 'FontSize', fontSize);
plot(t, Iin*nA);
xlabel('Time (s)');
ylabel('Current (nA)');
xlim(x_lim);

subplot(3,1,2, 'FontSize', fontSize);
plot(t, Vm*mV);
xlabel('Time (s)');
ylabel('V_m (mV)');
xlim(x_lim);


subplot(3,1,3, 'FontSize', fontSize);
plot(t(1:end-1), Im*nA);
xlabel('Time (s)');
ylabel('I_m (nA)');
xlim(x_lim);


figure('Position', [800 0 1000 1000]);
leg{1} = 'Original';
plot(t, Vm*mV, 'LineWidth', 1);
hold all;

leg_it = 2;
for cutoff = [1000 3000 5000]
    Fs = 1/dt;
    filt_N = 4;
    [b, a] = butter(filt_N, cutoff/(Fs/2), 'low');
    Vm_filt = filtfilt(b, a, Vm);
    plot(t, Vm_filt*mV, 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('V_m (mV)');
    leg{leg_it} = sprintf('cutoff: %d Hz', cutoff);
    leg_it = leg_it + 1;
    
end
xlim(x_lim);
legend(leg);


% figure;
% [Y f NFFT] = fourierTrans(Vm - mean(Vm), dt);
% Y_abs = 2*abs(Y(1:NFFT/2+1));
% plot(f, Y_abs);

figure;
plot(Vm(1:end-1)*mV, Im*pA, '.');
xlabel('V_m (mV)');
ylabel('I_m (pA)');


figure;
delta_V = 1e-3;
V_fnd = -50e-3;
V_ind = find(abs(Vm(1:end-1) - V_fnd) < delta_V);
Im_V = Im(V_ind); 
hist(Im_V*pA, 40);