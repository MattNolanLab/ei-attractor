from brian.library.electrophysiology import *
from matplotlib.pyplot import *

from scipy.io import loadmat
from numpy import diff

# This is a script which uses Brette et al. Active Electrode Compensation module
# implemented in brian (done offline)

###############################################################################
def pickVPreSpike(Vm, t_after, V_th=0):
    '''
    Pick only voltage traces that are more than t_after time units after the
    last spike. By setting V_th, one can fine_tune spike threshold.
    '''
    spike_times = spike_peaks(Vm, V_th)
    time_ids = np.arange(len(Vm))
    t_mask = np.ndarray(len(Vm), dtype=bool)
    t_mask[:] = True
    for t_spike in spike_times:
        t_mask = np.logical_and(t_mask, np.logical_or(time_ids <= t_spike,
            time_ids > t_spike+t_after))

    return time_ids[t_mask]
###############################################################################

dir = "../data/C_neutralisation/2012_02_16/"
file = "cell4 008 Copy Export"
inFile = dir + file + '.mat'

mV = 1e3
pA = 1e12
ms = 1e3
figSize = (12, 8)

data = loadmat(inFile)
times = data['c001_Time'][0]
V = data['c002_Membrane_Voltage_2'][0]
I = data['c003_Current_2'][0]
dt = times[1] - times[0]    # Assuming constant time step


# Full kernel estimation procedure
tstart = 0
I_T = 20
ksize = int(10e-3/dt)
Vk = V[tstart/dt:I_T/dt]
Ik = I[tstart/dt:I_T/dt]
K_full, V_mean = full_kernel(Vk, Ik, ksize, full_output=True)

# Electrode kernel
start_tail = int(3e-3/dt)
Ke, Km = electrode_kernel(K_full, start_tail, full_output=True)

figure(figsize=figSize)
subplot(311)
plot(times[0:ksize], K_full)
ylabel('Kernel($\Omega$)')
subplot(312)
plot(times[0:start_tail], Ke)
ylabel('Electrode kernel ($\Omega$)')
subplot(313)
plot(times[0:ksize], Km)
ylabel('Membrane kernel ($\Omega$)')
xlabel('Time (s)')


# Compensation
Vcorr = AEC_compensate(V, I, Ke)
figure(figsize=figSize)
subplot(311)
plot(times, V*mV)
ylabel('Uncomp. $V_m$ (mV)')
subplot(312)
plot(times, I*pA)
ylabel('Injected current $I_{in}$ (pA)')
subplot(313)
plot(times*ms, Vcorr*mV)
ylabel('Compensated $V_m$')
xlabel('Time (s)')


# IV curve
IV_tstart = 20/dt
C = 220e-12 # TODO: C estimation
t_after = int(200e-3/dt)
time_id = pickVPreSpike(Vcorr, t_after)
time_id = time_id[0:len(time_id)-1]
time_id = time_id[time_id > IV_tstart]

dV_pre = np.diff(Vcorr)[time_id]
Vm_pre = Vcorr[time_id]
Iin_pre = I[time_id]
Im_pre = C*dV_pre/dt - Iin_pre

figure(figsize=figSize)
plot(times, Vcorr*mV)
hold(True)
plot(times[time_id], Vm_pre*mV)
xlabel('Time (s)')
ylabel('$V_m$ (mV)')

f = figure(figsize=figSize)
plot(Vm_pre*mV, Im_pre*pA, '.')
xlabel('Membrane voltage (mV)')
ylabel('Membrane current (pA)')
xlim([-80, -40])
ylim([-1000, 2000])
title('Stellate cell I-V relationship')

f.savefig(dir + file + '_IV_curve.png')



show()


