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

def voltage_bins(V, I, binStart, binEnd, nbins):
    '''
    Average currents in each voltage bin. V and I must be the same length.
    Samples outside binStart and binEnd are discarded
    '''




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
hold(True)
traceLens = [0.1, 0.3, 0.5, 1.0]
for traceLen in traceLens:
    tstart = 20 - traceLen
    I_T = 20
    ksize = int(10e-3/dt)
    Vk = V[tstart/dt:I_T/dt]
    Ik = I[tstart/dt:I_T/dt]
    K_full, V_mean = full_kernel(Vk, Ik, ksize, full_output=True)
    
    # Electrode kernel
    start_tail = int(3e-3/dt)
    Ke, Km = electrode_kernel(K_full, start_tail, full_output=True)
    
    figure(1)
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
    figure(2)
    plot(times, Vcorr*mV)
    xlabel('Time (s)')
    ylabel('$V_m$ (mV)')
    xlim([30, 31])


show()


