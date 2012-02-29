from brian.library.electrophysiology import *
from matplotlib.pyplot import *

from scipy.io import loadmat
from numpy import diff
from numpy.random import random
from scipy.optimize import fmin

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

def voltage_bins(V, binStart, binEnd, nbins):
    '''
    Picks indices of V that are within bin ranges.
    Samples outside binStart and binEnd are discarded
    '''
    if binStart > binEnd:
        raise Exception("binStart must be less than binEnd")
    binCenters = np.linspace(binStart, binEnd, nbins+1)
    binWidth = (binEnd - binStart)/nbins
    Ibin = []
    for c in binCenters:
        Ibin.append((np.logical_and(V > c-binWidth/2, V <=
            c+binWidth/2)).nonzero()[0]) 

    return (binCenters, Ibin)


def findCapacitance(Iin, dVdt, C0):
    '''
    Estimate membrane capacitance, by using input current and dV/dt at a
    specific voltage value (set by the user, ideally somewhere near resting Vm)
    '''
    return fmin(lambda Ce: np.var(Iin/Ce - dVdt), C0, xtol=0.00001)[0]


###############################################################################

dir = "../data/C_neutralisation/2012_02_16/"
file = "cell4 008 Copy Export"
inFile = dir + file + '.mat'

mV = 1e3
pA = 1e12
pF = 1e12
ms = 1e-3
MOhm = 1e6
figSize = (8, 6)

data = loadmat(inFile)
times = data['c001_Time'][0]
V = data['c002_Membrane_Voltage_2'][0]
I = data['c003_Current_2'][0]
dt = times[1] - times[0]    # Assuming constant time step


# Full kernel estimation procedure
tstart = 0
I_T = 20
ksize = int(15e-3/dt)
Vk = V[tstart/dt:I_T/dt]
Ik = I[tstart/dt:I_T/dt]
K_full, V_mean = full_kernel(Vk, Ik, ksize, full_output=True)

# Electrode kernel
start_tail = int(3e-3/dt)
Ke, Km = electrode_kernel(K_full, start_tail, full_output=True)

figure(figsize=figSize)
subplot(311)
plot(times[0:ksize]/ms, K_full/MOhm)
ylabel('Kernel(M$\Omega$)')
subplot(312)
plot(times[0:start_tail]/ms, Ke/MOhm)
ylabel('Electrode kernel (M$\Omega$)')
subplot(313)
plot(times[0:ksize]/ms, Km/MOhm)
ylabel('Membrane kernel ($\Omega$)')
xlabel('Time (ms)')
savefig(dir + file + '_kernels.pdf')


# Compensation
Vcorr = AEC_compensate(V, I, Ke)

comp_xlim = [20, 30]
figure(figsize=figSize)
subplot(311)
plot(times, V*mV)
ylabel('Uncomp. $V_m$ (mV)')
xlim(comp_xlim)
subplot(312)
plot(times, I*pA)
ylabel('Injected current $I_{in}$ (pA)')
xlim(comp_xlim)
subplot(313)
plot(times, Vcorr*mV)
ylabel('Compensated $V_m$')
xlabel('Time (s)')
xlim(comp_xlim)
savefig(dir + file + '_compensation.pdf')


# Compensation detail
comp_detail_x0 = 30
comp_detail_x1 = 30.5
figure()
subplot(211)
plot(times, V*mV)
ylabel('Uncomp. $V_m$ (mV)')
xlim([comp_detail_x0, comp_detail_x1])
subplot(212)
plot(times, Vcorr*mV)
ylabel('Compensated $V_m$')
xlabel('Time (s)')
xlim([comp_detail_x0, comp_detail_x1])
savefig(dir + file + '_compensation_detail.pdf')



# IV curve
IV_tstart = 20/dt
C0 = 200e-12 # Initial capacitance estimation
#Cmax = 300e-12
t_after = int(200e-3/dt)
time_id = pickVPreSpike(Vcorr, t_after)
time_id = time_id[0:len(time_id)-1]
time_id = time_id[time_id > IV_tstart]

dV_pre = np.diff(Vcorr)[time_id]
Vm_pre = Vcorr[time_id]
Iin_pre = I[time_id]

# Sort voltage indices as a function of voltage (binned)
binStart = -80e-3
binEnd = -40e-3
nbins = 100
binCenters, binIds = voltage_bins(Vm_pre, binStart, binEnd, nbins)

Cest_bin_id = 27
C = findCapacitance(Iin_pre[binIds[Cest_bin_id]], dV_pre[binIds[Cest_bin_id]]/dt,
        C0)
print("Estimated capacitance: " + str(C*pF) + " pF.")
Im_pre = C*dV_pre/dt - Iin_pre

#figure(figsize=figSize)
#plot(times, Vcorr*mV)
#hold(True)
#plot(times[time_id], Vm_pre*mV)
#xlabel('Time (s)')
#ylabel('$V_m$ (mV)')

f = figure(figsize=figSize)
plot(Vm_pre*mV, Im_pre*pA, '.', zorder=0)
xlabel('Membrane voltage (mV)')
ylabel('Membrane current (pA)')
xlim([-80, -40])
ylim([-1000, 2000])
title('Stellate cell I-V relationship')
grid()
hold(True)

Imean = np.zeros(len(binCenters))
Istd  = np.zeros(len(binCenters))
IN    = np.zeros(len(binCenters))
for i in xrange(len(binCenters)):
    Imean[i] = np.mean(Im_pre[binIds[i]])
    Istd[i]  = np.std(Im_pre[binIds[i]])
    IN[i]    = len(binIds[i])
errorbar(binCenters*mV, Imean*pA, Istd/np.sqrt(IN)*pA, None, 'ro', zorder=1)
hold(False)

f.savefig(dir + file + '_IV_curve.png')


distrib_it = 27
figure(figsize=figSize)
h = hist(Im_pre[binIds[distrib_it]]*pA, 100, normed=True)
title('Current distribution at $V_m$ = ' + str(binCenters[distrib_it]*mV) + ' mV.')
xlabel('Membrane current (pA)')
ylabel('Frequency')

savefig(dir + file + '_current_distrib.pdf')


show()


