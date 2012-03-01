from brian.library.electrophysiology import *
from matplotlib.pyplot import *
from IV_curve import *

from scipy.io import loadmat

# This is a script which uses Brette et al. Active Electrode Compensation module
# implemented in brian (done offline)


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
dt = times[1] - times[0]


# Compensation
tstart = 0
I_T = 20
ksize = int(15e-3/dt)
tail_start = int(3e-3/dt)
Vk = V[tstart/dt:I_T/dt]
Ik = I[tstart/dt:I_T/dt]

Vcomp, kers, Vmean = compensate_voltage(I, V, Ik, Vk, ksize, tail_start)


figure(figsize=figSize)
subplot(311)
plot(times[0:ksize]/ms, kers.Kfull/MOhm)
ylabel('Kernel(M$\Omega$)')
subplot(312)
plot(times[0:tail_start]/ms, kers.Ke/MOhm)
ylabel('Electrode kernel (M$\Omega$)')
subplot(313)
plot(times[0:ksize]/ms, kers.Km/MOhm)
ylabel('Membrane kernel ($\Omega$)')
xlabel('Time (ms)')
savefig(dir + file + '_kernels.pdf')


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
plot(times, Vcomp*mV)
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
plot(times, Vcomp*mV)
ylabel('Compensated $V_m$')
xlabel('Time (s)')
xlim([comp_detail_x0, comp_detail_x1])
savefig(dir + file + '_compensation_detail.pdf')



# IV curve
IV_tstart = 20/dt
t_after = int(200e-3/dt)
binStart = -80e-3
binEnd = -40e-3
nbins = 100

ivc = DynamicIVCurveAfter(I[IV_tstart:], Vcomp[IV_tstart:], times[IV_tstart:],
        binStart, binEnd, nbins, t_after)
#region = (5, 100)
#ivc = DynamicIVCurveAfterRegion(I[IV_tstart:], Vcomp[IV_tstart:], times[IV_tstart:],
        binStart, binEnd, nbins, region)

figure(figsize=figSize)
plot(ivc.times_all, ivc.Vm_all*mV)
hold(True)
plot(ivc.times, ivc.Vm*mV)
xlabel('Time (s)')
ylabel('$V_m$ (mV)')

f = figure(figsize=figSize)
plot(ivc.Vm*mV, ivc.Im*pA, '.', zorder=0)
xlabel('Membrane voltage (mV)')
ylabel('Membrane current (pA)')
xlim([-80, -40])
ylim([-1000, 2000])
title('Stellate cell I-V relationship')
grid()
hold(True)
errorbar(ivc.binCenters*mV, ivc.Imean*pA, ivc.Istd/np.sqrt(ivc.IN)*pA, None, 'ro', zorder=1)
hold(False)

f.savefig(dir + file + '_IV_curve.png')


distrib_it = 27
figure(figsize=figSize)
h = hist(ivc.Im[ivc.binIds[distrib_it]]*pA, 100, normed=True)
title('Current distribution at $V_m$ = ' + str(ivc.binCenters[distrib_it]*mV) + ' mV.')
xlabel('Membrane current (pA)')
ylabel('Frequency')

savefig(dir + file + '_current_distrib.pdf')


show()


