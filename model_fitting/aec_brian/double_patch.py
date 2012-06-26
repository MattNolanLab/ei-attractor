from brian.library.electrophysiology import *
from matplotlib.pyplot import *
from IV_curve import *

from scipy.io import loadmat

# This is a script which uses Brette et al. Active Electrode Compensation module
# implemented in brian (done offline)

dir = '../../../data/IV_curves/OrnsteinUhlenbeck/Stellate/2012_04_11_IV_curves_recordings 2/IV_double_patch/'
file = "20120410h 002 Copy Export"
inFile = dir + file + '.mat'

mV = 1e-3
pA = 1e-12
pF = 1e-12
ms = 1e-3
MOhm = 1e6
figSize = (8, 6)

data = loadmat(inFile)
times = data['c001_Time'][0]
V = data['c004_Membrane_Voltage_2'][0]
I = data['c005_Current_2'][0]
V_double = data['c002_Membrane_Voltage_1'][0]
I_double = data['c003_Current_1'][0]
dt = times[1] - times[0]


# Compensation
tstart = 0
I_T = 20
ksize = int(15e-3/dt)
tail_start = int(3e-3/dt)
Vk = V[tstart/dt:I_T/dt]
Ik = I[tstart/dt:I_T/dt]

Vcomp, kers, Vmean = compensate_voltage(I, V, Ik, Vk, ksize, tail_start)


# IV curve
IV_tstart = 20/dt
t_after = int(100e-3/dt)
binStart = -70e-3
binEnd = -42e-3
nbins = 100

ivc_after = DynamicIVCurveAfter(I[IV_tstart:], Vcomp[IV_tstart:], times[IV_tstart:],
        binStart, binEnd, nbins, t_after)
region = (5e-3, 10e-3)
#ivc = DynamicIVCurveAfterRegion(I[IV_tstart:], Vcomp[IV_tstart:], times[IV_tstart:],
#        binStart, binEnd, nbins, region, ivc_after.Cest)
ivc = ivc_after


#
# Plot voltages and kernel
#
volt_xlim = [30, 30.100]
volt_ylim = [-80, -20]
volt_figSize = (14, 5)
figure(figsize=volt_figSize)
gcf().subplots_adjust(wspace=0.5)
gridSpec = (1, 5)
ax = subplot2grid(gridSpec, (0, 0), colspan=3)
hold('on')
plot(times[IV_tstart:], V[IV_tstart:]/mV, color='0.8')
plot(ivc.times_all, ivc.Vm_all/mV, linewidth=2)
plot(times[IV_tstart:], V_double[IV_tstart:]/mV, linewidth=2)
xlabel('Time (s)')
ylabel('$V_m$ (mV)')
xlim(volt_xlim)
ylim(volt_ylim)

legend(['Uncomp.', 'Comp.', 'Paired'],loc='best')

ax = subplot2grid(gridSpec, (0, 3), colspan=2)
plot(times[0:tail_start]/ms, kers.Ke/MOhm)
xlabel('Time (ms)')
ylabel('Electrode kernel (M$\Omega$)')

savefig(dir + file + '_voltages.pdf')


f = figure(figsize=figSize)
plot(ivc.Vm/mV, ivc.Im/pA, '.', zorder=0)
xlabel('Membrane voltage (mV)')
ylabel('Transmembrane current (pA)')
xlim([-80, -40])
ylim([-5000, 3000])
title('Stellate cell I-V relationship')
grid()
hold(True)
#plot(ivc.binCenters/mV, ivc.Imean/pA, 'r', zorder=1)
#plot(ivc.binCenters/mV, ivc.Imean/pA, 'ro', zorder=2)
errorbar(ivc.binCenters/mV, ivc.Imean/pA, ivc.Istd/np.sqrt(ivc.IN)/pA, None, 'ro', zorder=1)
hold(False)

f.savefig(dir + file + '_IV_curve.png')

expFit = ivc.fitExponentialIaF()
figure(figsize=figSize)
plot(ivc.cur.V/mV, -ivc.cur.I/ivc.cur.C, 'o')
hold(True)
plot(ivc.cur.V/mV, expFit.getCurve(ivc.cur.V))
grid(True)
savefig(dir + file + '_expFit.pdf')

print "Exp. IaF parameters:"
print expFit



show()


