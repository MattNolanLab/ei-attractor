from brian.library.electrophysiology import *
from matplotlib.pyplot import *
from IV_curve import *

from scipy.io import loadmat

# Compare kernels and compensated voltages with different levels of capacitance
# neutralisation

mV = 1e-3
pA = 1e-12
pF = 1e-12
ms = 1e-3
MOhm = 1e6
figSize = (10, 12)

dirName = "../data/C_neutralisation/C_values/"
files = [
        ('cell4_007', '-8 pF'),
        ('cell4_006', '-4 pF'),
        ('cell4_008', '0 pF'),
        ('cell4_004', 'small'),
        ('cell4_005', 'normal'),
        ('cell4_002', 'disabled')
        ]

fileIds = [1, 2, 3, 4]
leg = []

for file_it in fileIds:
    
    fileName = files[file_it][0]
    inFile = dirName + fileName + '.mat'
    
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

    figure(1, figsize=figSize)
    subplot(211)
    hold(True)
    plot(times[0:tail_start]/ms, kers.Ke/MOhm)
    xlabel('Time (ms)')
    ylabel('Electrode kernel (MOhm)')

    vm_tstart = 30/dt
    vm_tend = 30.5/dt
    subplot(212)
    plot(times[vm_tstart:vm_tend], Vcomp[vm_tstart:vm_tend])
    
    leg.append(files[file_it][1])

    # Separate kernel plot
    figure(figsize=(8, 6))
    plot(times[0:tail_start]/ms, kers.Ke/MOhm)
    xlabel('Time (ms)')
    ylabel('Electrode kernel (MOhm)')
    savefig(dirName + fileName + '_Ke.pdf')


    # IV curve
    IV_tstart = 20/dt
    t_after = int(200e-3/dt)
    binStart = -80e-3
    binEnd = -43.5e-3
    nbins = 50

    ivc = DynamicIVCurveAfter(I[IV_tstart:], Vcomp[IV_tstart:], times[IV_tstart:],
            binStart, binEnd, nbins, t_after)
    region = (5e-3, 10e-3)
    
    f = figure(figsize=figSize)
    plot(ivc.Vm/mV, ivc.Im/pA, '.', zorder=0)
    xlabel('Membrane voltage (mV)')
    ylabel('Transmembrane current (pA)')
    xlim([-80, -40])
    ylim([-3000, 1000])
    title('Stellate cell I-V relationship')
    grid()
    hold(True)
    errorbar(ivc.binCenters/mV, ivc.Imean/pA, ivc.Istd/np.sqrt(ivc.IN)/pA, None, 'ro', zorder=1)
    hold(False)
    
    f.savefig(dirName + fileName + '_IV_curve.png')
    
    expFit = ivc.fitExponentialIaF()
    figure(figsize=figSize)
    plot(ivc.cur.V/mV, -ivc.cur.I/ivc.cur.C, 'o')
    hold(True)
    plot(ivc.cur.V/mV, expFit.getCurve(ivc.cur.V))
    grid(True)

    savefig(dirName + fileName + '_IV_curve_fit.pdf')
    
#legend(leg)
figure(1)
savefig(dirName + fileName + '_c_values_kernels.pdf')
#show()
    
    
