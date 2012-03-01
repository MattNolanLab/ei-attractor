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
figSize = (8, 6)

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

    hold(True)
    figure(1, figsize=figSize)
    plot(times[0:tail_start]/ms, kers.Ke/MOhm)
    xlabel('Time (ms)')
    ylabel('Electrode kernel (MOhm)')

    vm_tstart = 30/dt
    vm_tend = 30.5/dt
    figure(2, figsize=figSize)
    plot(times[vm_tstart:vm_tend], Vcomp[vm_tstart:vm_tend])
    
    leg.append(files[file_it][1])
    
    
figure(1)
legend(leg)
figure(2)
legend(leg)
show()
    
    
