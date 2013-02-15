#
#   fig_spike_maps.py
#
#   Spike and rate maps of grid cells. Data analysis
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np

from scipy.io           import loadmat
from scipy.io           import savemat
from matplotlib.pyplot  import *
from tables             import *

from grid_cell_analysis import *

mV = 1e3
pA = 1e12
nA = 1e9


jobNum = 4400
trialNum = 0
dumpNum = 9


rcParams['font.size'] = 14


arenaDiam = 180.0     # cm
h = 3.0

# Neuron to extract spikes from
neuronNum = 10
spikeType = 'excitatory'


dirName = "output_local/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"

gridnessScores = []



print 'jobNum: ' + str(jobNum)

fileName = fileNameTemp
fileName = fileName.format(dirName, fileNamePrefix, jobNum,
                            trialNum, dumpNum)
try:
    data = loadmat(fileName +  '_output.mat')
except:
    print "warning: could not open: " + fileName
    exit(1)

simT       = data['options']['time'][0][0][0][0]
tend       = simT/1e3
tstart     = tend - 60

times      = data['stateMon_times'].ravel()
dt         = times[1] - times[0]
times      = times[tstart/dt:tend/dt]
Vm_e_0     = data['stateMon_e_values'][tstart/dt:tend/dt, 0]
Vm_e_1     = data['stateMon_e_values'][tstart/dt:tend/dt, 1]
Vm_i_0     = data['stateMon_i_values'][tstart/dt:tend/dt, 0]
Vm_i_1     = data['stateMon_i_values'][tstart/dt:tend/dt, 1]
Iclamp_e_0 = data['stateMon_Iclamp_e_values'][tstart/dt:tend/dt, 0]
Iclamp_e_1 = data['stateMon_Iclamp_e_values'][tstart/dt:tend/dt, 1]
Iclamp_i_0 = data['stateMon_Iclamp_i_values'][tstart/dt:tend/dt, 0]
Iclamp_i_1 = data['stateMon_Iclamp_i_values'][tstart/dt:tend/dt, 1]


h5file = openFile(fileName + '_voltage_igor.h5', mode = "w", title =
        "Voltage and current clamp export")

h5file.createArray(h5file.root, 'times'     , times)
h5file.createArray(h5file.root, 'Vm_e_0'    , Vm_e_0*mV)
h5file.createArray(h5file.root, 'Vm_e_1'    , Vm_e_1*mV)
h5file.createArray(h5file.root, 'Vm_i_0'    , Vm_i_0*mV)
h5file.createArray(h5file.root, 'Vm_i_1'    , Vm_i_1*mV)
h5file.createArray(h5file.root, 'Iclamp_e_0', Iclamp_e_0*nA)
h5file.createArray(h5file.root, 'Iclamp_e_1', Iclamp_e_1*nA)
h5file.createArray(h5file.root, 'Iclamp_i_0', Iclamp_i_0*nA)
h5file.createArray(h5file.root, 'Iclamp_i_1', Iclamp_i_1*nA)

h5file.close()

