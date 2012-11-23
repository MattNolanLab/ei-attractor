#
#   simulate_one_gc.py
#
#   Simulate one grid cell receiving inputs from a place cell sheet
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

from scipy.io import loadmat

import brian
from brian import *

from PlaceSpikingGroup import *

def loadRatPos(fname):
    ratData    = loadmat(fname)
    rat_dt     = ratData['dt'][0][0]

    rat_pos_x  = ratData['pos_x'].ravel()
    rat_pos_y  = ratData['pos_y'].ravel()
    pos = np.vstack((rat_pos_x, rat_pos_y)).T

    return pos, rat_dt

# Neuron parameters
taum      = 9.3 * ms
gL        = 22.73 * nS
C         = taum * gL
EL        = -68.5 * mV
Esyn_e    = 0 * mV
syn_tau_e = 2 * ms
Iext      = 0 * pA


eqs = Equations('''
    dvm/dt      = 1/C*Im                    : volt
    Im          = gL*(EL-vm) + Isyn + Iext  : amp
    Isyn        = ge*(Esyn_e - vm)          : amp
    dge/dt      = -ge/syn_tau_e             : siemens
    ''')


boxSize = (400, 400)
N = (20, 20)
totalSz = N[0]*N[1]
maxRates = 15
widths = 40

T = 60*second
dt = 0.2*ms
#traj_dt = 20*ms

#tmp_x = np.linspace(-boxSize[0]/4.0, boxSize[0]/4.0, T/traj_dt+1)
#traj = np.zeros((T/traj_dt + 1, 2))
#traj[:, 0] = tmp_x
ratPosFName = '../../../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat'
traj, traj_dt = loadRatPos(ratPosFName)

totalWeight = 100*nS

simulationClock = Clock(dt)

PC = UniformBoxPlaceCells(boxSize, N, maxRates, widths, random=True)
SG = PlaceSpikingGroup(PC, traj, traj_dt, simulationClock)
NrnG = NeuronGroup(1, model=eqs, threshold=np.infty*mV, reset="",
        clock=simulationClock)
NrnG.vm = EL

# All-to-1 connections from place cells onto a grid cell
conn = Connection(SG, NrnG, 'ge')
conn.connect_full(SG, NrnG, weight=totalWeight/PC.N)

net = Network(NrnG, SG, conn)

spikeMon = SpikeMonitor(SG)
stateMon = StateMonitor(NrnG, "vm", record=True)
net.add(spikeMon, stateMon)

net.run(T)

raster_plot(spikeMon)

figure()
plot(stateMon.times, stateMon.values.T/mV)
xlabel("Time (s)")
ylabel("$V_m$ (mV)")

figure()
plot(PC.centers[:, 0], PC.centers[:, 1], 'x')
hold('on')
traj_used = T//traj_dt+1
plot(traj[0:traj_used, 0], traj[0:traj_used, 1], 'r')

show()

