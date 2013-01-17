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
import numpy.ma as ma
from tables import *

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


def computeMeanSpatialField(sig, sig_t, pos, pos_dt, boxSize, dx):
    '''
    Compute mean spatial field of a quantity in sig, no smoothing is done here.
    Position is taken as 2D

    sig     Signal
    sig_t   Signal times
    pos     Position array of shape (len(pos), 2)
    pos_dt  dt of the pos array
    boxSize A tuple containing dimensions of the box: (X, Y)
    '''

    xedges = np.linspace(-boxSize[0]/2.0, boxSize[0]/2.0, double(boxSize[0])/dx+1)
    yedges = np.linspace(-boxSize[1]/2.0, boxSize[1]/2.0, double(boxSize[1])/dx+1)

    # Extract positions of each point in Vm in fact
    sigPos_i = np.array(sig_t//pos_dt, dtype=int)
    sigPos_x   = pos[sigPos_i, 0]
    sigPos_y   = pos[sigPos_i, 1]

    sigField = np.zeros((len(yedges), len(xedges)))
    sigMask = np.ones(sigField.shape, dtype=bool)

    for y_i in xrange(len(yedges)):
        for x_i in xrange(len(xedges)):
            x = xedges[x_i]
            y = xedges[y_i]
            sig_at_pos_x = np.logical_and(sigPos_x >= x-dx/2.0, sigPos_x <= x+dx/2.0)
            sig_at_pos_y = np.logical_and(sigPos_y >= y-dx/2.0, sigPos_y <= y+dx/2.0)
            sig_fin = sig[np.logical_and(sig_at_pos_x, sig_at_pos_y)]
            sigField[y_i, x_i] = np.mean(sig_fin)
            sigMask[y_i, x_i]  = not (len(sig_fin) > 0)

    return ma.masked_array(sigField, sigMask), (xedges, yedges)



# Neuron parameters
taum      = 9.3 * ms
gL        = 22.73 * nS
C         = taum * gL
EL        = -68.5 * mV
Esyn_e    = 0 * mV
syn_tau_e = 3 * ms
Iext      = 0 * pA


eqs = Equations('''
    dvm/dt      = 1/C*Im                    : volt
    Im          = gL*(EL-vm) + Isyn + Iext  : amp
    Isyn        = ge*(Esyn_e - vm)          : amp
    dge/dt      = -ge/syn_tau_e             : siemens
    ''')


boxSize = (200, 200)
N = (20, 20)
N_total = N[0]*N[1]
totalSz = N[0]*N[1]
maxRates = 15
widths = 40

T = 600*second
dt = 0.2*ms
#traj_dt = 20*ms

#tmp_x = np.linspace(-boxSize[0]/4.0, boxSize[0]/4.0, T/traj_dt+1)
#traj = np.zeros((T/traj_dt + 1, 2))
#traj[:, 0] = tmp_x
ratPosFName = '../../../../data/hafting_et_al_2005/rat_trajectory_original.mat'
traj, traj_dt = loadRatPos(ratPosFName)

totalWeight = 100*nS

simulationClock = Clock(dt)
placeClock = Clock(0.4*msecond)

PC = UniformBoxPlaceCells(boxSize, N, maxRates, widths, random=True)
SG = PlaceSpikingGroup(PC, traj, traj_dt, placeClock)
NrnG = NeuronGroup(1, model=eqs, threshold=np.infty*mV, reset="",
        clock=simulationClock)
NrnG.vm = EL

# All-to-1 connections from place cells onto a grid cell
conn = Connection(SG, NrnG, 'ge')
conn.connect_full(SG, NrnG, weight=totalWeight/PC.N)

net = Network(NrnG, SG, conn)

spikeMon = SpikeMonitor(SG)
stateMon = StateMonitor(NrnG, "vm", record=True, clock=simulationClock)
net.add(spikeMon, stateMon)

net.run(T, report='stdout')

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

figure()
boxSize_dx = 5.0    # cm
sigField, edges = computeMeanSpatialField(stateMon.values.T/mV, stateMon.times,
        traj, traj_dt, boxSize, boxSize_dx)
pcolor(sigField)
colorbar()
axis('equal')


###############################################################################
#                           Export to Igor Pro
###############################################################################
igorExport = True
if igorExport:
    h5file = openFile('output_local/figure_examples_export.h5', mode = "w", title = "Place cell remapping examples")
    
    # Save a raster plot of place cell spikes
    raster_start = 0
    raster_end = N_total
    raster_x = np.ndarray((0))
    raster_y = np.ndarray((0))
    for n_it in xrange(raster_start, raster_end):
        raster_x = np.hstack((raster_x, spikeMon[n_it]))
        raster_y = np.hstack((raster_y, np.zeros((len(spikeMon[n_it]))) + n_it -
            raster_start))
    
    h5file.createArray(h5file.root, 'PC_raster_x', raster_x)
    h5file.createArray(h5file.root, 'PC_raster_y', raster_y)
    
    # Save place cell centers
    h5file.createArray(h5file.root, 'PC_centers_x', PC.centers[:, 0])
    h5file.createArray(h5file.root, 'PC_centers_y', PC.centers[:, 1])
    
    # PC trajectory
    h5file.createArray(h5file.root, 'PC_traj_x', traj[0:traj_used, 0])
    h5file.createArray(h5file.root, 'PC_traj_y', traj[0:traj_used, 1])
    
    # Firing field example
    field_n_id = 0
    field_dx = 2    # cm
    field, positions = PC.getSingleCellFiringField(field_n_id, field_dx)
    h5file.createArray(h5file.root, 'PC_firing_field', field)
    h5file.createArray(h5file.root, 'PC_firing_field_X', positions[0])
    h5file.createArray(h5file.root, 'PC_firing_field_Y', positions[1])
    
    
    # Vm of the GC
    h5file.createArray(h5file.root, 'GC_Vm_time', stateMon.times)
    h5file.createArray(h5file.root, 'GC_Vm_values', stateMon.values.T/mV)
    
    
    # Color plot of average spatial Vm field
    h5file.createArray(h5file.root, 'GC_Vm_spatial', sigField)

    
    
    h5file.close()

show()
