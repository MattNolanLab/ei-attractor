#
#   simulation.py
#
#   Runs the gamma network simulation
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

import pylab    as pl
import numpy    as np
import logging  as lg
import nest
import time

from numpy      import exp
from nest       import raster_plot

from network    import *
from parameters import getOptParser

lg.basicConfig(level=lg.DEBUG)


# Parse command line parameters
parser = getOptParser()
parser.add_option("--print_time",   action  = "store_true",   help="Print time-progress of the simulation")
parser.add_option("--numThreads",   type    = "int",          help="Number of threads in the simulation")
(o, args) = parser.parse_args()

# Other parameters
N_rec       = 50 # record from 50 neurons
record_from = ['V_m', 'g_AMPA', 'g_GABA_A', 'I_stim', 'g_AHP', 'I_clamp_AMPA', 'I_clamp_GABA_A']


startbuild= time.time()
print "Setting up network."

np.random.seed(1234)

net = GammaNetwork(o)
net.connRandom()
#net.setConstantCurrent()
net.setRampCurrent()

# Recordings
net.monitorEState([0], record_from)
net.monitorIState([0], record_from)
net.monitorESpikes(range(N_rec))
net.monitorISpikes(range(N_rec))

endbuild=time.time()


print "Simulating."

nest.Simulate(o.time)

endsimulate= time.time()

build_time = endbuild-startbuild
sim_time   = endsimulate-endbuild

nest.raster_plot.from_device(net._espikes, hist=True)
nest.raster_plot.from_device(net._ispikes, hist=True)

# obtain and display data
meter_e = net._meter_e
meter_i = net._meter_i

for EI in ['E', 'I']:
    pl.figure(figsize=(12, 12))
    if EI == 'E':
        pl.title('Excitatory cell')
        events = nest.GetStatus(meter_e)[0]['events']
    else:
        pl.title('Inhibitory cell')
        events = nest.GetStatus(meter_i)[0]['events']

    t = events['times'];
    ax = pl.subplot(511)
    pl.plot(t, events['V_m'])
    pl.ylabel('$V_m$ [mV]')
    
    pl.subplot(512, sharex=ax)
    pl.plot(t, events['g_AMPA'], t, events['g_GABA_A'])
    pl.ylabel('Syn. cond. [nS]')
    pl.legend(('g_AMPA', 'g_GABA_A'), loc='best')

    pl.subplot(513, sharex=ax)
    pl.plot(t, events['I_clamp_AMPA'], t, events['I_clamp_GABA_A'])
    pl.ylabel('$I_{syn}$ (pA)')
    pl.legend(('I_clamp_AMPA', 'I_clamp_GABA_A'), loc='best')

    pl.subplot(514, sharex=ax)
    pl.plot(t, events['g_AHP'])
    pl.ylabel('$g_{AHP}$ (nS)')

    
    pl.subplot(515, sharex=ax)
    pl.plot(t, events['I_stim'])
    pl.xlabel('Time [ms]')
    pl.ylabel('$I_{ext}$ (pA)')


#nest.raster_plot.show()
pl.show()
