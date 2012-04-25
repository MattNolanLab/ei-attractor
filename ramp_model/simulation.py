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
# This version uses numpy to draw the random connections.

import pylab as pl
import numpy as np
import nest
import time

from numpy      import exp
from nest       import raster_plot

from network    import *
from parameters import getOptParser


# Parse command line parameters
parser = getOptParser()
(o, args) = parser.parse_args()

# Other parameters
N_rec     = 50 # record from 50 neurons


startbuild= time.time()

net = GammaNetwork(o)


espikes = nest.Create("spike_detector")
ispikes = nest.Create("spike_detector")
meter_e = nest.Create('multimeter',
        params = {'withtime': True,
                  'interval': 0.1,
                  'record_from': ['V_m', 'g_AMPA', 'g_GABA_A', 'g_NMDA']})
meter_i = nest.Create('multimeter',
        params = {'withtime': True,
                  'interval': 0.1,
                  'record_from': ['V_m', 'g_AMPA', 'g_GABA_A', 'g_NMDA']})

nest.SetStatus([espikes],[{"label": "ramp_model_ex",
                   "withtime": True,
                   "withgid": True}])

nest.SetStatus([ispikes],[{"label": "ramp_model_in",
                   "withtime": True,
                   "withgid": True}])

print "Connecting devices."

 
nest.ConvergentConnect(range(1,N_rec+1),          espikes, model="static_synapse")
nest.ConvergentConnect(range(o.Ne, o.Ne+1+N_rec), ispikes, model="static_synapse")

nest.Connect(meter_e, [1])
nest.Connect(meter_i, [o.Ne + 1])

print "Connecting network."

numpy.random.seed(1234)
net.connRandom()

endbuild=time.time()

print "Simulating."

nest.Simulate(simtime)

endsimulate= time.time()

build_time = endbuild-startbuild
sim_time   = endsimulate-endbuild

#nest.raster_plot.from_device(espikes, hist=True)
#nest.raster_plot.from_device(ispikes, hist=True)

# obtain and display data
events = nest.GetStatus(meter_e)[0]['events']
t = events['times'];

pl.figure()
pl.title('Excitatory cell')
ax = pl.subplot(211)
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')

pl.subplot(212, sharex=ax)
pl.plot(t, events['g_AMPA'], t, events['g_GABA_A'], t, events['g_NMDA'])
pl.xlabel('Time [ms]')
pl.ylabel('Synaptic conductance [nS]')
pl.legend(('g_AMPA', 'g_GABA_A', 'g_NMDA'))

# obtain and display data
events = nest.GetStatus(meter_i)[0]['events']
t = events['times'];

pl.figure()
pl.title('Inhibitory cell')
pl.subplot(211)
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')

pl.subplot(212)
pl.plot(t, events['g_AMPA'], t, events['g_GABA_A'], t, events['g_NMDA'])
pl.xlabel('Time [ms]')
pl.ylabel('Synaptic conductance [nS]')
pl.legend(('g_AMPA', 'g_GABA_A', 'g_NMDA'))

#nest.raster_plot.show()
pl.show()
