#
#   simulate_on_gc_nest.py
#
#   Script to test nest functionality
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

from scipy.optimize import fsolve

import nest
import nest.raster_plot

import numpy
from numpy import exp
 
import time

nest.Install('gridcellsmodule')
nest.ResetKernel()

startbuild= time.time()

dt      = 0.2    # the resolution in ms
simtime = 2e3    # Simulation time in ms
delay   = 2.0    # synaptic delay in ms

# Place cell parameters
PC_N    = 1
PC_rate = 1000.0   # Hz


# Initialize the parameters of the integrate and fire neuron
tau_m           = 10.0  # ms
E_L             = -70.0
C_m             = 250.0
t_ref           = 2.0
V_th            = 1000.0
V_reset         = E_L
E_ex            = 0.0
E_in            = -70.0
g_L             = C_m / tau_m
tau_syn_ex      = 3.0
tau_syn_in      = 20.0
I_e             = 0.0

GC_N    = 4000


J_ex  = 0.75/4
J_in  = 0

numThreads = 1
nest.SetKernelStatus({"resolution": dt, "print_time": True})
nest.SetKernelStatus({"local_num_threads": numThreads})

print "Building network"

neuron_params = {"V_m"              : E_L,
                 "E_L"              : E_L,
                 "C_m"              : C_m,
                 "t_ref"            : t_ref,
                 "V_th"             : V_th,
                 "V_reset"          : V_reset,
                 "E_ex"             : E_ex,
                 "E_in"             : E_in,
                 "g_L"              : g_L,
                 "tau_syn_ex"       : tau_syn_ex,
                 "tau_syn_in"       : tau_syn_in,
                 "I_e"              : I_e}

# Place cells
nest.SetDefaults("poisson_generator",{"rate": PC_rate})
PC = nest.Create("poisson_generator", PC_N)

PC_spikes = nest.Create("spike_detector")
nest.SetStatus([PC_spikes],[{"label": "Place cells",
                   "withtime": True,
                   "withgid": True}])


# Grid cells
GC_model = "iaf_cond_exp"
GC = nest.Create(GC_model, GC_N, params = neuron_params)

GC_meter = nest.Create('multimeter', 1, params = {'withtime': True, 'interval': 0.1, 'record_from': ['V_m']})

print "Connecting devices."

#nest.CopyModel("static_synapse","ex_AMPA",{"weight":J_ex, "delay":delay,
#    "receptor_type": receptors["AMPA"]})

nest.DivergentConnect(PC, GC, model="static_synapse", weight=J_ex, delay=delay)
#nest.ConvergentConnect(PC, PC_spikes, model="static_synapse")

nest.Connect(GC_meter, [GC[0]])

endbuild=time.time()


# Simulation
print "Simulating."

nest.Simulate(simtime)

endsimulate= time.time()


build_time = endbuild-startbuild
sim_time   = endsimulate-endbuild

print 'sim_time: ', sim_time

#nest.raster_plot.from_device(PC_spikes, hist=True)

# obtain and display data
events = nest.GetStatus(GC_meter)[0]['events']
t = events['times'];

pl.figure()
ax = pl.subplot(2, 1, 1)
pl.title('Grid cell')
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')

#events = nest.GetStatus(GC_meter)[1]['events']
#t = events['times'];
#pl.subplot(2, 1, 2, sharex=ax, sharey=ax)
#pl.plot(t, events['V_m'])


#nest.raster_plot.show()
pl.show()
