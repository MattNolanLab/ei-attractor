# This version uses numpy to draw the random connections.

import pylab as pl
from scipy.optimize import fsolve

import nest
import nest.raster_plot

import numpy
from numpy import exp
 
import time



nest.ResetKernel()

startbuild= time.time()

dt      = 0.1    # the resolution in ms
simtime = 10e3 # Simulation time in ms
delay   = 2.0    # synaptic delay in ms

# Parameters for asynchronous irregular firing
epsilon = 0.4    # connection probability

order     = 1000
NE        = 4*order
NI        = 1*order
N_neurons = NE+NI
N_rec     = 50 # record from 50 neurons

CE    = epsilon*NE   # number of excitatory synapses per neuron
CI    = epsilon*NI   # number of inhibitory synapses per neuron  
C_tot = int(CI+CE)  # total number of synapses per neuron

# Initialize the parameters of the integrate and fire neuron
tauMem = 10.0
E_L = -60.0
CMem = 250.0
t_ref = 2.0
V_th = -45.0
V_reset = E_L
E_ex = 0.0
E_in = -75.0
g_L = CMem / tauMem
tau_syn_ex = 2.0
tau_syn_in = 5.0

#I_e_e = 500.0
I_e_i = 0.
p_rate = 1000.0

# nS
J_ex  = 0.75/4
J_in  = -0.3/4
J_ext_ex  = 4.0
J_ext_in  = 0.0


numThreads = 8
nest.SetKernelStatus({"resolution": dt, "print_time": False})
nest.SetKernelStatus({"local_num_threads": numThreads})

print "Building network"

neuron_params= {"V_m"       : E_L,
                "E_L"       : E_L,
                "C_m"       : CMem,
                "t_ref"     : t_ref,
                "V_th"      : V_th,
                "V_reset"   : V_reset,
                "E_ex"      : E_ex,
                "E_in"      : E_in,
                "g_L"       : g_L,
                "tau_syn_ex": tau_syn_ex,
                "tau_syn_in": tau_syn_in}

nodes_ex=nest.Create("iaf_cond_exp",NE, params = neuron_params)
nodes_in=nest.Create("iaf_cond_exp",NI, params = neuron_params)

nest.SetDefaults("poisson_generator",{"rate": p_rate})
noise=nest.Create("poisson_generator")


espikes=nest.Create("spike_detector")
ispikes=nest.Create("spike_detector")
meter_e = nest.Create('multimeter', params = {'withtime': True, 'interval': 0.1, 'record_from': ['V_m', 'g_ex', 'g_in']})
meter_i = nest.Create('multimeter', params = {'withtime': True, 'interval': 0.1, 'record_from': ['V_m', 'g_ex', 'g_in']})

nest.SetStatus([espikes],[{"label": "brunel-py-ex",
                   "withtime": True,
                   "withgid": True}])

nest.SetStatus([ispikes],[{"label": "brunel-py-in",
                   "withtime": True,
                   "withgid": True}])

print "Connecting devices."

nest.CopyModel("static_synapse","excitatory",{"weight":J_ex, "delay":delay})
nest.CopyModel("static_synapse","inhibitory",{"weight":J_in, "delay":delay})

nest.DivergentConnect(noise,nodes_ex, weight=J_ext_ex, delay=delay, model="static_synapse")
nest.DivergentConnect(noise,nodes_in, weight=J_ext_in, delay=delay, model="static_synapse")

 
nest.ConvergentConnect(range(1,N_rec+1),espikes,model="excitatory")
nest.ConvergentConnect(range(NE+1,NE+1+N_rec),ispikes,model="excitatory")

nest.Connect(meter_e, [1])
nest.Connect(meter_i, [NE+1])

print "Connecting network."

# Here, we create the connections from the excitatory neurons to all other
# neurons. We exploit that the neurons have consecutive IDs, running from
# 1,...,NE           for the excitatory neurons and from
# (NE+1),...,(NE+NI) for the inhibitory neurons.

numpy.random.seed(1234)

nest.RandomConvergentConnect(nodes_ex, nodes_in, int(CE), model="excitatory")
nest.RandomConvergentConnect(nodes_in, nodes_ex, int(CI), model="inhibitory")

endbuild=time.time()

print "Simulating."

nest.Simulate(simtime)

endsimulate= time.time()

events_ex = nest.GetStatus(espikes,"n_events")[0]
rate_ex   = events_ex/simtime*1000.0/N_rec
events_in = nest.GetStatus(ispikes,"n_events")[0]
rate_in   = events_in/simtime*1000.0/N_rec

num_synapses = nest.GetDefaults("excitatory")["num_connections"]+\
nest.GetDefaults("inhibitory")["num_connections"]

build_time = endbuild-startbuild
sim_time   = endsimulate-endbuild

print "Brunel network simulation (Python)"
print "Number of neurons :", N_neurons
print "Number of synapses:", num_synapses
print "       Exitatory  :", int(CE*N_neurons)+N_neurons
print "       Inhibitory :", int(CI*N_neurons)
print "Excitatory rate   : %.2f Hz" % rate_ex
print "Inhibitory rate   : %.2f Hz" % rate_in
print "Building time     : %.2f s" % build_time
print "Simulation time   : %.2f s" % sim_time
    
nest.raster_plot.from_device(espikes, hist=True)
nest.raster_plot.from_device(ispikes, hist=True)

# obtain and display data
events = nest.GetStatus(meter_e)[0]['events']
t = events['times'];

pl.figure()
pl.title('Excitatory cell')
pl.subplot(211)
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')

pl.subplot(212)
pl.plot(t, events['g_ex'], t, events['g_in'])
pl.xlabel('Time [ms]')
pl.ylabel('Synaptic conductance [nS]')
pl.legend(('g_exc', 'g_inh'))

# obtain and display data
events = nest.GetStatus(meter_i)[0]['events']
t = events['times'];

pl.figure()
pl.title('Inhibitory cell')
pl.subplot(211)
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')

pl.subplot(212)
pl.plot(t, events['g_ex'], t, events['g_in'])
pl.xlabel('Time [ms]')
pl.ylabel('Synaptic conductance [nS]')
pl.legend(('g_exc', 'g_inh'))

nest.raster_plot.show()
