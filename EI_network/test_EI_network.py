#import brian_no_units
from brian import *

from brian import *
from brian.library.IF import *
from brian.library.synapses import *

from scipy import linspace
from scipy.io import loadmat
from scipy.io import savemat
from optparse import OptionParser
from datetime import datetime

import time
import math
import sys
import numpy as np
import logging as lg

from EI_network import *
from EI_network_sim_mod import *
from custombrian import *

lg.basicConfig(level=lg.DEBUG)


parser = getOptParser()
(options, args) = parser.parse_args()
options = setOptionDictionary(parser, options)

# Clock definitions
sim_dt = options.sim_dt*second
simulationClock = Clock(dt=sim_dt)
stimClock = Clock(50*msecond)

# Other
figSize = (12,8)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()

options.ndim = 2
ei_net = EI_Network(options, simulationClock)

# Mexican hat properties and AMPA/GABA connections
pAMPA_mu = 0.5
pAMPA_sigma = 0.25/6
pGABA_sigma = 0.5/6
ei_net.connMexicanHat(pAMPA_mu, pAMPA_sigma, pGABA_sigma)


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

stim_start = int(0.4*ei_net.o.Ne)
stim_range = int(0.2*ei_net.o.Ne)
stim_current = 1200*pA

@network_operation(stimClock)
def stimulateSubPopulation():
    if simulationClock.t > 500*msecond and simulationClock.t < 650*msecond:
        tmp = ei_net.E_pop.Iext.reshape((options.Ne, options.Ne))
        tmp[stim_start:stim_start+stim_range, stim_start:stim_start+stim_range] =\
            linspace(stim_current, stim_current, stim_range**2).reshape((stim_range, stim_range))
        ei_net.E_pop.Iext = tmp.ravel()
        print "Stimulation..."
    else:
        ei_net.E_pop.Iext = [ei_net.E_pop.Iext[0]] * len(ei_net.E_pop)


state_record_e = [465]
state_record_i = [465]

spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)
stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
stateMon_Isyn_e = StateMonitor(ei_net.E_pop, 'Isyn', record = state_record_e, clock=simulationClock)
stateMon_Isyn_i = StateMonitor(ei_net.I_pop, 'Isyn', record = state_record_i, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i, stateMon_Isyn_e,
        stateMon_Isyn_i)
ei_net.net.add(stimulateSubPopulation)


################################################################################
#                              Main cycle
################################################################################
for trial_it in range(ei_net.o.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."
    print "Simulation running..."
    start_time=time.time()
    
    ei_net.net.run(options.time*second, report='stdout',
            report_period=options.update_interval*second)
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)
    
    
    F_tstart = 0
    F_tend = options.time
    F_dt = 0.2
    F_winLen = 1.
    f = figure(figsize=figSize)
    title('Network size: ' + str(ei_net.o.Ne) + '(E), ' + str(ei_net.o.Ni) + '(I)')
    subplot2grid((5,1), (0, 0), rowspan=3)
    Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
    X, Y = meshgrid(Fe_t, np.arange(len(Fe)))
    pcolor(X, Y, Fe)
    colorbar()
    ylabel('Neuron no. (E)')
    xlabel('')
    subplot2grid((5,1), (3, 0), rowspan=2)
    Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)
    X, Y = meshgrid(Fi_t, np.arange(len(Fi)))
    pcolor(X, Y, Fi)
    colorbar()
    ylabel('Neuron no. (I)')
    xlabel('Time (s)')
    
    savefig(output_fname + '_firing_rates.png')
    
    
    # Print a plot of bump position
    (pos, times) = spikeMon_e.torusPopulationVector(ei_net.o.Ne, F_tstart, F_tend, F_dt,
            F_winLen)
    figure(figsize=figSize)
    plot(times, pos)
    xlabel('Time (s)')
    ylabel('Bump position (neurons)')
    
    savefig(output_fname + '_bump_position.pdf')
    
    
    outData = {}
    outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    if spikeMon_e is not None:
        outData['spikeCell_e'] = spikeMon_e.aspikes
    if spikeMon_i is not None:
        outData['spikeCell_i'] = spikeMon_i.aspikes
    outData['options'] = options._einet_optdict

    outData['bumpPos'] = pos
    outData['bumpPos_times'] = times

    outData['Fe'] = Fe
    outData['Fe_t'] = Fe_t
    outData['Fi'] = Fi
    outData['Fi_t'] = Fi_t
    outData['F_dt'] = F_dt



    savemat(output_fname + '_output.mat', outData, do_compression=True)


    print "End of trial no. " + str(trial_it) + "..."
    print 

    ei_net.reinit()
#                            End main cycle
################################################################################
