import brian_no_units
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

from tools import *

lg.basicConfig(level=lg.DEBUG)


parser = getOptParser()

parser.add_option("--Ivel", type="float", help="Velocity input (pA)")
parser.add_option("--pAMPA_sigma", type="float", help="AMPA profile spread (normalised)")
parser.add_option("--Iext_e_min", type=float,
        help="Minimal external current onto E cells (theta stim.) (A)")
parser.add_option("--Iext_i_min", type=float,
        help="Minimal external current onto I cells (theta stim.) (I)")
parser.add_option("--g_extraGABA_total", type=float,
        help="Uniform inhibition (E-->I only) total conductance (S)")
parser.add_option("--extraGABA_density", type=float,
        help="Uniform inhibition (E-->I only) connection density")

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
total_start_t = time.time()

options.ndim = 2
ei_net = EI_Network(options, simulationClock)

# Mexican hat properties and AMPA/GABA connections
pAMPA_size= 1.0
pGABA_sigma = 0.5/6
ei_net.connMexicanHat(pAMPA_size, pGABA_sigma)
ei_net.randomInhibition(options.g_extraGABA_total, options.extraGABA_density)

print('pAMPA_sigma = ' + str(options.pAMPA_sigma) + '/6')


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

stim_freq = 8*Hz
stim_omega = 2*np.pi*stim_freq
stim_e_A  = (options.Iext_e - options.Iext_e_min)/2*amp
stim_e_DC = (options.Iext_e + options.Iext_e_min)/2*amp
stim_i_A  = (options.Iext_i - options.Iext_i_min)/2*amp
stim_i_DC = (options.Iext_i + options.Iext_i_min)/2*amp


stim_start = int(0.45*ei_net.o.Ne)
stim_range = int(0.2*ei_net.o.Ne)
stim_current = 900*pA
#stim_current = options.Iext_e

@network_operation(stimClock)
def stimulateSubPopulation():
    if simulationClock.t >= 0*msecond and simulationClock.t < 100*msecond:
        #ei_net.E_pop.Iext = 0
        tmp = ei_net.E_pop.Iext.reshape((options.Ne, options.Ne))
        tmp[stim_start:stim_start+stim_range, stim_start:stim_start+stim_range] =\
            linspace(stim_current, stim_current, stim_range**2).reshape((stim_range, stim_range))
        ei_net.E_pop.Iext = tmp.ravel()
        print "Stimulation..."
    #elif simulationClock.t >= 0.5*second and simulationClock.t < options.time*second:
    #    v = np.array([[-1, -1]]).T
    #    Ivel = np.dot(ei_net.prefDirs, v) * options.Ivel*pA
    #    ei_net.E_pop.Iext = ei_net.o.Iext_e + Ivel.T
    else:
        ei_net.E_pop.Iext = [ei_net.E_pop.Iext[0]] * len(ei_net.E_pop)
    pass

@network_operation(simulationClock)
def thetaStimulation():
    if simulationClock.t >= 0.5*second and simulationClock.t < options.time*second:
        ph = stim_omega*simulationClock.t
        ei_net.E_pop.Iext = stim_e_DC + stim_e_A*np.sin(ph - np.pi/2)
        ei_net.I_pop.Iext = stim_i_DC + stim_i_A*np.sin(ph - np.pi/2)


state_record_e = [31, 2015]
state_record_i = [15, 527]

spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)

stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=simulationClock)
stateMon_Iext_e = StateMonitor(ei_net.E_pop, 'Iext', record = state_record_e, clock=simulationClock)
stateMon_Iext_i = StateMonitor(ei_net.I_pop, 'Iext', record = state_record_i, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i)
ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
ei_net.net.add(stateMon_Iext_e, stateMon_Iext_i)
ei_net.net.add(stimulateSubPopulation)
ei_net.net.add(thetaStimulation)


## Export connectivity matrices
#print "Exporting connections..."
#connOut = {}
#connOut['AMPA_conn'] = ei_net.AMPA_conn.W
#connOut['GABA_conn'] = ei_net.GABA_conn1.W
#connOut['options'] = options._einet_optdict
#
#conn_fname = "{0}/{1}job{2:04}_connections.mat".format(options.output_dir,
#        options.fileNamePrefix, options.job_num)
#savemat(conn_fname, connOut, do_compression=False)
#print "Finished exporting connections"


x_lim = [options.time-1, options.time]
#x_lim = [0, options.time]

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
    Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
    Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)

    # plot firing rates
    figure(figsize=figSize)
    subplot(211)
    T, FR = np.meshgrid(Fe_t, np.arange(ei_net.net_Ne))
    pcolormesh(T, FR, Fe)
    ylabel('E Neuron no.')
    colorbar()
    subplot(212)
    T, FR = np.meshgrid(Fi_t, np.arange(ei_net.net_Ni))
    pcolormesh(T, FR, Fi)
    xlabel('Time (s)')
    ylabel('I Neuron no.')
    colorbar()
    savefig(output_fname + '_firing_rate_e.png')

    figure()
    ax = subplot(211)
    plot(stateMon_e.times, stateMon_e.values[0:2].T/mV)
    ylabel('E membrane potential (mV)')
    subplot(212, sharex=ax)
    plot(stateMon_i.times, stateMon_i.values[0:2].T/mV)
    xlabel('Time (s)')
    ylabel('I membrane potential (mV)')
    xlim(x_lim)
    savefig(output_fname + '_Vm.pdf')
    
    
    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[0:2].T/pA)
    ylabel('E synaptic current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[0:2].T/pA)
    xlabel('Time (s)')
    ylabel('I synaptic current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Isyn.pdf')
    
    figure()
    ax = subplot(211)
    plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[0].T/pA)
    ylabel('E external current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[0].T/pA)
    xlabel('Time (s)')
    ylabel('I external current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Iext.pdf')
    
    # High pass filter these signals
    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, butterHighPass(stateMon_Iclamp_e.values[0].T/pA, options.sim_dt, 40))
    plot(stateMon_Iext_e.times, stateMon_Iext_e.values[0]/pA - stim_e_DC/pA)
    ylabel('E current (pA)')
    ylim([-500, 500])
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, butterHighPass(stateMon_Iclamp_i.values[0].T/pA, options.sim_dt, 40))
    #plot(stateMon_Iclamp_i.times, stateMon_Iext_i.values[0]/pA)
    xlabel('Time (s)')
    ylabel('I current (pA)')
    xlim(x_lim)
    ylim([-500, 500])
    savefig(output_fname + '_Isyn_filt.pdf')
    
    
    
    Ne = options.Ne
    figure()
    pcolormesh(np.reshape(ei_net.AMPA_conn.W.todense()[15, :], (options.Ni,
        options.Ni)));
    xlabel('I neuron no.')
    ylabel('I neuron no.')
    colorbar()
    savefig(output_fname + '_E2I_conn.png')

    Ni = options.Ni
    figure()
    pcolormesh(np.reshape(ei_net.GABA_conn1.W.todense()[15, :], (options.Ne,
        options.Ne)));
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    savefig(output_fname + '_I2E_conn.png')

    figure()
    Ne = options.Ne
    pcolormesh(np.reshape(np.dot(ei_net.AMPA_conn.W.todense(),
        ei_net.GABA_conn1.W.todense())[15, :], (Ne, Ne)));
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    savefig(output_fname + '_E2E_conn.png')


    figure()
    pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (options.Ne, options.Ne)))
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    savefig(output_fname + '_firing_snapshot_e.png')


    
    # Print a plot of bump position
    (pos, times) = spikeMon_e.torusPopulationVector(ei_net.o.Ne, F_tstart, F_tend, F_dt,
            F_winLen)
    figure(figsize=figSize)
    plot(times, pos)
    xlabel('Time (s)')
    ylabel('Bump position (neurons)')
    
    savefig(output_fname + '_bump_position.pdf')
    
    
    #outData = {}
    #outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    #if spikeMon_e is not None:
    #    outData['spikeCell_e'] = spikeMon_e.aspikes
    #if spikeMon_i is not None:
    #    outData['spikeCell_i'] = spikeMon_i.aspikes
    #outData['options'] = options._einet_optdict

    #outData['bumpPos'] = pos
    #outData['bumpPos_times'] = times

    #outData['Fe'] = Fe
    #outData['Fe_t'] = Fe_t
    #outData['Fi'] = Fi
    #outData['Fi_t'] = Fi_t


    #savemat(output_fname + '_output.mat', outData, do_compression=True)


    print "End of trial no. " + str(trial_it) + "..."
    print 

    #ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

