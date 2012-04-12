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

from scipy.signal import *

import time
import math
import sys
import numpy as np
import logging as lg

from EI_network import *
from EI_network_sim_mod import *
from custombrian import *

lg.basicConfig(level=lg.DEBUG)


def butterHighPass(sig, dt, f_pass):
    nyq_f = 1./dt/2
    norm_f_pass = f_pass/nyq_f

    # Low pass filter
    b, a = butter(3, norm_f_pass, btype='high')
    return filtfilt(b, a, sig)


parser = getOptParser()

parser.add_option("--theta_freq", type="float",
        help="Theta frequency stimulation (Hz)")

(options, args) = parser.parse_args()
options = setOptionDictionary(parser, options)

# Clock definitions
sim_dt = options.sim_dt*second
simulationClock = Clock(dt=sim_dt)
stimClock = Clock(50*msecond)

theta_omega = 2*np.pi*options.theta_freq/second

# Other
figSize = (12,8)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()

options.ndim = -1
ei_net = EI_Network(options, simulationClock)

ei_net.connRandom(options.AMPA_density, options.GABA_density)


duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

@network_operation(simulationClock)
def thetaStimulation():
    ph = theta_omega*simulationClock.t
    ei_net.E_pop.Iext = options.Iext_e/2*np.sin(ph + np.pi/2) + options.Iext_e/2
    ei_net.I_pop.Iext = options.Iext_i/2*np.sin(ph + np.pi/2) + options.Iext_i/2
    #pass


state_record_e = [10, 20]
state_record_i = [10, 20]

Iext_e_rec = [10]
Iext_i_rec = [10]

spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)
stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=simulationClock)
stateMon_Iext_e = StateMonitor(ei_net.E_pop, 'Iext', record=Iext_e_rec,
        clock=simulationClock)
stateMon_Iext_i = StateMonitor(ei_net.I_pop, 'Iext', record=Iext_i_rec,
        clock=simulationClock)
stateMon_g_ad_e = StateMonitor(ei_net.E_pop, 'g_ad', record=state_record_e,
        clock=simulationClock)
stateMon_g_ad_i = StateMonitor(ei_net.I_pop, 'g_ad', record=state_record_i,
        clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i, stateMon_Iclamp_e,
        stateMon_Iclamp_i)
ei_net.net.add(stateMon_Iext_e)
ei_net.net.add(stateMon_Iext_i)
ei_net.net.add(stateMon_g_ad_e)
ei_net.net.add(stateMon_g_ad_i)

ei_net.net.add(thetaStimulation)


################################################################################
#                              Main cycle
################################################################################
for trial_it in [1]: #range(ei_net.o.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."
    print "Simulation running..."
    start_time=time.time()
    
    ei_net.net.run(options.time*second, report='stdout',
            report_period=options.update_interval*second)
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)
    


    figure()
    ax = subplot(211)
    plot(stateMon_e.times, stateMon_e.values.T/mV)
    ylabel('E membrane potential (mV)')
    subplot(212, sharex=ax)
    plot(stateMon_i.times, stateMon_i.values.T/mV)
    xlabel('Time (s)')
    ylabel('I membrane potential (mV)')
    xlim([1.5, 2.5])
    

    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[0].T/pA)
    ylabel('E synaptic current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[0].T/pA)
    xlabel('Time (s)')
    ylabel('I synaptic current (pA)')
    xlim([1.5, 2.5])
    

    # Band pass filter these signals

    figure()
    hold(True)
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, butterHighPass(stateMon_Iclamp_e.values[0].T/pA +
            stateMon_Iext_e.values[0].T/pA, options.sim_dt, 40))
    plot(stateMon_Iclamp_e.times, stateMon_Iext_e.values[0]/pA)
    ylabel('E current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, butterHighPass(stateMon_Iclamp_i.values[0].T/pA +
            stateMon_Iext_i.values[0].T/pA, options.sim_dt, 40))
    plot(stateMon_Iclamp_i.times, stateMon_Iext_i.values[0]/pA)
    xlabel('Time (s)')
    ylabel('I current (pA)')
    xlim([1.5, 2.5])
    
    # Firing rate
    Favg_e = spikeMon_e.getNSpikes()/options.time
    mean_e = np.mean(Favg_e)
    #Favg_i = spikeMon_i.getNSpikes()/options.time
    figure()
    h = hist(Favg_e, 20)
    xlabel('E f. rate (Hz)')
    ylabel('Count')
    #text((h[1][-1] + h[1][0])/2, np.max(h[0])*1.2, 'Average: ' + str(mean_e) + ' Hz')
    title('Average: ' + str(mean_e) + ' Hz')

    
#    outData = {}
#    outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
#    if spikeMon_e is not None:
#        outData['spikeCell_e'] = spikeMon_e.aspikes
#    if spikeMon_i is not None:
#        outData['spikeCell_i'] = spikeMon_i.aspikes
#    outData['options'] = options._einet_optdict
#
#    outData['bumpPos'] = pos
#    outData['bumpPos_times'] = times
#
#    outData['Fe'] = Fe
#    outData['Fe_t'] = Fe_t
#    outData['Fi'] = Fi
#    outData['Fi_t'] = Fi_t
#    outData['F_dt'] = F_dt
#
#
#
#    savemat(output_fname + '_output.mat', outData, do_compression=True)

    show()

    print "End of trial no. " + str(trial_it) + "..."
    print 

    #ei_net.reinit()
#                            End main cycle
################################################################################
