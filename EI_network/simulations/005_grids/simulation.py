import brian_no_units
from brian import *

from brian import *
from brian.library.IF import *
from brian.library.synapses import *

from matplotlib.backends.backend_pdf import PdfPages

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
from plotting import *

lg.basicConfig(level=lg.DEBUG)


parser = getOptParser()

#parser.add_option("--Ivel", type="float", help="Velocity input (pA)")
parser.add_option("--pAMPA_sigma", type="float", help="AMPA profile spread (normalised)")
parser.add_option("--Iext_e_min", type=float,
        help="Minimal external current onto E cells (theta stim.) (A)")
parser.add_option("--Iext_i_min", type=float,
        help="Minimal external current onto I cells (theta stim.) (I)")
parser.add_option("--g_extraGABA_total", type=float,
        help="Uniform inhibition (E-->I only) total conductance (S)")
parser.add_option("--extraGABA_density", type=float,
        help="Uniform inhibition (E-->I only) connection density")
parser.add_option("--prefDirC", type=float,
        help="Preferred directtion multiplier")
parser.add_option("--Ivel_max", type=float,
        help="Rat velocity current will be between 0 and Ivel_max (pA)")

(options, args) = parser.parse_args()
options = setOptionDictionary(parser, options)


vel_fname = '../../../../data/hafting_et_al_2005/rat_trajectory_lowpass.mat'
#vel_fname = '../../../../data/hafting_et_al_2005/rat_data_original.mat'
ratData = loadmat(vel_fname)
rat_dt = ratData['dt'][0][0]
rat_vel_x = np.diff(ratData['pos_x'].ravel())/rat_dt
rat_vel_y = np.diff(ratData['pos_y'].ravel())/rat_dt

# Map velocities so that maximum is Ivel_max
rat_Ivel_x = rat_vel_x * options.Ivel_max/np.max(np.abs(rat_vel_x)) * amp
rat_Ivel_y = rat_vel_y * options.Ivel_max/np.max(np.abs(rat_vel_y)) * amp

print "mean rat_Ivel_x: " + str(np.mean(np.abs(rat_Ivel_x))/pA) + " pA"
print "mean rat_Ivel_y: " + str(np.mean(np.abs(rat_Ivel_y))/pA) + " pA"
print "max rat_Ivel_x: " + str(np.max(np.abs(rat_Ivel_x))/pA) + " pA"
print "max rat_Ivel_y: " + str(np.max(np.abs(rat_Ivel_y))/pA) + " pA"


# Clock definitions
sim_dt = options.sim_dt*second
simulationClock = Clock(dt=sim_dt)
stimClock = Clock(50*msecond)
ratClock = Clock(rat_dt*second)

# Other
figSize = (12,8)



################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

options.ndim = 'twisted_torus'
ei_net = EI_Network(options, simulationClock)

# Mexican hat properties and AMPA/GABA connections
y_dim = np.sqrt(3)/2.
pAMPA_mu = y_dim/2.
pAMPA_sigma = 0.5/6
pGABA_sigma = 0.5/6
ei_net.connMexicanHat(pAMPA_mu, pAMPA_sigma, pGABA_sigma)
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

theta_start_t = 0.5*second
theta_start_mon_t = 1.0*second


stim_start = 0.45
stim_start_x = int(stim_start*ei_net.Ne_x)
stim_start_y = int(stim_start*ei_net.Ne_y)
stim_range = 0.2
stim_range_x = int(stim_range*ei_net.Ne_x)
stim_range_y = int(stim_range*ei_net.Ne_y)
stim_current = 900*pA
#stim_current = options.Iext_e

@network_operation(stimClock)
def stimulateSubPopulation():
    if simulationClock.t >= 0*msecond and simulationClock.t < 100*msecond:
        #ei_net.E_pop.Iext = 0
        tmp = ei_net.E_pop.Iext.reshape((ei_net.Ne_y, ei_net.Ne_x))
        tmp[stim_start_y:stim_start_y+stim_range_y, stim_start_x:stim_start_x+stim_range_x] = linspace(stim_current, stim_current, stim_range_x*stim_range_y).reshape((stim_range_y,
                        stim_range_x))
        ei_net.E_pop.Iext = tmp.ravel()
        print "Stimulation..."
    elif simulationClock.t < theta_start_t:
        ei_net.E_pop.Iext = [ei_net.E_pop.Iext[0]] * len(ei_net.E_pop)
    else:
        pass


v = np.array([[0, 0]]).T
Ivel = np.dot(ei_net.prefDirs_e, v)
Ivel_it = 0
@network_operation(ratClock)
def velocityChange():
    global Ivel_it
    global Ivel
    if simulationClock.t >= theta_start_mon_t and simulationClock.t < options.time*second:
        v = np.array([[rat_Ivel_x[Ivel_it], rat_Ivel_y[Ivel_it]]]).T
        Ivel = np.dot(ei_net.prefDirs_e, v)
        Ivel_it += 1

@network_operation(simulationClock)
def thetaStimulation():
    if simulationClock.t >= theta_start_t and simulationClock.t < options.time*second:
        ph = stim_omega*simulationClock.t
        ei_net.E_pop.Iext = stim_e_DC + stim_e_A*np.sin(ph - np.pi/2)
        ei_net.I_pop.Iext = stim_i_DC + stim_i_A*np.sin(ph - np.pi/2)

    # Velocity inputs
    if simulationClock.t >= theta_start_mon_t and simulationClock.t < options.time*second:
        ei_net.E_pop.Iext += Ivel.ravel()


state_record_e = [31, 2015]
state_record_i = [15, 527]

theta_n_it_range = 1
theta_state_record_e = range(state_record_e[1] - theta_n_it_range/2,
        state_record_e[1] + theta_n_it_range/2 + 1)
theta_state_record_i = range(state_record_i[0] - theta_n_it_range/2,
        state_record_i[0] + theta_n_it_range/2 + 1)
theta_spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
theta_stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = theta_state_record_e, clock=simulationClock)
theta_stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = theta_state_record_i, clock=simulationClock)


ei_net.net.add(stimulateSubPopulation)
ei_net.net.add(thetaStimulation)
ei_net.net.add(velocityChange)

ei_net.net.add(theta_spikeMon_e, theta_stateMon_Iclamp_e,
        theta_stateMon_Iclamp_i)



#x_lim = [options.time-0.5, options.time]
x_lim = [0, options.time]

################################################################################
#                              Main cycle
################################################################################
for trial_it in range(ei_net.o.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."
    print "Simulation running..."
    start_time=time.time()
    
    print "  Network initialisation..."
    ei_net.net.run(theta_start_mon_t, report='stdout',
            report_period=options.update_interval*second)

    theta_spikeMon_e.reinit()
    theta_stateMon_Iclamp_e.reinit()
    theta_stateMon_Iclamp_i.reinit()

    print "  Theta stimulation..."
    ei_net.net.run(options.time*second, report='stdout',
            report_period=options.update_interval*second)
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)
    
    
    #F_tstart = 0
    #F_tend = options.time
    #F_dt = 0.02
    #F_winLen = 0.5
    #Fe, Fe_t = theta_spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 

    ## plot firing rates
    #figure(figsize=figSize)
    #T, FR = np.meshgrid(Fe_t, np.arange(ei_net.net_Ne))
    #pcolormesh(T, FR, Fe)
    #ylabel('E Neuron no.')
    #xlabel('Time (s)')
    #colorbar()
    #savefig(output_fname + '_firing_rate_e.png')

    #figure()
    #pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
    #xlabel('E neuron no.')
    #ylabel('E neuron no.')
    #colorbar()
    #axis('equal')
    #savefig(output_fname + '_firing_snapshot_e.png')

    
    ## Print a plot of bump position
    #(pos, bumpPos_times) = theta_spikeMon_e.torusPopulationVector(ei_net.o.Ne,
    #        theta_start_t, options.time, F_dt, F_winLen)
    #figure(figsize=figSize)
    #plot(bumpPos_times, pos)
    #xlabel('Time (s)')
    #ylabel('Bump position (neurons)')
    #
    #savefig(output_fname + '_bump_position.pdf')

    
    outData = ratData;
    #outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")

    #outData['bumpPos'] = pos
    #outData['bumpPos_times'] = bumpPos_times

    #outData['Fe'] = Fe
    #outData['Fe_t'] = Fe_t

    outData['theta_spikeCell_e'] = theta_spikeMon_e.aspikes[0:200]
    outData['options'] = options._einet_optdict

    #outData['theta_stateMon_Iclamp_e_times'] = theta_stateMon_Iclamp_e.times
    #outData['theta_stateMon_Iclamp_e_values'] = theta_stateMon_Iclamp_e.values
    #outData['theta_stateMon_Iclamp_i_times'] = theta_stateMon_Iclamp_i.times
    #outData['theta_stateMon_Iclamp_i_values'] = theta_stateMon_Iclamp_i.values
    
    savemat(output_fname + '_output.mat', outData, do_compression=False)


    print "End of trial no. " + str(trial_it) + "..."
    print 

    #ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

