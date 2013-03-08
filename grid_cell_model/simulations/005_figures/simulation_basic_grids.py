#
#   simulation_basic_grids.py
#
#   Main simulation run: grid fields with theta input and all the inhibition
#   (for gamma) and place input.
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

from scipy.io   import savemat
from optparse   import OptionParser

import nest
import nest.raster_plot

from parameters              import *
from models.gc_net_nest      import *
from analysis.spikes         import slidingFiringRateTuple, torusPopulationVector
from analysis.image          import Position2D, fitGaussianBumpTT

import time
import logging as lg


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--ndumps",  type="int",    help="Number of data output dumps during the simulation")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

# Other
figSize = (12,8)


x_lim = [options.time - 2e3, options.time]

################################################################################
#                              Main cycle
################################################################################
for trial_it in range(options.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."

    ################################################################################
    #                              Network setup
    ################################################################################
    print "Starting network and connections initialization..."
    start_time=time.time()
    total_start_t = time.time()

    ei_net = NestGridCellNetwork(options, simulationOpts=None)
    
    const_v = [options.Ivel, 0.0]
    ei_net.setVelocityCurrentInput_e()
    ei_net.setPlaceCells()
    
    

    rec_all_spikes = True
    if rec_all_spikes:
        nrecSpike_e = ei_net.Ne_x*ei_net.Ne_y
        nrecSpike_i = ei_net.Ni_x*ei_net.Ni_y
    else:
        nrecSpike_e = 200
        nrecSpike_i = 50
    
    state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
    state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]
    
    spikeMon_e = ei_net.getSpikeDetector("E")
    spikeMon_i = ei_net.getSpikeDetector("I")
    pc_spikemon = ei_net.getGenericSpikeDetector(ei_net.PC, "Place cells")
    
    #stateMon_params = {
    #        'withtime' : True,
    #        'interval' : 0.1,
    #        'record_from' : ['V_m', 'I_clamp_AMPA', 'I_clamp_NMDA',
    #            'I_clamp_GABA_A', 'I_stim']
    #}
    #stateMon_e = ei_net.getStateMonitor("E", state_record_e, stateMon_params)
    #stateMon_i = ei_net.getStateMonitor("I", state_record_i, stateMon_params)

    
    duration=time.time()-start_time
    print "Network setup time:",duration,"seconds"
    #                            End Network setup
    ################################################################################
    print "Simulation running..."
    start_time=time.time()

    
    for dump_it in range(options.ndumps):
        ei_net.simulate(options.time/options.ndumps, printTime=False)
        duration=time.time()-start_time
        print "Simulation time:",duration,"seconds"
        
        
        output_fname = "{0}/{1}job{2:05}_trial{3:04}_dump{4:03}".format(options.output_dir,
                options.fileNamePrefix, options.job_num, trial_it, dump_it)


        #events_e = nest.GetStatus(stateMon_e)[1]['events']
        #events_i = nest.GetStatus(stateMon_i)[0]['events']

        #F_tstart = 0.0
        #F_tend = options.time
        #F_dt = 20.0
        #F_winLen = 250.0
        senders_e      = nest.GetStatus(spikeMon_e)[0]['events']['senders'] - ei_net.E_pop[0]
        spikeTimes_e   = nest.GetStatus(spikeMon_e)[0]['events']['times']
        #Fe, Fe_t = slidingFiringRateTuple((senders_e, spikeTimes_e), ei_net.net_Ne,
        #        F_tstart, F_tend, F_dt, F_winLen)

        #if (len(ei_net.PC) != 0):
        #    senders_pc     = nest.GetStatus(pc_spikemon)[0]['events']['senders'] - ei_net.PC[0]
        #    spikeTimes_pc  = nest.GetStatus(pc_spikemon)[0]['events']['times']
        #    Fpc, Fpc_t = slidingFiringRateTuple((senders_pc, spikeTimes_pc),
        #            ei_net.N_pc_created, F_tstart, F_tend, F_dt, F_winLen)


        #bumpT = ei_net.no.time - 2*F_winLen
        #bumpI = bumpT / F_dt
        #bump_e = np.reshape(Fe[:, bumpI], (ei_net.Ne_y, ei_net.Ne_x))


        ## Flattened firing rate of E/I cells
        #figure()
        #T, N_id = np.meshgrid(Fe_t, np.arange(ei_net.net_Ne))
        #pcolormesh(T, N_id,  Fe)
        #xlabel("Time (ms)")
        #ylabel("Neuron #")
        #axis('tight')
        #colorbar()
        #title('Firing rate of E cells')
        #savefig(output_fname + '_firing_rate_flat.png')


        ## Flattened firing rate of Place cells
        #if (len(ei_net.PC) != 0):
        #    figure()
        #    T, N_id = np.meshgrid(Fpc_t, np.arange(ei_net.N_pc_created))
        #    xlim([0, 1000])
        #    pcolormesh(T, N_id, Fpc)
        #    xlabel("Time (s)")
        #    ylabel("Neuron #")
        #    axis('tight')
        #    colorbar()
        #    title('Firing rate of Place cells')
        #    savefig(output_fname + '_firing_rate_flat_pc.png')

        #
        #
        ## External currents
        #figure()
        #ax = subplot(211)
        #plot(events_e['times'], events_e['I_stim'])
        #ylabel('E cell $I_{stim}$')
        #axis('tight')
        #subplot(212)
        #plot(events_i['times'], events_i['I_stim'])
        #ylabel('I cell $I_{stim}$')
        #xlabel('Time (ms)')
        #savefig(output_fname + '_Iext.png')
        #
        #
        #
        ## E/I I_syn
        #figure()
        #ax = subplot(211)
        #plot(events_e['times'], events_e['I_clamp_GABA_A'])
        #ylabel('E synaptic current (pA)')
        #subplot(212, sharex=ax)
        #plot(events_i['times'], events_i['I_clamp_AMPA'] + events_i['I_clamp_NMDA'])
        #xlabel('Time (ms)')
        #ylabel('I synaptic current (pA)')
        #xlim(x_lim)
        #savefig(output_fname + '_Isyn.pdf')
        #
        #
        ## Firing rate of E cells on the twisted torus
        #figure()
        #pcolormesh(bump_e)
        #xlabel('E neuron no.')
        #ylabel('E neuron no.')
        #colorbar()
        #axis('equal')
        #title('Firing rates (torus) of E cells')
        #savefig(output_fname + '_firing_snapshot_e.png')
        #
        #
        #
        ## Print a plot of bump position
        #(pos, bumpPos_times) = torusPopulationVector(
        #        (senders_e, spikeTimes_e), [ei_net.Ne_x, ei_net.Ne_y],
        #        tstart = ei_net.no.theta_start_t,
        #        tend   = ei_net.no.time,
        #        dt     = F_dt,
        #        winLen = F_winLen)
        #figure(figsize=figSize)
        #plot(bumpPos_times, pos)
        #xlabel('Time (ms)')
        #ylabel('Bump position (neurons)')
        #legend(['X', 'Y'])
        #ylim([-ei_net.Ne_x/2 -5, ei_net.Ne_y/2 + 5])
        #savefig(output_fname + '_bump_position.pdf')



        
        outData = ei_net.getRatData()
        #outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")

        outData['senders_e']     = senders_e
        outData['spikeTimes_e']  = spikeTimes_e

        outData['options']       = options._einet_optdict
        outData['velocityStart'] = options.theta_start_t

        savemat(output_fname + '_output.mat', outData, do_compression=False)

        #print "Dump after" + str(simulationClock.t)


    print "End of trial no. " + str(trial_it) + "..."
    print 
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time:", total_time, "seconds"

