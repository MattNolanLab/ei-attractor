#
#   simulation_fig_model.py
#
#   Main simulation run: model description figures
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
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser
from brian    import *

from parameters              import getOptParser, setOptionDictionary
from grid_cell_network_brian import BrianGridCellNetwork
from custombrian             import ExtendedSpikeMonitor
from tools                   import butterHighPass, spikePhaseTrialRaster, \
        phaseCWT
from plotting                import phaseFigTemplate, rasterPhasePlot, \
        raster_bin_size

import time
import numpy as np
import logging as lg


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--theta_start_mon_t",  type="float",    help="theta start monitoring time")
(options, args) = parser.parse_args()

################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

options.ndim = 'twisted_torus'
ei_net = BrianGridCellNetwork(options, simulationOpts=None)
ei_net.uniformInhibition()
ei_net.uniformExcitation()
ei_net.setConstantCurrent()
ei_net.setStartCurrent()
ei_net.setThetaCurrentStimulation()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()

nrecSpike_e = ei_net.Ne_x*ei_net.Ne_y
nrecSpike_i = ei_net.Ni_x*ei_net.Ni_y

state_record_e = [ei_net.Ne_x/2 - 1, ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop[0:nrecSpike_e])
spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop[0:nrecSpike_i])

stateMon_e          = RecentStateMonitor(ei_net.E_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_i          = RecentStateMonitor(ei_net.I_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_ge_e       = RecentStateMonitor(ei_net.E_pop, 'ge',     duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iclamp_e   = RecentStateMonitor(ei_net.E_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i   = RecentStateMonitor(ei_net.I_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iext_e     = RecentStateMonitor(ei_net.E_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iext_i     = RecentStateMonitor(ei_net.I_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)


theta_n_it_range = 2
theta_state_record_e = range(state_record_e[1] - theta_n_it_range/2,
        state_record_e[1] + theta_n_it_range/2 + 1)
theta_state_record_i = range(state_record_i[1] - theta_n_it_range/2,
        state_record_i[1] + theta_n_it_range/2 + 1)
theta_spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
theta_spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)
theta_stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = theta_state_record_e, clock=simulationClock)
theta_stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = theta_state_record_e, clock=simulationClock)
theta_stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = theta_state_record_i, clock=simulationClock)
theta_stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = theta_state_record_i, clock=simulationClock)
theta_stateMon_Iext_e   = StateMonitor(ei_net.E_pop, 'Iext',   record = theta_state_record_e, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i)
ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
ei_net.net.add(stateMon_ge_e)
ei_net.net.add(stateMon_Iext_e, stateMon_Iext_i)

ei_net.net.add(theta_spikeMon_e, theta_spikeMon_i, theta_stateMon_Iclamp_e,
        theta_stateMon_Iclamp_i, theta_stateMon_Iext_e, theta_stateMon_e, theta_stateMon_i)


x_lim = [options.time/1e3 - 1, options.time/1e3]

################################################################################
#                              Main cycle
################################################################################
print "Simulation running..."
start_time=time.time()

print "  Network initialisation..."
ei_net.net.run(options.theta_start_mon_t*msecond, report='stdout')

theta_spikeMon_e.reinit()
theta_spikeMon_i.reinit()
theta_stateMon_e.reinit()
theta_stateMon_i.reinit()
theta_stateMon_Iclamp_e.reinit()
theta_stateMon_Iclamp_i.reinit()
theta_stateMon_Iext_e.reinit()

print "  Theta stimulation..."
ei_net.net.run((options.time - options.theta_start_mon_t)*msecond, report='stdout')
duration=time.time()-start_time
print "Simulation time:",duration,"seconds"


output_fname = "{0}/{1}job{2:04}".format(options.output_dir,
        options.fileNamePrefix, options.job_num)


F_tstart = 0
F_tend = options.time*1e-3
F_dt = 0.05
F_winLen = 0.25
Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)


# Plot membrane potentials of a cell in the middle of the sheet and at the
# bottom-center
figure()
ax = subplot(211)
plot(stateMon_e.times, stateMon_e.values[:, 0:2]/mV)
ylabel('E cell $V_m$ (mV)')
subplot(212, sharex=ax)
plot(stateMon_i.times, stateMon_i.values[:, 0:2]/mV)
xlabel('Time (s)')
ylabel('I cell $V_m$ (mV)')
xlim(x_lim)
tight_layout()
savefig(output_fname + '_Vm.pdf')


# Plot post-synaptic currents of the same cells
figure()
ax = subplot(211)
plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[:, 0:2]/pA)
ylabel('E cell $I_{syn}$ (pA)')
#ylim([0, 3000])
subplot(212, sharex=ax)
plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[:, 0:2]/pA)
xlabel('Time (s)')
ylabel('I cell $I_{syn}$ (pA)')
xlim(x_lim)
tight_layout()
savefig(output_fname + '_Isyn.pdf')


# Plot external current input
figure()
ax = subplot(211)
plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[:, 1]/pA)
ylabel('E cell $I_{ext}$ (pA)')
subplot(212, sharex=ax)
plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[:, 0]/pA)
xlabel('Time (s)')
ylabel('I cell $I_{ext}$ (pA)')
xlim(x_lim)
tight_layout()
savefig(output_fname + '_Iext.png')


# Plot a snapshot of the population firing rate on the twisted torus, E
# cells
figure()
pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
xlabel('E neuron #')
ylabel('E neuron #')
colorbar()
axis('equal')
tight_layout()
savefig(output_fname + '_firing_snapshot_e.png')


# Plot a snapshot of the population firing rate on the twisted torus, I
# cells
figure()
pcolormesh(np.reshape(Fi[:, len(Fi_t)/2], (ei_net.Ni_y, ei_net.Ni_x)))
xlabel('I neuron #')
ylabel('I neuron #')
colorbar()
axis('equal')
tight_layout()
savefig(output_fname + '_firing_snapshot_i.png')


###################################################################### 
#                        Wavelet and raster analysiS
###################################################################### 
print "Wavelet analysis..."
wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_e.pdf')
high_pass_freq = 40.
maxF = 200
for ei_it in [0]:
    if ei_it == 0:
        print '  E neurons...'
        wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_e.pdf')
        wavelet_sig_fname = output_fname + '_phase_wavelet_e'
        sig_epochs_fname = output_fname + '_sig_epochs_e'
        tmp_stateMon = theta_stateMon_Iclamp_e
        range_n_it = theta_state_record_e
    else:
        print '  I neurons...'
        wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_i.pdf')
        wavelet_sig_fname = output_fname + '_phase_wavelet_i'
        sig_epochs_fname = output_fname + '_sig_epochs_i'
        tmp_stateMon = theta_stateMon_Iclamp_i
        range_n_it = theta_state_record_i

    for n_it in range(len(range_n_it)):
        neuron_no = range_n_it[n_it]
        print('    Neuron no. ' + str(neuron_no))
        cwt_phases, sig_cwt, freq, sig_ph = \
            phaseCWT(butterHighPass(tmp_stateMon[neuron_no].T/pA,
                options.sim_dt*msecond, high_pass_freq), 1./options.theta_freq, options.sim_dt*1e-3, maxF)

        # Wavelet plot
        f = phaseFigTemplate()
        PH, F = np.meshgrid(cwt_phases, freq)
        pcolormesh(PH, F, sig_cwt, edgecolors='None', cmap=get_cmap('jet'))
        ylabel('F (Hz)')
        ylim([0, maxF])
        savefig(wavelet_sig_fname + '{0}.png'.format(n_it),
                dpi=300)
        close()


        # pcolor of signals over theta epochs
        figure(figsize=(12, 6))
        PH, T = np.meshgrid(cwt_phases, np.arange(1, len(sig_ph)+1))
        pcolormesh(PH, T, sig_ph, edgecolor='None')
        xlabel('Theta phase')
        ylabel('Theta epoch')
        xlim([-np.pi, np.pi])
        ylim([1, len(sig_ph)])
        xticks([-np.pi, 0, np.pi], ('$-\pi$', '',  '$\pi$'), fontsize=25)
        yticks([1, len(sig_ph)])
        savefig(sig_epochs_fname + '{0}.png'.format(n_it),
                dpi=300)

        # Average signal plot
        f = phaseFigTemplate()
        mn = np.mean(sig_ph, 0)
        st = np.std(sig_ph, 0)
        gca().fill_between(cwt_phases, mn+st, mn-st, facecolor='black', alpha=0.1, zorder=0)
        plot(cwt_phases, mn, 'k')
        ylabel('I (pA)')
        wavelet_sig_pp.savefig()
        close()

    wavelet_sig_pp.close()
print "Done"


# Raster plots
for ei_it in [0]:
    if ei_it == 0:
        raster_pp = PdfPages(output_fname + '_phase_raster_e.pdf')
        range_n_it = theta_state_record_e
        tmp_spikeMon = theta_spikeMon_e
    else:
        raster_pp = PdfPages(output_fname + '_phase_raster_i.pdf')
        range_n_it = theta_state_record_i
        tmp_spikeMon = theta_spikeMon_i

    for n_it in range(len(range_n_it)):
        neuron_no = range_n_it[n_it]
        print('Saving rasters for neuron no. ' + str(neuron_no))
        phases, times, trials = spikePhaseTrialRaster(tmp_spikeMon[neuron_no],
                options.theta_freq, options.theta_start_mon_t*msecond)
        phases -= np.pi

        # Raster plots (single cell over 'theta' epochs)
        ntrials = np.ceil((options.time - options.theta_start_mon_t) * msecond * options.theta_freq)
        f = rasterPhasePlot(phases, trials, ntrials)
        raster_pp.savefig()
        close()
    
    raster_pp.close()

#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

