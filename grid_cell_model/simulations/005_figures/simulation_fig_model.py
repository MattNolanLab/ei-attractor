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

from matplotlib.backends.backend_pdf import PdfPages

from scipy      import linspace
from scipy.io   import loadmat
from scipy.io   import savemat
from optparse   import OptionParser

from parameters              import *
from grid_cell_network_brian import *
from custombrian             import *
from tools                   import *
from plotting                import *

import time
import math
import sys
import numpy as np
import logging as lg


lg.basicConfig(level=lg.DEBUG)

parser          = getOptParser()
parser.add_option("--theta_start_mon_t",  type="float",    help="theta start monitoring time")

(options, args) = parser.parse_args()
options         = setOptionDictionary(parser, options)

# Other
figSize = (12,8)


################################################################################
#                              Network setup
################################################################################
print "Starting network and connections initialization..."
start_time=time.time()
total_start_t = time.time()

options.ndim = 'twisted_torus'
ei_net = BrianGridCellNetwork(options, simulationOpts=None)
ei_net.uniformInhibition()
ei_net.setConstantCurrent()
ei_net.setStartCurrent()
ei_net.setThetaCurrentStimulation()

duration=time.time()-start_time
print "Network setup time:",duration,"seconds"
#                            End Network setup
################################################################################

simulationClock = ei_net._getSimulationClock()

rec_all_spikes = True
if rec_all_spikes:
    nrecSpike_e = ei_net.Ne_x*ei_net.Ne_y
    nrecSpike_i = ei_net.Ni_x*ei_net.Ni_y
else:
    nrecSpike_e = 200
    nrecSpike_i = 50

state_record_e = [ei_net.Ne_x/2 - 1, ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

spikeMon_e          = ExtendedSpikeMonitor(ei_net.E_pop[0:nrecSpike_e])
spikeMon_i          = ExtendedSpikeMonitor(ei_net.I_pop[0:nrecSpike_i])

stateMon_e          = RecentStateMonitor(ei_net.E_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_i          = RecentStateMonitor(ei_net.I_pop, 'vm',     duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iclamp_e   = RecentStateMonitor(ei_net.E_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iclamp_i   = RecentStateMonitor(ei_net.I_pop, 'Iclamp', duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)
stateMon_Iext_e     = RecentStateMonitor(ei_net.E_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_e, clock=simulationClock)
stateMon_Iext_i     = RecentStateMonitor(ei_net.I_pop, 'Iext',   duration=options.stateMonDur*ms,   record = state_record_i, clock=simulationClock)


theta_n_it_range = 2
theta_state_record_e = range(state_record_e[1] - theta_n_it_range/2,
        state_record_e[1] + theta_n_it_range/2 + 1)
theta_state_record_i = range(state_record_i[0] - theta_n_it_range/2,
        state_record_i[0] + theta_n_it_range/2 + 1)
theta_spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
theta_spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)
theta_stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = theta_state_record_e, clock=simulationClock)
theta_stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = theta_state_record_i, clock=simulationClock)

ei_net.net.add(spikeMon_e, spikeMon_i)
ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
ei_net.net.add(stateMon_Iext_e, stateMon_Iext_i)

ei_net.net.add(theta_spikeMon_e, theta_spikeMon_i, theta_stateMon_Iclamp_e,
        theta_stateMon_Iclamp_i)


#x_lim = [options.time-0.5, options.time]
x_lim = [options.time/1e3 - 1, options.time/1e3]

################################################################################
#                              Main cycle
################################################################################
for trial_it in range(ei_net.no.ntrials):
    print "Starting trial no. " + str(trial_it) + "..."
    print "Simulation running..."
    start_time=time.time()


    print "  Network initialisation..."
    ei_net.net.run(options.theta_start_mon_t*msecond, report='stdout')

    theta_spikeMon_e.reinit()
    theta_spikeMon_i.reinit()
    theta_stateMon_Iclamp_e.reinit()
    theta_stateMon_Iclamp_i.reinit()

    print "  Theta stimulation..."
    ei_net.net.run((options.time - options.theta_start_mon_t)*msecond, report='stdout')
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)


    F_tstart = 0
    F_tend = options.time*1e-3
    F_dt = 0.02
    F_winLen = 0.25
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
    plot(stateMon_e.times, stateMon_e.values[:, 0:2]/mV)
    ylabel('E membrane potential (mV)')
    subplot(212, sharex=ax)
    plot(stateMon_i.times, stateMon_i.values[:, 0:2]/mV)
    xlabel('Time (s)')
    ylabel('I membrane potential (mV)')
    xlim(x_lim)
    savefig(output_fname + '_Vm.pdf')
    
    
    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[:, 0:2]/pA)
    ylabel('E synaptic current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[:, 0:2]/pA)
    xlabel('Time (s)')
    ylabel('I synaptic current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Isyn.pdf')
    
    figure()
    ax = subplot(211)
    plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[:, 1]/pA)
    ylabel('E external current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[:, 0]/pA)
    xlabel('Time (s)')
    ylabel('I external current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Iext.pdf')
    
    figure()
    pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_firing_snapshot_e.png')


    # Print a plot of bump position
    F_dt = 0.02
    F_winLen = 0.25
    (pos, bumpPos_times) = spikeMon_e.torusPopulationVector([ei_net.Ne_x,
        ei_net.Ne_y], options.theta_start_t*1e-3, options.time*1e-3, F_dt, F_winLen)
    figure(figsize=figSize)
    plot(bumpPos_times, pos)
    xlabel('Time (s)')
    ylabel('Bump position (neurons)')
    ylim([-ei_net.Ne_x/2 -5, ei_net.Ne_y/2 + 5])
    
    savefig(output_fname + '_bump_position.pdf')


    ###################################################################### 
    #                        Wavelet and raster analysiS
    ###################################################################### 
    print "Wavelet analysis..."
    wavelet_list_e = []
    wavelet_list_i = []
    sig_phase_list_e = []
    sig_phase_list_i = []
    wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_e.pdf')
    high_pass_freq = 40.
    maxF = 200
    for ei_it in [0, 1]:
        if ei_it == 0:
            print '  E neurons...'
            wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_e.pdf')
            wavelet_sig_fname = output_fname + '_phase_wavelet_e'
            sig_epochs_fname = output_fname + '_sig_epochs_e'
            max_fname = output_fname + '_Fmax_scatter_e.pdf'
            tmp_stateMon = theta_stateMon_Iclamp_e
            #w_Fmax_vec = w_Fmax_e_vec
            #w_phmax_vec = w_phmax_e_vec
            range_n_it = theta_state_record_e
            wavelet_list = wavelet_list_e
            sig_phase_list = sig_phase_list_e
        else:
            print '  I neurons...'
            wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_i.pdf')
            wavelet_sig_fname = output_fname + '_phase_wavelet_i'
            sig_epochs_fname = output_fname + '_sig_epochs_i'
            max_fname = output_fname + '_Fmax_scatter_i.pdf'
            tmp_stateMon = theta_stateMon_Iclamp_i
            #w_Fmax_vec = w_Fmax_i_vec
            #w_phmax_vec = w_phmax_i_vec
            range_n_it = theta_state_record_i
            wavelet_list = wavelet_list_i
            sig_phase_list = sig_phase_list_i
    
        w_Fmax = np.ndarray(len(range_n_it))
        w_phmax= np.ndarray(len(range_n_it))
        for n_it in range(len(range_n_it)):
            neuron_no = range_n_it[n_it]
            print('    Neuron no. ' + str(neuron_no))
            cwt_phases, sig_cwt, freq, sig_ph = \
                phaseCWT(butterHighPass(tmp_stateMon[neuron_no].T/pA,
                    options.sim_dt*msecond, high_pass_freq), 1./options.theta_freq, options.sim_dt*1e-3, maxF)

            expt_freq = np.linspace(freq[0], freq[-1], len(freq)+1)
            expt_phases = np.linspace(cwt_phases[0], cwt_phases[-1],
                    len(cwt_phases)+1)
    
            w_max = sig_cwt.argmax()
            w_Fmax[n_it] = freq[w_max//len(sig_cwt[0])]
            w_phmax[n_it]= cwt_phases[np.mod(w_max, len(sig_cwt[0]))]

            sig_phase_list.append((sig_ph, cwt_phases))
                    
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

            # Append everything to lists for export
            wavelet_list.append(sig_cwt)

        wavelet_sig_pp.close()
    
        ## Insert this to avg and std vectors and make a plot
        #w_Fmax_vec.append(w_Fmax)
        #w_phmax_vec.append(w_phmax)
        #f = phaseFigTemplate()
        #plot(w_phmax, w_Fmax, 'ko', markersize=10, alpha=0.25)
        #errorbar(np.mean(w_phmax), np.mean(w_Fmax), np.std(w_Fmax),
        #        np.std(w_phmax), 'ko', markersize=10)
        #ylim([50, 120])
        #savefig(max_fname)
    print "Done"


    # Raster plots
    raster_list_e = []
    raster_list_i = []
    hist_list_e = []
    hist_list_i = []
    hists_mean_e = []
    hists_mean_i = []
    hists_std_e = []
    hists_std_i = []
    for ei_it in [0, 1]:
        if ei_it == 0:
            raster_pp = PdfPages(output_fname + '_phase_raster_e.pdf')
            avg_fname = output_fname + '_phase_pspike_avg_e.pdf'
            range_n_it = theta_state_record_e
            tmp_spikeMon = theta_spikeMon_e
            middle_it = options.Ne**2 / 2 + options.Ne/2
            raster_list = raster_list_e
            hists = hist_list_e
            hists_mean = hists_mean_e
            hists_std = hists_std_e
        else:
            raster_pp = PdfPages(output_fname + '_phase_raster_i.pdf')
            avg_fname = output_fname + '_phase_pspike_avg_i.pdf'
            range_n_it = theta_state_record_i
            tmp_spikeMon = theta_spikeMon_i
            middle_it = options.Ni**2 / 2 + options.Ni/2
            raster_list = raster_list_i
            hists = hist_list_i
            hists_mean = hists_mean_i
            hists_std = hists_std_i

        hist_nbins = 1./options.theta_freq/raster_bin_size
       
        for n_it in range(len(range_n_it)):
            neuron_no = range_n_it[n_it]
            print('Saving rasters for neuron no. ' + str(neuron_no))
            phases, times, trials = spikePhaseTrialRaster(tmp_spikeMon[neuron_no],
                    options.theta_freq, options.theta_start_mon_t*msecond)

            phases -= np.pi

            raster_list.append((phases, trials+1))
        
            # Raster plots (single cell over 'theta' epochs)
            ntrials = np.ceil((options.time - options.theta_start_mon_t) * msecond * options.theta_freq)
            f = rasterPhasePlot(phases, trials, ntrials)
            raster_pp.savefig()
            close()
        
            ## Histograms
            #figure(figsize=small_plot_figsize)
            if (len(phases) != 0):
                h = hist(phases, hist_nbins, [-np.pi, np.pi],
                        normed=False)
                hists.append(h[0]/double(ntrials))
                hist_ph = h[1][0:len(h[1])-1]
            else:
                hists.append([])
            delaxes(gca())
        raster_pp.close()

        # Average histogram + one neuron
        f = phaseFigTemplate()
        hists_mean.append(np.mean(hists, 0))
        hists_std.append(np.std(hists, 0))
        plot(hist_ph, hists_mean[0], 'k', linewidth=2., zorder=1)
        gca().fill_between(hist_ph, hists_mean[0]+hists_std[0],
                hists_mean[0]-hists_std[0],
            facecolor='black', alpha=0.1)
        plot(hist_ph, hists[0], 'k--', dashes=(5, 2.5), zorder=2)
        ylabel('p(spike)')
        ylim([-0.01, 0.7])
        yticks([0, 0.7])
        savefig(avg_fname)
        
        
    saveIgor = True
    if saveIgor:
        from tables import *
        h5file = openFile(output_fname + '_igor_export.h5', mode = "w", title =
                "Attractor export figures")

        # Save slice of the raster so it is 1D
        raster_start = options.Ne**2/2
        raster_end = raster_start + options.Ne
        raster_x = np.ndarray((0))
        raster_y = np.ndarray((0))
        for n_it in xrange(raster_start, raster_end):
            raster_x = np.hstack((raster_x, spikeMon_e[n_it]))
            raster_y = np.hstack((raster_y, np.zeros((len(spikeMon_e[n_it]))) + n_it -
                raster_start))

        h5file.createArray(h5file.root, 'bump_raster_x', raster_x)
        h5file.createArray(h5file.root, 'bump_raster_y', raster_y)

        # Bump snapshot
        snapTime = 1.5*second
        firingSnapshot = np.reshape(Fe[:, snapTime/F_dt], (ei_net.Ne_y, ei_net.Ne_x))
        h5file.createArray(h5file.root, 'bump_snapshot', firingSnapshot)

        # Unfiltered input current, rasters and wavelets
        for n_it in range(len(theta_state_record_e)):
            nm_e_t = 'gamma_inhib_{0}_t'.format(n_it)
            nm_e_values = 'gamma_inhib_{0}_values'.format(n_it)
            nm_ras_e_ph = 'gamma_ras_{0}_phases'.format(n_it)
            nm_ras_e_trials = 'gamma_ras_{0}_trials'.format(n_it)
            nm_h_e = 'gamma_hist_{0}'.format(n_it)
            nm_wave_e = 'gamma_wavelet_e_{0}'.format(n_it)
            nm_sig_ph_e = 'gamma_sig_phase_e_{0}'.format(n_it)
            nm_sig_ph_mn_e = 'gamma_sig_phase_mean_e_{0}'.format(n_it)
            nm_sig_ph_std_e = 'gamma_sig_phase_std_e_{0}'.format(n_it)
            it = theta_state_record_e[n_it]
            h5file.createArray(h5file.root, nm_e_t, theta_stateMon_Iclamp_e.times)
            h5file.createArray(h5file.root, nm_e_values,
                    theta_stateMon_Iclamp_e[it]/pA)

            # Rasters
            h5file.createArray(h5file.root, nm_ras_e_ph, raster_list_e[n_it][0])
            h5file.createArray(h5file.root, nm_ras_e_trials, raster_list_e[n_it][1])

            # Hists
            h5file.createArray(h5file.root, nm_h_e, hist_list_e[n_it])


            # Wavelets and signals
            h5file.createArray(h5file.root, nm_wave_e, wavelet_list_e[n_it])
            h5file.createArray(h5file.root, nm_sig_ph_e,
                    sig_phase_list_e[n_it][0])
            h5file.createArray(h5file.root, nm_sig_ph_mn_e,
                    np.mean(sig_phase_list_e[n_it][0], 0))
            h5file.createArray(h5file.root, nm_sig_ph_std_e,
                    np.std(sig_phase_list_e[n_it][0], 0))

        nm_h_mean_e = 'gamma_hists_mean_e'
        nm_h_std_e = 'gamma_hists_std_e'
        h5file.createArray(h5file.root, nm_h_mean_e, hists_mean[0])
        h5file.createArray(h5file.root, nm_h_std_e, hists_std[0])

        for n_it in range(len(theta_state_record_i)):
            nm_i_t = 'gamma_exc_{0}_t'.format(n_it)
            nm_i_values = 'gamma_exc_{0}_values'.format(n_it)
            it = theta_state_record_i[n_it]
            h5file.createArray(h5file.root, nm_i_t, theta_stateMon_Iclamp_i.times)
            h5file.createArray(h5file.root, nm_i_values,
                    theta_stateMon_Iclamp_i[it]/pA)


        nm_h_ph = 'gamma_hist_ph'
        nm_sig_ph_ph = 'gamma_sig_pase_phase'
        h5file.createArray(h5file.root, nm_h_ph , hist_ph)
        h5file.createArray(h5file.root, nm_sig_ph_ph , sig_phase_list_e[0][1])


        # Export frequency and phase labels
        h5file.createArray(h5file.root, 'freq_label' , expt_freq)
        h5file.createArray(h5file.root, 'phase_label' , expt_phases)


        # Export connection profiles in the middle of the sheet
        n_e_it = ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1
        n_i_it = ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1
        AMPA_export = np.reshape(ei_net.AMPA_conn.W.todense()[n_e_it, :],
                (ei_net.Ni_y, ei_net.Ni_x))
        GABA_export = np.reshape(ei_net.GABA_conn1.W.todense()[n_i_it, :] +
            ei_net.extraGABA_conn1.W.todense()[n_i_it, :], (ei_net.Ne_y,
                ei_net.Ne_x))

        h5file.createArray(h5file.root, 'AMPA_profile', AMPA_export)
        h5file.createArray(h5file.root, 'GABA_profile', GABA_export)



        h5file.close()

        print "Dump after " + str(simulationClock.t)


    print "End of trial no. " + str(trial_it) + "..."
    print 

    ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

