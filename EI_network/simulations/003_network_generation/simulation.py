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
parser.add_option("--prefDirC", type=float,
        help="Preferred directtion multiplier")

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
    else:
        ei_net.E_pop.Iext = [ei_net.E_pop.Iext[0]] * len(ei_net.E_pop)
    pass


v = np.array([[0, 1]]).T
Ivel = np.dot(ei_net.prefDirs_e, v) * options.Ivel*pA

@network_operation(simulationClock)
def thetaStimulation():
    if simulationClock.t >= theta_start_t and simulationClock.t < options.time*second:
        ph = stim_omega*simulationClock.t
        ei_net.E_pop.Iext = stim_e_DC + stim_e_A*np.sin(ph - np.pi/2)
        ei_net.I_pop.Iext = stim_i_DC + stim_i_A*np.sin(ph - np.pi/2)

    # Velocity inputs
    if simulationClock.t >= theta_start_t and simulationClock.t < options.time*second:
        #ei_net.I_pop.Iext += Ivel.ravel()
        ei_net.E_pop.Iext += np.max((Ivel.ravel(), np.zeros(len(Ivel))), 0)


state_record_e = [ei_net.Ne_x/2 -1 , ei_net.Ne_y/2*ei_net.Ne_x + ei_net.Ne_x/2 - 1]
state_record_i = [ei_net.Ni_x/2 - 1, ei_net.Ni_y/2*ei_net.Ni_x + ei_net.Ni_x/2 - 1]

#spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
#spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)

#stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
#stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
#stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=simulationClock)
#stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=simulationClock)
#stateMon_Iext_e = StateMonitor(ei_net.E_pop, 'Iext', record = state_record_e, clock=simulationClock)
#stateMon_Iext_i = StateMonitor(ei_net.I_pop, 'Iext', record = state_record_i, clock=simulationClock)

theta_n_it_range = 2
theta_state_record_e = range(state_record_e[1] - theta_n_it_range/2,
        state_record_e[1] + theta_n_it_range/2 + 1)
theta_state_record_i = range(state_record_i[0] - theta_n_it_range/2,
        state_record_i[0] + theta_n_it_range/2 + 1)
theta_spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
theta_spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)
theta_stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = theta_state_record_e, clock=simulationClock)
theta_stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = theta_state_record_i, clock=simulationClock)


#ei_net.net.add(spikeMon_e, spikeMon_i)
#ei_net.net.add(stateMon_e, stateMon_i, stateMon_Iclamp_e, stateMon_Iclamp_i)
#ei_net.net.add(stateMon_Iext_e, stateMon_Iext_i)
ei_net.net.add(stimulateSubPopulation)
ei_net.net.add(thetaStimulation)

ei_net.net.add(theta_spikeMon_e, theta_spikeMon_i, theta_stateMon_Iclamp_e,
        theta_stateMon_Iclamp_i)


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


x_lim = [options.time-0.5, options.time]
#x_lim = [0, options.time]

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
    theta_spikeMon_i.reinit()
    theta_stateMon_Iclamp_e.reinit()
    theta_stateMon_Iclamp_i.reinit()

    print "  Theta stimulation..."
    ei_net.net.run(options.time*second - theta_start_mon_t, report='stdout',
            report_period=options.update_interval*second)
    duration=time.time()-start_time
    print "Simulation time:",duration,"seconds"
    
    
    output_fname = "{0}/{1}job{2:04}_trial{3:04}".format(options.output_dir,
            options.fileNamePrefix, options.job_num, trial_it)
    
    
    F_tstart = 0
    F_tend = options.time
    F_dt = 0.02
    F_winLen = 0.25
    #Fe, Fe_t = spikeMon_e.getFiringRate(F_tstart, F_tend, F_dt, F_winLen) 
    #Fi, Fi_t = spikeMon_i.getFiringRate(F_tstart, F_tend, F_dt, F_winLen)

    ## plot firing rates
    #figure(figsize=figSize)
    #subplot(211)
    #T, FR = np.meshgrid(Fe_t, np.arange(ei_net.net_Ne))
    #pcolormesh(T, FR, Fe)
    #ylabel('E Neuron no.')
    #colorbar()
    #subplot(212)
    #T, FR = np.meshgrid(Fi_t, np.arange(ei_net.net_Ni))
    #pcolormesh(T, FR, Fi)
    #xlabel('Time (s)')
    #ylabel('I Neuron no.')
    #colorbar()
    #savefig(output_fname + '_firing_rate_e.png')

    #figure()
    #ax = subplot(211)
    #plot(stateMon_e.times, stateMon_e.values[0:2].T/mV)
    #ylabel('E membrane potential (mV)')
    #subplot(212, sharex=ax)
    #plot(stateMon_i.times, stateMon_i.values[0:2].T/mV)
    #xlabel('Time (s)')
    #ylabel('I membrane potential (mV)')
    #xlim(x_lim)
    #savefig(output_fname + '_Vm.pdf')
    #
    #
    #figure()
    #ax = subplot(211)
    #plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[0:2].T/pA)
    #ylabel('E synaptic current (pA)')
    #subplot(212, sharex=ax)
    #plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[0:2].T/pA)
    #xlabel('Time (s)')
    #ylabel('I synaptic current (pA)')
    #xlim(x_lim)
    #savefig(output_fname + '_Isyn.pdf')
    #
    #figure()
    #ax = subplot(211)
    #plot(stateMon_Iext_e.times, -stateMon_Iext_e.values[0].T/pA)
    #ylabel('E external current (pA)')
    #subplot(212, sharex=ax)
    #plot(stateMon_Iext_i.times, -stateMon_Iext_i.values[0].T/pA)
    #xlabel('Time (s)')
    #ylabel('I external current (pA)')
    #xlim(x_lim)
    #savefig(output_fname + '_Iext.pdf')
    #
    ## High pass filter these signals
    #figure()
    #ax = subplot(211)
    #plot(stateMon_Iclamp_e.times, butterHighPass(stateMon_Iclamp_e.values[1].T/pA, options.sim_dt, 40))
    #plot(stateMon_Iext_e.times, -(stateMon_Iext_e.values[0]/pA - stim_e_DC/pA))
    #ylabel('E current (pA)')
    #ylim([-500, 500])
    #subplot(212, sharex=ax)
    #plot(stateMon_Iclamp_i.times, butterHighPass(stateMon_Iclamp_i.values[0].T/pA, options.sim_dt, 40))
    ##plot(stateMon_Iclamp_i.times, stateMon_Iext_i.values[0]/pA)
    #xlabel('Time (s)')
    #ylabel('I current (pA)')
    #xlim(x_lim)
    #ylim([-500, 500])
    #savefig(output_fname + '_Isyn_filt.pdf')
    
    
    
    Ne = options.Ne
    figure()
    pcolormesh(np.reshape(ei_net.AMPA_conn.W.todense()[57*68 + 33, :], (ei_net.Ni_y,
        ei_net.Ni_x)));
    xlabel('I neuron no.')
    ylabel('I neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_E2I_conn.png')

    Ni = options.Ni
    figure()
    pcolormesh(np.reshape(ei_net.GABA_conn1.W.todense()[0, :], (ei_net.Ne_y,
        ei_net.Ne_x)));
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_I2E_conn.png')

    #figure()
    #pcolormesh(np.reshape(np.dot(ei_net.AMPA_conn.W.todense(),
    #    ei_net.GABA_conn1.W.todense())[15, :], (ei_net.Ne_y, ei_net.Ne_x)));
    #xlabel('E neuron no.')
    #ylabel('E neuron no.')
    #colorbar()
    #savefig(output_fname + '_E2E_conn.png')


    figure()
    pcolormesh(np.reshape(Fe[:, len(Fe_t)/2], (ei_net.Ne_y, ei_net.Ne_x)))
    xlabel('E neuron no.')
    ylabel('E neuron no.')
    colorbar()
    axis('equal')
    savefig(output_fname + '_firing_snapshot_e.png')


    #
    ## Print a plot of bump position
    #(pos, bumpPos_times) = theta_spikeMon_e.torusPopulationVector(ei_net.o.Ne,
    #        theta_start_t, options.time, F_dt, F_winLen)
    #figure(figsize=figSize)
    #plot(bumpPos_times, pos)
    #xlabel('Time (s)')
    #ylabel('Bump position (neurons)')
    #axis('equal')
    #
    #savefig(output_fname + '_bump_position.pdf')

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
                options.sim_dt, high_pass_freq), 1./stim_freq, options.sim_dt, maxF)

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

        hist_nbins = 1./stim_freq/raster_bin_size
       
        for n_it in range(len(range_n_it)):
            neuron_no = range_n_it[n_it]
            print('Saving rasters for neuron no. ' + str(neuron_no))
            phases, times, trials = spikePhaseTrialRaster(tmp_spikeMon[neuron_no],
                    stim_freq, theta_start_mon_t)

            phases -= np.pi

            raster_list.append((phases, trials+1))
        
            # Raster plots (single cell over 'theta' epochs)
            ntrials = np.ceil((options.time - theta_start_mon_t) * stim_freq)
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
    ##
    ##
    ##outData = {}
    ###outData['timeSnapshot'] = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    ###if spikeMon_e is not None:
    ###    outData['spikeCell_e'] = spikeMon_e.aspikes
    ###if spikeMon_i is not None:
    ###    outData['spikeCell_i'] = spikeMon_i.aspikes
    ###outData['options'] = options._einet_optdict

    ##outData['bumpPos'] = pos
    ##outData['bumpPos_times'] = bumpPos_times

    ###outData['Fe'] = Fe
    ###outData['Fe_t'] = Fe_t
    ###outData['Fi'] = Fi
    ###outData['Fi_t'] = Fi_t

    ##outData['theta_stateMon_Iclamp_e_times'] = theta_stateMon_Iclamp_e.times
    ##outData['theta_stateMon_Iclamp_e_values'] = theta_stateMon_Iclamp_e.values
    ##outData['theta_stateMon_Iclamp_i_times'] = theta_stateMon_Iclamp_i.times
    ##outData['theta_stateMon_Iclamp_i_values'] = theta_stateMon_Iclamp_i.values
    ##
    ##savemat(output_fname + '_output.mat', outData, do_compression=True)

    #saveIgor = False
    #if saveIgor:
    #    from tables import *
    #    h5file = openFile(output_fname + 'igor_export.h5', mode = "w", title =
    #            "Attractor export figures")

    #    # Save slice of the raster so it is 1D
    #    raster_start = options.Ne**2/2
    #    raster_end = raster_start + options.Ne
    #    raster_x = np.ndarray((0))
    #    raster_y = np.ndarray((0))
    #    for n_it in xrange(raster_start, raster_end):
    #        raster_x = np.hstack((raster_x, spikeMon_e[n_it]))
    #        raster_y = np.hstack((raster_y, np.zeros((len(spikeMon_e[n_it]))) + n_it -
    #            raster_start))

    #    h5file.createArray(h5file.root, 'bump_raster_x', raster_x)
    #    h5file.createArray(h5file.root, 'bump_raster_y', raster_y)

    #    # Bump snapshot
    #    snapTime = 0.5*second
    #    firingSnapshot = np.reshape(Fe[:, snapTime/F_dt], (options.Ne, options.Ne))
    #    h5file.createArray(h5file.root, 'bump_snapshot', firingSnapshot)

    #    # Unfiltered input current, rasters and wavelets
    #    for n_it in range(len(theta_state_record_e)):
    #        nm_e_t = 'gamma_inhib_{0}_t'.format(n_it)
    #        nm_e_values = 'gamma_inhib_{0}_values'.format(n_it)
    #        nm_ras_e_ph = 'gamma_ras_{0}_phases'.format(n_it)
    #        nm_ras_e_trials = 'gamma_ras_{0}_trials'.format(n_it)
    #        nm_h_e = 'gamma_hist_{0}'.format(n_it)
    #        nm_wave_e = 'gamma_wavelet_e_{0}'.format(n_it)
    #        nm_sig_ph_e = 'gamma_sig_phase_e_{0}'.format(n_it)
    #        nm_sig_ph_mn_e = 'gamma_sig_phase_mean_e_{0}'.format(n_it)
    #        nm_sig_ph_std_e = 'gamma_sig_phase_std_e_{0}'.format(n_it)
    #        it = theta_state_record_e[n_it]
    #        h5file.createArray(h5file.root, nm_e_t, theta_stateMon_Iclamp_e.times)
    #        h5file.createArray(h5file.root, nm_e_values,
    #                theta_stateMon_Iclamp_e[it]/pA)

    #        # Rasters
    #        h5file.createArray(h5file.root, nm_ras_e_ph, raster_list_e[n_it][0])
    #        h5file.createArray(h5file.root, nm_ras_e_trials, raster_list_e[n_it][1])

    #        # Hists
    #        h5file.createArray(h5file.root, nm_h_e, hist_list_e[n_it])


    #        # Wavelets and signals
    #        h5file.createArray(h5file.root, nm_wave_e, wavelet_list_e[n_it])
    #        h5file.createArray(h5file.root, nm_sig_ph_e,
    #                sig_phase_list_e[n_it][0])
    #        h5file.createArray(h5file.root, nm_sig_ph_mn_e,
    #                np.mean(sig_phase_list_e[n_it][0], 0))
    #        h5file.createArray(h5file.root, nm_sig_ph_std_e,
    #                np.std(sig_phase_list_e[n_it][0], 0))

    #    nm_h_mean_e = 'gamma_hists_mean_e'
    #    nm_h_std_e = 'gamma_hists_std_e'
    #    h5file.createArray(h5file.root, nm_h_mean_e, hists_mean[0])
    #    h5file.createArray(h5file.root, nm_h_std_e, hists_std[0])

    #    for n_it in range(len(theta_state_record_i)):
    #        nm_i_t = 'gamma_exc_{0}_t'.format(n_it)
    #        nm_i_values = 'gamma_exc_{0}_values'.format(n_it)
    #        it = theta_state_record_i[n_it]
    #        h5file.createArray(h5file.root, nm_i_t, theta_stateMon_Iclamp_i.times)
    #        h5file.createArray(h5file.root, nm_i_values,
    #                theta_stateMon_Iclamp_i[it]/pA)


    #    nm_h_ph = 'gamma_hist_ph'
    #    nm_sig_ph_ph = 'gamma_sig_pase_phase'
    #    h5file.createArray(h5file.root, nm_h_ph , hist_ph)
    #    h5file.createArray(h5file.root, nm_sig_ph_ph , sig_phase_list_e[0][1])


    #    # Export frequency and phase labels
    #    h5file.createArray(h5file.root, 'freq_label' , expt_freq)
    #    h5file.createArray(h5file.root, 'phase_label' , expt_phases)

    #    h5file.close()



    print "End of trial no. " + str(trial_it) + "..."
    print 

    #ei_net.reinit()
#                            End main cycle
################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"

