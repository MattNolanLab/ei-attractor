import brian_no_units
from brian import *

from matplotlib.backends.backend_pdf import PdfPages

from scipy.io import savemat
from optparse import OptionParser

import time
import numpy as np
import logging as lg

from EI_network import *
from EI_network_sim_mod import *
from custombrian import *
from plotting import *
from tools import *

lg.basicConfig(level=lg.DEBUG)


parser = getOptParser()

parser.add_option("--net_generations", type="int",
        help="Number of repetitions of network generation. Each generation has"
            "ntrials*3 runs of the same network")

(options, args) = parser.parse_args()
options = setOptionDictionary(parser, options)
print options

# Clock definitions
sim_dt = options.sim_dt*second
simulationClock = Clock(dt=sim_dt)

stim_omega = nan

# Other
x_lim = [options.time-1, options.time]


################################################################################
#                      Stimulation frequency list (Hz)

stim_freq_list = numpy.array([2, 8, 16])

################################################################################


total_start_t = time.time()


outputDir = createJobDir(options)
output_fname_prefix = '{0}/{1}'.format(outputDir,
        options.fileNamePrefix)
output_fname_job = output_fname_prefix + 'job{0:04}'.format(options.job_num)


for net_it in xrange(options.net_generations):
    output_fname_gen = output_fname_job + '_gen{0:02}'.format(net_it)

    print 'Network generation no. ' + str(net_it)
    ################################################################################
    #                              Network setup
    ################################################################################
    print "Starting network and connections initialization..."
    start_time=time.time()
    
    options.ndim = -1
    ei_net = EI_Network(options, simulationClock)
    
    ei_net.connRandom(options.AMPA_density, options.GABA_density)
    
    # After network generation, save connections
    print "Saving connections..."
    connItems = [
            ['AMPA_conn_W', ei_net.AMPA_conn.W],
            ['GABA_conn1_W', ei_net.GABA_conn1.W],
            ['options', options._einet_optdict]]
    matFormatSaver(output_fname_gen + "_connections.mat", connItems,
            do_compression=True)
    print "Done"
    
    
    duration=time.time()-start_time
    print "Network setup time:",duration,"seconds"
    #                            End Network setup
    ################################################################################
    
    @network_operation(simulationClock)
    def thetaStimulation():
        ph = stim_omega*simulationClock.t
        ei_net.E_pop.Iext = options.Iext_e/2*np.sin(ph - np.pi/2) + options.Iext_e/2
        ei_net.I_pop.Iext = options.Iext_i/2*np.sin(ph - np.pi/2) + options.Iext_i/2
        #pass
    
    
    state_record_e = range(10)
    state_record_i = range(10)
    
    spikeMon_e = ExtendedSpikeMonitor(ei_net.E_pop)
    spikeMon_i = ExtendedSpikeMonitor(ei_net.I_pop)
    stateMon_e = StateMonitor(ei_net.E_pop, 'vm', record = state_record_e, clock=simulationClock)
    stateMon_i = StateMonitor(ei_net.I_pop, 'vm', record = state_record_i, clock=simulationClock)
    stateMon_Iclamp_e = StateMonitor(ei_net.E_pop, 'Iclamp', record = state_record_e, clock=simulationClock)
    stateMon_Iclamp_i = StateMonitor(ei_net.I_pop, 'Iclamp', record = state_record_i, clock=simulationClock)
    stateMon_Iext_e = StateMonitor(ei_net.E_pop, 'Iext', record=state_record_e,
            clock=simulationClock)
    stateMon_Iext_i = StateMonitor(ei_net.I_pop, 'Iext', record=state_record_i,
            clock=simulationClock)
    
    ei_net.net.add(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i, stateMon_Iclamp_e,
            stateMon_Iclamp_i)
    ei_net.net.add(stateMon_Iext_e)
    ei_net.net.add(stateMon_Iext_i)
    
    ei_net.net.add(thetaStimulation)
    
    
    ################################################################################
    #                              Main cycle
    ################################################################################
    for trial_it in range(ei_net.o.ntrials):
        print "Starting trial no. " + str(trial_it) + "..."
        close('all')

        F_mean_e_vec = np.ndarray(len(stim_freq_list))
        F_mean_i_vec = np.ndarray(len(stim_freq_list))
        F_std_e_vec = np.ndarray(len(stim_freq_list))
        F_std_i_vec = np.ndarray(len(stim_freq_list))

        # Frequencies and phases of power maxima extracted from wavelets
        w_Fmax_e_vec = []
        w_Fmax_i_vec = []
        w_phmax_e_vec= []
        w_phmax_i_vec= []

        for stim_freq_it in xrange(len(stim_freq_list)):
            stim_freq = stim_freq_list[stim_freq_it]

            print "Stimulation frequency: " + str(stim_freq) + ' Hz'
            stim_omega = 2*np.pi*stim_freq*Hz
    
            print "Simulation running..."
            start_time=time.time()
            
            ei_net.net.run(options.time*second, report='stdout',
                    report_period=options.update_interval*second)
            duration=time.time()-start_time
            print "Simulation time:",duration,"seconds"
            
            
            output_fname_trial = output_fname_gen + '_trial{0:04}'.format(trial_it)
            output_fname = output_fname_trial + "_stim{0}".format(stim_freq)

            ## Save current and voltage traces
            #printAndSaveTraces(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i,
            #    stateMon_Iclamp_e, stateMon_Iclamp_i, stateMon_Iext_e, stateMon_Iext_i,
            #    options, output_fname, x_lim)
    
            ## Firing rate
            #Favg_e = spikeMon_e.getNSpikes()/options.time
            #mean_e = np.mean(Favg_e)
            #F_mean_e_vec[stim_freq_it] = mean_e
            #F_std_e_vec[stim_freq_it] = np.std(Favg_e)
            #Favg_i = spikeMon_i.getNSpikes()/options.time
            #mean_i = np.mean(Favg_i)
            #F_mean_i_vec[stim_freq_it] = mean_i
            #F_std_i_vec[stim_freq_it] = np.std(Favg_i)

            #printFiringRatesBar(Favg_e, Favg_i, mean_e, mean_i, output_fname)

            print "Wavelet analysis..."
            wav_n_range = 10
            wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_e.pdf')
            high_pass_freq = 40.
            maxF = 200
            for ei_it in [0, 1]:
                if ei_it == 0:
                    print '  E neurons...'
                    wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_e.pdf')
                    wavelet_sig_fname = output_fname + '_phase_wavelet_e'
                    max_fname = output_fname + '_Fmax_scatter_e.pdf'
                    tmp_stateMon = stateMon_Iclamp_e
                    w_Fmax_vec = w_Fmax_e_vec
                    w_phmax_vec = w_phmax_e_vec
                else:
                    print '  I neurons...'
                    wavelet_sig_pp = PdfPages(output_fname + '_phase_sig_i.pdf')
                    wavelet_sig_fname = output_fname + '_phase_wavelet_i'
                    max_fname = output_fname + '_Fmax_scatter_i.pdf'
                    tmp_stateMon = stateMon_Iclamp_i
                    w_Fmax_vec = w_Fmax_i_vec
                    w_phmax_vec = w_phmax_i_vec

                w_Fmax = np.ndarray(wav_n_range)
                w_phmax= np.ndarray(wav_n_range)
                for n_it in xrange(wav_n_range):
                    print('    Neuron no. ' + str(n_it))
                    cwt_phases, sig_cwt, freq, sig_ph = \
                        phaseCWT(butterHighPass(tmp_stateMon.values[n_it].T/pA,
                        options.sim_dt, high_pass_freq), 1./stim_freq, options.sim_dt, maxF)

                    w_max = sig_cwt.argmax()
                    w_Fmax[n_it] = freq[w_max//len(sig_cwt[0])]
                    w_phmax[n_it]= cwt_phases[np.mod(w_max, len(sig_cwt[0]))]
                            

                    # Wavelet plot
                    f = phaseFigTemplate()
                    PH, F = np.meshgrid(cwt_phases, freq)
                    pcolormesh(PH, F, sig_cwt, edgecolors='None', cmap=get_cmap('jet'))
                    ylabel('F (Hz)')
                    ylim([0, maxF])
                    savefig(wavelet_sig_fname + '{0}.png'.format(n_it),
                            dpi=300)
                    close()

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

                # Insert this to avg and std vectors and make a plot
                w_Fmax_vec.append(w_Fmax)
                w_phmax_vec.append(w_phmax)
                f = phaseFigTemplate()
                plot(w_phmax, w_Fmax, 'ko', markersize=10, alpha=0.25)
                errorbar(np.mean(w_phmax), np.mean(w_Fmax), np.std(w_Fmax),
                        np.std(w_phmax), 'ko', markersize=10)
                ylim([50, 120])
                savefig(max_fname)
            print "Done"
    
            
            for ei_it in [0, 1]:
                if ei_it == 0:
                    raster_pp = PdfPages(output_fname + '_phase_raster_e.pdf')
                    avg_fname = output_fname + '_phase_pspike_avg_e.pdf'
                    n_range = 20
                    tmp_spikeMon = spikeMon_e
                else:
                    raster_pp = PdfPages(output_fname + '_phase_raster_i.pdf')
                    avg_fname = output_fname + '_phase_pspike_avg_i.pdf'
                    n_range = 10
                    tmp_spikeMon = spikeMon_i

                hist_nbins = 1./stim_freq/raster_bin_size
                hists = [] #np.ndarray((n_range, hist_nbins)) #    hist_ph = np.ndarray(hist_nbins)
    
                for n_it in xrange(n_range):
                    print('Saving rasters for neuron no. ' + str(n_it))
                    phases, times, trials = spikePhaseTrialRaster(tmp_spikeMon[n_it],
                            stim_freq)
    
                    # Raster plots (single cell over 'theta' epochs)
                    ntrials = np.ceil(options.time * stim_freq)
                    f = rasterPhasePlot(phases - np.pi, trials)
                    raster_pp.savefig()
                    close()
    
                    ## Histograms
                    #figure(figsize=small_plot_figsize)
                    if (len(phases) != 0):
                        h = hist(phases, hist_nbins, [0, 2*np.pi],
                                normed=False)
                        hists.append(h[0]/double(len(phases)))
                    hist_ph = h[1][0:len(h[1])-1] - np.pi
                    delaxes(gca())
                raster_pp.close()

                # Average histogram + one neuron
                f = phaseFigTemplate()
                h_avg = np.mean(hists, 0)
                h_std = np.std(hists, 0)
                plot(hist_ph, h_avg, 'k', linewidth=2., zorder=1)
                gca().fill_between(hist_ph, h_avg+h_std, h_avg-h_std,
                    facecolor='black', alpha=0.1)
                plot(hist_ph, hists[0], 'k--', dashes=(5, 2.5), zorder=2)
                ylabel('p(spike)')
                ylim([-0.01, 0.7])
                yticks([0, 0.7])
                savefig(avg_fname)

    
            
            ei_net.reinit()
            print "done"


        # Bar plot of mean firing rates for different stimulation frequencies
        f = firingRateBarPlot(stim_freq_list, F_mean_e_vec, F_std_e_vec)
        savefig(output_fname_trial + '_Fmean_freq_bar_e.pdf')
    
        f = firingRateBarPlot(stim_freq_list, F_mean_i_vec, F_std_i_vec)
        savefig(output_fname_trial + '_Fmean_freq_bar_i.pdf')

        # Bar plot of mean maximal power for different stimulation frequencies
    
        print "End of trial no. " + str(trial_it) + "..."
        print 



    #                            End main cycle
    ################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"


