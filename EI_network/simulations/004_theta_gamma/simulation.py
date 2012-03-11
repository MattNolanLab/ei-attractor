import brian_no_units
from os import system
from brian import *

from matplotlib.backends.backend_pdf import PdfPages

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

def spikePhaseTrialRaster(spikeTimes, f):
    '''Here assuming that phase(t=0) = 0'''
    trials = np.floor(f*spikeTimes)
    phases = np.mod(2*np.pi*f*spikeTimes, 2*np.pi)
    times  = np.mod(spikeTimes, 1./f)
    return (phases, times, trials)

def set_axis_params(ax):
    ax.tick_params(direction='out', length=6, zorder=0)
    ax.tick_params(bottom=True, top=False, left=True, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.margins(0.05, tight=False)

def createJobDir(options):
    # Create job directory in options.output_dir/jobXXXX
    outputDir = options.output_dir + '/job{0:04}'.format(options.job_num) +'/'
    ec = system('mkdir ' + outputDir)
    if ec != 0:
        print "Could not create output directory: " + outputDir
        ec = system('ls ' + outputDir)
        if ec == 0:
            print "But it can be listed --> continuing!"
        else:
            print "And it cannot be listed. Check your permissions and rerun!"
            exit(1)
    return outputDir



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
stimClock = Clock(50*msecond)

stim_omega = nan

# Other
figSize = (12,8)
x_lim = [options.time-1, options.time]

small_plot_figsize = (3.5, 2.75)
small_plot_axsize = [0.25, 0.15, 0.70, 0.80]
small_plot_fontsize = 16
small_plot_texsize = 25
raster_bin_size = 2e-3

rcParams['font.size'] = small_plot_fontsize

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
    
    
    state_record_e = range(int(len(ei_net.E_pop)/20.))
    state_record_i = range(int(len(ei_net.I_pop)/10.))
    
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

        F_mean_e_vec = np.ndarray(len(stim_freq_list))
        F_mean_i_vec = np.ndarray(len(stim_freq_list))
        F_std_e_vec = np.ndarray(len(stim_freq_list))
        F_std_i_vec = np.ndarray(len(stim_freq_list))

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
            close()
            
    
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
            close()
            
    
            # Band pass filter these signals
            figure()
            ax = subplot(211)
            plot(stateMon_Iclamp_e.times, butterHighPass(stateMon_Iclamp_e.values[0].T/pA +
                    stateMon_Iext_e.values[0].T/pA, options.sim_dt, 40))
            #plot(stateMon_Iclamp_e.times, stateMon_Iext_e.values[0]/pA)
            ylabel('E current (pA)')
            ylim([-500, 500])
            subplot(212, sharex=ax)
            plot(stateMon_Iclamp_i.times, butterHighPass(stateMon_Iclamp_i.values[0].T/pA +
                    stateMon_Iext_i.values[0].T/pA, options.sim_dt, 40))
            #plot(stateMon_Iclamp_i.times, stateMon_Iext_i.values[0]/pA)
            xlabel('Time (s)')
            ylabel('I current (pA)')
            xlim(x_lim)
            ylim([-500, 500])
            savefig(output_fname + '_Isyn_filt.pdf')
            close()
            
    
            # Firing rate
            Favg_e = spikeMon_e.getNSpikes()/options.time
            mean_e = np.mean(Favg_e)
            F_mean_e_vec[stim_freq_it] = mean_e
            F_std_e_vec[stim_freq_it] = np.std(Favg_e)
            Favg_i = spikeMon_i.getNSpikes()/options.time
            mean_i = np.mean(Favg_i)
            F_mean_i_vec[stim_freq_it] = mean_i
            F_std_i_vec[stim_freq_it] = np.std(Favg_i)
            figure()
            subplot(121)
            h = hist(Favg_e, 20)
            xlabel('E f. rate (Hz)')
            ylabel('Count')
            title('Average: ' + str(mean_e) + ' Hz')
            subplot(122)
            hist(Favg_i, 20)
            xlabel('I f. rate (Hz)')
            title('Average: ' + str(mean_i) + ' Hz')
            savefig(output_fname + '_Fhist.pdf')
            close()
    
    
            
            for ei_it in [0, 1]:
                if ei_it == 0:
                    raster_pp = PdfPages(output_fname + '_phase_raster_e.pdf')
                    #pspike_pp = PdfPages(output_fname + '_phase_pspike_e.pdf')
                    avg_fname = output_fname + '_phase_pspike_avg_e.pdf'
                    n_range = int(len(ei_net.E_pop)/10)
                    tmp_spikeMon = spikeMon_e
                else:
                    raster_pp = PdfPages(output_fname + '_phase_raster_i.pdf')
                    #pspike_pp = PdfPages(output_fname + '_phase_pspike_i.pdf')
                    avg_fname = output_fname + '_phase_pspike_avg_i.pdf'
                    n_range = int(len(ei_net.I_pop)/10)
                    tmp_spikeMon = spikeMon_i

                hist_nbins = 1./stim_freq/raster_bin_size
                hists = [] #np.ndarray((n_range, hist_nbins))
                hist_ph = np.ndarray(hist_nbins)
    
                for n_it in xrange(n_range):
                    print('Saving rasters for neuron no. ' + str(n_it))
                    phases, times, trials = spikePhaseTrialRaster(tmp_spikeMon[n_it],
                            stim_freq)
    
                    # Raster plots (single cell over 'theta' epochs)
                    ntrials = np.ceil(options.time * stim_freq)
                    figure(figsize=small_plot_figsize)
                    axes(small_plot_axsize)
                    plot(phases - np.pi, trials, 'k|', markeredgewidth=3)
                    set_axis_params(gca())
                    ylabel('Trial')
                    xlim([-np.pi, np.pi])
                    ylim([-1, ntrials+1])
                    xticks([-np.pi, 0, np.pi], ('$-\pi$', '',  '$\pi$'), fontsize=25)
                    yticks([0, ntrials])
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
                    #axes(small_plot_axsize)
                    #plot(hist_ph, hists[n_it, :], 'k--', dashes=(5,
                    #    2.5))
                    #set_axis_params(gca())
                    #ylabel('p(spike)')
                    #xlim([-np.pi, np.pi])
                    #ylim([-0.01, 0.7])
                    #yticks([0, 0.7])
                    #xticks([-np.pi, 0, np.pi], ('$-\pi$', '',  '$\pi$'), fontsize=small_plot_texsize)
                    #pspike_pp.savefig()
                    #close()
                raster_pp.close()
                #pspike_pp.close()

                # Average histogram + one neuron
                figure(figsize=small_plot_figsize)
                ax = axes(small_plot_axsize)
                h_avg = np.mean(hists, 0)
                h_std = np.std(hists, 0)
                plot(hist_ph, h_avg, 'k', linewidth=2., zorder=1)
                ax.fill_between(hist_ph, h_avg+h_std, h_avg-h_std,
                    facecolor='black', alpha=0.1)
                plot(hist_ph, hists[0], 'k--', dashes=(5, 2.5), zorder=2)
                set_axis_params(gca())
                ylabel('p(spike)')
                xlim([-np.pi, np.pi])
                ylim([-0.01, 0.7])
                yticks([0, 0.7])
                xticks([-np.pi, 0, np.pi], ('$-\pi$', '',  '$\pi$'), fontsize=small_plot_texsize)
                savefig(avg_fname)
                close()

    
            
            print "Saving output data to matlab file..."
            saveItems = [
                    ['spikeMon_e', spikeMon_e],
                    ['spikeMon_i', spikeMon_i],
                    ['stateMon_e', stateMon_e],
                    ['stateMon_i', stateMon_i],
                    ['stateMon_Iclamp_e', stateMon_Iclamp_e],
                    ['stateMon_Iclamp_i', stateMon_Iclamp_i],
                    ['stateMon_Iext_e', stateMon_Iext_e],
                    ['stateMon_Iext_i', stateMon_Iext_i],
                    ['options', options._einet_optdict]]
            matFormatSaver(output_fname + '_output.mat', saveItems, do_compression=True)
    
            ei_net.reinit()
            print "done"

        # Bar plot of mean firing rates for different stimulation frequencies
        figure(figsize=(2.5, 4))
        ax = axes(small_plot_axsize)
        bar(range(len(stim_freq_list)), F_mean_e_vec, color='k',
                yerr=F_std_e_vec, ecolor='k', align='center', width=0.8)
        xticks(range(len(stim_freq_list)), stim_freq_list)
        gca().tick_params(bottom=True, top=False, left=True, right=False)
        gca().spines['top'].set_visible(False)
        gca().spines['right'].set_visible(False)
        xlabel('Stim. freq. (Hz)')
        ylabel('F. rate (Hz)')
        ylim([0, max(F_std_e_vec+F_mean_e_vec)+10])
        savefig(output_fname_trial + '_Fmean_freq_bar_e.pdf')
        close()
    
        figure(figsize=(2.5, 4))
        ax = axes(small_plot_axsize)
        bar(range(len(stim_freq_list)), F_mean_i_vec, color='k',
                yerr=F_std_i_vec, ecolor='k', align='center', width=0.8)
        xticks(range(len(stim_freq_list)), stim_freq_list)
        gca().tick_params(bottom=True, top=False, left=True, right=False)
        gca().spines['top'].set_visible(False)
        gca().spines['right'].set_visible(False)
        xlabel('Stim. freq (Hz)')
        ylabel('F. rate (Hz)')
        ylim([0, max(F_std_i_vec+F_mean_i_vec)+10])
        savefig(output_fname_trial + '_Fmean_freq_bar_i.pdf')
        close()
    
        print "End of trial no. " + str(trial_it) + "..."
        print 



    #                            End main cycle
    ################################################################################

total_time = time.time()-total_start_t
print "Overall time: ", total_time, " seconds"


