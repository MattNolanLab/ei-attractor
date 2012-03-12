from matplotlib.pyplot import *
from brian import *

from tools import *

small_plot_figsize = (3.75, 2.75)
small_plot_axsize = [0.3, 0.15, 0.65, 0.80]
small_plot_fontsize = 16
small_plot_texsize = 25
raster_bin_size = 2e-3

rcParams['font.size'] = small_plot_fontsize

def set_axis_params(ax):
    ax.tick_params(direction='out', length=6, zorder=0)
    ax.tick_params(bottom=True, top=False, left=True, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.margins(0.05, tight=False)


def phaseFigTemplate():
    f = figure(figsize=small_plot_figsize)
    axes(small_plot_axsize)
    set_axis_params(gca())
    xlim([-np.pi, np.pi])
    xticks([-np.pi, 0, np.pi], ('$-\pi$', '',  '$\pi$'), fontsize=25)
    return f


def rasterPhasePlot(phases, trials):
    f = phaseFigTemplate()
    plot(phases, trials, 'k|', markeredgewidth=3)
    ylabel('Trial')
    ylim([-1, len(trials)+1])
    yticks([0, len(trials)])
    return f

def firingRateBarPlot(stim_freq_list, F_mean_vec, F_std_vec):
    f= figure(figsize=(2.5, 4))
    ax = axes(small_plot_axsize)
    bar(range(len(stim_freq_list)), F_mean_vec, color='k',
            yerr=F_std_vec, ecolor='k', align='center', width=0.8)
    xticks(range(len(stim_freq_list)), stim_freq_list)
    gca().tick_params(bottom=True, top=False, left=True, right=False)
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    xlabel('Stim. freq. (Hz)')
    ylabel('F. rate (Hz)')
    ylim([0, max(F_std_vec+F_mean_vec)+10])
    return f
    

def printAndSaveTraces(spikeMon_e, spikeMon_i, stateMon_e, stateMon_i,
        stateMon_Iclamp_e, stateMon_Iclamp_i, stateMon_Iext_e, stateMon_Iext_i,
        options, output_fname, x_lim):
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
    plot(stateMon_Iclamp_e.times, stateMon_Iclamp_e.values[0:2].T/pA + \
            stateMon_Iext_e.values[0:2].T/pA)
    ylabel('E synaptic current (pA)')
    subplot(212, sharex=ax)
    plot(stateMon_Iclamp_i.times, stateMon_Iclamp_i.values[0:2].T/pA + \
            stateMon_Iext_i.values[0:2].T/pA)
    xlabel('Time (s)')
    ylabel('I synaptic current (pA)')
    xlim(x_lim)
    savefig(output_fname + '_Isyn.pdf')
    
    
    # High pass filter these signals
    figure()
    ax = subplot(211)
    plot(stateMon_Iclamp_e.times, butterHighPass(stateMon_Iclamp_e.values[0].T/pA, options.sim_dt, 40))
    #plot(stateMon_Iclamp_e.times, stateMon_Iext_e.values[0]/pA)
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
    

def printFiringRatesBar(Favg_e, Favg_i, mean_e, mean_i, output_fname):
    f = figure()
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
    
