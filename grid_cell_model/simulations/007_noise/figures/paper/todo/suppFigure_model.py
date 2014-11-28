#!/usr/bin/env python
'''
Supplementary figure: model description and output.
'''
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter

from fig_conn_func        import plotWeights
from grid_cell_model.data_storage         import DataStorage
from grid_cell_model.data_storage.sim_models.ei import extractSummedSignals
from EI_plotting.base     import plotStateSignal, plotThetaSignal, \
        getOption, thetaLim
from grid_cell_model.plotting.grids       import plotGridRateMap, plotAutoCorrelation, plotSpikes2D
from grid_cell_model.plotting.global_defs import globalAxesSettings
from grid_cell_model.plotting.low_level   import xScaleBar
from grid_cell_model.analysis.visitors import AutoCorrelationVisitor
from grid_cell_model.parameters           import DictDataSet
from grid_cell_model.submitting import flagparse

parser = flagparse.FlagParser()
parser.add_flag('--examples')
parser.add_flag('--grids')
parser.add_flag('--gamma')
args = parser.parse_args()

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')
plt.rcParams['font.size'] = 11

outputDir = "panels"

jobNum = 573
dataRootDir = 'output_local'
root0   = "{0}/single_neuron".format(dataRootDir)
root150 = "{0}/single_neuron".format(dataRootDir)
root300 = "{0}/single_neuron".format(dataRootDir)
gridRootDir = '{0}/grids'.format(dataRootDir)
fileTemplate = "noise_sigma{0}_output.h5"


##############################################################################

def openJob(rootDir, noise_sigma):
    fileName = rootDir + '/' + fileTemplate.format(noise_sigma)
    return DataStorage.open(fileName, 'r')

def openGridJob(rootDir, noise_sigma, jobNum):
    fileName = rootDir + '/' + \
            'EI_param_sweep_{0}pA/job{1:05}_output.h5'.format(noise_sigma, jobNum)
    return DataStorage.open(fileName, 'a')


def plotHistogram(ax, sig, color='black', labelx="", labely="",
        labelyPos=-0.5, powerLimits=(0, 3)):
    hist(sig, bins=100, normed=True, histtype='step', align='mid', color=color)

    # y label manually
    if (labely is None):
        labely = 'Count'
    ax.text(labelyPos, 0.5, labely,
        verticalalignment='center', horizontalalignment='right',
        transform=ax.transAxes,
        rotation=90)
    xlabel(labelx)
    
    ax.minorticks_on()
    ax.xaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.set_major_locator(LinearLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    f = ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits(powerLimits)
    ax.yaxis.set_major_formatter(f)
    ax.tick_params(
            which='both',
            direction='out'
    )
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def plotGamma(gs, data, gsRow, gsCol, plotTStart, plotTEnd, yLabelOn=True,
        scaleBar=None):
    if (yLabelOn):
        IsynText = "I (nA)"
    else:
        IsynText = ""

    d = data['trials'][0]
    mon_e = d['stateMon_e']

    # E cell Isyn
    labelYPos = -0.175
    ax0 = subplot(gs[gsRow,gsCol]) 
    t, IsynMiddle = extractSummedSignals(mon_e, ['I_clamp_GABA_A'],
            plotTStart, plotTEnd)
    plotStateSignal(ax0, t, IsynMiddle*1e-3, labely=IsynText,
            labelyPos=labelYPos, color='red', scaleBar=scaleBar, scaleX=0.85,
            scaleY=-0.15, scaleText=None)

    # Autocorrelation of the 10s signal sample
    acStateList = ['I_clamp_GABA_A']
    acTEnd = 10e3 # ms
    v = AutoCorrelationVisitor('stateMon_e', acStateList, tEnd=acTEnd,
            forceUpdate=False)
    v.visitDictDataSet(DictDataSet(d))
    a = d['analysis']
    acVec = a['acVec'][0, :]
    freq = a['freq'][0]
    freq_T = 1. /  freq * 1e3
    acVal  = a['acVal'][0]
    dt    = a['ac_dt']
    times = np.arange(len(acVec))*dt
    ax1 = subplot(gs[gsRow + 1,gsCol]) 
    globalAxesSettings(ax1)
    ax1.plot(times, acVec, color='black')
    ax1.set_xlim([0, times[-1]])
    ax1.set_ylim([-1, 1])
    ax1.xaxis.set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.set_major_locator(LinearLocator(3))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.text(labelYPos, 0.5, 'Correlation',
        verticalalignment='center', horizontalalignment='right',
        transform=ax1.transAxes,
        rotation=90)
    xScaleBar(scaleBar, x=0.75, y=-0.15, ax=ax1, size='small',
            unitsText='ms')

    # Frequency annotation
    ann_x = freq_T
    ann_y = acVal
    txt = "Frequency: {0:.1f} Hz ({1:.1f} ms)\nCorrelation: {2:.1f}"
    txt = txt.format(freq, freq_T, acVal)
    ax1.annotate(txt,
        xy=(ann_x, ann_y ), xycoords='data',
        xytext=(0.2, 1.2), textcoords='axes fraction',
        va='center', ha='left',
        arrowprops=dict(arrowstyle="-|>",
                        connectionstyle="arc3",
                        relpos=(0, 0.5)),
        )



def drawSignals(gs, data, colStart, noise_sigma, yLabelOn=True, letter='',
        letterPos=None, scaleBar=None):
    if (yLabelOn):
        VmText = "V (mV)"
        IsynText = "I (nA)"
        countText = None
    else:
        VmText = ""
        IsynText = ""
        countText = ""
    histLabelX = "V (mV)"

    ncols = 4
    plotTStart = 5e3
    plotTEnd   = 5.25e3

    stateYlim = [-80, -40]

    theta_start_t = getOption(data, 'theta_start_t')
    #theta_start_t = 1e3
    simTime = getOption(data, 'time')

    mon_e = data['stateMon_e']
    mon_i = data['stateMon_i']

    ax0 = subplot(gs[0, colStart:colStart+ncols])
    t, IStim_e = extractSummedSignals(mon_e, ['I_stim'], plotTStart,
            plotTEnd)
    plotThetaSignal(ax0, t, IStim_e, noise_sigma, yLabelOn, thetaLim,
            color='red')
    t, IStim_i = extractSummedSignals(mon_i, ['I_stim'], plotTStart,
            plotTEnd)
    plotThetaSignal(ax0, t, IStim_i, noise_sigma, yLabelOn, thetaLim,
            color='blue')

    # E cell Vm
    ax1 = subplot(gs[1, colStart:colStart+ncols])
    t, VmMiddle = extractSummedSignals(mon_e, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax1, t, VmMiddle, labely=VmText, color='red')
    ylim(stateYlim)

    # I cell Vm
    ax2 = subplot(gs[2, colStart:colStart+ncols])
    t, VmMiddle = extractSummedSignals(mon_i, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax2, t, VmMiddle, labely=VmText, color='blue',
            scaleBar=scaleBar)
    ylim(stateYlim)

    ## E cell Vm histogram
    #ax3 = subplot(gs[2, colStart:colStart+2])
    #t, VmMiddle = extractSummedSignals(mon_e, ['V_m'], theta_start_t,
    #        simTime)
    #plotHistogram(ax3, VmMiddle, labelx = histLabelX, labely=countText, color='red')

    ## I cell Vm histogram
    #ax4 = subplot(gs[2, colStart+2:colStart+4])
    #t, VmMiddle = extractSummedSignals(mon_i, ['V_m'], theta_start_t,
    #        simTime)
    #plotHistogram(ax4, VmMiddle, labelx = histLabelX, labely="", color='blue')


    #if (yLabelOn):
    #    ax1.legend(['E cell', 'I cell'], fontsize='small', frameon=False,
    #            loc=[0.0, 1.1], ncol=2)


def plotGrids(gs, data, gsRowStart=0, gsCol=0):
    a = data['trials'][0]['analysis']
    arenaDiam = data['trials'][0]['options']['arenaSize']
    rateMap = a['rateMap_e']

    # Spikes
    ax0 = subplot(gs[gsRowStart, gsCol])
    plotSpikes2D(a['spikes_e'], a['rat_pos_x'], a['rat_pos_y'], a['rat_dt'],
            scaleBar=50, scaleText=False, spikeDotSize=2)

    # Grid field
    ax1 = subplot(gs[gsRowStart+1, gsCol])
    X = a['rateMap_e_X']
    Y = a['rateMap_e_Y']
    plotGridRateMap(rateMap, X, Y, diam=arenaDiam, scaleBar=50,
            scaleText=False, rasterized=True, maxRate=False)

    # Grid field autocorrelation
    ax2 = subplot(gs[gsRowStart+2, gsCol])
    X = a['corr_X']
    Y = a['corr_Y']
    ac = a['corr']
    plotAutoCorrelation(ac, X, Y, diam=arenaDiam, scaleBar=50, rasterized=True)






hr = 0.75
vh = 1.  # Vm height
th = 0.75 # top plot height
height_ratios = [th, vh, vh]

top = 0.87
bottom = 0.075
margin = 0.075
div = 0.06
width = 0.26
hspace = 0.3
wspace = 1.2

letter_top=0.95
letter_div = 0.02
letter_left=0.01
letter_va='bottom'
letter_ha='left'

gs_rows = 3
gs_cols = 4


if args.examples or args.all:
    fig_examples = figure(figsize=(9.6, 2.6))

    # noise_sigm = 0 pA
    left = margin
    right = left + width
    ds = openJob(root0, noise_sigma=0)
    gs = GridSpec(gs_rows, gs_cols, height_ratios=height_ratios, hspace=hspace,
            wspace=wspace)
    # do not update left and right
    gs.update(left=left, right=right, bottom=bottom, top=top)
    drawSignals(gs, ds, colStart=0, noise_sigma=0)
    #fig.text(left, top+letter_div, "B", va=letter_va, ha=letter_ha,
    #        fontsize=19, fontweight='bold')


    # noise_sigma = 150 pA
    ds = openJob(root150, noise_sigma=150)
    gs = GridSpec(gs_rows, gs_cols, height_ratios=height_ratios, hspace=hspace,
            wspace=wspace)
    left = right + div
    right = left + width
    gs.update(left=left, right=right, bottom=bottom, top=top)
    drawSignals(gs, ds, colStart=0, yLabelOn=False, noise_sigma=150, letterPos=-0.2)



    # noise_sigma = 300 pA
    ds = openJob(root300, noise_sigma=300)
    gs = GridSpec(gs_rows, gs_cols, height_ratios=height_ratios, hspace=hspace,
            wspace=wspace)
    left = right + div
    right = left + width
    gs.update(left=left, right=right, bottom=bottom, top=top)
    drawSignals(gs, ds, colStart=0, yLabelOn=False, noise_sigma=300,
            letterPos=-0.2, scaleBar=50)

    fname = outputDir + "/suppFigure_model_examples.pdf"
    savefig(fname, dpi=300)


grids_ds = openGridJob(gridRootDir, noise_sigma=150, jobNum=340)
if args.grids or args.all:
    ## Grid fields and gamma
    fig_grids = figure(figsize=(1.3, 2.6))
    gs = GridSpec(3, 1)
    gs.update(left=0, right=1, bottom=0.1, top=1, hspace=0)
    plotGrids(gs, grids_ds) 
    #fig.text(g_left-0.2*div, letter_top, "B", va=letter_va, ha=letter_ha, fontsize=19,
    #        fontweight='bold')
    fname = outputDir + "/suppFigure_model_grids.pdf"
    savefig(fname, dpi=300, transparent=True)

if args.gamma or args.all:
    fig_grids = figure(figsize=(4, 2.6))
    gs = GridSpec(2, 1, height_ratios=[1, 0.8])
    gs.update(left=0, right=1, bottom=0.1, top=1)
    plotGamma(gs, grids_ds, 0, 0, plotTStart=582e3, plotTEnd=582.25e3,
            scaleBar=25)
    gs.tight_layout(fig_grids, rect=(0.05, 0.05, 1, 1), h_pad=3)
    fname = outputDir + "/suppFigure_model_gamma.pdf"
    savefig(fname, dpi=300, transparent=True)
