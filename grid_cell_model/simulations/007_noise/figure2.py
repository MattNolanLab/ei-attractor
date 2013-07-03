#!/usr/bin/env python
#
#   figure1.py
#
#   Noise publication Figure 1.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
from matplotlib.pyplot   import figure, subplot, savefig
from matplotlib.gridspec import GridSpec

from data_storage          import DataStorage
from figures_shared        import plotThetaSignal, thetaLim, extractStateVars,\
        plotStateSignal, plotEIRaster

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
rc('pdf', fonttype=42)

outputDir = "."
figSize = (13, 8)


trialNum = 0
jobNum = 573
dataRootDir = 'output_local'
root0   = "{0}/2013-04-24T15-27-30_EI_param_sweep_0pA_big".format(dataRootDir)
root150 = "{0}/2013-04-24T21-37-47_EI_param_sweep_150pA_big".format(dataRootDir)
root300 = "{0}/2013-04-24T21-43-32_EI_param_sweep_300pA_big".format(dataRootDir)
fileTemplate = "job{0:05}_output.h5"

##############################################################################

def openJob(rootDir, jobNum, trialNum):
    fileTemplate = "job{0:05}_output.h5"
    fileName = rootDir + '/' + fileTemplate.format(jobNum)
    return DataStorage.open(fileName, 'r')['trials'][trialNum]



def plotBump(ax, rateMap):
    rateMap = np.zeros((10, 10))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.pcolormesh(rateMap)
    axis("scaled")


def drawSignals(gs, data, colStart, noise_sigma, yLabelOn=True):
    if (yLabelOn):
        VmText = "V (mV)"
        IsynText = "I (nA)"
        rasterText = "neuron #"
    else:
        VmText = ""
        IsynText = ""
        rasterText = ""

    ncols = 4
    plotTStart = 2e3
    plotTEnd   = 2.25e3

    mon_e = data['stateMon_e']
    mon_i = data['stateMon_i']

    ax0 = subplot(gs[0, colStart:colStart+ncols])
    t, IStim = extractStateVars(mon_e, ['I_stim'], plotTStart,
            plotTEnd)
    plotThetaSignal(ax0, t, IStim, noise_sigma, yLabelOn, thetaLim)

    # E cell Vm
    ax1 = subplot(gs[1, colStart:colStart+ncols])
    t, VmMiddle = extractStateVars(mon_e, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax1, t, VmMiddle, labely=VmText, color='red')

    # E cell Isyn
    ax2 = subplot(gs[2, colStart:colStart+ncols])
    t, IsynMiddle = extractStateVars(mon_e, ['I_clamp_GABA_A'],
            plotTStart, plotTEnd)
    plotStateSignal(ax2, t, IsynMiddle*1e-3, labely=IsynText, color='blue')

    # I cell Vm
    ax3 = subplot(gs[3, colStart:colStart+ncols])
    t, VmMiddle = extractStateVars(mon_i, ['V_m'], plotTStart,
            plotTEnd)
    plotStateSignal(ax3, t, VmMiddle, labely=VmText, color='red')

    # I cell Isyn
    ax4 = subplot(gs[4, colStart:colStart+ncols])
    t, IsynMiddle = extractStateVars(mon_i, ['I_clamp_AMPA',
        'I_clamp_NMDA'], plotTStart, plotTEnd)
    plotStateSignal(ax4, t, IsynMiddle*1e-3, labely=IsynText, color='blue')

    ax5 = subplot(gs[5, colStart:colStart+ncols])
    plotEIRaster(ax5, data, 'spikeMon_e', 'spikeMon_i', (plotTStart, plotTEnd),
            labely=rasterText)


    #ax5 = subplot(gs[5, colStart])
    #plotBump(ax5, None)

    #ax6 = subplot(gs[5, colStart+1])
    #plotBump(ax6, None)

    #ax7 = subplot(gs[5, colStart+2])
    #plotBump(ax7, None)

    #ax8 = subplot(gs[5, colStart+3])
    #plotBump(ax8, None)

    if (yLabelOn):
        ax1.text(-0.32, -0.15, "Excitatory layer",
                verticalalignment='center', horizontalalignment='right',
                transform=ax1.transAxes,
                rotation=90,
                fontsize=16)
        ax3.text(-0.32, -0.15, "Inhibitory layer",
                verticalalignment='center', horizontalalignment='right',
                transform=ax3.transAxes,
                rotation=90,
                fontsize=16)



fig = figure(figsize=figSize)

hr = 1
th = 0.75 # top plot height
rr = 1.5  # raster height
top = 0.93
bottom = 0.02
margin = 0.1
div = 0.085
width = 0.23
hspace = 0.3
gsHeightRatios = [th, hr, hr, hr, hr, rr]

letter_top=0.97
letter_div = 0.02
letter_left=0.01
letter_va='bottom'
letter_ha='left'

gsRows = 6
gsCols = 4

# noise_sigm = 0 pA
left = margin
right = left + width
ds = openJob(root0, jobNum, trialNum)
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
# do not update left and right
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, ds, colStart=0, noise_sigma=0)
fig.text(letter_left, top+letter_div, "A", va=letter_va, ha=letter_ha,
        fontsize=19, fontweight='bold')


# noise_sigma = 150 pA
ds = openJob(root150, jobNum, trialNum)
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, ds, colStart=0, yLabelOn=False, noise_sigma=150)
fig.text(letter_left+margin+width+0.5*div, top+letter_div, "B", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')



# noise_sigma = 300 pA
ds = openJob(root300, jobNum, trialNum)
gs = GridSpec(gsRows, gsCols, height_ratios=gsHeightRatios, hspace=hspace)
left = right + div
right = left + width
gs.update(left=left, right=right, bottom=bottom, top=top)
drawSignals(gs, ds, colStart=0, yLabelOn=False, noise_sigma=300)
fig.text(letter_left+margin+2*width+1.5*div, top+letter_div, "C", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')

fname = outputDir + "/figure2.png"
savefig(fname)

