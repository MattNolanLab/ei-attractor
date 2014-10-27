#!/usr/bin/env python
#
#   figure2.py
#
#   Noise publication Figure 2.
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
import matplotlib.pyplot as plt
from matplotlib.pyplot   import figure, subplot, savefig
from matplotlib.gridspec import GridSpec

from grid_cell_model.data_storage    import DataStorage
from grid_cell_model.data_storage.sim_models.ei import extractSummedSignals
from EI_plotting.base import plotThetaSignal, thetaLim, plotStateSignal


outputDir = "."
figSize = (12, 8)


trialNum = 0
jobNum = 573
dataRootDir = 'output_local/one_to_one'
root0   = "{0}/EI_param_sweep_0pA".format(dataRootDir)
root150 = "{0}/EI_param_sweep_150pA".format(dataRootDir)
root300 = "{0}/EI_param_sweep_300pA".format(dataRootDir)
fileTemplate = "job{0:05}_output.h5"

##############################################################################

def openJob(rootDir, jobNum, trialNum):
    fileTemplate = "job{0:05}_output.h5"
    fileName = rootDir + '/' + fileTemplate.format(jobNum)
    return DataStorage.open(fileName, 'r')['trials'][trialNum]



def drawSignals(gs, data, colStart, noise_sigma, yLabelOn=True, scaleBar=None):
    if (yLabelOn):
        VmText = "V (mV)"
        IsynText = "I (nA)"
        rasterText = "neuron #"
    else:
        VmText = ""
        IsynText = ""
        rasterText = ""

    ncols = 4
    plotTStart = 1e3
    plotTEnd   = 1.25e3

    mon_e = data['stateMon_e']
    mon_i = data['stateMon_i']

    ax0 = subplot(gs[0, colStart:colStart+ncols])
    t, IStim = extractSummedSignals(mon_e, ['I_stim'], plotTStart, plotTEnd)
    plotThetaSignal(ax0, t, IStim, noise_sigma, yLabelOn, thetaLim)

    # E cell Vm
    ax1 = subplot(gs[1, colStart:colStart+ncols])
    t, VmMiddle = extractSummedSignals(mon_e, ['V_m'], plotTStart, plotTEnd)
    plotStateSignal(ax1, t, VmMiddle, labely=VmText, color='red',
            scaleBar=scaleBar)

    # E cell Isyn
    ax2 = subplot(gs[2, colStart:colStart+ncols])
    t, IsynMiddle = extractSummedSignals(mon_e, ['I_clamp_GABA_A'], plotTStart,
            plotTEnd)
    plotStateSignal(ax2, t, IsynMiddle*1e-3, labely=IsynText, color='red')

    # I cell Vm
    ax3 = subplot(gs[3, colStart:colStart+ncols])
    t, VmMiddle = extractSummedSignals(mon_i, ['V_m'], plotTStart, plotTEnd)
    plotStateSignal(ax3, t, VmMiddle, labely=VmText, color='blue')

    # I cell Isyn
    ax4 = subplot(gs[4, colStart:colStart+ncols])
    t, IsynMiddle = extractSummedSignals(mon_i, ['I_clamp_AMPA',
        'I_clamp_NMDA'], plotTStart, plotTEnd)
    plotStateSignal(ax4, t, IsynMiddle*1e-3, labely=IsynText, color='blue')


    # TODO: use EI_plotting.plotEIRaster
    #ax5 = subplot(gs[5, colStart:colStart+ncols])
    #plotEIRaster(ax5, data, 'spikeMon_e', 'spikeMon_i', (plotTStart, plotTEnd),
    #        labely=rasterText)


    if (yLabelOn):
        labelPos = -0.32
        ax1.text(labelPos, -0.15, "Excitatory cell",
                verticalalignment='center', horizontalalignment='right',
                transform=ax1.transAxes,
                rotation=90,
                fontsize='large')
        ax3.text(labelPos, -0.15, "Inhibitory cell",
                verticalalignment='center', horizontalalignment='right',
                transform=ax3.transAxes,
                rotation=90,
                fontsize='large')



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
drawSignals(gs, ds, colStart=0, yLabelOn=False, noise_sigma=300, scaleBar=50)
fig.text(letter_left+margin+2*width+1.5*div, top+letter_div, "C", va=letter_va,
        ha=letter_ha, fontsize=19, fontweight='bold')

fname = outputDir + "/figure_gamma_examples.png"
savefig(fname, dpi=300)

