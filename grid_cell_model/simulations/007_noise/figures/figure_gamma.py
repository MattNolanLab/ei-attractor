#!/usr/bin/env python
#
#   figure2.py
#
#   Theta/gamma analysis using a custom "peak" method - E/I coupling parameter
#   sweep.
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
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from copy import deepcopy


from EI_plotting          import sweeps, examples, details, scatter
from EI_plotting          import aggregate as aggr
from parameters           import JobTrialSpace2D, DataSpace
from plotting.global_defs import globalAxesSettings, prepareLims
from EI_plotting.base     import plotOneHist, NoiseDataSpaces

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

###############################################################################
cFreq = 'blue'
cAC = 'green'
cCount = 'red'

outputDir = "panels"
NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
exampleIdx    = [(0, 0), (0, 0), (0, 0)] # (row, col)
bumpDataRoot  = 'output_local/even_spacing/gamma_bump'
velDataRoot   = None
gridsDataRoot = 'output_local/even_spacing/grids'
shape    = (31, 31)

gammaSweep      = 0
threshold       = 0
freqHist        = 0
detailed_noise  = 0
examplesFlag    = 0
scatterPlot     = 1
scatterPlot_all = 0

###############################################################################

def extractACExample(sp, r, c, trialNum):
    data = sp[r][c][trialNum].data
    ac = data['analysis']['acVec'][0]
    dt = data['stateMonF_e'][0]['interval']
    freq = data['analysis']['freq'][0]
    acVal = data['analysis']['acVal'][0]
    noise_sigma = data['options']['noise_sigma']
    return ac, dt, freq, acVal, noise_sigma


def aggregateBar2(spList, varLists, trialNumList, func=(None, None)):
    vars = ([], [])
    noise_sigma = []
    for idx in xrange(len(spList)):
        for varIdx in range(len(varLists)):
            f = func[varIdx]
            if f is None:
                f = lambda x: x
            vars[varIdx].append(f(aggr.aggregate2DTrial(spList[idx], varLists[varIdx],
                trialNumList).flatten()))
        noise_sigma.append(spList[idx][0][0][0].data['options']['noise_sigma'])

    noise_sigma = np.array(noise_sigma, dtype=int)
    return vars, noise_sigma
 


def getACFreqThreshold(spList, trialNumList, ACThr):
    varLists = [['acVal'], ['freq']]
    vars, noise_sigma = aggregateBar2(spList, varLists, trialNumList)
    AC = vars[0]
    freq = vars[1]
    ACMean   = []
    freqMean = []
    ACStd    = []
    freqStd  = []
    thrCount = []
    for spIdx in range(len(spList)):
        thrIdx = np.logical_and(AC[spIdx] >= ACThr,
                np.logical_not(np.isnan(AC[spIdx])))
        ac_filt = AC[spIdx][thrIdx]
        ACMean.append(np.mean(ac_filt))
        ACStd.append(np.std(ac_filt))
        freq_filt = freq[spIdx][thrIdx]
        freqMean.append(np.mean(freq_filt))
        freqStd.append(np.std(freq_filt))
        thrCount.append(float(len(AC[spIdx][thrIdx])) / len (AC[spIdx]))

    return (ACMean, ACStd), (freqMean, freqStd), thrCount, noise_sigma


def plotThresholdComparison(spList, trialNumList, ACThrList):
    counts = []
    noise_sigma = None
    for ACThr in ACThrList:
        _, _, thrCount, noise_sigma = getACFreqThreshold(spList, trialNumList,
                ACThr)
        counts.append(thrCount)
    counts = np.array(counts)

    print ACThrList, counts

    ax = plt.gca()
    globalAxesSettings(ax)
    plt.plot(ACThrList, counts, 'o-')
    plt.plot([0], [1], linestyle='None', marker='None')
    ax.set_xlabel('Correlation threshold', labelpad=5)
    ax.set_ylabel('Count', labelpad=5)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='small')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.3))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.margins(0.025)
    

def plotFreqHistogram(spList, trialNumList, ylabelPos=-0.2, CThreshold=0.1):
    FVarList = ['freq']
    CVarList = ['acVal']
    noise_sigma = [0, 150, 300]
    colors = ['red', 'green', 'blue']

    ax = plt.gca()
    plt.hold('on')
    globalAxesSettings(ax)

    for idx, sp in enumerate(spList):
        F = aggr.aggregate2DTrial(sp, FVarList, trialNumList).flatten()
        C = aggr.aggregate2DTrial(sp, CVarList, trialNumList).flatten()
        filtIdx = np.logical_and(np.logical_not(np.isnan(F)), C > CThreshold)
        plotOneHist(F[filtIdx], bins=20, normed=True)
    leg = []
    for s in noise_sigma:
        leg.append("{0}".format(int(s)))
    l = ax.legend(leg, loc=(0.8, 0.5), title='$\sigma$ (pA)', frameon=False,
            fontsize='x-small', ncol=1)
    plt.setp(l.get_title(), fontsize='x-small')

    ax.set_xlabel("Oscillation frequency (Hz)")
    #ax.text(ylabelPos, 0.5, 'p(F)', rotation=90, transform=ax.transAxes,
    #        va='center', ha='right')
    ax.set_ylabel('p(Frequency)')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(ti.MultipleLocator(20))
    ax.yaxis.set_major_locator(ti.MaxNLocator(4))
    ax.xaxis.set_minor_locator(ti.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(2))
    f = ti.ScalarFormatter(useMathText=True)
    f.set_scientific(True)
    f.set_powerlimits([0, 3])
    ax.yaxis.set_major_formatter(f)
    #ax.margins(0.01, 0.00)

    thStr = 'Frequencies with C > {0}'.format(CThreshold)
    ax.text(0.99, 1.1, thStr, transform=ax.transAxes, va='bottom',
            ha='right')
    





###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)

# gamma example rows and columns
exampleRC = ( (5, 15), (15, 5) )


sweepFigSize = (2.8, 2)
sweepLeft   = 0.15
sweepBottom = 0.2
sweepRight  = 0.87
sweepTop    = 0.85
transparent  = True

AC_vmin = -0.09
AC_vmax = 0.675
F_vmin  = 30
F_vmax  = 120

ACVarList = ['acVal']
FVarList  = ['freq']

AC_cbar_kw = dict(
        orientation='vertical',
        ticks=ti.MultipleLocator(0.3),
        fraction=0.25,
        shrink=0.8,
        pad=0.05,
        labelpad=8,
        label='$1^{st}$ autocorrelation\npeak')
F_cbar_kw = dict(
        orientation='vertical',
        ticks=ti.MultipleLocator(30),
        fraction=0.25,
        shrink=0.8,
        pad=0.05,
        labelpad=8,
        label='Oscillation\nfrequency (Hz)',
        extend='max', extendfrac=0.1)


ann_color = 'white'
ann0 = dict(
        txt='b',
        rc=exampleRC[0],
        xytext_offset=(1.5, 0),
        color=ann_color)
ann1 = dict(
        txt='a',
        rc=exampleRC[1],
        xytext_offset=(1.5, 1),
        color=ann_color)
ann = [ann0, ann1]
annF = [deepcopy(ann0), deepcopy(ann1)]


if (gammaSweep):
    # noise_sigma = 0 pA
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotACTrial(ps.bumpGamma[0], ACVarList, iterList,
            noise_sigma=ps.noise_sigmas[0],
            ax=ax,
            xlabel='', xticks=False,
            trialNumList=xrange(NTrials),
            cbar=False, cbar_kw=AC_cbar_kw,
            vmin=AC_vmin, vmax=AC_vmax,
            annotations=ann)
    fname = outputDir + "/gamma_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotACTrial(ps.bumpGamma[0], FVarList, iterList,
            noise_sigma=ps.noise_sigmas[0],
            ax=ax,
            trialNumList=xrange(NTrials),
            sigmaTitle=False,
            cbar=False, cbar_kw=F_cbar_kw,
            vmin=F_vmin, vmax=F_vmax,
            annotations=annF)
    fname = outputDir + "/gamma_freq_sweeps0.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        

    # noise_sigma = 150 pA
    for a in ann:
        a['color'] = 'black'
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotACTrial(ps.bumpGamma[1], ACVarList, iterList,
            noise_sigma=ps.noise_sigmas[1],
            ax=ax,
            xlabel='', xticks=False,
            trialNumList=xrange(NTrials),
            ylabel='', yticks=False,
            cbar=False, cbar_kw=AC_cbar_kw,
            vmin=AC_vmin, vmax=AC_vmax,
            annotations=ann)
    fname = outputDir + "/gamma_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotACTrial(ps.bumpGamma[1], FVarList, iterList,
            noise_sigma=ps.noise_sigmas[1],
            ax=ax,
            trialNumList=xrange(NTrials),
            ylabel='', yticks=False,
            sigmaTitle=False,
            cbar=False, cbar_kw=F_cbar_kw,
            vmin=F_vmin, vmax=F_vmax,
            annotations=annF)
    fname = outputDir + "/gamma_freq_sweeps150.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        

    # noise_sigma = 300 pA
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotACTrial(ps.bumpGamma[2], ACVarList, iterList,
            noise_sigma=ps.noise_sigmas[2],
            ax=ax,
            xlabel='', xticks=False,
            trialNumList=xrange(NTrials),
            ylabel='', yticks=False,
            cbar=True, cbar_kw=AC_cbar_kw,
            vmin=AC_vmin, vmax=AC_vmax,
            annotations=ann)
    fname = outputDir + "/gamma_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        
    fig = plt.figure(figsize=sweepFigSize)
    ax = fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))
    sweeps.plotACTrial(ps.bumpGamma[2], FVarList, iterList,
            noise_sigma=ps.noise_sigmas[2],
            ax=ax,
            sigmaTitle=False,
            trialNumList=xrange(NTrials),
            ylabel='', yticks=False,
            cbar=True, cbar_kw=F_cbar_kw,
            vmin=F_vmin, vmax=F_vmax,
            annotations=annF)
    fname = outputDir + "/gamma_freq_sweeps300.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)
        


if (threshold):
    ###############################################################################
    plt.figure(figsize=(3.5, 2))
    plotThresholdComparison(ps.bumpGamma,
            trialNumList=range(NTrials),
            ACThrList=np.arange(0, 0.65, 0.05))
    plt.tight_layout()
    fname = outputDir + '/gamma_AC_threshold_comparison.pdf'
    plt.savefig(fname, transparent=True, dpi=300)


if (freqHist):
    ylabelPos = -0.16
    fig = plt.figure(figsize=(3.7, 2.5))
    plotFreqHistogram(ps.bumpGamma, range(NTrials), ylabelPos=ylabelPos)
    plt.tight_layout()
    fname = outputDir + "/gamma_freq_histograms.pdf"
    plt.savefig(fname, dpi=300, transparent=True)


##############################################################################
EI13Root  = 'output_local/detailed_noise/gamma_bump/EI-1_3'
EI31Root  = 'output_local/detailed_noise/gamma_bump/EI-3_1'
detailedShape = (31, 9)

EI13PS = JobTrialSpace2D(detailedShape, EI13Root)
EI31PS = JobTrialSpace2D(detailedShape, EI31Root)
detailedNTrials = 5

detailFigSize = (4.25, 1.8)
detailLeft   = 0.2
detailBottom = 0.3
detailRight  = 0.98
detailTop    = 0.9
if (detailed_noise):
    ylabelPos = -0.17

    # 1st autocorrelation peak (gamma power)
    types = ('gamma', 'acVal')
    fig = plt.figure(figsize=detailFigSize)
    ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
        detailTop))
    _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
            ylabelPos=ylabelPos,
            xlabel='', xticks=False,
            color='red', markerfacecolor='red', zorder=10)
    _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
            xlabel='', xticks=False,
            ylabel='$1^{st}$\nautocorrelation\npeak', ylabelPos=ylabelPos,
            color='#505050')
    ax.xaxis.set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.6))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(6))
    ax.set_ylim(prepareLims((0, 0.6), margin=0.03))
    leg = ['a', 'b']
    l = ax.legend([p31, p13], leg, loc=(0.85, 0.7), fontsize='small', frameon=False,
            numpoints=1, handletextpad=0.05)
    plt.setp(l.get_title(), fontsize='small')


    fname = outputDir + "/gamma_detailed_noise_power.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()

    # Gamma frequency
    types = ('gamma', 'freq')
    fig = plt.figure(figsize=detailFigSize)
    ax = fig.add_axes(Bbox.from_extents(detailLeft, detailBottom, detailRight,
        detailTop))
    _, p13, l13 = details.plotDetailedNoise(EI13PS, detailedNTrials, types, ax=ax,
            ylabelPos=ylabelPos,
            xlabel='',
            color='red', markerfacecolor='red', zorder=10)
    _, p31, l31 = details.plotDetailedNoise(EI31PS, detailedNTrials, types, ax=ax,
            ylabel='Oscillation\nfrequency (Hz)', ylabelPos=ylabelPos,
            color='#505050')
    ax.yaxis.set_major_locator(ti.MultipleLocator(30))
    ax.yaxis.set_minor_locator(ti.AutoMinorLocator(3))
    ax.set_ylim(prepareLims((30, 90), margin=0.03))

    fname = outputDir + "/gamma_detailed_noise_freq.pdf"
    plt.savefig(fname, dpi=300, transparent=True)
    plt.close()



##############################################################################
exampleFName = outputDir + "/gamma_example{0}_{1}.pdf"
exampleTrialNum = 0
exampleFigSize = (2, 1.1)
exampleLeft   = 0.08
exampleBottom = 0.2
exampleRight  = 0.99
exampleTop    = 0.85
example_xscale_kw = dict(
        scaleLen=50,
        x=0.75, y=-0.1,
        size='x-small')

if (examplesFlag):
    for nsIdx, ns in enumerate(ps.noise_sigmas):
        for idx, rc in enumerate(exampleRC):
            fname = exampleFName.format(ns, idx)
            fig = plt.figure(figsize=exampleFigSize)
            ax = fig.add_axes(Bbox.from_extents(exampleLeft, exampleBottom,
                exampleRight, exampleTop))
            nsAnn = None
            xscale_kw = None
            if (idx == 0):
                nsAnn = ns
                if (nsIdx == len(ps.noise_sigmas)-1):
                    xscale_kw = example_xscale_kw
            examples.plotGammaExample(ps.bumpGamma[nsIdx], ax=ax,
                    r=exampleRC[idx][0], c=exampleRC[idx][1],
                    trialNum=exampleTrialNum,
                    tStart = 2e3, tEnd=2.25e3,
                    noise_sigma=nsAnn, noise_sigma_xy=(0.95, 1),
                    xscale_kw=xscale_kw)
            plt.savefig(fname, dpi=300, transparent=True)
            plt.close()


##############################################################################
# Separate scatter plot of gridness score vs. gamma power
scatterFigSize = (3.8, 3.2)
scatterLeft   = 0.2
scatterBottom = 0.32
scatterRight  = 0.98
scatterTop    = 0.87

scatterColorFigSize = (0.75, 0.75)

ignoreNaNs = True

if (scatterPlot):
    NTrialsGamma = 5
    NTrialsGrids = 3
    typesGamma = ['gamma', 'acVal']
    typesGrids = ['grids', 'gridnessScore']

    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        fig = plt.figure(figsize=scatterFigSize)
        ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
            scatterTop))

        if (ns_idx != 0):
            ylabel = ''
        else:
            ylabel = 'Gridness score'
        if (ns_idx == 1):
            xlabel = '$1^{st}$ autocorrelation peak'
        else:
            xlabel = ''

        scatterPlot = scatter.ScatterPlot(
                ps.bumpGamma[ns_idx], ps.grids[ns_idx], typesGamma,
                typesGrids, iterList, NTrialsGamma, NTrialsGrids,
                s=15,
                linewidth=0.3,
                color2D=True,
                xlabel=xlabel,
                ylabel=ylabel,
                sigmaTitle=True,
                noise_sigma=noise_sigma,
                ignoreNaNs=ignoreNaNs)
        scatterPlot.plot()
        ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
        ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
        ax.set_ylim(prepareLims((-0.5, 1.2), margin=0.02))


        fname = outputDir + "/gamma_scatter_gamma_grids{0}.pdf"
        fig.savefig(fname.format(int(noise_sigma)), dpi=300,
                transparent=transparent)

    fig = plt.figure(figsize=scatterColorFigSize)
    ax = fig.gca()
    scatterPlot.plotColorbar(ax)
    fig.tight_layout(pad=0)
    fname = outputDir + "/gamma_scatter_gamma_grids_colorbar.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)



##############################################################################
# Scatter plot of gridness score vs. gamma power 
# All in one plot
if (scatterPlot_all):
    fig = plt.figure(figsize=scatterFigSize)
    ax = fig.add_axes(Bbox.from_extents(scatterLeft, scatterBottom, scatterRight,
        scatterTop))

    NTrialsGamma = 5
    NTrialsGrids = 3
    typesGamma = ['gamma', 'acVal']
    typesGrids = ['grids', 'gridnessScore']
    ax.hold('on')
    scatterColors = ['blue', 'green', 'red']

    for ns_idx, noise_sigma in enumerate(ps.noise_sigmas):
        color = scatterColors[ns_idx]
        scatterPlot = scatter.ScatterPlot(
                ps.bumpGamma[ns_idx], ps.grids[ns_idx], typesGamma,
                typesGrids, iterList, NTrialsGamma, NTrialsGrids,
                c=color,
                s=15,
                linewidth=0.3,
                xlabel='$1^{st}$ autocorrelation peak',
                ylabel='Gridness score')
        scatterPlot.plot()
    ax.xaxis.set_major_locator(ti.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(ti.MultipleLocator(0.5))
    leg = ['0', '150', '300']
    l = ax.legend(leg, loc=(0.9, 0.4), fontsize='small', frameon=False,
            numpoints=1, title='$\sigma$ (pA)')
    plt.setp(l.get_title(), size='small')

    fname = outputDir + "/gamma_scatter_gamma_grids_all.pdf"
    fig.savefig(fname, dpi=300, transparent=transparent)

