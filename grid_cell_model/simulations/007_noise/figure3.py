#!/usr/bin/env python
#
#   analysis_EI.py
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
import matplotlib
#matplotlib.use('cairo')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator, LinearLocator, MaxNLocator

import numpy.ma as ma

from parameters  import JobTrialSpace2D, DataSpace
from EI_plotting import plot2DTrial
from plotting.global_defs import globalAxesSettings
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
rc('pdf', fonttype=42)

plt.rcParams['font.size'] = 25

###############################################################################

def aggregate2DTrial(sp, varList, trialNumList):
    varList = ['analysis'] + varList
    retVar = sp.aggregateData(varList, trialNumList, funReduce=np.mean,
            saveData=True)
    return np.mean(retVar, 2)



def computeYX(sp, iterList):
    E, I = sp.getIteratedParameters(iterList)
    Ne = DataSpace.getNetParam(sp[0][0][0].data, 'net_Ne')
    Ni = DataSpace.getNetParam(sp[0][0][0].data, 'net_Ni')
    return E/Ne, I/Ni


def plotACTrial(sp, varList, iterList, trialNumList=[0], xlabel="", ylabel="",
        colorBar=True, clBarLabel="", vmin=None, vmax=None, title="", clbarNTicks=2,
        xticks=True, yticks=True):
    C = aggregate2DTrial(sp, varList, trialNumList)
    C = ma.MaskedArray(C, mask=np.isnan(C))
    Y, X = computeYX(sp, iterList)
    plot2DTrial(X, Y, C, xlabel, ylabel, colorBar, clBarLabel, vmin, vmax,
            title, clbarNTicks, xticks, yticks)

 
def drawColorbar(drawAx, label):
    pos = drawAx.get_position()
    pos.y0 -= 0.12
    pos.y1 -= 0.12
    pos.y1 = pos.y0 + 0.1*(pos.y1 - pos.y0)
    w = pos.x1 - pos.x0
    pos.x0 += 0.1*w
    pos.x1 -= 0.1*w
    clba = plt.gcf().add_axes(pos)
    globalAxesSettings(clba)
    clba.minorticks_on()
    cb = plt.colorbar(cax=clba, orientation='horizontal',
            ticks=LinearLocator(2))
    cb.set_label(label)


def plot2DNoiseACFreq(spList, iterList, trialNumList=[0]):
    NSp = len(spList)
    Sz2D = 1.
    clBarSz = 0.1
    hr = [Sz2D]*NSp + [clBarSz]
    gs = GridSpec(len(spList), 2, wspace=0.15, hspace=0.15, height_ratios=hr,
            bottom=0.22)
    for spIdx in range(NSp):
        if (spIdx == NSp - 1):
            xlabel = "I (nS)"
            xticks = True
            clbar = True
        else:
            xlabel = ""
            xticks = False
            clbar = False

        ax_AC = plt.subplot(gs[spIdx, 0])
        acVal = plotACTrial(spList[spIdx], ['acVal'], iterList,
                trialNumList=trialNumList,
                xlabel=xlabel,
                ylabel='E (nS)',
                colorBar=False,
                clBarLabel = "Correlation",
                clbarNTicks=3,
                xticks=xticks,
                vmin=0,
                vmax=0.7)
        if (clbar):
            drawColorbar(ax_AC, 'Correlation')

        ax_Freq = plt.subplot(gs[spIdx, 1])
        freq = plotACTrial(spList[spIdx], ['freq'], iterList,
                trialNumList=range(NTrials),
                xlabel=xlabel,
                colorBar=False,
                clBarLabel = "Frequency (Hz)",
                clbarNTicks=3,
                xticks=xticks,
                yticks=False,
                vmin=0,
                vmax=150)
        if (clbar):
            drawColorbar(ax_Freq, 'Frequency (Hz)')




def extractACExample(sp, r, c, trialNum):
    data = sp[r][c][trialNum].data
    ac = data['analysis']['acVec'][0]
    dt = data['stateMonF_e'][0]['interval']
    freq = data['analysis']['freq'][0]
    acVal = data['analysis']['acVal'][0]
    noise_sigma = data['options']['noise_sigma']
    return ac, dt, freq, acVal, noise_sigma

def plotACExamples(spList, r, c, trialNum=0):
    gs = GridSpec(len(spList), 1, hspace=0.2)
    #gs.update(left=left, right=right, bottom=bottom, top=top)
    plt.hold('on')
    idx = 0
    for sp in spList:
        ax = plt.subplot(gs[idx, :])
        globalAxesSettings(ax)
        ax.xaxis.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_major_locator(LinearLocator(2))
        ax.yaxis.set_major_locator(LinearLocator(3))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ac, dt, freq, acVal, noise_sigma = extractACExample(sp, r, c, trialNum)
        t = np.arange(len(ac))*dt
        plt.plot(t, ac)
        ax.set_ylim([-1, 1])
        ax.set_xlim([0, 125])
        if (idx == 1):
            ax.set_ylabel('Correlation')
        ax.text(0.9, 1.0, '$\sigma$ = ' + str(noise_sigma) + ' pA',
                verticalalignment='center', horizontalalignment='right',
                transform=ax.transAxes, fontsize='small')
        ann_x = 1./freq*1e3
        ann_y = acVal
        ax.annotate("",
            xy=(ann_x, ann_y ), xycoords='data',
            xytext=(ann_x, ann_y+0.9), textcoords='data',
            arrowprops=dict(arrowstyle="-|>",
                            connectionstyle="arc3"),
            )
        idx += 1
    ax.xaxis.set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.set_xlabel('Lag (ms)')
    gs.tight_layout(plt.gcf(), h_pad=0.2)


def getStatData(sp, r, c, trialNum):
    data = sp[r][c][trialNum].data
    freq = data['analysis']['freq']
    acVal = data['analysis']['acVal']
    noise_sigma = data['options']['noise_sigma']
    return freq, acVal, noise_sigma

def plotFreqACStat(spList, r, c, trialNum=0):
    ax_freq = plt.gca()
    ax_ac   = ax_freq.twinx()
    globalAxesSettings(ax_freq)
    globalAxesSettings(ax_ac, setTickPos=False)
    freqAvg = []
    freqStd = []
    acAvg   = []
    acStd   = []
    nsa     = []
    for sp in spList:
        freq, acVal, noise_sigma = getStatData(sp, r, c, trialNum)
        freqAvg.append(np.mean(freq))
        freqStd.append(np.std(freq))
        acAvg.append(np.mean(acVal))
        acStd.append(np.std(acVal))
        nsa.append(noise_sigma)

    cFreq = 'blue'
    cAC = 'green'
    pFreq = ax_freq.errorbar(nsa, freqAvg, yerr=freqStd, fmt='-o', color=cFreq)
    ax_ac.errorbar(nsa, acAvg, yerr=acStd, fmt='-o', color=cAC)

    ax_freq.set_ylabel('Frequency (Hz)')
    ax_freq.set_xlabel('Noise std.dev. (pA)')
    ax_freq.set_xticks(nsa)
    noise_range = nsa[-1] - nsa[0]
    ax_freq.set_xlim([nsa[0]-0.1*noise_range, nsa[-1]+0.1*noise_range])
    ax_ac.set_ylabel('Correlation')
    ax_freq.yaxis.get_label().set_color(cFreq)
    ax_ac.yaxis.get_label().set_color(cAC)

    ax_ac.set_ylim([0, 0.7])
    ax_ac.yaxis.set_major_locator(LinearLocator(2))
    ax_ac.yaxis.set_minor_locator(AutoMinorLocator(6))


    plt.gcf().tight_layout()
        
            

###############################################################################

NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']
r = 14
c = 13

baseDir = 'output_local/one_to_one'
dirs = [ \
    ('EI_param_sweep_0pA',    (40, 40)),
    ('EI_param_sweep_150pA',  (40, 40)),
    ('EI_param_sweep_300pA',  (40, 40))
]

dataSpaces = []
for (dir, shape) in dirs:
    rootDir = "{0}/{1}".format(baseDir, dir)
    dataSpaces.append(JobTrialSpace2D(shape, rootDir))
    

#################################################################################
plt.figure(figsize=(7, 5))
plotACExamples(dataSpaces, r, c)
plt.savefig('{0}/analysis_AC_Examples.pdf'.format(baseDir))
#################################################################################
plt.figure(figsize=(7, 5))
plotFreqACStat(dataSpaces, r, c)
plt.savefig('{0}/analysis_Freq_AC_stat.pdf'.format(baseDir))


#################################################################################
plt.figure(figsize=(6.5, 8))
plot2DNoiseACFreq(dataSpaces, iterList)
#plt.tight_layout()
plt.savefig('{0}/aggregated_AC_Freq.png'.format(baseDir), transparent=True,
        dpi=300)
