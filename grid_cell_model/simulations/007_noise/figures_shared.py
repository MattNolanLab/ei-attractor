#
#   figures_shared.py
#
#   Shared routines between figure plotting. Noise paper.
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
import numpy as np
#from matplotlib.pyplot import plot,
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, LinearLocator
import matplotlib.transforms as transforms

from analysis.visitors.interface import extractStateVariable, sumAllVariables,\
        extractSpikes
from plotting.global_defs        import globalAxesSettings


theta_T = 250.0 # ms
theta_dt = 1.0 # ms
theta_f = 8.0  # Hz
theta_DC = 300.0 # pA
theta_A  = 375.0 # pA
thetaMax = 1500
thetaMin = -500
thetaLim = [thetaMin, thetaMax]

def getOption(data, optStr):
    return data['options'][optStr]


def generateThetaSignal(noise_sigma):
    t = np.arange(0, theta_T, theta_dt) * 1e-3
    phase = np.pi
    normSig = 0.5*(1. + np.cos(2*np.pi*theta_f*t + phase))
    noise = noise_sigma * np.random.randn(len(t))
    return t, theta_DC + theta_A*normSig + noise


def setSignalAxes(ax, leftSpineOn):
    globalAxesSettings(ax)
    ax.minorticks_on()
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if (leftSpineOn == False):
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_visible(False)

    ax.yaxis.set_major_locator(MaxNLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

def plotStateSignal(ax, t, sig, leftSpineOn=True, labely="", labelyPos=-0.3,
        color='black', scaleBar=None):
    setSignalAxes(ax, leftSpineOn)

    if (sig is not None):
        ax.plot(t, sig, color=color)

    #ax.set_ylabel(labely)
    ax.text(labelyPos, 0.5, labely,
        verticalalignment='center', horizontalalignment='right',
        transform=ax.transAxes,
        rotation=90)
    ax.set_xlim([t[0], t[-1]])

    if (scaleBar is not None):
        scaleHeight = 0.05  # axis coords
        scaleBottom = 0.05 # data coords
        scaleOffset = 20
        scaleCenter = t[-1] - 0.5*scaleBar - scaleOffset
        ax.axvspan(
                xmin = t[-1] - scaleBar - scaleOffset,
                xmax = t[-1] - scaleOffset,
                ymin = scaleBottom,
                ymax = scaleBottom + scaleHeight,
                color = 'black')
        trans = transforms.blended_transform_factory(ax.transData,
                ax.transAxes)
        ax.text(scaleCenter, scaleBottom - 0.05,
                "{0} ms".format(scaleBar),
                va='top', ha='center', transform=trans)



def plotThetaSignal(ax, t, theta, noise_sigma, yLabelOn, thetaLim):
    setSignalAxes(ax, leftSpineOn=False)
    ax.plot(t, theta, color="grey")
    ax.set_ylim(thetaLim)
    txt = '$\sigma = ' + str(noise_sigma) + '\ \mathrm{pA}$'
    ax.text(0.5, 1.1, txt,
            verticalalignment='bottom', horizontalalignment='center',
            transform=ax.transAxes,
            fontsize='large', fontweight='bold')
    ax.axhline(0.0, color='grey', linestyle=':', linewidth=0.5)
    if (yLabelOn):
        ax.text(t[-1] - 10, -50, "0 pA", ha="right", va='top', fontsize='small')
    

def plotBump(ax, rateMap):
    rateMap = np.zeros((10, 10))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.pcolormesh(rateMap)
    axis("scaled")

def plotSpikes(ax, t, trajectory, spikeTimes):
    pass



def sliceSignal(t, sig, tStart, tEnd):
    idx = np.logical_and(t >= tStart, t <= tEnd)
    return t[idx], sig[idx], idx


def extractStateVars(mon, varName, plotTStart, plotTEnd):
    '''
    Extract state variables from a pair of monitors. One in the centre, the
    other one at the edge of the neural sheet.
    '''
    nIdxMiddle = 0

    t, dt = extractStateVariable(mon, nIdxMiddle, 'times')
    sig, dt = sumAllVariables(mon, nIdxMiddle, varName)
    t, sigMiddle, idx = sliceSignal(t, sig, plotTStart, plotTEnd)
    return t, sigMiddle


def sliceSpikes(senders, times, tLimits):
    idx = np.logical_and(times >= tLimits[0], times <= tLimits[1])
    return senders[idx], times[idx]

def _getSpikeTrain(data, mon, NParam, tLimits):
    '''
    Return the senders and spike times from a monitor in the data
    dictionary

    Parameters
    ----------
    data : dict
        A data dictionary containing monName
    mon : str
        The spike monitor dictionary name
    NParam : string
        Name of the parameter that specifies the total number of neurons
    output: tuple
        A tuple containing (senders, times, N). 'N' is the total number of
        neurons in the population.
    '''
    N = data['net_attr'][NParam]
    senders, times = extractSpikes(data[mon])
    senders, times = sliceSpikes(senders, times, tLimits)
    return senders, times, N


def plotEIRaster(ax, data, mon_e, mon_i, tLimits, labely='', labelyPos=-0.3):
    ESenders, ETimes, EN = _getSpikeTrain(data, mon_e, 'net_Ne', tLimits)
    ISenders, ITimes, IN = _getSpikeTrain(data, mon_i, 'net_Ni', tLimits)
    ISenders += EN

    globalAxesSettings(ax)
    ax.minorticks_on()
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ms = 1.0
    ax.plot(ETimes, ESenders+1, '.', color='red', markersize=ms)
    ax.plot(ITimes, ISenders+1, '.', color='blue', markersize=ms)

    ax.set_xlim(tLimits)
    ax.set_ylim([0, EN+IN])
    ax.invert_yaxis()
    ax.text(labelyPos, 0.5, labely,
        verticalalignment='center', horizontalalignment='right',
        transform=ax.transAxes,
        rotation=90)

