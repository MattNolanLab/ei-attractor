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
import numpy       as np
import matplotlib.transforms as transforms

from matplotlib.pyplot   import gca, axis, colorbar
from matplotlib.ticker   import MaxNLocator, AutoMinorLocator, LinearLocator
from matplotlib.colorbar import make_axes

from parameters.param_space import JobTrialSpace2D
from plotting.global_defs import globalAxesSettings
from plotting.low_level   import xScaleBar



theta_T = 250.0 # ms
theta_dt = 1.0 # ms
theta_f = 8.0  # Hz
theta_DC = 300.0 # pA
theta_A  = 375.0 # pA
thetaMax = 1500
thetaMin = -500
thetaLim = [thetaMin, thetaMax]



def getNoiseRootDir(prefix, noise_sigma):
    return  "{0}/{1}pA".format(prefix, int(noise_sigma))



def getNoiseRoots(prefix, noise_sigmas):
    roots = []
    for s in noise_sigmas:
        roots.append(getNoiseRootDir(prefix, s))
    return roots


def getNoiseDataSpaces(dataRoot, noise_sigmas, shape):
    if (dataRoot is None):
        return None
    roots = getNoiseRoots(dataRoot, noise_sigmas)
    ds = []
    for root in roots:
        ds.append(JobTrialSpace2D(shape, root))
    return ds


class NoiseDataSpaces(object):
    class Roots(object):
        def __init__(self, *args):
            la = len(args)
            if (la == 1):
                self.bump  = args[0]
                self.v     = args[0]
                self.grids = args[0]
            elif (la == 3):
                self.bump  = args[0]
                self.v     = args[1]
                self.grids = args[2]
            else:
                raise IndexError("Roots class constructor needs either 1 or"+\
                        " three arguments")

    '''
    A container for the data spaces for all levels of noise.
    '''
    def __init__(self, roots, shape, noise_sigmas):
        self.bumpGamma    = getNoiseDataSpaces(roots.bump,  noise_sigmas,
                shape)
        self.v            = getNoiseDataSpaces(roots.v,     noise_sigmas,
                shape)
        self.grids        = getNoiseDataSpaces(roots.grids, noise_sigmas,
                shape)
        self.noise_sigmas = noise_sigmas




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

def plotStateSignal(ax, t, sig, leftSpineOn=True, labely="", labelyPos=-0.2,
        color='black', scaleBar=None, scaleX=0.75, scaleY=0.05, scaleText='ms'):
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
        xScaleBar(scaleBar, x=scaleX, y=scaleY, ax=ax, size='small',
                unitsText=scaleText)



def plotThetaSignal(ax, t, theta, noise_sigma, yLabelOn, thetaLim, color='grey'):
    setSignalAxes(ax, leftSpineOn=False)
    ax.plot(t, theta, color=color)
    ax.set_ylim(thetaLim)
    txt = '$\sigma$ = ' + str(noise_sigma) + ' pA'
    ax.text(0.5, 1.1, txt,
            verticalalignment='bottom', horizontalalignment='center',
            transform=ax.transAxes,
            fontsize='large', fontweight='normal')
    ax.axhline(0.0, color='grey', linestyle=':', linewidth=0.5)
    if (yLabelOn):
        ax.text(t[-1] - 10, -50, "0 pA", ha="right", va='top', fontsize='small')
    

def plotBump(ax, rateMap, cmap='jet', maxRate=True, **kw):
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    fs = kw.pop('fontsize', 'small')
    rx = kw.pop('rateXPos', 0.95)
    ry = kw.pop('rateYPos', 1.025)
    ax.pcolormesh(rateMap, cmap=cmap, **kw)
    axis("scaled")
    axis('off')
    if (maxRate):
        rStr = '{0:.1f} Hz'.format(np.max(rateMap.flatten()))
        ax.text(rx, ry, rStr, ha="right", va='bottom', fontsize=fs,
                transform=ax.transAxes)


def plotOneHist(data, bins=40, normed=False, **kw):
    ax = gca()
    globalAxesSettings(ax)
    ax.hist(data, bins=bins, normed=normed, histtype='step',
            align='mid', **kw)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))

###############################################################################
# Parameter sweep plotting
def createColorbar(ax, **kwargs):
    cbLabel     = kwargs.pop('label', '')
    orientation = kwargs.get('orientation', 'horizontal')
    rasterized  = kwargs.pop('rasterized', None)

    cax, kwargs = make_axes(ax, **kwargs)
    globalAxesSettings(cax)
    cb = colorbar(ax=ax, cax=cax, **kwargs)
    cb.set_label(cbLabel)
    cb.solids.set_rasterized(rasterized)

    if ('orientation' == 'horizontal'):
        cax.xaxis.set_minor_locator(AutoMinorLocator(2))
    else:
        cax.yaxis.set_minor_locator(AutoMinorLocator(2))
    return cax
