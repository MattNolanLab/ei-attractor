#
#   base.py
#
#   Basic routines for EI plotting.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
    #orientation = kwargs.get('orientation', 'horizontal')
    rasterized  = kwargs.pop('rasterized', None)
    labelpad    = kwargs.pop('labelpad', None)
    mappable    = kwargs.pop('mappable', None)

    cax, kwargs = make_axes(ax, **kwargs)
    globalAxesSettings(cax)
    fig = ax.figure
    cb = fig.colorbar(mappable, ax=ax, cax=cax, **kwargs)
    cb.set_label(cbLabel, labelpad=labelpad)
    cb.solids.set_rasterized(rasterized)

    return cax
