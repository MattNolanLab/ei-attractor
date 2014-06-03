'''Print synaptic connection profiles.'''
from __future__ import absolute_import, print_function

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from grid_cell_model.plotting.global_defs import globalAxesSettings, createColorbar

from .base import FigurePlotter

__all__ = [
    'ConnectionFunctionPlotter',
]

dx = 0.001
x0 = -0.5
x1 = 0.5

def plotWeights(ax, d, exc_profile, inh_profile, inh_const, linewidth):
    plt.hold('on')
    ax = plt.gca()
    globalAxesSettings(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ep, = plt.plot(d, exc_profile, linewidth=linewidth, color='red', label="E")
    ip, = plt.plot(d, inh_profile, linewidth=linewidth, color='blue', label="I")
    icp, = plt.plot(d, [inh_const]*len(d), ':', color='blue')
    plt.xlabel('Distance')
    plt.ylabel('G/G$_\mathrm{max}$')
    ax.yaxis.set_ticks([0, 1])
    ax.xaxis.set_ticks([x0, 0, x1])
    leg1 = ['E$\\rightarrow$I', 'I$\\rightarrow$E']
    leg2 = ['I$\\rightarrow$E uniform\nrandom']
    l1 = ax.legend([ep, ip], leg1, loc=(0.02, 1.0), frameon=False, fontsize='x-small',
            ncol=1)
    l2 = ax.legend([icp], leg2, loc=(0.45, 1.03), frameon=False, fontsize='x-small')
    plt.setp(l1.get_title(), fontsize='x-small')
    plt.setp(l2.get_title(), fontsize='x-small')
    ax.add_artist(l1)
    ax.margins(0.02)

    #arrow_clr='grey'
    #arrowprops = dict(
    #    arrowstyle = "->",
    #    linewidth=.5,
    #    color = arrow_clr,
    #    connectionstyle = "angle,angleA=0,angleB=90,rad=10")
    #
    #rnd_x, rnd_y = 0., inh_const
    #plt.annotate('Random uniform',
    #            (rnd_x, rnd_y), xytext=(0.4, 1.2), textcoords='axes fraction',
    #            arrowprops=arrowprops, ha='left', va='center',  size='small',
    #            color=arrow_clr, zorder=-1)
    

class ConnectionFunctionPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(ConnectionFunctionPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        d = np.arange(x0, x1+dx, dx)
        y_dim = np.sqrt(3)/2.0
        ES_pAMPA_mu = y_dim/2.0
        ES_pAMPA_sigma = 0.5/6
        ES_pGABA_sigma = 0.5/6
        ES_pGABA_const = 0.1
        shift = 0.1

        figsize = (3, 1.5)
        left    = 0.2
        bottom  = 0.25
        top     = 0.75
        right   = 0.95

        # Excitatory surround
        ES_exc_profile         = np.exp(-(np.abs(d) - ES_pAMPA_mu)**2/2/ES_pAMPA_sigma**2)
        ES_exc_profile_shifted = np.exp(-(np.abs(d - shift) - ES_pAMPA_mu)**2/2/ES_pAMPA_sigma**2)
        ES_inh_profile         = (1-ES_pGABA_const)*np.exp(-d**2/2./ES_pGABA_sigma**2) + ES_pGABA_const

        fig = self._get_final_fig(figsize)
        ax = fig.add_axes(Bbox.from_extents(left, bottom, right, top))
        plotWeights(ax, d, ES_exc_profile, ES_inh_profile, ES_pGABA_const,
                    linewidth=self.config['scale_factor'])
        ax.set_xticks([-.5, .5])
        ax.xaxis.set_minor_locator(ti.MultipleLocator(0.5))
        ax.xaxis.set_label_coords(x=0.5, y=-0.2)
        fileBase = 'fig_conn_func'
        fileName = "{0}/{1}_E_surr.pdf".format(self.config['output_dir'], fileBase)
        plt.savefig(fileName, transparent=True)

