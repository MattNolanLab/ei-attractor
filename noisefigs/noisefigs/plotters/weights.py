'''Print synaptic connection profiles.

.. currentmodule:: noisefigs.plotters.weights

Classes
-------

.. autosummary::

    ConnectionFunctionPlotter
'''
from __future__ import absolute_import, print_function

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ti
from matplotlib.transforms import Bbox
from grid_cell_model.plotting.global_defs import globalAxesSettings, createColorbar
from simtools.plotting.plotters import FigurePlotter

__all__ = [
    'ConnectionFunctionPlotter',
]

class ConnectionFunctionPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(ConnectionFunctionPlotter, self).__init__(*args, **kwargs)

    def plotWeights(self, ax, d, exc_profile, inh_profile, inh_const,
                    linewidth, x_range):
        x0, x1, _ = x_range

        plt.hold('on')
        ax = plt.gca()
        globalAxesSettings(ax)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ep, = plt.plot(d, exc_profile, linewidth=linewidth, color='red', label="E")
        ip, = plt.plot(d, inh_profile, linewidth=linewidth, color='blue', label="I")

        ax.set_xlabel(self.myc.get('xlabel', "'Distance'"))
        ax.set_ylabel('G (nS)')
        ax.yaxis.set_ticks([0, 1])
        ax.yaxis.set_ticklabels([0, '$g_{E/I}$'])

        leg1 = ['E$\\rightarrow$I', 'I$\\rightarrow$E']
        l1 = ax.legend([ep, ip], leg1, **self.myc['leg1_kwargs'])
        plt.setp(l1.get_title(), fontsize='x-small')

        # If we want the uniform random, draw it
        if self.myc.get('uniform_random', True):
            leg2 = ['I$\\rightarrow$E uniform\nrandom']
            icp, = plt.plot(d, [inh_const]*len(d), ':', color='blue')
            l2 = ax.legend([icp], leg2, **self.myc['leg2_kwargs'])
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


    def plot(self, *args, **kwargs):
        ylabel_coords = self.myc.get('ylabel_coords', None)
        x0, x1, dx = self.myc.get('x_range', (-.5, .5, .001))
        l, b, r, t = self.myc['bbox_rect']

        d = np.arange(x0, x1+dx, dx)
        y_dim = np.sqrt(3)/2.0
        ES_pAMPA_mu = y_dim/2.0
        ES_pAMPA_sigma = 0.5/6
        ES_pGABA_sigma = 0.5/6
        ES_pGABA_const = 0.013
        shift = 0.1

        # Excitatory surround
        ES_exc_profile         = np.exp(-(np.abs(d) - ES_pAMPA_mu)**2/2/ES_pAMPA_sigma**2)
        ES_exc_profile_shifted = np.exp(-(np.abs(d - shift) - ES_pAMPA_mu)**2/2/ES_pAMPA_sigma**2)
        ES_inh_profile         = (1-ES_pGABA_const)*np.exp(-d**2/2./ES_pGABA_sigma**2) + ES_pGABA_const

        fig = self._get_final_fig(self.myc['fig_size'])
        ax = fig.add_axes(Bbox.from_extents(l, b, r, t))
        self.plotWeights(ax, d, ES_exc_profile, ES_inh_profile, ES_pGABA_const,
                         linewidth=self.config['scale_factor'],
                         x_range=(x0, x1, dx))
        ax.set_xticks(self.myc.get('xticks', [-.5, .5]))
        ax.xaxis.set_minor_locator(ti.MultipleLocator(0.5))
        ax.xaxis.set_label_coords(x=0.5, y=-0.2)
        if ylabel_coords is not None:
            ax.yaxis.set_label_coords(ylabel_coords[0], ylabel_coords[1])
        fileName = self.get_fname("{base}_E_surr.pdf", base='fig_conn_func')
        plt.savefig(fileName, transparent=True)

