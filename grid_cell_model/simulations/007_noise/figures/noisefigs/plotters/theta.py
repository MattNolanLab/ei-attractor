from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
from grid_cell_model.plotting.signal import signalPlot

from .base import FigurePlotter

__all__ = [
    'ThetaSignalPlotter'
]

##############################################################################

class ThetaSignalPlotter(FigurePlotter):
    '''Theta signal plotter.'''
    dt = .1  # ms
    freq = 8. # Hz
    const = .4 # Fraction of max. theta

    def __init__(self, *args, **kwargs):
        super(ThetaSignalPlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        T = self.myc['T']
        l, b, r, top = self.myc['bbox']

        self.fig = self._get_final_fig(self.myc['fig_size'])
        self.ax = self.fig.add_subplot(111)

        t = np.arange(0, T+self.dt, self.dt) * 1e-3
        theta = self.const + .5 * (1. + 
                np.cos(2*np.pi*self.freq*t - np.pi)) * (1 - self.const)
        self.ax.fill_between(
                t, theta,
                edgecolor='None',
                color=self.myc['color'])
        self.ax.axis('off')
        self.ax.set_xlim([t[0], t[-1]]) 
        self.ax.set_ylim(-.02, 1.02)
        self.fig.subplots_adjust(left=l, bottom=b, right=r, top=top)

    def save(self, *args, **kwargs):
        fname = self.config['output_dir'] + "/theta_signal.pdf"
        self.fig.savefig(fname, dpi=300, transparent=True)
        plt.close(self.fig)

