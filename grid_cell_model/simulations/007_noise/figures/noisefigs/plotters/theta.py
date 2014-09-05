from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
from grid_cell_model.plotting.signal import signalPlot
from grid_cell_model.plotting.low_level import xScaleBar

from .base import FigurePlotter

__all__ = [
    'ThetaSignalPlotter',
    'PACExamplePlotter',
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
        for n_theta_cycles in range(1, 9):
            T = n_theta_cycles / self.freq * 1e3
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

            fname = self.config['output_dir'] + "/theta_signal_{0}.pdf"
            self.fig.savefig(fname.format(n_theta_cycles), dpi=300,
                             transparent=True)
            plt.close(self.fig)


##############################################################################
class PACExamplePlotter(FigurePlotter):
    '''Phase-amplitude coupling example plotter.'''
    dt = .1  # ms
    theta_freq = 8.   # Hz
    gamma_freq = 90.  # Hz
    const = 0 # Fraction of max. theta

    def __init__(self, *args, **kwargs):
        super(PACExamplePlotter, self).__init__(*args, **kwargs)

    def plot(self, *args, **kwargs):
        for n_theta_cycles in range(1, 9):
            T = n_theta_cycles / self.theta_freq * 1e3
            l, b, r, top = self.myc['bbox']

            self.fig = self._get_final_fig(self.myc['fig_size'])
            self.ax_theta = self.fig.add_subplot(3, 1, 1)

            # Only theta
            t = np.arange(0, T+self.dt, self.dt)
            theta = self.const + .5 * (1. + 
                    np.cos(2*np.pi*self.theta_freq*t*1e-3 - np.pi)) * (1 - self.const)
            signalPlot(
                    t, theta,
                    ax=self.ax_theta,
                    color=self.myc['theta_color'],
                    zeroLine=False)
            self.ax_theta.axis('off')
            self.ax_theta.set_xlim([t[0], t[-1]]) 
            self.ax_theta.set_ylim(-.2, 1.2)
            self.ax_theta.set_title('A', x=0, y=.95, size=14, weight='bold',
                                    ha='left', va='top')

            # PAC only gamma
            self.ax_pac = self.fig.add_subplot(3, 1, 2)
            gamma = np.cos(2*np.pi*self.gamma_freq*t*1e-3 - np.pi) * theta
            signalPlot(
                    t, gamma,
                    ax=self.ax_pac,
                    color=self.myc['gamma_color'],
                    zeroLine=False)

            self.ax_pac.axis('off')
            self.ax_pac.set_xlim([t[0], t[-1]]) 
            self.ax_pac.set_ylim(-1.02, 1.02)
            self.ax_pac.set_title('B', x=0, y=1, size=14, weight='bold',
                                  ha='left', va='top')

            # PAC Theta + gamma
            self.ax_pac_all = self.fig.add_subplot(3, 1, 3)
            self.ax_pac_all.hold('on')

            signalPlot(
                    t, theta,
                    ax=self.ax_pac_all,
                    color=self.myc['theta_color'],
                    zeroLine=False)

            gamma = .25 * np.cos(2*np.pi*self.gamma_freq*t*1e-3 - np.pi) * theta
            signalPlot(
                    t, theta + gamma,
                    ax=self.ax_pac_all,
                    color=self.myc['gamma_color'],
                    zeroLine=False)
            bar_x = 1 - 5. / 8 / n_theta_cycles
            xScaleBar(50, x=bar_x, y=.1, ax=self.ax_pac_all, size='small')

            self.ax_pac_all.axis('off')
            self.ax_pac_all.set_xlim([t[0], t[-1]]) 
            self.ax_pac_all.set_ylim(-.02, 1.5)
            self.ax_pac_all.set_title('C', x=0, y=1, size=14, weight='bold',
                                      ha='left', va='top')

            self.fig.subplots_adjust(left=l, bottom=b, right=r, top=top)

            fname = self.config['output_dir'] + "/pac_example_{0}.pdf"
            self.fig.savefig(fname.format(n_theta_cycles), dpi=300,
                             transparent=True)
            plt.close(self.fig)

