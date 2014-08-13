'''Base classes for figure plotters'''
from __future__ import absolute_import, division, print_function

import logging
from abc import ABCMeta, abstractmethod

import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

from grid_cell_model.plotting.global_defs import globalAxesSettings
import pyentropy
from minepy import MINE

logger = logging.getLogger(__name__)

class _FigureContextManager(object):
    def __init__(self, fname, plotter, transparent):
        self.fname = fname
        self.plotter = plotter
        self.transparent = transparent

    def __enter__(self):
        self.fig = self.plotter.get_fig()
        return self.fig, self.plotter.get_ax(self.fig)

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None and exc_value is None and traceback is None:
            self.fig.savefig(self.fname, dpi=300, transparent=self.transparent)
        plt.close(self.fig)
        return False


class FigurePlotter(object):
    __metaclass__ = ABCMeta
    '''Performs plotting without returning any value.'''
    def __init__(self, config, env):
        self._env = env

    @property
    def config(self):
        '''Global configuration.'''
        return self.env.config

    @property
    def myc(self):
        '''Class-specific configuration.'''
        return self._get_class_config()

    @property
    def env(self):
        return self._env

    def _get_class_config(self):
        logger.debug("Fetching config section for class <%s>",
                     self.__class__.__name__)
        try:
            return self.config[self.__class__.__name__]
        except KeyError:
            logger.warn('Config section [%s] missing.' %
                        self.__class__.__name__)
            return {}

    def pre_plot(self, *args, **kwargs):
        pass

    @abstractmethod
    def plot(self, *args, **kwargs):
        raise NotImplementedError()

    def post_plot(self, *args, **kwargs):
        pass

    def save(self, *args, **kwargs):
        pass

    def run_all(self, *args, **kwargs):
        self.pre_plot(*args, **kwargs)
        self.plot(*args, **kwargs)
        self.post_plot(*args, **kwargs)
        self.save(*args, **kwargs)

    def figure_and_axes(self, fname, config):
        return _FigureContextManager(fname, self, config['transparent'])

    def get_fig(self):
        raise NotImplementedError()

    def get_ax(self, fig):
        raise NotImplementedError()

    def _get_final_fig(self, nominal_fig_size):
        global_scale = self.config['scale_factor']
        class_scale  = self._get_class_config().get('scale_factor', 1.)
        final_scale = global_scale*class_scale
        fig = plt.figure(figsize=np.asarray(nominal_fig_size)*final_scale)
        return fig


class SweepPlotter(FigurePlotter):
    '''Parameter sweeps plotter'''
    def __init__(self, *args, **kwargs):
        super(SweepPlotter, self).__init__(*args, **kwargs)

    def _get_sweep_config(self):
        return self.config['sweeps']

    def get_fig(self):
        fig_size = np.asarray(self.config['sweeps']['fig_size'])
        return self._get_final_fig(fig_size)
    
    def get_ax(self, fig):
        color_bar_pos = self._get_class_config()['cbar_kw']['location']
        l, b, w, h = self.config['sweeps']['bbox']
        if color_bar_pos == 'right':
            left = l
        else:
            left = .12
    
        right = left + w
        top = b + h
        return fig.add_axes(Bbox.from_extents(left, b, right, top))


class ExampleSetting(object):
    '''A setting that specifies where an example is in the 2D sweep parameter
    space
    '''
    def __init__(self, r, c, trialNum, ps, noise_sigma):
        self.r = r
        self.c = c
        self.trialNum = trialNum
        self.ps = ps
        self.noise_sigma = noise_sigma


class ProbabilityPlotter(FigurePlotter):
    def __init__(self, *args, **kwargs):
        super(ProbabilityPlotter, self).__init__(*args, **kwargs)

    def plotDistribution(self, X, Y, ax, noise_sigma=None, **kw):
        xlabel = kw.get('xlabel', 'P(bump)') 
        ylabel = kw.get('ylabel', '$Power_\gamma$')
        yticks = kw.get('yticks', True)
        bins   = kw.get('bins', [40, 50])
        range  = kw.get('range', [[0, 1], [-.2, .8]])

        H, xedges, yedges = np.histogram2d(
                X.flatten(),
                Y.flatten(),
                bins=bins,
                range=range,
                normed=True)

        globalAxesSettings(ax)
        ax.pcolormesh(xedges, yedges, H.T, vmin=0, rasterized=True)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if noise_sigma is not None:
            ax.set_title("$\sigma$ = %d pA" % int(noise_sigma))
        else:
            ax.set_title("All noise levels")
        if not yticks:
            ax.yaxis.set_ticklabels([])

    def mutual_information(self, X, Y, title=None, nbins_X=50, nbins_Y=50,
            noise_sigma='all'): 
        #import pdb; pdb.set_trace()
        no_nans_idx = np.logical_not(np.logical_or(np.isnan(X), np.isnan(Y)))
        Xq, _, _ = pyentropy.quantise(X[no_nans_idx], nbins_X)
        Yq, _, _ = pyentropy.quantise(Y[no_nans_idx], nbins_Y)
        s = pyentropy.DiscreteSystem(Yq, (1, nbins_Y), Xq, (1, nbins_X))
        s.calculate_entropies()

        # MINE
        mine = MINE()
        mine.compute_score(X.flatten(), Y.flatten())

        # Linear regression
        slope, intercept, r, p, stderr = \
                scipy.stats.linregress(X[no_nans_idx], Y[no_nans_idx])

        #import pdb; pdb.set_trace()
        if title is not None:
            print(title)
        print(" MIC/MI/r^2/p/slope for %s:\t%.3f\t%.3f\t%s\t%s\t%s" %
                (noise_sigma, mine.mic(), s.I(), r**2, p, slope))
       

