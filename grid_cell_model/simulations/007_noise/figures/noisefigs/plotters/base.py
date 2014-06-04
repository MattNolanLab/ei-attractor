'''Base classes for figure plotters'''
from __future__ import absolute_import, division, print_function

import logging
from abc import ABCMeta, abstractmethod

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

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



