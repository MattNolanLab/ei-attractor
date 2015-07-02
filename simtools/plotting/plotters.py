'''Base classes for figure plotters'''
from __future__ import absolute_import, division, print_function

import logging

import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

class _FigureContextManager(object):
    '''Context manager for matplotlib figures that contain one axes.

    When the context is created, the manager creates the figure and axes. When
    leaving the context, the figure is saved.'''
    def __init__(self, fname, plotter, transparent):
        self.fig = None
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


class Computation(object):
    '''Performs computation without returning any value.'''
    def __init__(self, config, env):
        self._env = env
        self._config = config

    def _get_class_config(self):
        '''Get class-specific configuration.'''
        logger.debug("Fetching config section for class <%s>",
                     self.__class__.__name__)
        try:
            return self.config[self.__class__.__name__]
        except KeyError:
            logger.warn('Config section [%s] missing.',
                        self.__class__.__name__)
            return {}

    @property
    def config(self):
        '''Global configuration.'''
        return self._config

    @property
    def env(self):
        '''The parent environment.'''
        return self._env

    def get_fname(self, template, *args, **kwargs):
        '''Format the template of the file name into a ready-to-use file path
        for writing data.

        Parameters
        ----------
        template : str
            String of the file name template. This can be prefixed and suffixed
            with various paths. *args and **kwargs will be passed to
            template.format() method.

        Returns
        -------
        fname : str
            Fully qualifying path to save the data into.
        '''
        return "{output_dir}/{fname_prefix}{obj_prefix}{file_str}".format(
            output_dir=self.config['output_dir'],
            fname_prefix=self.config.get('fname_prefix', ''),
            obj_prefix=self.myc.get('fname_prefix', ''),
            file_str=template.format(*args, **kwargs))

    @property
    def myc(self):
        '''Class-specific configuration.'''
        return self._get_class_config()

    def run_all(self, *args, **kwargs):
        '''Run all computation.'''
        raise NotImplementedError()


class FigurePlotter(Computation):
    '''Performs plotting without returning any value.'''
    def pre_plot(self, *args, **kwargs):
        '''Do this before plotting.'''
        pass

    def plot(self, *args, **kwargs):
        '''Actual plotting.'''
        raise NotImplementedError()

    def post_plot(self, *args, **kwargs):
        '''Post-plotting work.'''
        pass

    def save(self, *args, **kwargs):
        '''Figure saving.'''
        pass

    def run_all(self, *args, **kwargs):
        self.pre_plot(*args, **kwargs)
        self.plot(*args, **kwargs)
        self.post_plot(*args, **kwargs)
        self.save(*args, **kwargs)

    def figure_and_axes(self, fname, config):
        '''Return a context manager with the appropriate figure and axes.'''
        return _FigureContextManager(fname, self, config['transparent'])

    def get_fig(self):
        '''Get the figure for this object.'''
        raise NotImplementedError()

    def get_ax(self, fig):
        '''Get the axes for this object.'''
        raise NotImplementedError()

    def _get_final_fig(self, nominal_fig_size):
        '''Create the figure that is appropriately scaled.

        The parameters taken into account are the root scale_factor, and
        per-object scale_factor, plus ``nominal_fig_size``.
        '''
        global_scale = self.config['scale_factor']
        class_scale  = self._get_class_config().get('scale_factor', 1.)
        final_scale = global_scale*class_scale
        fig = plt.figure(figsize=np.asarray(nominal_fig_size)*final_scale)
        return fig
