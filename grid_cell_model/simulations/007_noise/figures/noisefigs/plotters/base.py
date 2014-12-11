'''Base classes for figure plotters'''
from __future__ import absolute_import, division, print_function

import logging

import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.backends.backend_pdf import PdfPages

from grid_cell_model.plotting.global_defs import globalAxesSettings
from ..EI_plotting import sweeps
from ..EI_plotting import aggregate as aggr
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


class Computation(object):
    '''Performs computation without returning any value.'''
    def __init__(self, config, env):
        self._env = env
        self._config = config

    def _get_class_config(self):
        logger.debug("Fetching config section for class <%s>",
                     self.__class__.__name__)
        try:
            return self.config[self.__class__.__name__]
        except KeyError:
            logger.warn('Config section [%s] missing.' %
                        self.__class__.__name__)
            return {}

    @property
    def config(self):
        '''Global configuration.'''
        return self._config

    @property
    def env(self):
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
        if "ns" in kwargs:
            kwargs["ns"] = int(kwargs["ns"])
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
        raise NotImplementedError()


class FigurePlotter(Computation):
    '''Performs plotting without returning any value.'''
    def pre_plot(self, *args, **kwargs):
        pass

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
        elif color_bar_pos == 'left':
            left = .12
        else:
            left = .2

        right = left + w
        top = b + h
        return fig.add_axes(Bbox.from_extents(left, b, right, top))

    def plot_grid_contours(self, ns_idx, ax, ps_grids):
        '''Plot grid field contours if necessary.'''
        grids_example_idx = self.config['grids']['example_idx']
        iter_list = self.config['iter_list']
        if self.myc['plot_grid_contours'][ns_idx]:
            gridData = aggr.GridnessScore(ps_grids[ns_idx], iter_list,
                                            ignoreNaNs=True, normalizeTicks=True,
                                            r=grids_example_idx[ns_idx][0],
                                            c=grids_example_idx[ns_idx][1])
            contours = sweeps.Contours(gridData,
                    self.config['sweeps']['grid_contours'])
            contours.plot(
                    ax,
                    **self.config['sweeps']['contours_kwargs'])


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
        title_size = kw.get('title_size', 'medium')

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
            ax.set_title("%d pA" % int(noise_sigma), size=title_size)
        else:
            ax.set_title("All", size=title_size)
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


class DummyPlotter(FigurePlotter):
    '''Does not plot anything. This is only for semantic identification of the
    plotter.
    '''
    def __init__(self, *args, **kwargs):
        super(DummyPlotter, self).__init__(*args, **kwargs)

    def plot(self):
        pass


##############################################################################

class MultiFigureSaver(object):
    '''Matplotlib figure saver that processes several images.

    This is an abstract class that defines basic functionality. Use some of the
    derived class ideally combined with dependency injection or some kind of
    factory.

    The concept is that the client gets this interface and simply calls savefig
    multiple times. Depending on the concrete implementation, either the
    figures will be saved as `one` combined PDF, or separate PDFs.
    '''
    def __init__(self, file_name, ext='pdf', start_cnt=0):
        self.file_name = file_name
        self.ext = ext
        self._cnt = 0

    def set_file_name(self, file_name):
        self.file_name = file_name
        self.reset()

    def reset(self, cnt_val=None):
        if cnt_val is not None:
            self._cnt = cnt_val

    def set_backend_params(self, **kwargs):
        self._backend_params = kwargs

    def savefig(self, fig):
        raise NotImplementedError()

    def close(self):
        pass


class PdfOutputSaver(MultiFigureSaver):
    '''Create a single-PDF multipage document, by calling the savefig()
    method.
    '''
    def __init__(self, file_name, ext='pdf'):
        super(PdfOutputSaver, self).__init__(file_name, ext, start_cnt=0)
        self._saver = None
        self._reset()

    def reset(self, cnt_val=None):
        self._reset()

    def _reset(self):
        self.close()
        if self.file_name is not None:
            fname = "%s.%s" % (self.file_name, self.ext)
            self._saver = PdfPages(fname)

    def savefig(self, fig):
        self._saver.savefig(fig, **self._backend_params)

    def close(self):
        if self._saver is not None:
            self._saver.close()


class SeparateMultipageSaver(MultiFigureSaver):
    '''Create PDF documents, that are numbered by a counter. This counter is
    appended to the file name when it is being saved.
    '''
    def __init__(self, file_name, ext='pdf'):
        super(SeparateMultipageSaver, self).__init__(file_name, ext,
                                                     start_cnt=0)

    def savefig(self, fig):
        fname = "%s_%d.%s" % (self.file_name, self._cnt, self.ext)
        fig.savefig(fname, **self._backend_params)
        self._cnt += 1


