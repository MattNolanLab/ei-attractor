'''Plotters for network testing.

.. currentmodule:: noisefigs.plotters.nettests

Classes
-------

.. autosummary::

    EIDataSet
    PopulationActivityPlotter
'''
from __future__ import absolute_import, print_function

import os.path

from matplotlib.transforms import Bbox

from grid_cell_model.parameters.param_space import TrialSet
from grid_cell_model.parameters.data_sets import DictDataSet
from grid_cell_model.data_storage.sim_models import ei
from simtools.plotting.plotters import FigurePlotter
from noisefigs.EI_plotting.rasters import plotEIRaster

from ..EI_plotting import examples

__all__ = [
    'EIDataSet',
    'PopulationActivityPlotter',
]


class EIDataSet(DictDataSet):
    '''A data set for the excitatory-inhibitory network.'''
    def __init__(self, data_dict):
        super(EIDataSet, self).__init__(data_dict)

    def get_e_spikes(self):
        '''Retrieve excitatory spikes.'''
        return ei.MonitoredTorusSpikes(self.data, 'spikeMon_e', 'Ne_x', 'Ne_y')

    def get_i_spikes(self):
        '''Retrieve inhibitory spikes.'''
        return ei.MonitoredTorusSpikes(self.data, 'spikeMon_i', 'Ni_x', 'Ni_y')


class PopulationActivityPlotter(FigurePlotter):
    '''Simply plot population activity as a raster plot and 2D population rate
    snapshots.
    '''
    def __init__(self, *args, **kwargs):
        super(PopulationActivityPlotter, self).__init__(*args, **kwargs)
        self._trial_data = None

    @property
    def trial_data(self):
        '''Return data containing all the trials as a TrialSet.'''
        if self._trial_data is None:
            data_path = os.path.join(self.config['data_root'],
                                     self.config['data_file_name'])
            self._trial_data = TrialSet(data_path, 'r', data_set_cls=EIDataSet)
        return self._trial_data

    def plot(self, *args, **kwargs):
        # Plot all the trials into separate pages.
        data = self.trial_data
        rl, rb, rr, rt = self.myc['raster_rect']
        t_limits = self.myc['t_limits']
        max_e_rate = self.myc.get('max_e_rate', True)
        max_i_rate = self.myc.get('max_i_rate', True)
        scale_bar = self.myc.get('scale_bar', 500)
        scale_x = self.myc.get('scale_x', .85)
        scale_y = self.myc.get('scale_y', -.05)
        ann_ei = self.myc.get('ann_ei', True)
        y_label_pos = self.myc.get('y_label_pos', -0.22)
        reshape_senders = self.myc.get('reshape_senders', True)

        saver = self.myc['fig_saver']
        saver.set_file_name(self.get_fname('population_activity'))
        saver.ext = "pdf"
        saver.set_backend_params(dpi=300, transparent=True)

        for trial_no in range(len(data)):
            fig = self._get_final_fig(self.myc['fig_size'])
            e_spikes = data[trial_no].get_e_spikes()
            i_spikes = data[trial_no].get_i_spikes()

            # Raster plot
            ax_raster = fig.add_axes(Bbox.from_extents(rl, rb, rr, rt))
            plotEIRaster(e_spikes,
                         i_spikes,
                         t_limits,
                         ax=ax_raster,
                         sigmaTitle=False,
                         markersize=2 * self.config['scale_factor'],
                         rasterized=True,
                         ann_EI=ann_ei,
                         scaleBar=scale_bar,
                         scaleX=scale_x,
                         scaleY=scale_y,
                         scaleTextYOffset=.02,
                         reshape_senders=reshape_senders,
                         ylabelPos=y_label_pos)

            # E and I 2D population plots
            win_dt = 125.  # ms
            winLen = 250.  # ms
            tStart = t_limits[0]
            tEnd   = t_limits[1] - win_dt
            fr_e, frt_e = e_spikes.slidingFiringRate(tStart, tEnd, win_dt,
                                                     winLen)
            e_max = examples.plotBumpSnapshots(fr_e, frt_e, self.myc['snapshot_tstep'],
                                               fig=fig,
                                               axesCoords=self.myc['e_snapshots_rect'],
                                               bumpQuality=False,
                                               timeTitles=True,
                                               maxRate=False,
                                               fontsize='x-small',
                                               rateYPos=-.2)
            rateText = "%.0f Hz" % e_max
            rect = self.myc['e_snapshots_rect']
            fig.text(rect[2], rect[3], rateText, ha='left',
                     va='top', size='small', weight='bold',
                     clip_on=False)

            fr_i, frt_i = i_spikes.slidingFiringRate(tStart, tEnd, win_dt,
                                                     winLen)
            i_max = examples.plotBumpSnapshots(fr_i, frt_i, self.myc['snapshot_tstep'],
                                               fig=fig,
                                               axesCoords=self.myc['i_snapshots_rect'],
                                               bumpQuality=False,
                                               timeTitles=False,
                                               maxRate=False,
                                               fontsize='x-small')
            rateText = "%.0f Hz" % i_max
            rect = self.myc['i_snapshots_rect']
            fig.text(rect[2], rect[3], rateText, ha='left',
                     va='top', size='small', weight='bold',
                     clip_on=False)

            saver.savefig(fig)
        saver.close()
