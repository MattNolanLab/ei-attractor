'''Plotters for network testing.'''
from __future__ import absolute_import, print_function

import os.path

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

from grid_cell_model.plotting.global_defs import globalAxesSettings
from grid_cell_model.parameters.param_space import TrialSet
from grid_cell_model.parameters.data_sets import DictDataSet
from grid_cell_model.data_storage.sim_models import ei
from simtools.plotting.plotters import FigurePlotter
from noisefigs.EI_plotting.rasters import plotEIRaster

__all__ = [
    'PopulationActivityPlotter',
]


class EIDataSet(DictDataSet):
    '''A data set for the excitatory-inhibitory network.'''
    def __init__(self, data_dict):
        super(EIDataSet, self).__init__(data_dict)

    def get_e_spikes(self):
        '''Retrieve excitatory spikes.'''
        return ei.MonitoredSpikes(self.data, 'spikeMon_e', 'net_Ne')

    def get_i_spikes(self):
        '''Retrieve inhibitory spikes.'''
        return ei.MonitoredSpikes(self.data, 'spikeMon_i', 'net_Ni')


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

        saver = self.myc['fig_saver']
        saver.set_file_name(self.get_fname('population_activity'))
        saver.ext = "pdf"
        saver.set_backend_params(dpi=300, transparent=True)

        for trial_no in range(len(data)):
            fig = self._get_final_fig(self.myc['fig_size'])
            ax_raster = fig.add_axes(Bbox.from_extents(rl, rb, rr, rt))
            plotEIRaster(data[trial_no].get_e_spikes(),
                         data[trial_no].get_i_spikes(),
                         self.myc['t_limits'],
                         ax=ax_raster,
                         sigmaTitle=False,
                         markersize=2*self.config['scale_factor'],
                         rasterized=True,
                         ann_EI=True,
                         scaleBar=1000,
                         scaleX=.9,
                         scaleY=-.05,
                         scaleTextYOffset=.02)
            saver.savefig(fig)

        saver.close()

