'''Slope manipulations for the noise simulations.'''
from __future__ import absolute_import, print_function, division

import os.path

import numpy as np
from grid_cell_model.otherpkg.log import getClassLogger
from simtools.storage import DataStorage

baseLogger = getClassLogger('SlopeSelector', __name__)
probLogger = getClassLogger('ProbabilisticConnectionsSelector', __name__)


__all__ = [
    'DefaultSelector',
    'NoThetaSelector',
    'ProbabilisticConnectionsSelector',
]


class SlopeSelector(object):
    '''Extracts the slope data from bump slope data files.

    Parameters
    ----------
    data_root : str
        Path to the directory containing the slope data files.
    threshold : float
        Threshold below which the values will be ignored. This can avoid
        unnecessary simulations in which bump slope is close to zero. Set to
        -infinity if you want to ignore this threshold.
    '''
    def __init__(self, data_root, threshold, file_name_template):
        self._data_root = data_root
        self._threshold = threshold
        self._file_name_template = file_name_template

    @property
    def data_root(self):
        '''Root directory for the slope data.'''
        return self._data_root

    @property
    def threshold(self):
        '''Threshold for the slopes.'''
        return self._threshold

    def get_slopes(self, noise_sigma):
        '''Retrieve the slopes for the given noise parameter value.'''
        file_name = os.path.join(self.data_root, self._file_name_template)
        file_name = file_name.format(int(noise_sigma))
        baseLogger.info('Using the following file for bump slope data:\n  %s',
                        file_name)
        return self._retrieve_slopes(file_name)

    def _retrieve_slopes(self, file_name):
        '''Retrieve the slopes from a file ``file_name``.'''
        ds = DataStorage.open(file_name, 'r')
        slopes = ds['lineFitSlope'].flatten()
        ds.close()
        slopes[slopes < self.threshold] = np.nan
        return slopes


class DefaultSelector(SlopeSelector):
    '''A selector that loads the default slope data (E-surround, basic network
    without any additions).'''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_{0}pA.h5'
        super(DefaultSelector, self).__init__(data_root, threshold, template)


class NoThetaSelector(SlopeSelector):
    '''A slope selector for no-theta simulations.'''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_no_theta_{0}pA.h5'
        super(NoThetaSelector, self).__init__(data_root, threshold, template)


class ProbabilisticConnectionsSelector(SlopeSelector):
    '''A selector that retrieves data for simulations with probabilistic
    connections.'''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_probabilistic_connections_{0}pA.h5'
        super(ProbabilisticConnectionsSelector, self).__init__(data_root,
                                                               threshold,
                                                               template)
