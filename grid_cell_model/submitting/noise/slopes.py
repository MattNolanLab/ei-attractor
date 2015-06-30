'''Slope manipulations specific to [SOLANKA2015]_

.. currentmodule:: grid_cell_model.submitting.noise.slopes

These classes read bump slope data from the HDF5 files stored in the repository
in ``grid_cell_model/simulations/007_noise/bump_slope_data``. Each simulation
uses its specific slope extractor that accesses an appropriate file.

.. warning::

    In general, it is not good to manipulate these classes unless you
    **exactly** know what you are doing, although there should not be any
    problems in adding new slope extractors.

Classes
-------

.. autosummary::

    SlopeSelector
    DefaultSelector
    PickedDefaultSelector
    NoThetaSelector
    ProbabilisticConnectionsSelector
    IIConnectionsSelector
    EEConnectionsSelector
    ISurroundOrigSelector
    ISurroundPastollSelector
    PickedISurroundPastollSelector
'''
from __future__ import absolute_import, print_function, division

import os.path

import numpy as np
from grid_cell_model.otherpkg.log import getClassLogger
from simtools.storage import DataStorage

baseLogger = getClassLogger('SlopeSelector', __name__)
probLogger = getClassLogger('ProbabilisticConnectionsSelector', __name__)


__all__ = [
    'SlopeSelector',
    'DefaultSelector',
    'PickedDefaultSelector',
    'NoThetaSelector',
    'ProbabilisticConnectionsSelector',
    'IIConnectionsSelector',
    'EEConnectionsSelector',
    'ISurroundOrigSelector',
    'ISurroundPastollSelector',
    'PickedISurroundPastollSelector',
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

    def get_slopes(self, noise_sigma, flatten=True):
        '''Retrieve the slopes for the given noise parameter value.'''
        file_name = os.path.join(self.data_root, self._file_name_template)
        file_name = file_name.format(int(noise_sigma))
        baseLogger.info('Using the following file for bump slope data:\n  %s',
                        file_name)
        return self._retrieve_slopes(file_name, flatten)

    def _retrieve_slopes(self, file_name, flatten):
        '''Retrieve the slopes from a file ``file_name``.'''
        ds = DataStorage.open(file_name, 'r')
        slopes = ds['lineFitSlope']
        if flatten:
            slopes = slopes.flatten()
        ds.close()
        slopes[slopes < self.threshold] = np.nan
        return slopes


class DefaultSelector(SlopeSelector):
    '''A selector that loads the default slope data (E-surround, basic network
    without any additions).'''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_{0}pA.h5'
        super(DefaultSelector, self).__init__(data_root, threshold, template)


class PickedDefaultSelector(DefaultSelector):
    '''Bump slope selector that picks a value from the data obtained by
    DefaultSelector.

    The picked value is then copied to all simulation runs.

    Parameters
    ----------
    data_root, threshold
    row, col : int
        Indexes into data obtained by DefaultSelector
    dimensions : pair of int
        Dimensions of the resulting parameter sweep.
    '''
    def __init__(self, data_root, threshold, row, col, dimensions):
        super(PickedDefaultSelector, self).__init__(data_root, threshold)
        self._row = row
        self._col = col
        self._dimensions = dimensions

    def get_slopes(self, noise_sigma):
        slopes = super(PickedDefaultSelector, self).get_slopes(noise_sigma,
                                                               flatten=False)
        return (np.ones(self._dimensions[0] * self._dimensions[1]) *
                slopes[self._row, self._col])


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


class IIConnectionsSelector(SlopeSelector):
    '''A selector that retrieves data for simulations with I-->I
    connectivity.'''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_ii_connections_{0}pA.h5'
        super(IIConnectionsSelector, self).__init__(data_root, threshold,
                                                    template)


class EEConnectionsSelector(SlopeSelector):
    '''A selector that retrieves data for simulations with E-->E
    connectivity.'''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_ee_connections_{0}pA.h5'
        super(EEConnectionsSelector, self).__init__(data_root, threshold,
                                                    template)


class ISurroundOrigSelector(SlopeSelector):
    '''I-surround slope selector with other parameters default.

    A selector that retrieves data for the I-surround simulations (with all
    other config as in the E-surround).
    '''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_i_surround_original_{0}pA.h5'
        super(ISurroundOrigSelector, self).__init__(data_root, threshold,
                                                    template)


class ISurroundPastollSelector(SlopeSelector):
    '''A selector that retrieves data for the I-surround simulations (Pastoll
    et al. configuration).'''
    def __init__(self, data_root, threshold):
        template = 'bump_slope_i_surround_pastoll_{0}pA.h5'
        super(ISurroundPastollSelector, self).__init__(data_root, threshold,
                                                       template)


class PickedISurroundPastollSelector(ISurroundPastollSelector):
    '''Picks one value from :class:`~ISurroundPastollSelector`.

    Bump slope selector that picks a value from the data obtained by
    ISurroundPastollSelector and copies this slope value to all simulations.

    Parameters
    ----------
    data_root, threshold
    row, col : int
        Indexes into data obtained by DefaultSelector
    dimensions : pair of int
        Dimensions of the resulting parameter sweep.
    '''
    def __init__(self, data_root, threshold, row, col, dimensions):
        super(PickedISurroundPastollSelector, self).__init__(data_root,
                                                             threshold)
        self._row = row
        self._col = col
        self._dimensions = dimensions

    def get_slopes(self, noise_sigma):
        slopes = super(PickedISurroundPastollSelector,
                       self).get_slopes(noise_sigma, flatten=False)
        return (np.ones(self._dimensions[0] * self._dimensions[1]) *
                slopes[self._row, self._col])
