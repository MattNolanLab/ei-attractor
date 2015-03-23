'''Handle metadata of parameter spaces.'''
from __future__ import absolute_import, print_function, division

from .param_space import DataSpace


class MetaDataExtractor(object):
    '''Extracts parameter names and parameter iteration data.

    Parameters
    ----------
    space : JobTrialSpace2D
        A data space to extract the information from.
    iter_names : a tuple of strings
        A tuple that describes the names of the parameters during the iteration
        process. The format is (ROW_NAME, COLUMN_NAME), which should correspond
        to (Y, X) respectively.
    '''
    def __init__(self, space, iter_names):
        self._sp = space
        self._iter_names = iter_names

    @property
    def x_data(self):
        '''Return 2D X coordinates in the parameter space.'''
        raise NotImplementedError()

    @property
    def y_data(self):
        '''Return 2D Y coordinates in the parameter space.'''
        raise NotImplementedError()

    @property
    def xy_data(self):
        '''Return a tuple of 2D X and Y coordinates in the parameter space.'''
        return (self.x_data, self.y_data)

    @property
    def x_label(self):
        '''X data label.'''
        return self._iter_names[1]

    @property
    def y_label(self):
        '''Y data label.'''
        return self._iter_names[0]

    @property
    def xy_labels(self):
        '''X and Y data labels as a tuple.'''
        return (self.x_label, self.y_label)


class EISweepExtractor(MetaDataExtractor):
    '''An extractor for gE and gI iteration data.

    Parameters
    ----------
    space : JobTrialSpace2D
        A data space to extract the information from.
    r, c, trial_num : int
        Row and column of the parameter space files to extract additional
        information from, e.g. network size.
    normalize : bool
        If ``True``, the iteration data will be divided by network size.
    '''
    ROW_LABEL = 'g_AMPA_total'
    COL_LABEL = 'g_GABA_total'

    def __init__(self, space, r=0, c=0, trial_num=0, normalize=True):
        super(EISweepExtractor, self).__init__(space, (self.ROW_LABEL,
                                                       self.COL_LABEL))
        self._r = r
        self._c = c
        self._trial_num = trial_num
        self._normalize = normalize

    @property
    def x_data(self):
        I, _ = self._sp.getIteratedParameters(self.xy_labels)

        if self._normalize:
            Ni = DataSpace.getNetParam(
                self._sp[self._r][self._c][self._trial_num].data, 'net_Ni')
        else:
            Ni = 1.

        return I / Ni

    @property
    def y_data(self):
        _, E = self._sp.getIteratedParameters(self.xy_labels)

        if self._normalize:
            Ne = DataSpace.getNetParam(
                self._sp[self._r][self._c][self._trial_num].data, 'net_Ne')
        else:
            Ne = 1.

        return E / Ne


class GEProfileWidthExtractor(EISweepExtractor):
    '''Metadata extractor for GE on Y axis and pAMPA_sigma on X axis.'''
    ROW_LABEL = 'g_AMPA_total'
    COL_LABEL = 'pAMPA_sigma'

    def __init__(self, space, r=0, c=0, trial_num=0, normalize=True):
        super(GEProfileWidthExtractor, self).__init__(space, r, c, trial_num,
                                                      normalize)

    @property
    def x_data(self):
        sigma, _ = self._sp.getIteratedParameters(self.xy_labels)
        return sigma
