'''Synaptic weight constructors.'''
from __future__ import absolute_import, print_function, division

import numpy as np


class WeightConstructor(object):
    '''Generates synaptic weights from template synaptic profiles.'''
    def generate_weights(self, profile, weight):
        '''Generate synaptic weights from the synaptic profile function.

        Parameters
        ----------
        profile : np.ndarray
            Synaptic profile array. It should be normalized between <0, 1>.
        weight : float
            The actual synaptic weight.

        Returns
        -------
        weights : np.ndarray
            Synaptic weights with the same shape as ``profile``
        '''
        raise NotImplementedError()


class IsomorphicConstructor(WeightConstructor):
    '''Generates weights in an isomorhic way, i.e. the weights will equal the
    synaptic profile.
    '''
    def generate_weights(self, profile, weight):
        return profile * weight


class ProbabilisticConstructor(WeightConstructor):
    '''A constructor that treats the synaptic profile as a probability.

    The weights have a constant value but their presence is determined by the
    probability given in the synaptic profile.
    '''
    def generate_weights(self, profile, weight):
        return (np.random.rand(*profile.shape) < profile) * weight


