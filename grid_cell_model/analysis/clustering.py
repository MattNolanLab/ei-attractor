#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
'''
Simple clustering techniques.

.. currentmodule:: analysis.clustering


Functions
---------
.. autosummary::

'''
import numpy as np
import collections
from abc import ABCMeta, abstractmethod


class ClusteredData(collections.Sequence):
    '''
    Abstract base class to hold information about clustered data.
    '''
    __metaclass__ = ABCMeta

    def __init__(self, data):
        '''
        That that has been clustered.
        '''
        self._data = data


class ThresholdClusters(ClusteredData):

    def __init__(self, data, thresholds):
        '''
        **Parameters**

        data : numpy array
            A list containing data coordinates. There is one list for each
            dimension.
        thresholds : a list of sequences
            A list, containing a sequence of thresholds, one for each
            dimension. The thresholds must be in an ascending order.

            Cluster ``i`` will contain data that fit
            :math::
            
                thresholds[dim][i] <= data[dim] < thresholds[dim][i+1]

            Therefore, the number of clusters (segments) per each dimension
            will be one less than the number of threshold items in each
            dimension.
        '''
        self.checkSizes(data, thresholds)
        self.checkThresholds(thresholds)

        ClusteredData.__init__(self, data)
        self._thresholds = thresholds

        self._nClusters = self._numClusters()
        self._ndims = len(self._data)
        self._clusters = self._doClustering(True, 0)


    def checkSizes(self, data, thresholds):
        if len(data) != len(thresholds):
            raise TypeError("For each data dimension, you need one item in " +\
                    "the threshold list")

    def checkThresholds(self, tList):
        for thIdx, thresholds in enumerate(tList):
            if np.any(np.sort(thresholds) != thresholds):
                msg = ("Threshold list for dimension {0} is not in ascending"+\
                        "order: {1}").format(thIdx, thresholds)
                raise ValueError(msg)
 

    def _doClustering(self, currentRestriction, currentDim):
        currThresholds = self._thresholds[currentDim]
        result = []
        for tIdx in xrange(len(currThresholds) - 1):
            myRestriction = np.logical_and(
                    self._data[currentDim] >= currThresholds[tIdx],
                    self._data[currentDim] <  currThresholds[tIdx + 1])
            myRestriction = np.logical_and(myRestriction, currentRestriction)
            if currentDim < self._ndims - 1:
                result.extend(self._doClustering(myRestriction, currentDim + 1))
            else:
                result.append(np.nonzero(myRestriction)[0])

        return result
            

    def _numClusters(self):
        nClusters = 1
        for t in self._thresholds:
            nClusters *= len(t) - 1
        return nClusters

    
    def assignClusters(self):
        '''
        Assign clusters to data.

        Create an array of the same size as items from the ``data`` list, and
        fill it with cluster IDs, so that one can easily identify which cluster
        each data item belongs to.

        **Returns:**
            An array of the same size as all the items in the ``data`` list,
            filled with cluster numbers. Unclustered data points will be NaNs.

        .. todo::

            Check that ``len(data)`` > 0. In the constructor of course
        '''
        result = np.zeros(len(self._data[0])) * np.nan
        for cIdx, c in enumerate(self):
            result[c] = cIdx
        return result


    ##########################################################################
    # collections.Sequence method implementations
    def __getitem__(self, key):
        return self._clusters[key]

    def __len__(self):
        return self._nClusters

    def __repr__(self):
        return self._clusters.__repr__()

    def __str__(self):
        ret = ""
        for idx, c in enumerate(self):
            ret += "Cluster {0}: {1}\n".format(idx, c)
        return ret

