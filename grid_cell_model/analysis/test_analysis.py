#
#   test_analysis.py
#
#   Unit tests for this package.
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
.. currentmodule: analysis.test_analysis

The :mod:`analysis.test_analysis` module defines a set of classes for unit
testing of the most important functions and classes of this package.  The
module is currently based on the unittest module.

A list of currently supported tests:

.. autosummary::
    TestCorrelation
    TestPopulationSpikes
'''
import unittest as ut
import collections
import numpy as np
from unittest.util import safe_repr

import analysis.signal as asignal
import analysis.spikes as aspikes

notImplMsg = "Not implemented"


##############################################################################
# Continuous signal tests (analysis.signal)

class TestCorrelation(ut.TestCase):
    '''
    Test the analysis.signal.corr function (and effectively the core of the
    autoCorrelation) function.
    '''
    def setUp(self):
        self.places = 10
        self.maxN = 1000
        self.maxLoops = 1000

    def assertSequenceAlmostEqual(self, first, second, places=None, msg=None, delta=None):
        """
        Fail if the two objects are unequal as determined by the difference
        between all of their values, rounded to the given number of decimal
        places (default 7) and comparing to zero, or by comparing that the
        differences between any two items of the two objects is more than the
        given delta.

        The test will fail if the conditions for any of the elements are not
        met.

        Note that decimal places (from zero) are usually not the same
        as significant digits (measured from the most signficant digit).

        If the two objects compare equal then they will automatically
        compare almost equal.
        """
        if delta is not None and places is not None:
            raise TypeError("specify delta or places not both")

        if delta is not None:
            if np.all(np.abs(first - second) <= delta):
                return

            standardMsg = '%s != %s within %s delta' % \
                    (safe_repr(first),
                     safe_repr(second),
                     safe_repr(delta))
        else:
            if places is None:
                places = 7

            if np.all(np.round(np.abs(second-first), places) == 0):
                return

            standardMsg = '%s != %s within %r places' % \
                    (safe_repr(first),
                    safe_repr(second),
                    places)
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)


    def checkBlitzNumpyCorr(self, a1, a2):
        c_blitz = asignal.corr(a1, a2, mode='twosided')
        c_np    = np.correlate(a1, a2, mode='full')[::-1]
        self.assertSequenceAlmostEqual(c_blitz, c_np, places=self.places) 
        

    def test_blitz_numpy(self):
        for loop_it in xrange(self.maxLoops):
            N1 = np.random.randint(self.maxN) + 1
            N2 = np.random.randint(self.maxN) + 1
            if (N1 == 0 and N2 == 0):
                continue

            a1 = np.random.rand(N1)
            a2 = np.random.rand(N2)
            
            self.checkBlitzNumpyCorr(a1, a2)


    @ut.skip(notImplMsg)
    def test_onesided(self):
        pass


    @ut.skip(notImplMsg)
    def test_twosided(self):
        pass

    
    @ut.skip(notImplMsg)
    def test_range(self):
        pass

    def test_zero_len(self):
        a1 = np.array([])
        a2 = np.arange(10)

        # corr(a1, a2)
        lag_start = 0
        lag_end   = 0
        for mode in ("onesided", "twosided", "range"):
            self.assertRaises(TypeError, asignal.corr, a1, a2, mode, lag_start,
                lag_end)
            self.assertRaises(TypeError, asignal.corr, a2, a1, mode, lag_start,
                lag_end)
            self.assertRaises(TypeError, asignal.corr, a1, a1, mode, lag_start,
                lag_end)
    


##############################################################################
# Spike analysis tests (analysis.spikes)
def _computeOutputSequence(train1, train2):
    res = np.array([])
    for t1 in train1:
        res = np.hstack((res, train2 - t1))
    return res

def _createTestSequence(trainSize, N):
    senders = np.random.randint(N, size=trainSize * N)
    times   = np.random.rand(trainSize * N)
    sp      = aspikes.PopulationSpikes(N, senders, times)
    return senders, times, sp


class TestPopulationSpikes(ut.TestCase):
    '''
    Unit tests of :class:`analysis.spikes.PopulationSpikes`.
    '''

    def test_negative_N(self):
        self.assertRaises(ValueError, aspikes.PopulationSpikes, -10, [], [])

    def test_zero_N(self):
        aspikes.PopulationSpikes(0, [], [])
    

    @ut.skip(notImplMsg)
    def testAvgFiringRate(self):
        pass


    @ut.skip(notImplMsg)
    def test_slidingFiringRate(self):
        pass


    
    def test_lists(self):
        trainSize = 100
        N = 10
        senders = list(np.random.randint(N, size=trainSize * N))
        times   = list(np.random.rand(trainSize * N))
        sp      = aspikes.PopulationSpikes(N, senders, times)

        # try to retrieve spike trains
        for nIdx in xrange(N):
            train = sp[nIdx]

        # Try to run all the methods, None should raise an exception
        sp.avgFiringRate(0, 1)
        sp.slidingFiringRate(0, 1, 0.05, 0.1)
        sp.windowed((0, 1))
        sp.rasterData()
        sp.spikeTrainDifference(range(N))




class TestSpikeTrainDifference(ut.TestCase):

    def test_full(self):
        # full must be True
        N       = 10
        senders, times, sp = _createTestSequence(0, N)
        std     = sp.spikeTrainDifference
        self.assertRaises(NotImplementedError, std, 1, 10, False)

    def test_empty(self):
        # Empty spike trains
        N = 10
        senders, times, sp = _createTestSequence(0, N)
        std = sp.spikeTrainDifference
        res = std(1, 2, True)
        self.assertIsInstance(res, collections.Sequence)
        self.assertEqual(len(res), 1) #
        self.assertEqual(res[0].shape[0], 0)

    def test_out_of_bounds(self):
        # Out of bounds spike trains
        trainSize = 100
        N         = 100
        senders, times, sp = _createTestSequence(trainSize, N)
        std       = sp.spikeTrainDifference

        self.assertRaises(IndexError, std, N, 1)
        self.assertRaises(IndexError, std, 1, N)

        # boundaries
        std(N - 1, 1)
        std(1, N-1)
        std(-1, 0)
        std(0, -1)


    def test_result_length(self):
        trainSize = 100
        N         = 100
        senders, times, sp = _createTestSequence(trainSize, N)
        std       = sp.spikeTrainDifference

        # result length must be correct
        trainLengths = [np.count_nonzero(senders == x) for x in xrange(N)]
        res = std(range(N), None, True)
        self.assertEqual(len(res), N)
        for nIdx in xrange(N):
            self.assertEqual(len(res[nIdx]), N)

        for n1 in xrange(N):
            for n2 in xrange(N):
                expectedLen = trainLengths[n1] * trainLengths[n2]
                self.assertEqual(len(res[n1][n2]), expectedLen)

    def test_correct_values(self):
        trainSize = 100
        N         = 50
        senders, times, sp = _createTestSequence(trainSize, N)
        std       = sp.spikeTrainDifference
        res       = std(range(N), None, True)

        for n1 in xrange(N):
            train1 = times[senders == n1]
            for n2 in xrange(N):
                train2 = times[senders == n2]
                diff = res[n1][n2]
                expectedDiff = _computeOutputSequence(train1, train2)
                self.assertTrue(np.all(diff == expectedDiff))


class TestSpikeTrainXCorrelation(ut.TestCase):

    def test_bin_edges(self):
        trainSize = 100
        N = 50
        bins = 37
        senders, times, sp = _createTestSequence(trainSize, N)
        xcf = sp.spikeTrainXCorrelation

        trainLens = [np.count_nonzero(senders == x) for x in xrange(N)]
        res, bin_edges = xcf(range(N), None, (0, 1), bins)
        for n1 in xrange(N):
            for n2 in xrange(N):
                self.assertEqual(len(res[n1][n2]) + 1, len(bin_edges))


    @ut.skip(notImplMsg)
    def test_correct_values(self):
        '''
        Since we are running this on numpy.histogram, it should be ok for these
        purposes.
        '''
        pass

