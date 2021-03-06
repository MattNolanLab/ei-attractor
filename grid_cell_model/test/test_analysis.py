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
.. currentmodule: test.test_analysis

The :mod:`~test.test_analysis` module defines a set of classes for unit
testing of the most important functions and classes of this package.  The
module is currently based on the unittest module.

A list of currently supported tests:

.. autosummary::
    TestCorrelation
    TestPopulationSpikes
'''
import unittest
import collections
import numpy as np
import sys
import warnings
req_version = (2, 7)
if (sys.version_info >= req_version):
    from unittest.util import safe_repr
else:
    safe_repr = lambda x: x
    warnings.warn('It is better to use unittest from python >= 2.7. Consider upgrading.')


import analysis.signal as asignal
import analysis.spikes as aspikes


notImplMsg = "Not implemented"

##############################################################################
# Spike analysis tests (analysis.spikes)
def _computeOutputSequence(train1, train2):
    res = np.array([])
    for t1 in train1:
        res = np.hstack((res, train2 - t1))
    return res

def _createTestSequence(trainSize, N):
    '''
    Create a test sequence of ``N`` spike trains with exactly ``trainSize``
    number of spikes. Spike times are random, but **time sorted**.
    '''
    senders = np.repeat(np.arange(N), trainSize)
    np.random.shuffle(senders)
    times   = np.random.rand(trainSize * N)
    times.sort()
    sp      = aspikes.PopulationSpikes(N, senders, times)
    return senders, times, sp


class TestPopulationSpikes(unittest.TestCase):
    '''
    Unit tests of :class:`analysis.spikes.PopulationSpikes`.
    '''

    def test_negative_N(self):
        self.assertRaises(ValueError, aspikes.PopulationSpikes, -10, [], [])

    def test_zero_N(self):
        aspikes.PopulationSpikes(0, [], [])


    @unittest.skip(notImplMsg)
    def testAvgFiringRate(self):
        pass


    @unittest.skip(notImplMsg)
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


class TestISI(unittest.TestCase):

    def test_empty(self):
        # empty spike trains
        N       = 100
        for nSpikes in [0, 1]:
            senders, times, sp = _createTestSequence(nSpikes, N)
            res = sp.ISI()
            self.assertEqual(len(res), N)
            for ISIs in res:
                self.assertEqual(len(ISIs), 0)


    def test_results_length(self):
        trainSize = 101
        N       = 100
        senders, times, sp = _createTestSequence(trainSize, N)
        res = sp.ISI()
        for ISIs in res:
            self.assertEqual(len(ISIs), trainSize-1)


    def test_positive(self):
        trainSize = 1000
        N         = 100
        senders, times, sp = _createTestSequence(trainSize, N)
        res = sp.ISI()
        for ISIs in res:
            self.assertTrue(np.all(ISIs >=0))


    def test_constant_ISI(self):
        '''
        .. todo::

            this will only work if dt = 2^x. For now it should be enough.
        '''
        maxSize = 1011
        dt = 0.25
        for trainSize in xrange(2, maxSize):
            senders = [0] * trainSize
            times   = np.arange(trainSize, dtype=float) * dt
            sp = aspikes.PopulationSpikes(1, senders, times)
            res = sp.ISI(n=0)
            self.assertTrue(np.all(res[0] == dt))



class TestISICV(unittest.TestCase):

    def test_empty(self):
        # empty spike trains
        N       = 137
        for nSpikes in [0, 1]:
            senders, times, sp = _createTestSequence(nSpikes, N)
            res = sp.ISICV()
            self.assertEqual(len(res), N)


    def test_results_length(self):
        trainSize = 101
        N       = 137
        senders, times, sp = _createTestSequence(trainSize, N)
        res = sp.ISICV()
        self.assertEqual(len(res), N)
        for CV in res:
            self.assertTrue(isinstance(CV, int) or isinstance(CV, float))


    def test_positive(self):
        trainSize = 10
        N         = 137
        senders, times, sp = _createTestSequence(trainSize, N)
        res = sp.ISICV()
        self.assertTrue(np.all(np.asarray(res >= 0)))
