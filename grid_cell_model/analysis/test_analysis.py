#!/usr/bin/env python
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
import numpy as np
from unittest.util import safe_repr

import analysis.signal as asignal


class TestCorrelation(ut.TestCase):
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


    def checkBlitzNumpyCorr(self, a1, a2, mode):
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


    @ut.skip("Not implemented")
    def test_onesided(self):
        pass


    @ut.skip("Not implemented")
    def test_twosided(self):
        pass

    
    @ut.skip("Not implemented")
    def test_range(self):
        pass

    def test_zero_len(self):
        a1 = np.array([])
        a2 = np.arange(10)

        # corr(a1, a2)
        lag_start = 0
        lag_end   = 0
        for mode in ("onesided", "twosided", "range"):
            self.assertRaises(TypeError, asignal.corr, (a1, a2, mode, lag_start,
                lag_end))
            self.assertRaises(TypeError, asignal.corr, (a2, a1, mode, lag_start,
                lag_end))
    

