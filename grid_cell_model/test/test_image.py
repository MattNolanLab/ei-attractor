'''
.. currentmodule: test.test_image

The :mod:`~test.test_image` module defines a set of classes for unit testing of
the most important functions and classes of the :mod:`analysis.image` module.
The module is currently based on the unittest module.

A list of currently supported tests:

.. autosummary::
'''
import unittest
import numpy as np

import analysis.image as aimage


notImplMsg = "Not implemented"


##############################################################################
# Gaussian fitting on the twisted torus

def generateGaussianTT(A, mu_x, mu_y, sigma, X, Y):
    dim = aimage.Position2D(X.shape[1], X.shape[0])
    a = aimage.Position2D(mu_x, mu_y)
    others = aimage.Position2D(X.ravel(), Y.ravel())
    dist = aimage.remapTwistedTorus(a, others, dim)
    G = A * np.exp( -dist**2 / (2*sigma**2)  )
    return np.reshape(G, (dim.y, dim.x))


class Test_fitgaussianTT(unittest.TestCase):
    decimalAlmostEqual = 1
    gaussianAMax = 40
    nIter = 1000
    minSigma = 1.
    maxFailures = .02

    def assertNdarrayAlmostEqual(self, first, second, msg=None):
        #print first, second
        np.testing.assert_almost_equal(first, second,
                self.decimalAlmostEqual)


    def setUp(self):
        self.addTypeEqualityFunc(np.ndarray, self.assertNdarrayAlmostEqual)


    def test_bumpFitting(self):
        dim = aimage.Position2D(34, 30)
        X, Y = np.meshgrid(np.arange(dim.x), np.arange(dim.y))
        failures = 0
        for it in xrange(self.nIter):
            A     = np.random.rand() * self.gaussianAMax
            mu_x  = np.random.rand() * dim.x
            mu_y  = np.random.rand() * dim.y
            sigma = self.minSigma + np.random.rand() * (np.min((dim.x, dim.y)) / 4 - self.minSigma)

            ourG = generateGaussianTT(A, mu_x, mu_y, sigma, X, Y) + \
                    np.random.randn(dim.y, dim.x)*1e-3*A
            estParams, err = aimage.fitGaussianBumpTT(ourG, dim)
            estParams[0] = np.abs(estParams[0])
            estParams[3] = np.abs(estParams[3])

            try:
                first = np.array([A, mu_x, mu_y, sigma])
                second = np.array(estParams)
                self.assertNdarrayAlmostEqual(first, second)
            except AssertionError as e:
                failures += 1
                #print failures
                #print('%s != \n%s within %r places' % (first, second,
                #    self.decimalAlmostEqual))
                if failures / float(self.nIter) > self.maxFailures:
                    msg = '%.1f%%  fitting errors reached.' % (self.maxFailures*100)
                    raise self.failureException(msg)
                






