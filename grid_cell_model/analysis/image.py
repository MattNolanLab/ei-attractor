'''
Classes and functions for analysis of image-like structures.

.. currentmodule:: analysis.image

Classes
-------

.. autosummary::

    SingleBumpPopulation
    Position2D
    MLGaussianFit
    LikelihoodGaussianInitParams



Functions
---------

.. autosummary::

    remapTwistedTorus
    fitGaussianTT
    fitGaussianBumpTT


'''
'''
Image analysis fitting - rateMaps, generic image analysis for GridCells, etc.
'''
from abc import ABCMeta, abstractmethod
import collections
import logging

from scipy    import weave

import numpy as np
import scipy.optimize

from . import spikes


logger = logging.getLogger(__name__)


class Position2D(object):

    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y

        if not isinstance(x, float) and isinstance(y, float):
            if len(x) != len(y):
                msg = "Position2D: x and y parameters must have the same length: %d != %d."
                logger.warn(msg, len(x), len(y))

    def __str__(self):
        return "Position2D: x: " + str(self.x) + ", y: " + str(self.y)

    def __len__(self):
        return len(self.x)



class FittingParams(object):
    __meta__ = ABCMeta

    @abstractmethod
    def __init__(self):
        raise NotImplementedError()



class SymmetricGaussianParams(FittingParams):
    def __init__(self, A, mu_x, mu_y, sigma, err2):
        self.A = A
        self.mu_x = mu_x
        self.mu_y = mu_y
        self.sigma = sigma
        self.err2 = err2




##############################################################################
# Simple ML solutions and lists
class MLFit(FittingParams):
    def __init__(self, mu, sigma2, ln_L, err2):
        self.mu     = mu
        self.sigma2 = sigma2
        self.ln_L   = ln_L
        self.err2    = err2


class MLFitList(MLFit, collections.Sequence):
    def __init__(self, mu=None, sigma2=None, ln_L=None, err2=None, times=None):
        if mu     is None: mu = []
        if sigma2 is None: sigma2 = []
        if ln_L   is None: ln_L = []
        if err2   is None: err2 = []
        if times  is None: times = []
        super(MLFitList, self).__init__(mu, sigma2, ln_L, err2)
        self.times = times

        if not self._consistent():
            raise ValueError('All input arguments mus have same length')

    def _consistent(self):
        return len(self.mu) == len(self.sigma2) and \
               len(self.mu) == len(self.ln_L) and   \
               len(self.mu) == len(self.err2) and    \
               len(self.mu) == len(self.times)

    def __getitem__(self, key):
        return MLFit(self.mu[key],
                     self.sigma2[key],
                     self.ln_L[key],
                     self.err2), \
               times


    def __len__(self):
        return len.self.mu


    def _appendData(self, d, t):
        '''`d` must be an instance of :class:`MLFit`'''
        if not isinstance(d, MLFit):
            raise TypeError('ML data must be an instance of MLFit')
        self.mu.append(d.mu)
        self.sigma2.append(d.sigma2)
        self.ln_L.append(d.ln_L)
        self.err2.append(d.err2)
        self.times.append(t)



##############################################################################
# Symmetric Gaussian ML solutions and lists
class MLGaussianFit(SymmetricGaussianParams):
    def __init__(self, A, mu_x, mu_y, sigma, err2, ln_L, lh_precision):
        super(MLGaussianFit, self).__init__(A, mu_x, mu_y, sigma, err2)
        self.ln_L = ln_L
        self.lh_precision = lh_precision


class MLGaussianFitList(MLGaussianFit, collections.Sequence):
    def __init__(self, A=None, mu_x=None, mu_y=None, sigma=None, err2=None,
            ln_L=None, lh_precision=None, times=None):
        if A            is None: A     = []
        if mu_x         is None: mu_x  = []
        if mu_y         is None: mu_y  = []
        if sigma        is None: sigma = []
        if err2         is None: err2  = []
        if ln_L         is None: ln_L  = []
        if lh_precision is None: lh_precision = []
        if times        is None: times = []

        super(MLGaussianFitList, self).__init__(\
                A, mu_x, mu_y, sigma, err2, ln_L, lh_precision)
        self.times = times

        if not self._consistent():
            raise ValueError('All input arguments mus have same length')


    def _consistent(self):
        return \
            len(self.A) == len(self.mu_x) and         \
            len(self.A) == len(self.mu_y) and         \
            len(self.A) == len(self.sigma) and        \
            len(self.A) == len(self.err2) and          \
            len(self.A) == len(self.ln_L) and         \
            len(self.A) == len(self.lh_precision) and \
            len(self.A) == len(self.times)


    def _appendData(self, d, t):
        '''`d` must be an instance of :class:`MLGaussianFit`'''
        if not isinstance(d, MLGaussianFit):
            raise TypeError('Data must be an instance of MLGaussianFit')

        self.A.append(d.A)
        self.mu_x.append(d.mu_x)
        self.mu_y.append(d.mu_y)
        self.sigma.append(d.sigma)
        self.err2.append(d.err2)
        self.ln_L.append(d.ln_L)
        self.lh_precision.append(d.lh_precision)
        self.times.append(t)


    def __getitem__(self, key):
        return MLGaussianFit(self.A[key],
                             self.mu_x[key],
                             self.mu_y[key],
                             self.sigma[key],
                             self.err2[key],
                             self.ln_L,
                             self.lh_precision), \
                self._times[key]

    def __len__(self):
        return len(self.A) # All same length



def remapTwistedTorus(a, others, dim):
    ''' Calculate a distance between ``a`` and ``others`` on a twisted torus.
    
    Take ``a`` which is a 2D position and others, which is a vector of 2D
    positions and compute the distances between them based on the topology of
    the twisted torus.
    
    If you just want to remap a function of (X, Y), set a==[[0, 0]].

    Parameters
    ----------
    
    a : a Position2D instance
        Specifies the initial position. ``a.x`` and ``a.y`` must be convertible
        to floats
    others : Position2D instance
        Positions for which to compute the distance.
    dim : Position2D
        Dimensions of the torus. ``dim.x`` and ``dim.y`` must be convertible to
        floats.

    Returns
    -------
    An array of positions, always of the length of others
    '''
    

    a_x      = float(a.x)
    a_y      = float(a.y)
    others_x = np.asarray(others.x)
    others_y = np.asarray(others.y)
    szO      = others.x.shape[0]
    x_dim    = float(dim.x)
    y_dim    = float(dim.y)
    ret      = np.ndarray((szO,))

    # Remap the values modulo torus size.
    a_x = a_x % x_dim
    a_y = a_y % y_dim
    others_x = others_x % x_dim
    others_y = others_y % y_dim

    code = '''
    #define SQ(x) ((x) * (x))
    #define MIN(x1, x2) ((x1) < (x2) ? (x1) : (x2))

    for (int i = 0; i < szO; i++)
    {
        double o_x = others_x(i);
        double o_y = others_y(i);

        double d1 = sqrt(SQ(a_x - o_x            ) + SQ(a_y - o_y        ));
        double d2 = sqrt(SQ(a_x - o_x - x_dim    ) + SQ(a_y - o_y        ));
        double d3 = sqrt(SQ(a_x - o_x + x_dim    ) + SQ(a_y - o_y        ));
        double d4 = sqrt(SQ(a_x - o_x + 0.5*x_dim) + SQ(a_y - o_y - y_dim));
        double d5 = sqrt(SQ(a_x - o_x - 0.5*x_dim) + SQ(a_y - o_y - y_dim));
        double d6 = sqrt(SQ(a_x - o_x + 0.5*x_dim) + SQ(a_y - o_y + y_dim));
        double d7 = sqrt(SQ(a_x - o_x - 0.5*x_dim) + SQ(a_y - o_y + y_dim));

        ret(i) = MIN(d7, MIN(d6, MIN(d5, MIN(d4, MIN(d3, MIN(d2, d1))))));
    }
    '''
    
    weave.inline(code,
        ['others_x', 'others_y', 'szO', 'ret', 'a_x', 'a_y', 'x_dim', 'y_dim'],
        type_converters=weave.converters.blitz,
        compiler='gcc',
        extra_compile_args=['-O3'],
        verbose=2)

    return ret



##############################################################################
#                      Image analysis/manipulation functions
##############################################################################



def fitGaussianTT(sig_f, i):
    '''Fit a 2D circular Gaussian function to a 2D signal using a maximum
    likelihood estimator.

    The Gaussian is not generic: :math:`\sigma_x = \sigma_y = \sigma`, i.e.
    it is circular only.

    The function fitted looks like this:
    .. math::

        f(\mathbf{X}) = |A| \exp\left{\frac{-||X - mu||^2}{2*\sigma^2}\right},

    where :math:`||\cdot||` is a distance metric on the twisted torus.

    Paramters
    ---------
    sig_f : np.ndarray
        A 2D array that specified the signal to fit the Gaussian onto. The
        dimensions of the torus will be inferred from the shape of `sig_f`:
        (dim.y, dim.x) = `sig_f.shape`.
    i : SymmetricGaussianParams
        Guassian initialisation parameters. The `err2` field will be ignored.

    Returns
    -------
    :class:`MLGaussianFit`
        Estimated values, together with maximum likelihood value and precision
        (inverse variance of noise: *NOT* of the fitted Gaussian).
    '''
    # Fit the Gaussian using least squares
    f_flattened = sig_f.ravel()
    dim         = Position2D(sig_f.shape[1], sig_f.shape[0])
    X, Y        = np.meshgrid(np.arange(dim.x), np.arange(dim.y))
    others      = Position2D()
    others.x    = X.flatten()
    others.y    = Y.flatten()

    a = Position2D()
    def gaussDiff(x):
        a.x = x[1] # mu_x
        a.y = x[2] # mu_y
        dist = remapTwistedTorus(a, others, dim)
        return np.abs(x[0]) * np.exp( -dist**2/2./ x[3]**2 ) - f_flattened
#                       |                            |
#                       A                          sigma

    x0 = np.array([i.A, i.mu_x, i.mu_y, i.sigma])
    xest, ierr = scipy.optimize.leastsq(gaussDiff, x0)
    err2 = gaussDiff(xest)**2

    # Remap the values modulo torus size
    xest[1] = xest[1] % dim.x
    xest[2] = xest[2] % dim.y

    # Compute the log-likelihood
    N = dim.x * dim.y
    AIC_correction = 5 # Number of optimized parameters
    beta = 1.0 / ( np.mean(err2) )
    ln_L = -beta / 2. * np.sum(err2) +  \
            N / 2. * np.log(beta) -     \
            N / 2. * np.log(2*np.pi) -  \
            AIC_correction

    res = MLGaussianFit(xest[0], xest[1], xest[2], xest[3], err2, ln_L, beta)
    return res



def fitGaussianBumpTT(sig):
    '''Fit a 2D Gaussian onto a (potential) firing rate bump on the twisted torus.

    Parameters
    ----------
    sig : np.ndarray
        2D firing rate map to fit. Axis 0 is the Y position. This will be
        passed directly to :func:`~analysis.image.fitGaussianTT`.

    Returns
    -------
    :class:`analysis.image.MLGaussianFit`
        Estimated values of the fit.

    Notes
    -----
    The function initialises the Gaussian fitting parameters to a position at
    the maximum of `sig`.
    '''
    mu0_y, mu0_x = np.unravel_index(np.argmax(sig), sig.shape)
    A0           = sig[mu0_y, mu0_x]
    sigma0       = np.max(sig.shape) / 4.
    init = SymmetricGaussianParams(A0, mu0_x, mu0_y, sigma0, None)
    return fitGaussianTT(sig, init)



def fitMaximumLikelihood(sig):
    '''Fit a maximum likelihood solution under Gaussian noise.

    Parameters
    ----------
    sig : np.ndarray
        A vector containing the samples

    Returns
    fit : MLFit
        Maximum likelihood parameters
    '''
    sig = sig.flatten()
    mu = np.mean(sig)
    sigma2 = np.var(sig)

    N = len(sig)
    AIC_correction = 2
    ln_L = -.5 / sigma2 * np.sum((sig - mu)**2) - \
            .5 * N * np.log(sigma2) -             \
            .5 * N * np.log(2*np.pi) -            \
            AIC_correction
    err2 = (sig - mu)**2

    return MLFit(mu, sigma2, ln_L, err2)




class SingleBumpPopulation(spikes.TwistedTorusSpikes):
    '''
    A population of neurons that is supposed to form a bump on a twisted torus.

    This class contains methods for processing  the population activity over
    time.
    '''
    def __init__(self, senders, times, sheetSize):
        super(SingleBumpPopulation, self).__init__(senders, times, sheetSize)


    def _performFit(self, tStart, tEnd, dt, winLen, fitCallable, listCls, fullErr=True):
        F, Ft = self.slidingFiringRate(tStart, tEnd, dt, winLen)
        dims = Position2D(self.Nx, self.Ny)
        res = listCls()
        for tIdx in xrange(len(Ft)):
            logger.debug('%s:: fitting: %d/%d, %.3f/%.3f ',
                    fitCallable.__name__, tIdx+1, len(Ft), Ft[tIdx], Ft[-1])
            fitParams = fitCallable(F[:, :, tIdx])
            
            if fullErr == False:
                fitParams.err2 = np.sum(fitParams.err2)

            res._appendData(fitParams, Ft[tIdx])
        return res


    def bumpPosition(self, tStart, tEnd, dt, winLen, fullErr=True):
        '''Estimate bump positions during the simulation time.
        
            1. Use :py:meth:`~.slidingFiringRate`

            2. Apply the bump position estimation procedure to each of the
               population activity items.

        Parameters
        ----------
            tStart, tEnd, dt, winLen
                As in :py:meth:`~analysis.spikes.slidingFiringRate`.
            fullErr : bool
                If ``True``, save the full error of fit. Otherwise a sum only.

        Returns
        -------
        MLGaussianFitList
            A list of fitted Gaussian parameters

        Notes
        -----
        This method uses the Maximum likelihood estimator to fit the Gaussian
        function (:meth:`~analysis.image.fitGaussianBumpTT`)
        '''
        return self._performFit(tStart, tEnd, dt, winLen, fitGaussianBumpTT,
                MLGaussianFitList, fullErr=fullErr)
            

    def uniformFit(self, tStart, tEnd, dt, winLen, fullErr=True):
        '''Estimate the mean firing rate using maximum likelihood estimator
        (:func:`~analysis.image.fitMaximumLikelihood`)
        
            1. Use :py:meth:`~.slidingFiringRate`.

            2. Apply the estimator.

        Parameters
        ----------
            tStart, tEnd, dt, winLen
                As in :py:meth:`~analysis.spikes.slidingFiringRate`.
            fullErr : bool
                If ``True``, save the full error of fit. Otherwise a sum only.

        Returns
        -------
        MLFitList
            A list of fitted parameters.
        '''
        return self._performFit(tStart, tEnd, dt, winLen, fitMaximumLikelihood,
                MLFitList, fullErr=fullErr)




##############################################################################
#                           Tests
##############################################################################

def remapAndPlot(a, dims):
    X, Y = np.meshgrid(np.arange(dims.x), np.arange(dims.y))
    others = Position2D()
    others.x = X.flatten()
    others.y = Y.flatten()
    dist = remapTwistedTorus(a, others, dims)
    plt.figure()
    plt.pcolormesh(np.reshape(dist, (dims.y, dims.x)))
    plt.title('a:(x, y): ' + str(a.x) + ', ' + str(a.y) + ', x_dim:' + str(dims.x) +
            ', y_dim: ' + str(dims.y))
    plt.axis('equal')


if (__name__ == '__main__'):

    import matplotlib.pyplot as plt
    from matplotlib.pyplot import *

    dims = Position2D()
    dims.x = 34
    dims.y = 34
    a = Position2D()
    a.x = 0; a.y = 0
    remapAndPlot(a, dims)

    a.x = dims.x/2.
    a.y = 0
    remapAndPlot(a, dims)

    a.x = dims.x - 1;
    a.y = 0
    remapAndPlot(a, dims)

    a.x = 0
    a.y = dims.y/2.0
    remapAndPlot(a, dims)

    a.x = dims.x/2.0
    a.y = dims.y/2.0
    remapAndPlot(a, dims)

    a.x = dims.x - 1
    a.y = dims.y/2.0
    remapAndPlot(a, dims)

    a.x = 0
    a.y = dims.y - 1
    remapAndPlot(a, dims)

    a.x = dims.x/2.0
    a.y = dims.y - 1
    remapAndPlot(a, dims)

    a.x = dims.x - 1
    a.y = dims.y - 1
    remapAndPlot(a, dims)

    show()

