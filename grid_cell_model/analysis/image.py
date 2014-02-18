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

__all__ = ['Position2D', 'remapTwistedTorus', 'fitGaussianTT',
        'fitGaussianBumpTT']

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



class SymmetricGaussianFit(FittingParams):
    def __init__(self, A, mu_x, mu_y, sigma, err):
        self.A = A
        self.mu_x = mu_x
        self.mu_y = mu_y
        self.sigma = sigma
        self.err = err



class GaussianInitParams(object):
    def __init__(self):
        self.A0      = 10.0
        self.mu0_x   = 1.0
        self.mu0_y   = 1.0
        self.sigma0  = 1.0




def remapTwistedTorus(a, others, dim):
    '''
    Calculate a distance between ``a`` and ``others`` on a twisted torus.
    
    Take ``a`` which is a 2D position and others, which is a vector of 2D
    positions and compute the distances between them based on the topology of
    the twisted torus.
    
    If you just want to remap a function of (X, Y), set a==[[0, 0]].

    **Parameters**
    
    a : a Position2D instance
        Specifies the initial position. ``a.x`` and ``a.y`` must be convertible
        to floats
    others : Position2D instance
        Positions for which to compute the distance.
    dim : Position2D
        Dimensions of the torus. ``dim.x`` and ``dim.y`` must be convertible to
        floats.

    **Returns**

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



## Fit a 2D Gaussian function to a 2D signal using a least squares method.
#
# The Gaussian is not generic - sigmax = sigmay = sigma, i.e. it is circular
# only.
#
# The function fitted looks like this:
#     fun = A*exp(-||X - mu||^2 / (2*sigma^2))
#
# The initialisation parameters are A, mu_x, mu_y, sigma
#
#
# @param sig_f  Signal function value - this is fitted, use the Position2D class
# @param i      A struct with initialisation values
# @param dim    Dimensions of the twisted torus (Position2D)
# @return Estimated values in the order specified above for A...sigma
#
def fitGaussianTT(sig_f, i, dim):
    X, Y = np.meshgrid(np.arange(dim.x), np.arange(dim.y))
    others = Position2D()
    others.x = X.flatten()
    others.y = Y.flatten()

    a = Position2D()
    def gaussDiff(x):
        a.x = x[1] # mu_x
        a.y = x[2] # mu_y
        dist = remapTwistedTorus(a, others, dim)
        #dist = np.sqrt((others.x - a.x)**2 + (others.y - a.y)**2)
        #print "A:", x[0], a, "sigma:", x[3]
        return np.abs(x[0]) * np.exp( -dist**2/2./ x[3]**2 ) - sig_f
#                       |                            |
#                       A                          sigma

    x0 = np.array([i.A0, i.mu0_x, i.mu0_y, i.sigma0])

    xest,ierr = scipy.optimize.leastsq(gaussDiff, x0)
    # Remap the values module torus size
    xest[1] = xest[1] % dim.x
    xest[2] = xest[2] % dim.y
    err = gaussDiff(xest)
    return xest, err**2



## Fit a 2D Gaussian onto a (potential) firing rate bump. On Twisted torus
#
# @param sig    Firing rate map to fit - a 2D numpy array
# @param dim    Dimensions of the twisted torus (Position2D)
# @return Estimated values for the Gaussian (see Position2D)
def fitGaussianBumpTT(sig, dim):
    '''
    Fit a Gaussian to a rate map, using least squares method.
    '''
    init = GaussianInitParams()
    init.mu0_y, init.mu0_x  = np.unravel_index(np.argmax(sig), sig.shape)
    init.A0 = sig[init.mu0_y, init.mu0_x]
    init.sigma0  = np.max([dim.x, dim.y]) / 2.0
    
    return fitGaussianTT(sig.flatten(), init, dim)






class SingleBumpPopulation(spikes.TwistedTorusSpikes):
    '''
    A population of neurons that is supposed to form a bump on a twisted torus.

    This class contains methods for processing  the population activity over
    time.
    '''
    class BumpFitList(collections.Sequence):
        def __init__(self, AList=None, mu_xList=None, mu_yList=None,
                sigmaList=None, err=None, times=None):
            self._A     = AList     if AList is not None else []
            self._mu_x  = mu_xList  if mu_xList is not None else []
            self._mu_y  = mu_yList  if mu_yList is not None else []
            self._sigma = sigmaList if sigmaList is not None else []
            self._err   = err       if err is not None else []
            self._times = times     if times is not None else []

            if not self._consistent():
                raise ValueError('All input arguments mus have same length')


        def _consistent(self):
            return \
                len(self._A) == len(self._mu_x) and   \
                len(self._A) == len(self._mu_y) and   \
                len(self._A) == len(self._sigma) and  \
                len(self._A) == len(self._err) and    \
                len(self._A) == len(self._times)


        def getA(self): return self._A
        def getMu_x(self): return self._mu_x
        def getMu_y(self): return self._mu_y
        def getSigma(self): return self._sigma
        def getErr(self): return self._err
        def getT(self): return self._times
        A     = property(getA)
        mu_x  = property(getMu_x)
        mu_y  = property(getMu_y)
        sigma = property(getSigma)
        err   = property(getErr)
        t     = property(getT)


        def _appendData(self, A, mu_x, mu_y, sigma, err, t):
            self._A.append(A)
            self._mu_x.append(mu_x)
            self._mu_y.append(mu_y)
            self._sigma.append(sigma)
            self._err.append(err)
            self._times.append(t)


        def __getitem__(self, key):
            return SymmetricGaussianFit(self._A[key], self._mu_x[key],
                    self._mu_y[key], self._sigma[key], self._err[key]), \
                            self._times[key]

        def __len__(self):
            return len(self._A) # All same length


    def __init__(self, senders, times, sheetSize):
        super(SingleBumpPopulation, self).__init__(senders, times, sheetSize)



    def bumpPosition(self, tStart, tEnd, dt, winLen, fullErr=True):
        '''
        Estimate bump positions during the simulation time:
        
            1. Use :py:meth:`~.slidingFiringRate`

            2. Apply the bump position estimation procedure to each of the
               population activity items.

        **Parameters**

            tStart, tEnd, dt, winLen : as in :py:meth:`~.slidingFiringRate`.
            fullErr : bool
                If ``True``, save the full error of fit. Otherwise a sum only.

        **Output*
            
            A list of estimation result objects
        '''
        F, Ft = self.slidingFiringRate(tStart, tEnd, dt, winLen)
        dims = Position2D(self.Nx, self.Ny)
        res = self.BumpFitList()
        for tIdx in xrange(len(Ft)):
            logger.debug('Bump fitting: %s/%s', tIdx+1, len(Ft))
            (A, mu_x, mu_y, sigma), err2 = fitGaussianBumpTT(F[:, :, tIdx],
                    dims)
            
            err = np.sqrt(err2)
            if fullErr == False:
                err = np.sum(err)

            res._appendData(A, mu_x, mu_y, sigma, err, Ft[tIdx])
        return res
            





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


    ###################################################################### 
    # Gaussian bump fitting on a twisted torus
    dims = Position2D()
    dims.x = 100
    dims.y = 100

    X, Y = np.meshgrid(np.arange(dims.x), np.arange(dims.y))

    pos = Position2D()
    pos.x = X.flatten()
    pos.y = Y.flatten()

    A       = 20.0
    mu_x    = 99
    mu_y    = 50.0
    sigma   = 20

    a0 = Position2D()
    a0.x = mu_x
    a0.y = mu_y
    dist = remapTwistedTorus(a0, pos, dims)
    rateMap = A*np.exp(- dist**2 / (2*sigma**2))
    rateMap += A*0.2*np.random.randn(rateMap.shape[0])
    rateMap[rateMap < 0] = 0
    rateMap = np.reshape(rateMap, (dims.y, dims.x))
  
    figure()
    pcolormesh(np.reshape(dist, (dims.y, dims.x))); colorbar()
    figure()
    pcolormesh(X, Y, rateMap); colorbar()
  
    param_est = fitGaussianBumpTT(rateMap, dims)
    Aest, muxest, muyest, sigmaest = param_est

    print param_est
    
    show()
