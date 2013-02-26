#
#   data_analysis.py
#
#   Some useful data analysis functions
#
#       Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
from scipy.io import loadmat
from scipy    import weave
from matplotlib.pyplot import *

import numpy as np
import scipy.optimize

class Position2D(object):

    def __init__(self):
        self.x = None
        self.y = None

    def __str__(self):
        return "Position2D: x: " + str(self.x) + ", y: " + str(self.y)



class GaussianInitParams(object):
    def __init__(self):
        self.A0      = 10.0
        self.mu0_x   = 1.0
        self.mu0_y   = 1.0
        self.sigma0  = 1.0



def circularVariance(x, range):
    ''' Returns circular moments of a vector of circular variable x, defined in
    the 'range' (this will be mapped to a circle)'''
    c = np.exp(1j*2*np.pi*x/range)
    avg = np.mean(c)
    theta_avg = np.angle(avg)
    theta_var = 1 - np.abs(avg)
    return (theta_avg, theta_var)
    

def bumpVariancePos(dirName, fileNamePrefix, jobRange, trialRange, t_start):
    '''t_start - in discrete time'''
    ret = Position2D()

    fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}"
    F_dt = 0.2
    it_start = t_start/F_dt
    
    jobN = jobRange[1] - jobRange[0]+1
    trialN = trialRange[1] - trialRange[0] + 1
    pos_var = []
    pos_N = []
    pos_var_mean = np.ndarray((jobN, 2))
    pos_var_std = np.ndarray((jobN, 2))
    pos_x_vec = np.ndarray((jobN, trialN), dtype=object)
    pos_y_vec = np.ndarray((jobN, trialN), dtype=object)
    o_vec = np.ndarray((jobN, trialN), dtype=object)
    
    for job_it in range(jobN):
        jobNum = job_it + jobRange[0]
        pos_var.append([])
        for trial_it in range(trialN):
            trialNum = trialRange[0] + trial_it
            try:
                fileName = fileNameTemp.format(dirName, fileNamePrefix, jobNum,
                    trialNum)
                print "Processing file: " + fileName + "_output.mat"
                data = loadmat(fileName + "_output.mat")
            except:
                print "Warning: could not open: " + fileName
                continue
    
            o = data['options']
            o_vec[job_it, trial_it] = o
            pos_x = data['bumpPos'][:, 0]
            pos_y = data['bumpPos'][:, 1]
            pos_x_vec[job_it, trial_it] = pos_x
            pos_y_vec[job_it, trial_it] = pos_y
            pos_x = pos_x[it_start:]
            pos_y = pos_y[it_start:]
            
            Ne = o['Ne'][0][0][0][0] # Ugly but these bastards cannot export to
                                     # matlab
            mean_x, var_x = circularVariance(pos_x, Ne)
            mean_y, var_y = circularVariance(pos_y, Ne)
            pos_var[job_it].append([var_x, var_y])
    
        pos_var_mean[job_it, :] = np.mean(pos_var[job_it], 0)
        pos_var_std[job_it, :] = np.std(pos_var[job_it], 0)
        pos_N.append(Ne)
    
    ret.pos_x_vec = pos_x_vec
    ret.pos_y_vec = pos_y_vec
    ret.o_vec = o_vec
    ret.pos_var = pos_var
    ret.pos_var_mean = pos_var_mean
    ret.pos_var_std = pos_var_std
    ret.pos_N = pos_N
    return ret




## Remap a distance between 'a' and others on a twisted torus
#
# Take 'a' which is a 2D position and others, which is a vector of 2D positions
# and compute the distances between them based on the topology of the twisted
# torus, mainly its dimensions.
# 
# If you just want to remap a function of (X, Y), set a==[[0, 0]].
#
# @param a  An array of shape (1, 2) that specifies a position on the twisted
#           torus.
# @param others An array of shape (N, 2) that gives the positions on the twisted
#               torus to compute distance from.
# @param x_dim  X dimension of the torus
# @param y_dim  Y dimension of the torus
# @return       An array of shape (N, ) with all the distances
#
def remapTwistedTorus(a, others, dim):

    a_x = float(a.x)
    a_y = float(a.y)
    others_x = others.x
    others_y = others.y
    szO = others.x.shape[0]
    x_dim = dim.x
    y_dim = dim.y

    ret = np.ndarray((szO,))

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
    X, Y = np.meshgrid(np.arange(dims.x), np.arange(dims.y))
    others = Position2D()
    others.x = X.flatten()
    others.y = Y.flatten()

    a = Position2D()
    def gaussDiff(x):
        a.x = x[1] # mu_x
        a.y = x[2] # mu_y
        dist = remapTwistedTorus(a, others, dim)
        #dist = np.sqrt((others.x - a.x)**2 + (others.y - a.y)**2)
        print "C:", x[0], a, "sigma:", x[3]
        return np.abs(x[0]) * np.exp( -dist**2/2./ x[3]**2 ) - sig_f
#                       |                            |
#                       A                          sigma

    x0 = np.array([i.A0, i.mu0_x, i.mu0_y, i.sigma0])

    xest,ierr = scipy.optimize.leastsq(gaussDiff, x0)
    return xest



## Fit a 2D Gaussian onto a (potential) firing rate bump. On Twisted torus
#
# @param sig    Firing rate map to fit - a 2D numpy array
# @param sim    Dimensions of the twisted torus (Position2D)
# @return Estimated values for the Gaussian (see Position2D)
def fitGaussianBumpTT(sig, dim):
    '''
    Fit a Gaussian to a rate map, using least squares method.
    '''
    init = GaussianInitParams()
    init.mu0_y, init.mu0_x  = np.unravel_index(np.argmax(sig), sig.shape)
    init.A0 = sig[init.mu0_y, init.mu0_x]
    init.sigma0  = np.max([dim.x, dim.y]) / 2.0
    
    return fitGaussianTT(rateMap.flatten(), init, dim)




##############################################################################
#                           Tests
##############################################################################

def remapAndPlot(a, dims):
    X, Y = np.meshgrid(np.arange(dims.x), np.arange(dims.y))
    others = Position2D()
    others.x = X.flatten()
    others.y = Y.flatten()
    dist = remapTwistedTorus(a, others, dims)
    figure()
    pcolormesh(np.reshape(dist, (dims.y, dims.x)))
    title('a:(x, y): ' + str(a.x) + ', ' + str(a.y) + ', x_dim:' + str(dims.x) +
            ', y_dim: ' + str(dims.y))
    axis('equal')


if (__name__ == '__main__'):

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
