#!/usr/bin/env python
'''
Create a 3D plot of the excitatory population firing rate snapshot, overlaid
with a fitted Gaussian.
'''

import numpy as np
import mayavi.mlab as mlab

from EI_plotting import aggregate as aggr
from EI_plotting.base     import NoiseDataSpaces

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


outputDir = "panels/"

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
exampleIdx   = [(0, 0), (0, 0), (0, 0)] # (row, col)
gridsDataRoot= None
bumpDataRoot= 'output_local/even_spacing/gamma_bump'
velDataRoot = None
shape = (31, 31)


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


exampleRC = ( (5, 15), (15, 5) )


class GaussianGenerator(object):
    def __init__(self, A, mu, sigma):
        self.A = A
        self.sigma = sigma
        self.mu = mu

    def __call__(self, x, y):
        mu_x = self.mu[0]
        mu_y = self.mu[1]
        return A * np.exp(- ((x-mu_x)**2 + (y-mu_y)**2) / (2. * self.sigma**2))


##############################################################################
# Do it!
rateMapTypes = ['bump', 'rateMap_e']
resolution = 8
scale_factor = 0.25

# Plot the firing rate snapshot (E cells)
#for ns_idx in [1]:
#    rateMaps, _, _ = aggr.aggregateType(ps.bumpGamma[ns_idx], iterList,
#            rateMapTypes, NTrials)
#    r = exampleRC[0][0]
#    c = exampleRC[0][1]
#    rateMap_e = rateMaps[r][c]
#    [nY, nX] = rateMap_e.shape
#
#    X, Y = np.meshgrid(np.arange(nX), np.arange(nY))
#    X = X.ravel()
#    Y = Y.ravel()
#    mlab.points3d(X, Y, rateMap_e.ravel(), scale_factor=scale_factor,
#            resolution=resolution)

# Plot the fitted Gaussian
mu = [0, 0]
sigma = 2.5
A = 10
f = GaussianGenerator(A, mu, sigma)
x, y = np.mgrid[-17:17:0.5, -15:15:0.5]
mlab.surf(x, y, f, opacity=0.5)

#fname = outputDir + "/bump_fitting_3d.png"
#mlab.savefig(fname)
