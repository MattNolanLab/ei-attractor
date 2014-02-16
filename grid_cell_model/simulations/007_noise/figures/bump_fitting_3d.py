#!/usr/bin/env python
'''
Create a 3D plot of the excitatory population firing rate snapshot, overlaid
with a fitted Gaussian.
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from EI_plotting import aggregate as aggr
from EI_plotting.base     import NoiseDataSpaces

from submitting import flagparse
parser = flagparse.FlagParser()
args = parser.parse_args()

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
trialIdx = 0

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the firing rate snapshot (E cells)
for ns_idx in [1]:
    rateMaps, _, _ = aggr.aggregateType(ps.bumpGamma[ns_idx], iterList,
            rateMapTypes, NTrials)
    r = exampleRC[0][0]
    c = exampleRC[0][1]
    rateMap_e = rateMaps[r][c][trialIdx]
    [nY, nX] = rateMap_e.shape

    X, Y = np.meshgrid(np.arange(nX), np.arange(nY))
    X = X.ravel()
    Y = Y.ravel()
    rateMap = rateMap_e.ravel()
    ax.plot_wireframe(X, Y, rateMap)


## Plot the fitted Gaussian
#mu = [0, 0]
#sigma = 2.5
#A = 10
#x, y = np.mgrid[-17:17:0.25, -15:15:0.25]
#f = GaussianGenerator(A, mu, sigma)(x, y)
#ax.plot_surface(x, y, f, rstride=1, cstride=1, linewidth=0, antialiased=False,
#        cmap='coolwarm', shade=True, alpha=0.5)

plt.show()
#fname = outputDir + "/bump_fitting_3d.png"
#fig.savefig(fname)
