#! /usr/bin/env python
'''
Create a 3D scatterplot of the firing rates of neurons. On top of that,
superimpose a surface plot of a Gaussian with parameters from the fitting
procedure.

.. TODO::

    Currently, the surface data is taken from ``/analysis/bump_e/positions``,
    however the scatter plot is taken from a single bump fitting time at the
    end of the simulations. In order for both of these plots to be properly
    aligned the times should also be aligned (which should be the case now).
'''
import logging

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import scipy.ndimage as ndi
import numpy as np

from EI_plotting import aggregate as aggr
from EI_plotting.base     import NoiseDataSpaces
import grid_cell_model.analysis.image as image

from grid_cell_model.submitting import flagparse
parser = flagparse.FlagParser()
args = parser.parse_args()

logger = logging.getLogger(__name__)
outputDir = "panels/"

NTrials=10
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas = [0, 150, 300]
gridsDataRoot= None
bumpDataRoot= 'output_local/even_spacing/gamma_bump'
velDataRoot = None
shape = (31, 31)


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)




class GaussianGenerator(object):
    def __init__(self, A, mu, sigma):
        self.A = A
        self.sigma = sigma
        self.mu = mu

    def __call__(self, x, y):
        mu_x = self.mu[0]
        mu_y = self.mu[1]
        return A * np.exp(- ((x-mu_x)**2 + (y-mu_y)**2) / (2. * self.sigma**2))


## Create a GL View widget to display data
app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: GLSurfacePlot')
w.resize(400, 400)
w.setCameraPosition(distance=70)
w.qglClearColor(pg.mkColor('w'))


##############################################################################
# Do it!
rateMapTypes = ['bump', 'rateMap_e']
gaussianVarList = ['analysis', 'bump_e', 'positions']
exampleRC = ( (23, 1), (15, 5) )
trialIdx = 0
gaussT = 9.5e3 # ms

# Plot the firing rate snapshot (E cells)
for ns_idx in [0]:
    r = exampleRC[0][0]
    c = exampleRC[0][1]

    # Plot rate maps
    rateMaps, _, _ = aggr.aggregateType(ps.bumpGamma[ns_idx], iterList,
            rateMapTypes, NTrials)
    rateMap_e = rateMaps[r][c][trialIdx]
    [nY, nX] = rateMap_e.shape

    X2D, Y2D = np.meshgrid(np.arange(nX), np.arange(nY))
    X = X2D.ravel()
    Y = Y2D.ravel()
    rateMap = rateMap_e.ravel()
    pos = np.vstack((X, Y, rateMap)).T
    pRates = gl.GLScatterPlotItem(pos=pos, color=(1., 0, 0, 1), size=5,
            pxMode=True)

    # Plot the fitted Gaussian
    A = ps.bumpGamma[ns_idx].aggregateData(gaussianVarList + ['A'], [trialIdx],
            output_dtype='list', loadData=True, saveData=False)
    mu_x = ps.bumpGamma[ns_idx].aggregateData(gaussianVarList + ['mu_x'],
            [trialIdx], output_dtype='list', loadData=True, saveData=False)
    mu_y = ps.bumpGamma[ns_idx].aggregateData(gaussianVarList + ['mu_y'],
            [trialIdx], output_dtype='list', loadData=True, saveData=False)
    sigma = ps.bumpGamma[ns_idx].aggregateData(gaussianVarList + ['sigma'],
            [trialIdx], output_dtype='list', loadData=True, saveData=False)
    t = ps.bumpGamma[ns_idx].aggregateData(gaussianVarList + ['t'], [trialIdx],
            output_dtype='list', loadData=True, saveData=False)

    timeIdx = np.nonzero(t[r][c][trialIdx] <= gaussT)[0][-1]
    A     = A[r][c][trialIdx][timeIdx]
    mu_x  = mu_x[r][c][trialIdx][timeIdx]
    mu_y  = mu_y[r][c][trialIdx][timeIdx]
    sigma = sigma[r][c][trialIdx][timeIdx]
    print("A: {0}, mu_x: {1}, mu_y: {2}, sigma: {3}".format(A, mu_x, mu_y,
        sigma))

    mu = image.Position2D(mu_x, mu_y)
    xy = image.Position2D(X, Y)
    torusSize = image.Position2D(nX, nY)
    distances = image.remapTwistedTorus(mu, xy, torusSize)
    gaussianZ = A * np.exp( - distances**2 / 2. / sigma**2 )
    gaussPos = np.vstack((X, Y, gaussianZ)).T
    Z2D = np.reshape(gaussianZ, (nY, nX)).T
    pGaussian = gl.GLSurfacePlotItem(
            x=np.arange(nX),
            y=np.arange(nY),
            z=Z2D,
            color=(0.9, 0.9, 0.9, 0.1),
            shader='shaded',
            smooth=True)
    w.addItem(pGaussian)
    w.addItem(pRates)

    

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()

    #w.grabFrameBuffer().save('panels/bump_fitting_3d.png')

