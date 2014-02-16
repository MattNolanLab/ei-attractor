from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import scipy.ndimage as ndi
import numpy as np

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


## Create a GL View widget to display data
app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: GLSurfacePlot')
w.setCameraPosition(distance=50)
w.qglClearColor(pg.mkColor('w'))


##############################################################################
# Do it!
rateMapTypes = ['bump', 'rateMap_e']
trialIdx = 0

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
    pos = np.vstack((X, Y, rateMap)).T
    p1 = gl.GLScatterPlotItem(pos=pos, color=(0.5, 0.5, 1, 1), size=5)
    p1.setGLOptions('translucent')
    w.addItem(p1)





### Add a grid to the view
#g = gl.GLGridItem()
#g.scale(2,2,1)
#g.setDepthValue(10)  # draw grid after surfaces since they may be translucent
#w.addItem(g)



### Manually specified colors
#z = ndi.gaussian_filter(np.random.normal(size=(50,50)), (1,1))
#x = np.linspace(-12, 12, 50)
#y = np.linspace(-12, 12, 50)
#colors = np.ones((50,50,4), dtype=float)
#colors[...,0] = np.clip(np.cos(((x.reshape(50,1) ** 2) + (y.reshape(1,50) ** 2)) ** 0.5), 0, 1)
#colors[...,1] = colors[...,0]
#
#p3 = gl.GLSurfacePlotItem(z=z, colors=colors.reshape(50*50,4), shader='shaded', smooth=False)
#p3.scale(16./49., 16./49., 1.0)
#p3.translate(2, -18, 0)
#w.addItem(p3)
#


## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()

    w.grabFrameBuffer().save('bump3d.png')

