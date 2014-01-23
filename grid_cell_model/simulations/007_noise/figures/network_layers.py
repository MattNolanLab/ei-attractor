import numpy as np
import mayavi.mlab as mlab

nX = 34
nY = 30
ZOffset = 20.
resolution = 8
scale_factor = 0.25
outputDir = 'panels'
arrowOffset = 0.25
ao = arrowOffset

X, Y = np.meshgrid(np.arange(nX), np.arange(nY))
X = X.ravel()
Y = Y.ravel()

ZE = np.zeros(nX * nY)
ZI = np.zeros(nX * nY) + ZOffset

mlab.figure(bgcolor=(1, 1, 1), size=(610*2, 582*2))

# Layers
mlab.points3d(X, Y, ZE, scale_factor=scale_factor, color=(1, 0, 0),
        resolution=resolution)
mlab.points3d(X, Y, ZI, scale_factor=scale_factor, color=(0, 0, 1),
        resolution=resolution)

# gE, gI arrows
x = [ao * nX, (0.9-ao)*nX]
y = [ao * nY, (0.9-ao)*nY]
zStart = [0.2 * ZOffset, 0.8 * ZOffset]
zComp = np.array([0.6 * ZOffset, -0.6 * ZOffset])
mlab.quiver3d(x, y, zStart, [0, 0], [0, 0], zComp, color=(0, 0, 0), line_width=20,
        scale_factor=1, mode='arrow', resolution=resolution)

textScale = 2.0
mlab.text3d(x[0] - 0.075*nX, y[0] - 0.075*nY, ZOffset/2, "gE", orientation=[0, 0, 0],
        orient_to_camera=True, color=(0, 0, 0), line_width=3, scale=textScale)

mlab.text3d(x[1] + 0.075*nX, y[1] + 0.075*nY, ZOffset/2, "gI", orientation=[0, 0, 0],
        orient_to_camera=True, color=(0, 0, 0), line_width=3, scale=textScale)

# Neuron labels
startLabelOffset = 1
labelScale = 1.5
mlab.text3d(-startLabelOffset, -startLabelOffset, ZOffset - 1,
        str(nX), orientation=[0, 0, 0], orient_to_camera=True, color=(0, 0, 0),
        line_width=3, scale=labelScale)
mlab.text3d(nX, -startLabelOffset, ZOffset - 1,
        "1", orientation=[0, 0, 0], orient_to_camera=True, color=(0, 0, 0),
        line_width=3, scale=labelScale)
mlab.text3d(nX - 1, nY + 1, ZOffset + 0.5,
        str(nY), orientation=[0, 0, 0], orient_to_camera=True, color=(0, 0, 0),
        line_width=3, scale=labelScale)

roll = 177.9584710619396
view = (-20.96663248113742,
         107.30927449790735,
          94.066015884153884,
           np.array([ 17.14404891,  16.30124532,   9.85753332]))
mlab.view(*view)
mlab.roll(roll)


fname = outputDir + "/network_layers.png"
mlab.savefig(fname)
