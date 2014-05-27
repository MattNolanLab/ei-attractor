'''
Default figure plotting and analysis settings.
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

from EI_plotting.base     import NoiseDataSpaces

plt.rcParams['font.size'] = 30

# Math fonts are all sans-serif
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

scaleFactor = 2.5

linewidth = 1. * scaleFactor
rc('lines', linewidth=linewidth, markersize=6.*scaleFactor)
rc('axes', linewidth=linewidth)

iterList  = ['g_AMPA_total', 'g_GABA_total']
figOutputDir = 'panels/'
evenShape = (31, 31)
noise_sigmas = [0, 150, 300]


sweepFigSize = np.array([3.7*scaleFactor, 2.6*scaleFactor])
sweepW = .72
sweepH = .65
sweepLeft   = 0.08
sweepBottom = 0.2


def getDefaultSweepFig(scale=1.0, colorBarPos='right'):
    fig = plt.figure(figsize=sweepFigSize*scale)
    return fig, getDefaultSweepAxes(fig, colorBarPos)

def getDefaultSweepAxes(fig, colorBarPos='right'):
    if colorBarPos == 'right':
        left = sweepLeft
    else:
        left = .12

    right = left + sweepW
    top = sweepBottom + sweepH
    return fig.add_axes(Bbox.from_extents(left, sweepBottom, right, top))

###############################################################################
# Parameter spaces roots
gridsDataRoot    = 'simulation_data/submission/grids'
bumpDataRoot     = 'simulation_data/submission/gamma_bump'
velDataRoot      = 'simulation_data/submission/velocity'
constPosDataRoot = 'simulation_data/submission/const_position'

def getDefaultParamSpaces():
    roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot,
            constPos=constPosDataRoot)
    return NoiseDataSpaces(roots, evenShape, noise_sigmas)
