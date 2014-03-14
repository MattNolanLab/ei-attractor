'''
Default figure plotting and analysis settings.
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

from EI_plotting.base     import NoiseDataSpaces

plt.rcParams['font.size'] = 11

# Math fonts are all sans-serif
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

iterList  = ['g_AMPA_total', 'g_GABA_total']
figOutputDir = 'panels/'
evenShape = (31, 31)
noise_sigmas = [0, 150, 300]


sweepFigSize = np.array([3.7, 2.6])
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
gridsDataRoot    = 'output_local/even_spacing/grids'
bumpDataRoot     = 'output_local/even_spacing/gamma_bump'
velDataRoot      = 'output_local/even_spacing/velocity_vertical'
constPosDataRoot = 'output_local/even_spacing/const_position'

def getDefaultParamSpaces():
    roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot,
            constPos=constPosDataRoot)
    return NoiseDataSpaces(roots, evenShape, noise_sigmas)
