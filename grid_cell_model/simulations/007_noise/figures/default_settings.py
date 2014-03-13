'''
Default figure plotting and analysis settings.
'''
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


sweepFigSize = (3.7, 2.6)
sweepLeft   = 0.08
sweepBottom = 0.2
sweepRight  = 0.8
sweepTop    = 0.85
sweepTransparent  = True


def getDefaultSweepFig():
    fig = plt.figure(figsize=sweepFigSize)
    return fig, getDefaultSweepAxes(fig)

def getDefaultSweepAxes(fig):
    return fig.add_axes(Bbox.from_extents(sweepLeft, sweepBottom, sweepRight,
        sweepTop))

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
