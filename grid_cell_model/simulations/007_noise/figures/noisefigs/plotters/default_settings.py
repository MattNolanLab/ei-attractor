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
