#!/usr/bin/env python
#
#   figure4.py
#
#   Final figure: network mechanisms
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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot   import figure, plot, savefig
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, MaxNLocator, \
        ScalarFormatter
from matplotlib.transforms import Bbox

import EI_plotting as EI
from plotting.global_defs import globalAxesSettings
from figures_shared import getNoiseRoots, NoiseDataSpaces
from data_storage.sim_models.ei import MonitoredSpikes

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 11

###############################################################################
cFreq = 'blue'
cAC = 'green'
cCount = 'red'

outputDir = "."
NTrials = 5
iterList  = ['g_AMPA_total', 'g_GABA_total']

noise_sigmas  = [0, 150, 300]
rasterRC      = [(18, 15), (18, 15), (18, 15)] # (row, col)
bumpDataRoot  = 'output_local/even_spacing/gamma_bump'
velDataRoot   = None
gridsDataRoot = None
shape = (31, 31)

raster0   = 1
raster150 = 1
raster300 = 1

###############################################################################

def rasterPlot(dataSpaces, noise_sigma_idx, r, c, trialNum=0, **kw):
    title = kw.pop('title', True)
    tLimits = [2e3, 2.5e3] # ms

    noise_sigma = dataSpaces.noise_sigmas[noise_sigma_idx]
    data = dataSpaces.bumpGamma[noise_sigma_idx][r][c][trialNum].data
    ESpikes = MonitoredSpikes(data, 'spikeMon_e', 'net_Ne')
    ISpikes = MonitoredSpikes(data, 'spikeMon_i', 'net_Ni')
    ax = EI.plotEIRaster(ESpikes, ISpikes, tLimits, **kw) 
    ax.set_title('$\sigma$ = {0} pA'.format(int(noise_sigma)), y=1.02,
            va='bottom')
    return ax


###############################################################################
roots = NoiseDataSpaces.Roots(bumpDataRoot, velDataRoot, gridsDataRoot)
ps    = NoiseDataSpaces(roots, shape, noise_sigmas)


rasterFigSize = (3.5, 1.75)
transparent   = True
rasterLeft    = 0.2
rasterBottom  = 0.1
rasterRight   = 0.99
rasterTop     = 0.82


if (raster0):
    # noise_sigma = 0 pA
    fig = figure("rasters0", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasterPlot(ps, 
            noise_sigma_idx=0,
            r=rasterRC[0][0], c=rasterRC[0][1])
    fname = outputDir + "/figure4_raster0.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

if (raster150):
    # noise_sigma = 150 pA
    fig = figure("rasters150", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasterPlot(ps, 
            noise_sigma_idx=1,
            r=rasterRC[1][0], c=rasterRC[1][1],
            ylabel='', yticks=False)
    fname = outputDir + "/figure4_raster150.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        

if (raster300):
    # noise_sigma = 300 pA
    fig = figure("rasters300", figsize=rasterFigSize)
    ax = fig.add_axes(Bbox.from_extents(rasterLeft, rasterBottom, rasterRight,
        rasterTop))
    rasterPlot(ps, 
            noise_sigma_idx=2,
            r=rasterRC[2][0], c=rasterRC[2][1],
            ylabel='', yticks=False)
    fname = outputDir + "/figure4_raster300.png"
    fig.savefig(fname, dpi=300, transparent=transparent)
    plt.close()
        
