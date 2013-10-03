#!/usr/bin/env python
#
#   suppFigure_velocity.py
#
#   Supplementary figure that illustrates bump speed tracking simulations.
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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

from parameters  import JobTrialSpace2D, DataSpace
from analysis.spikes import TorusPopulationSpikes
from analysis.visitors.interface import extractSpikes
from figures_shared import plotBump, getNoiseRoots


import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)

from matplotlib import rc
rc('pdf', fonttype=42)
rc('mathtext', default='regular')

plt.rcParams['font.size'] = 10


DS = DataSpace


noise_sigmas = [0, 150, 300]
velDataRoot = 'output_local/even_spacing/velocity'
velShape  = (31, 31)

# Flags
bump_movement = 1

###############################################################################
velRoots = getNoiseRoots(velDataRoot, noise_sigmas)
velDataSpace0   = JobTrialSpace2D(velShape, velRoots[0])
velDataSpace150 = JobTrialSpace2D(velShape, velRoots[1])
velDataSpace300 = JobTrialSpace2D(velShape, velRoots[2])



class BumpTiming(object):
    def __init__(self, tStart=0, nRows=2, nCols=5, dt=500, winLen=500):
        self.tStart = tStart
        self.nRows  = nRows
        self.nCols  = nCols
        self.dt     = dt     # ms
        self.winLen = winLen # ms

    @property
    def tEnd(self):
        return (self.nRows*self.nCols - 1) * self.dt



def _getSpikeTrain(data, monName, dimList):
    N_x = DS.getNetParam(data, dimList[0])
    N_y = DS.getNetParam(data, dimList[1])
    senders, times = extractSpikes(data[monName])
    return senders, times, (N_x, N_y)


def drawBumpMovement(dataSpace, r, c, fileName, **kw):
    # Keyword args
    figSize  = kw.pop('figSize', (5, 2))
    trialNum = kw.pop('trialNum', 0)
    IvelIdx  = kw.pop('IvelIdx', 3)
    bt       = kw.pop('bumpTiming', BumpTiming())
    h_pad    = kw.pop('h_pad', 2)

    ds     = dataSpace[r][c].getAllTrialsAsDataSet()
    data   = ds.data
    trials = data['trials']
    Ivel = trials[trialNum]['IvelVec'][IvelIdx]
    iData  = trials[trialNum]['IvelData'][IvelIdx]

    # Extract spikes
    senders, times, sheetSize =  _getSpikeTrain(iData, 'spikeMon_e', ['Ne_x',
        'Ne_y'])
    pop = TorusPopulationSpikes(senders, times, sheetSize)
    Fe, Ft = pop.slidingFiringRate(bt.tStart, bt.tEnd, bt.dt, bt.winLen)

    # Create the example plot
    fig = plt.figure(figsize=figSize)
    idx = 0
    for r in range(bt.nRows):
        for c in range(bt.nCols):
            ax = plt.subplot(bt.nRows, bt.nCols, idx+1)
            plotBump(ax, Fe[:, :, idx], maxRate=False)
            plt.title("{0} s".format(Ft[idx]*1e-3), fontsize='small')
            idx += 1

    print("Ivel: {0} pA".format(Ivel))
    fig.tight_layout(h_pad=h_pad)
    plt.savefig(fileName, dpi=300, transparent=False)
    plt.close()


# Bump Movement
if (bump_movement):
    drawBumpMovement(velDataSpace150,
            r = 6, c = 12,
            fileName = "suppFigure_bump_moving.png")

