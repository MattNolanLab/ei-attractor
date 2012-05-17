#
#   grid_cell_analysis.py
#
#   Grid cell analysis module. Use this to analyse spikes/membrane potentials
#   of grid cell models.
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

from scipy.integrate    import trapz
from matplotlib.pyplot  import *


def gaussianFilter(X, sigma):
    '''Simple Gaussian function'''
    return np.exp(-X**2/ 2.0 / sigma**2)



def extractSpikePositions2D(spikeTimes, rat_pos_x, rat_pos_y, dt):
    '''
    Extract spike positions from the rat tracking data and cell spike times.
    Both positions and spikes must be aligned!
    '''
    neuronPos_i = np.array(spikeTimes/dt, dtype=int)
    neuronPos_x = rat_pos_x[neuronPos_i]
    neuronPos_y = rat_pos_y[neuronPos_i]

    return (neuronPos_x, neuronPos_y, np.max(neuronPos_i))


def plotSpikes2D(spikeTimes, rat_pos_x, rat_pos_y, dt):
    '''
    Plot spike positions into the figure. Both positions and spikes must be aligned!
    '''
    neuronPos_x, neuronPos_y, m_i = extractSpikePositions2D(spikeTimes, rat_pos_x, rat_pos_y, dt)

    plot(rat_pos_x[0:m_i], rat_pos_y[0:m_i])
    hold('on')
    plot(neuronPos_x, neuronPos_y, 'or', markersize=5)
    hold('off')
    axis('off')
    axis('equal')


def SNSpatialRate2D(spikeTimes, rat_pos_x, rat_pos_y, dt, arenaDiam, h):
    '''
    Preprocess neuron spike times into a spatial rate map, given arena parameters.
    Both spike times and rat tracking data must be aligned in time!
    '''
    precision = arenaDiam/h
    xedges = np.linspace(-arenaDiam/2, arenaDiam/2, precision+1)
    yedges = np.linspace(-arenaDiam/2, arenaDiam/2, precision+1)

    rateMap = np.zeros((len(xedges), len(yedges))) * np.nan

    for x_i in xrange(len(xedges)):
        for y_i in xrange(len(yedges)):
            x = xedges(x_i)
            y = yedges(y_i)
            isNearTrack = np.count_nonzero(np.sqrt((rat_pos_x - x)**2 + (rat_pos_y - y)**2) <= h) > 0

            if isNearTrack:
                normConst = trapz(gaussianFilter(np.sqrt((rat_pos_x - x)**2 + (rat_pos_y - y)**2), sigma=h), dt)
                neuronPos_x, neuronPos_y, m_i = extractSpikePositions2D(spikeTimes, rat_pos_x, rat_pos_y, dt)
                spikes = np.sum(gaussianFilter(np.sqrt((neuronPos_x - x)**2 + (neuronPos_y - y)**2), sigma=h))
                rateMap[x_i, y_i] = spikes/normConst

    return  rateMap, xedges, yedges


def plotSNSpatialRate2D(spikeTimes, rat_pos_x, rat_pos_y, dt, arenaDiam, h):
    '''
    Create a 2D rate map from spike times and rat tracking data and plot this.
    '''

    rateMap, xedges, yedges = SNSpatialRate2D(spikeTimes, rat_pos_x, rat_pos_y, dt, arenaDiam, h)

    X, Y = np.meshgrid(xedges, yedges)
    figure(fig)
    pcolormesh(X, Y, rateMap)
    colormap('jet')


