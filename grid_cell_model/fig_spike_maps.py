#
#   fig_spike_maps.py
#
#   Spike and rate maps of grid cells. Data analysis
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

from scipy.io           import loadmat
from scipy.io           import savemat
from matplotlib.pyplot  import *
from numpy.fft          import fft2

from grid_cell_analysis import *


jobRange = [0, 0]
dumpNum = 0

jobN = jobRange[1] - jobRange[0] + 1

rcParams['font.size'] = 14


arenaDiam = 180.0     # cm
h = 3.0

# Neuron to extract spikes from
neuronNum = 10
spikeType = 'excitatory'


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_dump{3:03}"

gridnessScores = []

for job_it in range(jobN):
    jobNum = job_it + jobRange[0]
    print 'jobNum: ' + str(jobNum)

    fileName = fileNameTemp
    fileName = fileName.format(dirName, fileNamePrefix, jobNum, dumpNum)
    try:
        data = loadmat(fileName +  '_output.mat')
    except:
        print "warning: could not open: " + fileName
        continue

    pos_x           = data['pos_x'].ravel()
    pos_y           = data['pos_y'].ravel()
    rat_dt          = data['dt'][0][0]
    velocityStart   = data['velocityStart'][0][0]
    if spikeType == 'excitatory':
        spikeTimes  = data['spikeCell_e'].ravel()
    if spikeType == 'inhibitory':
        spikeTimes  = data['spikeCell_i'].ravel()

    gridSep         = data['options']['gridSep'][0][0][0][0]
    corr_cutRmin    = gridSep / 2

    spikes = spikeTimes[neuronNum] - velocityStart*1e-3
    spikes = np.delete(spikes, np.nonzero(spikes < 0)[0])

    figure()
    plotSpikes2D(spikes, pos_x, pos_y, rat_dt)
    savefig(fileName + '_spikePlot_' + spikeType + '.png')

    figure()
    rateMap, xedges, yedges = SNSpatialRate2D(spikes, pos_x, pos_y, rat_dt, arenaDiam, h)
    X, Y = np.meshgrid(xedges, yedges)
    pcolormesh(X, Y, rateMap)
    colorbar()
    axis('equal')
    axis('off')
    savefig(fileName + '_rateMap_' + spikeType + '.png')

    
    figure()
    FT_size = 256
    Fs = 1.0/(h/100.0) # units: 1/m
    rateMap_pad = np.ndarray((FT_size, FT_size))
    rateMap_pad[:, :] = 0
    rateMap_pad[0:rateMap.shape[0], 0:rateMap.shape[0]] = rateMap - np.mean(rateMap.flatten())
    FT = fft2(rateMap_pad)
    fxy = np.linspace(-1.0, 1.0, FT_size)
    fxy_igor = Fs/2.0*np.linspace(-1.0, 1.0, FT_size+1)
    FX, FY = np.meshgrid(fxy, fxy)
    FX *= Fs/2.0
    FY *= Fs/2.0
    PSD_centered = np.abs(np.fft.fftshift(FT))**2
    pcolormesh(FX, FY, PSD_centered)
    #axis('equal')
    xlim([-10, 10])
    ylim([-10, 10])
    savefig(fileName + '_fft2' + spikeType + '.png')


    figure()
    corr, xedges_corr, yedges_corr = SNAutoCorr(rateMap, arenaDiam, h)
    X, Y = np.meshgrid(xedges_corr, yedges_corr)
    pcolormesh(X, Y, corr)
    axis('equal')
    axis('off')
    savefig(fileName + '_rateCorr_' + spikeType + '.png')


    figure()
    G, crossCorr, angles = cellGridnessScore(rateMap, arenaDiam, h, corr_cutRmin)
    gridnessScores.append(G)
    plot(angles, crossCorr)
    xlabel('Angle (deg.)')
    ylabel('Corr. coefficient')
    savefig(fileName + '_gridnessCorr_' + spikeType + '.png')

    close('all')


print "Gridness scores:"
print gridnessScores

