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
from matplotlib.pyplot  import *

from grid_cell_analysis import *


jobRange = [1200, 1249]
trialNum = 0
dumpNum = 9

jobN = jobRange[1] - jobRange[0] + 1

rcParams['font.size'] = 14


arenaDiam = 180.0     # cm
h = 3.0

# Neuron to extract spikes from
neuronNum = 10


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"


for job_it in range(jobN):
    jobNum = job_it + jobRange[0]
    print 'jobNum: ' + str(jobNum)

    fileName = fileNameTemp +  '_output.mat'
    fileName = fileName.format(dirName, fileNamePrefix, jobNum,
                                trialNum, dumpNum)
    try:
        data = loadmat(fileName)
    except:
        print "warning: could not open: " + fileName
        continue

    pos_x           = data['pos_x'].ravel()
    pos_y           = data['pos_y'].ravel()
    spikeTimes_e    = data['spikeCell_e'].ravel()
    rat_dt          = data['dt'][0][0]

    figure()
    plotSpikes2D(spikeTimes_e[neuronNum], pos_x, pos_y, rat_dt)
    savefig(fileName + '_spikePlot.pdf')
    close()


