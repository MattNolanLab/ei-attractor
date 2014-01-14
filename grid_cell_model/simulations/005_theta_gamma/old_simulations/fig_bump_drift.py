#
#   fig_bump_drift.py
#
#   Bump drift over time: determining whether the source is random or systematic
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
from tables             import *

from grid_cell_analysis import *


jobNums = range(7000, 7020)
genNum = 0
trialNums = range(0, 20)

rcParams['font.size'] = 16
nTicks_y = 3


arenaDiam = 180.0     # cm
h = 3.0
Ne = 34


# Neuron to extract spikes from
neuronNum = 10
spikeType = 'excitatory'


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_gen{3:04}_trial{4:04}"

for jobNum in jobNums:
    bumpPos_all_x = []
    bumpPos_all_y = []
    
    for trial_it in trialNums:
        print 'trial_it: ' + str(trial_it)
    
        fileName = fileNameTemp
        fileName = fileName.format(dirName, fileNamePrefix, jobNum, genNum, trial_it)
        try:
            data = loadmat(fileName +  '_output.mat')
        except:
            print "warning: could not open: " + fileName
            continue
    
        
        bumpPos_x = np.unwrap(data['bumpPos'][:, 0].ravel(), discont=Ne)
        bumpPos_y = data['bumpPos'][:, 1].ravel()
    
        times     = data['bumpPos_times'].ravel()
        bumpPos_all_x.append(bumpPos_x - bumpPos_x[0])
        bumpPos_all_y.append(bumpPos_y - bumpPos_y[0])
    
    bumpPos_all_x = np.array(bumpPos_all_x)
    bumpPos_all_y = np.array(bumpPos_all_y)
    
    figure()
    subplot(2, 1, 1)
    plot(times, bumpPos_all_x.T, color="grey")
    hold('on')
    plot(times, np.mean(bumpPos_all_x, axis=0), "b", linewidth=2)
    grid('on')
    gca().yaxis.set_major_locator(MaxNLocator(nTicks_y))
    ylabel('X drift')
    ylim([-30, 30])
    
    subplot(2, 1, 2)
    plot(times, bumpPos_all_y.T, color="grey")
    plot(times, np.mean(bumpPos_all_y, axis=0), "b", linewidth=2)
    grid('on')
    gca().yaxis.set_major_locator(MaxNLocator(nTicks_y))
    xlabel('Time (s)')
    ylabel('Y drift')
    ylim([-30, 30])

    savefig(fileName + '_bump_drifts.pdf')

    close()

