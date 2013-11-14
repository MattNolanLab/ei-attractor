#
#   fig_bump_stabilitypy
#
#   Bump drift over time
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


jobRange = [1000, 1024]
genNum = 0

jobN = jobRange[1] - jobRange[0] + 1

rcParams['font.size'] = 14


arenaDiam = 180.0     # cm
h = 3.0
Ne = 34

# Neuron to extract spikes from
neuronNum = 10
spikeType = 'excitatory'


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_gen{3:04}"

h5file = openFile('bump_stability.h5', mode = "w", title =
        "Bump stability export figures")
    
bumpPos_all = []

for job_it in range(jobN):
    jobNum = job_it + jobRange[0]
    print 'jobNum: ' + str(jobNum)

    fileName = fileNameTemp +  '_output.mat'
    fileName = fileName.format(dirName, fileNamePrefix, jobNum, genNum)
    try:
        data = loadmat(fileName)
    except:
        print "warning: could not open: " + fileName
        continue

    
    bumpPos = np.unwrap(data['bumpPos'][:, 0].ravel(), discont=Ne)
    #bumpPos = data['bumpPos'][:, 0].ravel()
    times   = data['bumpPos_times'].ravel()
    bumpPos_all.append(bumpPos)

    h5file.createArray(h5file.root, 'bumpPos_' + str(job_it), bumpPos)

h5file.createArray(h5file.root, 'bumpPos_all', bumpPos_all)
h5file.createArray(h5file.root, 'bumpPos_times', times)
h5file.createArray(h5file.root, 'bumpPos_var', np.var(bumpPos_all, 0))

h5file.close()

