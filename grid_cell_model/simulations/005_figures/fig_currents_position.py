#
#   fig_position.py
#
#   Average/maximal syn. currents as a function of the distance from the
#   centre of the grid field.
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
from tools              import *


jobRange = [1200, 1200]
trialNum = 0
dumpNum = 0

jobN = jobRange[1] - jobRange[0] + 1

rcParams['font.size'] = 14


arenaDiam = 180.0     # cm
h = 3.0
Ne = 34

# Neuron to extract spikes from
neuronNum = 10
spikeType = 'excitatory'

pA = 1e-12
ms = 1e-3

f_start = 40
f_stop   = 200

theta_freq = 8  # Hz

dirName = "tmp/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"


def sliceArrayToThetaCycles(sig, times, theta_freq, removeFirst=0):
    '''
    Slice an array to theta cycles. Take a 1D array and return a 2D array,
    in which each row is one theta cycle. The array must be properly aligned
    '''
    dt = times[1] - times[0]
    thetaPoints = 1. / theta_freq / dt
    numCycles = np.floor(len(sig) / thetaPoints)
    sig_sliced =  np.reshape(sig[0:thetaPoints*numCycles], (numCycles, thetaPoints))
    time_centers = (np.arange(0, thetaPoints*numCycles, thetaPoints) + thetaPoints/2)*dt
    return (sig_sliced[removeFirst:, :], time_centers[removeFirst:])



def findThetaMaxMin(sig, times, theta_freq):
    '''
    For each theta cycle (given by theta_freq), apply the extremum function
    and return an array of times of theta cycles and corresponding values
    for each cycle.
    Assuming 'times' is uniform
    '''
    dt = times[1] - times[0]
    sig_sliced, time_centers = sliceArrayToThetaCycles(syn_current, times, theta_freq, removeFirst=5)
    sig_max = np.max(sig_sliced, axis=1)
    sig_min = np.min(sig_sliced, axis=1)
    return (sig_max, sig_min, time_centers)
     


for job_it in range(jobN):
    jobNum = job_it + jobRange[0]
    print 'jobNum: ' + str(jobNum)

    fileName = fileNameTemp
    fileName = fileName.format(dirName, fileNamePrefix, jobNum, trialNum, dumpNum)
    try:
        data = loadmat(fileName +  '_output.mat')
    except:
        print "warning: could not open: " + fileName
        continue

    h5file = openFile(fileName + '_export.h5', mode = "w", title =
            "Bump stability export figures")
    

    pos_x           = data['pos_x'].ravel()
    pos_y           = data['pos_y'].ravel()
    rat_dt          = data['dt'][0][0]

    times           = data['stateMon_times'].ravel()
    rec_dt          = data['options']['stateRec_dt'][0][0][0][0] * ms
    velocityStart   = data['velocityStart'].ravel()[0] * ms

    syn_current       = data['stateMon_Iclamp_e_values'][0, :].ravel()/pA
    syn_current_bpass = syn_current
    #syn_current_bpass = butterBandPass(syn_current, rec_dt, f_start, f_stop)

    sig_max, sig_min, time_centers = findThetaMaxMin(syn_current, times, theta_freq)
    # Recompute time_centers into positions
    time_centers -= velocityStart # relative to the start of velocity modulation
    centers_pos_i = time_centers / rat_dt
    centers_pos_x = pos_x[centers_pos_i]
    centers_pos_y = pos_y[centers_pos_i]


#h5file.createArray(h5file.root, 'bumpPos_all', bumpPos_all)
#h5file.createArray(h5file.root, 'bumpPos_times', times)
#h5file.createArray(h5file.root, 'bumpPos_var', np.var(bumpPos_all, 0))
#
#h5file.close()

