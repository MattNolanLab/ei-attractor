#
#   fig_currents_position.py
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

from scipy.io                   import loadmat
from scipy.ndimage.measurements import label
from matplotlib.pyplot          import *
from tables                     import *

from grid_cell_analysis         import *
from tools                      import *


jobRange = [3101, 3101]
trialNum = 0
dumpNum = 7

jobN = jobRange[1] - jobRange[0] + 1

rcParams['font.size'] = 14


arenaDiam = 180.0     # cm
h = 3.0
Ne = 34
rateThreshold_e = 2    # fields are > this
rateThreshold_i = 20   # fields are < this
minFieldArea = 30

# Neuron to extract spikes from
spikeType = 'excitatory'

pA = 1e-12
ms = 1e-3

f_start = 40
f_stop   = 200

theta_freq = 8  # Hz

dirName = "output/"
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
     

def getCurrentChargeTheta(sig, times, theta_freq):
    '''
    For each theta cycle, find the total charge of synaptic current.
    '''
    dt = times[1] - times[0]
    sig_sliced, time_centers = sliceArrayToThetaCycles(sig, times, theta_freq, removeFirst=5)
    return np.trapz(sig_sliced, dx=dt, axis=1), time_centers


def getRateMapThreshold_e(rateMap, threshold):
    return np.where(rateMap > threshold, 1, 0)

def getRateMapThreshold_i(rateMap, threshold):
    return np.where(rateMap < threshold, 1, 0)
    

def findGridCenters(rateMap, xedges, yedges, thresholdFun, threshold, minFieldArea):
    '''
    Segment the single neuron rate map into grid fields and find their centers.
    Only contiguous areas which are larger than 'minFieldArea' are considered
    grid fields (assuming pixels of area 1)
    '''
    # Segmentation
    dx = np.abs(xedges[1] - xedges[0])
    rateMapTh = thresholdFun(rateMap, threshold)
    segments, nseg = label(rateMapTh, structure=None) # 4-neighborhood

    field_ctr_x = []
    field_ctr_y = []

    # Filter smaller segments and compute center of mass
    for seg_it in range(1, nseg+1):
        seg_y, seg_x = np.nonzero(segments == seg_it)
        if len(seg_x) < minFieldArea:
            continue

        #totalRates = np.sum(rateMap[seg_y, seg_x])
        rates = rateMap[seg_y, seg_x]
        max_i = np.argmax(rates)
        field_ctr_x.append(xedges[seg_x[max_i]])
        field_ctr_y.append(yedges[seg_y[max_i]])

        #field_ctr_x.append(1./totalRates * np.dot(xedges[seg_x], rates))
        #field_ctr_y.append(1./totalRates * np.dot(yedges[seg_y], rates))
        #field_ctr_x.append(np.mean(xedges[seg_x]))
        #field_ctr_y.append(np.mean(yedges[seg_y]))

    return field_ctr_x, field_ctr_y, segments


def timeCentersToPositions(time_centers, pos_x, pos_y, field_ctr_x, field_ctr_y,
        velocityStart, gridSep):
    # Recompute time_centers into positions
    time_centers -= velocityStart # relative to the start of velocity modulation
    centers_pos_i = np.array(time_centers // rat_dt, dtype=int)
    centers_pos_x = pos_x[centers_pos_i]
    centers_pos_y = pos_y[centers_pos_i]
    dist_min_grids = []
    for c_it in xrange(len(centers_pos_i)):
        dist =  np.min(np.sqrt((centers_pos_x[c_it] - field_ctr_x)**2 +
                (centers_pos_y[c_it] - field_ctr_y)**2))
        dist_min_grids.append(dist)
    dist_min_grids = np.array(dist_min_grids)
    dist_th = np.where(dist_min_grids > gridSep/2)
    dist_min_grids[dist_th] = gridSep - dist_min_grids[dist_th]
    return dist_min_grids, centers_pos_x, centers_pos_y


def distanceMeanMedian(dDistBins, dist_min_grids, sig_areas):
    distBins = np.arange(-dDistBins, np.max(dist_min_grids) + dDistBins,
            dDistBins)
    sig_areas_dist_binned_i = np.digitize(dist_min_grids, distBins)
    sig_areas_dist_mean   = np.ndarray(len(distBins))
    sig_areas_dist_median = np.ndarray(len(distBins))
    sig_areas_dist_std    = np.ndarray(len(distBins))
    for it in range(len(distBins)):
        sig_bin = sig_areas[np.where(sig_areas_dist_binned_i == it)]
        sig_areas_dist_mean[it]   = np.mean(sig_bin)
        sig_areas_dist_median[it] = np.median(sig_bin)
        sig_areas_dist_std[it]    = np.std(sig_bin)/np.sqrt(len(sig_bin))
    return distBins[2:]-dDistBins, sig_areas_dist_mean[2:], sig_areas_dist_median[2:], sig_areas_dist_std[2:]


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

    print "Data loaded..."

    Ne_x            = data['options']['Ne'][0][0][0][0]
    Ni_x            = data['options']['Ni'][0][0][0][0]
    gridSep         = data['options']['gridSep'][0][0][0][0]
    neuronNum_e     = Ne_x/2 - 1 
    neuronNum_i     = Ni_x/2 - 1
    pos_x           = data['pos_x'].ravel()
    pos_y           = data['pos_y'].ravel()
    rat_dt          = data['dt'][0][0]
    velocityStart   = data['velocityStart'].ravel()[0] * ms

    spikeTimes_e    = data['spikeCell_e'][neuronNum_e][0].ravel()
    spikes_e        = spikeTimes_e - velocityStart
    spikes_e        = np.delete(spikes_e, np.nonzero(spikes_e < 0)[0])

    spikeTimes_i    = data['spikeCell_i'][neuronNum_i][0].ravel()
    spikes_i        = spikeTimes_i - velocityStart
    spikes_i        = np.delete(spikes_i, np.nonzero(spikes_i < 0)[0])

    times           = data['stateMon_times'].ravel()
    rec_dt          = data['options']['stateRec_dt'][0][0][0][0] * ms

    syn_current_e   = data['stateMon_Iclamp_e_values'][0, :].ravel()/pA
    syn_current_i   = data['stateMon_Iclamp_i_values'][0, :].ravel()/pA


    # Segment the rate map 
    # E neuron
    rateMap_e, xedges, yedges  = SNSpatialRate2D(spikes_e, pos_x, pos_y, rat_dt, arenaDiam, h)
    field_ctr_x_e, field_ctr_y_e, segments_e = findGridCenters(rateMap_e, xedges,
            yedges, getRateMapThreshold_e, rateThreshold_e, minFieldArea)
    sig_areas_e, time_centers_e =  getCurrentChargeTheta(syn_current_e, times, theta_freq)
    dist_min_grids_e, centers_pos_x_e, centers_pos_y_e =  \
        timeCentersToPositions(time_centers_e, pos_x, pos_y, field_ctr_x_e,
                    field_ctr_y_e, velocityStart, gridSep)


    # Charge dependent on distance from grid field
    figure()
    plot(dist_min_grids_e, sig_areas_e, 'o')
    xlabel('Distance from grid field centre (cm)')
    ylabel('Theta cycle synaptic charge (pC)')
    # mean, median
    dDistBins = 2
    distBins_e, sig_areas_dist_mean_e, sig_areas_dist_median_e, sig_areas_dist_std_e = \
            distanceMeanMedian(dDistBins, dist_min_grids_e, sig_areas_e)

    figure()
    subplot(2, 1, 1)
    plot(distBins_e, sig_areas_dist_mean_e)
    ylabel('Mean charge (pC)')
    subplot(2, 1, 2)
    plot(distBins_e, sig_areas_dist_median_e)
    ylabel('Median charge (pC)')
    xlabel('Distance from the centre of grid field (cm)')


    # I neurons - distance from E grid field centers
    sig_areas_i, time_centers_i = getCurrentChargeTheta(syn_current_i, times,
            theta_freq)
    # !!! dist. from E centers here !!!
    dist_min_grids_i, centers_pos_x_i, centers_pos_y_i = \
            timeCentersToPositions(time_centers_i, pos_x, pos_y, field_ctr_x_e,
                    field_ctr_y_e, velocityStart, gridSep)

    dDistBins = 2
    distBins_i, sig_areas_dist_mean_i, sig_areas_dist_median_i, sig_areas_dist_std_i = \
            distanceMeanMedian(dDistBins, dist_min_grids_i, sig_areas_i)

    figure()
    subplot(2, 1, 1)
    plot(distBins_i, sig_areas_dist_mean_i)
    ylabel('Mean charge (pC)')
    subplot(2, 1, 2)
    plot(distBins_i, sig_areas_dist_median_i)
    ylabel('Median charge (pC)')
    xlabel('Distance from the centre of grid field (cm)')


    ## I neuron, distance from I 'grid' field center
    #rateMap_i, xedges, yedges = SNSpatialRate2D(spikes_i, pos_x, pos_y, rat_dt,
    #        arenaDiam, h)
    #field_ctr_x_i, field_ctr_y_i, segments_i = findGridCenters(rateMap_i, xedges,
    #        yedges, getRateMapThreshold_i, rateThreshold_i, minFieldArea)


    # HDF5 export here

    h5file = openFile(fileName + '_currents_position.h5', mode = "w")

    h5file.createArray(h5file.root, 'dist_min_grids_e', dist_min_grids_e)
    h5file.createArray(h5file.root, 'sig_areas_e', sig_areas_e)
    h5file.createArray(h5file.root, 'distBins_e', distBins_e)
    h5file.createArray(h5file.root, 'sig_areas_dist_mean_e',
            sig_areas_dist_mean_e)
    h5file.createArray(h5file.root, 'sig_areas_dist_median_e',
            sig_areas_dist_median_e)
    h5file.createArray(h5file.root, 'sig_areas_dist_std_e', sig_areas_dist_std_e)
    

    h5file.createArray(h5file.root, 'dist_min_grids_i', dist_min_grids_i)
    h5file.createArray(h5file.root, 'sig_areas_i', sig_areas_i)
    h5file.createArray(h5file.root, 'distBins_i', distBins_i)
    h5file.createArray(h5file.root, 'sig_areas_dist_mean_i',
            sig_areas_dist_mean_i)
    h5file.createArray(h5file.root, 'sig_areas_dist_median_i',
            sig_areas_dist_median_i)
    h5file.createArray(h5file.root, 'sig_areas_dist_std_i', sig_areas_dist_std_i)
    
    h5file.close()




    #figure()
    #X, Y = np.meshgrid(xedges, yedges)
    #pcolormesh(X, Y, segments_i)
    #hold('on')
    #plot(field_ctr_x_i, field_ctr_y_i, 'or')

    #figure()
    #pcolormesh(X, Y, rateMap_i)
    #plot(field_ctr_x_i, field_ctr_y_i, 'or')



#    # Charge dependent on average firing rate
#    dx_edge = xedges[1] - xedges[0]
#    center_x_i = np.array(centers_pos_x//dx_edge, dtype=int) + len(xedges)/2
#    center_y_i = np.array(centers_pos_y//dx_edge, dtype=int) + len(xedges)/2
#    if (any(center_x_i < 0) or any(center_y_i < 0)):
#        raise Exception()
#    centers_rates = rateMap[center_y_i, center_x_i]
#    # mean, median
#    figure()
#    dRateBins = 1
#    rateBins = np.arange(-dRateBins, np.max(centers_rates) + dRateBins, dRateBins)
#    sig_areas_binned_i = np.digitize(centers_rates, rateBins)
#    sig_areas_mean   = np.ndarray(len(rateBins), dtype=object)
#    sig_areas_median = np.ndarray(len(rateBins), dtype=object)
#    for it in range(len(rateBins)):
#        sig_bin = sig_areas[np.where(sig_areas_binned_i == it)]
#        sig_areas_mean[it]   = np.mean(sig_bin)
#        sig_areas_median[it] = np.median(sig_bin)
#    
#    subplot(2, 1, 1)
#    plot(rateBins, sig_areas_mean)
#    ylabel('Mean charge (pC)')
#    subplot(2, 1, 2)
#    plot(rateBins, sig_areas_median)
#    ylabel('Median charge (pC)')
#    xlabel('Firing rate (Hz)')
#
#
#    figure()
#    plot(centers_rates, sig_areas, 'o')
#    xlabel('Firing rate (Hz)')
#    ylabel('Theta cycle synaptic charge (pC)')
#
#
#    # Heat map
#    figure()
#    hist = np.ndarray((len(xedges), len(yedges)), dtype=object)
#    avgHist = np.ndarray((len(xedges), len(yedges)))
#    for it in xrange(len(sig_areas)):
#        x_bin = center_x_i[it]
#        y_bin = center_y_i[it]
#        if hist[y_bin, x_bin] == None:
#            hist[y_bin, x_bin] = []
#
#        hist[y_bin, x_bin].append(sig_areas[it])
#
#    for x_it in range(len(xedges)):
#        for y_it in range(len(yedges)):
#            if hist[y_it, x_it] == None:
#                avgHist[y_it, x_it] = 0
#            else:
#                avgHist[y_it, x_it] = np.mean(hist[y_it, x_it])
#
#pcolormesh(X, Y, avgHist)
#colorbar()
#hold('on')
#plot(pos_x, pos_y, color='white')



