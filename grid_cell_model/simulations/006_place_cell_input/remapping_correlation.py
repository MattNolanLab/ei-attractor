#
#   remapping_correlation.py
#
#   Remapping grid cells and their population vector correlations
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
import numpy.ma as ma
from tables import *

from scipy.io import loadmat

import brian
from brian import *

from PlaceSpikingGroup import *


def skewness(X):
    mean = np.mean(X)
    std  = np.std(X)
    return 1.0/len(X)*np.sum((X - mean)**3) / std**3

def fieldPVCC(PC1, PC2):
    '''
    Compute population vector correlation coefficient for two place cell environments
    '''
    ff1, pos = PC1.getFiringFields(None, dx)
    ff2, pos = PC2.getFiringFields(None, dx)
    
    cc = np.ndarray(ff1.shape[1]*ff1.shape[2])
    it = 0
    for it_y in xrange(ff1.shape[1]):
        for it_x in xrange(ff1.shape[2]):
            cc[it] = np.corrcoef(ff1[:, it_y, it_x], ff2[:, it_y, it_x])[0, 1]
            it += 1

    return cc



boxSize = (200, 200)
N = (20, 20)
N_total = N[0]*N[1]
totalSz = N[0]*N[1]
maxRates = 15
widths = 20

dx = 5  # cm

numRepeat = 500


h5file = openFile('output_local/remapping_correlation.h5', mode = "w", title = "Place cell remapping examples")
cc_grp = h5file.createGroup(where='/', name='CrossCorrelations')

doExample = True
doPFSize  = True
doNetSize = True


if doExample:
    cc_it = 0
    cc_mean = []
    cc_skewness = []
    for rep in range(numRepeat):
        
        PC1 = UniformBoxPlaceCells(boxSize, N, maxRates, widths, random=True)
        PC2 = UniformBoxPlaceCells(boxSize, N, maxRates, widths, random=True)
        #N_same = 40
        #PC2.centers[0:N_same] = PC1.centers[0:N_same]
        
        # Histogram of corr. coeffiecnts (CC) for each spatial location in the box
        cc = fieldPVCC(PC1, PC2)
        
        #figure()
        #hist(cc, 50)
    
        
        print "Mean cross-correlation: ", np.mean(cc)
        print "XC distribution skewness: ", skewness(cc)
    
        cc_mean.append(np.mean(cc))
        cc_skewness.append(skewness(cc))
        
        h5file.createArray(cc_grp, "cc_{0}".format(cc_it), cc)
    
        cc_it += 1
    
    
    h5file.createArray(cc_grp, "cc_skewness", cc_skewness)
    h5file.createArray(cc_grp, "cc_mean", cc_mean)



################################################################################
#                             Place field size
################################################################################
if doPFSize:
    pfSizes = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 40, 60, 80, 100, 150, 200])
    pfSizeNumRepeat = 500
    
    cc_size = []    # A container for each size (2D array)
    pf_sizes = []
    for pfSize in pfSizes:
        print "pfSize: ", pfSize
        cc_PV_mean = [] # Mean for the whole population (one value for each remap trial)
        for rep_it in range(pfSizeNumRepeat):
            PC1 = UniformBoxPlaceCells(boxSize, N, maxRates, pfSize, random=True)
            PC2 = UniformBoxPlaceCells(boxSize, N, maxRates, pfSize, random=True)
        
            cc_PV_mean.append(np.mean(fieldPVCC(PC1, PC2)))
    
        cc_size.append(cc_PV_mean)
        pf_sizes.append(np.zeros(len(cc_PV_mean))+pfSize)
        
    
    h5file.createArray(h5file.root, "cc_size", cc_size)
    h5file.createArray(h5file.root, "pf_sizes", pf_sizes)
    h5file.createArray(h5file.root, "pfSizes", pfSizes)




################################################################################
#                             Network size
################################################################################
if doNetSize:
    netSizes = np.array([6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40])
    netSizeNumRepeat = 500

    cc_net_size = []
    net_size_all = []
    for netSize in netSizes:
        print "netSize: ", netSize
        cc_PV_mean = []
        netSize_N = (netSize, netSize)
        for rep_it in range(netSizeNumRepeat):
            PC1 = UniformBoxPlaceCells(boxSize, netSize_N, maxRates, widths, random=True)
            PC2 = UniformBoxPlaceCells(boxSize, netSize_N, maxRates, widths, random=True)

            print "PC1.N: ", PC1.N

            cc_PV_mean.append(np.mean(fieldPVCC(PC1, PC2)))

        cc_net_size.append(cc_PV_mean)
        net_size_all.append(np.zeros(len(cc_PV_mean)) + netSize)

    h5file.createArray(h5file.root, "cc_net_size", cc_net_size)
    h5file.createArray(h5file.root, "net_size_all", net_size_all)
    h5file.createArray(h5file.root, "netSizes", netSizes)


h5file.close()

