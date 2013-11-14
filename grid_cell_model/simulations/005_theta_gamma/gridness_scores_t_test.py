#
#   gridness_scores_t_test.py
#
#   T-test of gridness score of two populations
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
import scipy as sp

from scipy.io           import loadmat
from scipy.io           import savemat
from matplotlib.pyplot  import *

from scipy.stats        import ttest_ind


jobNum1 = 3100
spikeType1 = 'excitatory'
jobNum2 = 200
spikeType2 = 'excitatory'
fileName1 = 'job{0:04}_gridness_scores_'.format(jobNum1) + '_' + spikeType1 + '.mat'
fileName2 = 'job{0:04}_gridness_scores_'.format(jobNum2) + '_' + spikeType2 + '.mat'

data1 = loadmat(fileName1)
data2 = loadmat(fileName2)


GS1 = data1['gridnessScores'].flatten()
GS2 = data2['gridnessScores'].flatten()

t, p = ttest_ind(GS1, GS2)
print 'GS stellate:'
print GS1
print 'GS interneuron:'
print GS2
print fileName1, ': mean GS: ', np.mean(GS1)
print fileName2, ': mean GS: ', np.mean(GS2)
print fileName1, ': sem GS: ', sp.stats.sem(GS1)
print fileName2, ': sem GS: ', sp.stats.sem(GS2)
print 'p: ', p

