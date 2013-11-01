#
#   fig_gridness_field_size.py
#
#   Figure of gridness score depending on spacing between grid fields
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
from scipy              import stats
from scipy.stats        import sem
from matplotlib.pyplot  import *
from tables             import *

from grid_cell_analysis import *


gridSeps = [40, 50, 60, 70]
Ntrials = 10

jobs_e_start = [3600, 3800, 3810, 3100]
jobs_i_start = [3700, 3900, 3910, 3200]

rcParams['font.size'] = 14


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/job{1:04}_gridness_scores__excitatory"

figure(figsize=(5, 6))
subplot(2, 1, 1)
hold('on')

gridSeps_e_ind = []
gridness_e_ind = []
gridSeps_i_ind = []
gridness_i_ind = []
avg_gridness_e = []
std_gridness_e = []
avg_gridness_i = []
std_gridness_i = []
    
# Modulation onto E cells
for gridSep_it in range(len(gridSeps)):
    data = loadmat(fileNameTemp.format(dirName, jobs_e_start[gridSep_it]) +  '.mat')
    GS = data['gridnessScores'].ravel()
    
    gSep = [gridSeps[gridSep_it]]*len(GS) 
    plot(gSep, GS, '.', color=[0.1, 0.1, 0.1])

    gridSeps_e_ind += gSep
    gridness_e_ind += list(GS)
    avg = np.mean(GS)
    std = sem(GS)
    avg_gridness_e.append(avg)
    std_gridness_e.append(std)
    errorbar(gridSeps[gridSep_it], avg, std, fmt='ob')
    xlim([30, 80])
    ylim([0.7, 1.3])
    ylabel('Gridness score')
    gca().spines['right'].set_visible(False)
    gca().spines['top'].set_visible(False)
    #xlabel('Grid field spacing (cm)')
    #title('Vel. modulation onto E cells')
    

subplot(2, 1, 2)
hold('on')
# Modulation onto I cells
for gridSep_it in range(len(gridSeps)):
    data = loadmat(fileNameTemp.format(dirName, jobs_i_start[gridSep_it]) +  '.mat')
    GS = data['gridnessScores'].ravel()
    
    gSep = [gridSeps[gridSep_it]]*len(GS) 
    plot(gSep, GS, '.', color=[0.1, 0.1, 0.1])


    gridSeps_i_ind += gSep
    gridness_i_ind += list(GS)
    avg = np.mean(GS)
    std = sem(GS)
    avg_gridness_i.append(avg)
    std_gridness_i.append(std)
    errorbar(gridSeps[gridSep_it], avg, std, fmt='ob')
    xlim([30, 80])
    ylim([0.7, 1.3])
    xlabel('Grid field spacing (cm)')
    ylabel('Gridness score')
    gca().spines['right'].set_visible(False)
    gca().spines['top'].set_visible(False)
    #title('Vel. modulation onto I cells')
    
gcf().subplots_adjust(hspace=0.3)
gcf().subplots_adjust(left=0.2)

savefig('fig_gridness_field_size.png')


h5file = openFile('fig_gridness_field_size.h5', mode = "w")

h5file.createArray(h5file.root, 'gridSeps', gridSeps)
h5file.createArray(h5file.root, 'gridSeps_e_ind', gridSeps_e_ind)
h5file.createArray(h5file.root, 'gridSeps_i_ind', gridSeps_i_ind)
h5file.createArray(h5file.root, 'gridness_e_ind', gridness_e_ind)
h5file.createArray(h5file.root, 'gridness_i_ind', gridness_i_ind)
h5file.createArray(h5file.root, 'avg_gridness_e', avg_gridness_e)
h5file.createArray(h5file.root, 'std_gridness_e', std_gridness_e)
h5file.createArray(h5file.root, 'avg_gridness_i', avg_gridness_i)
h5file.createArray(h5file.root, 'std_gridness_i', std_gridness_i)

h5file.close()
