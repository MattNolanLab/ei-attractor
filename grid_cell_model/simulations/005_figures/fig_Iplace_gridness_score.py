#
#   fig_Iplac_gridness_score.py
#
#   Figure of gridness score depending on Iplace
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
from tables             import *

from grid_cell_analysis import *


job_e = 3400
job_e_no_place = 200
job_i = 3500

Iplace = [50, 100, 150, 200, 250]
Ntrials = 10

rcParams['font.size'] = 14


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/job{1:04}_gridness_scores__excitatory"

# vel. modulation of E cells
data_e = loadmat(fileNameTemp.format(dirName, job_e) +  '.mat')
gridnessScores_e = np.reshape(data_e['gridnessScores'].ravel(), (len(Iplace),
    Ntrials))
data_e_no_place = loadmat(fileNameTemp.format(dirName, job_e_no_place) +  '.mat')
gridnessScores_e_no_place = data_e_no_place['gridnessScores'].ravel()


data_i = loadmat(fileNameTemp.format(dirName, job_i) +  '.mat')
gridnessScores_i = np.reshape(data_i['gridnessScores'].ravel(), (len(Iplace),
    Ntrials))

Iplace_individual = []
avg_gridness_e = []
avg_gridness_i = []
std_gridness_e = []
std_gridness_i = []

figure()
subplot(2, 1, 1)
hold('on')
for it in range(len(Iplace)):
    plot([Iplace[it]]*Ntrials, gridnessScores_e[it, :], '.', color=[0.1, 0.1, 0.1])
    avg = np.mean(gridnessScores_e[it, :])
    std = np.std(gridnessScores_e[it, :])
    avg_gridness_e.append(avg)
    std_gridness_e.append(std)
    errorbar(Iplace[it], avg, std, fmt='ob')

plot([0]*len(gridnessScores_e_no_place), gridnessScores_e_no_place, '.', color=[0.1, 0.1, 0.1])
errorbar(0, np.mean(gridnessScores_e_no_place),
        np.std(gridnessScores_e_no_place), fmt='ob')
xlim([-25, Iplace[-1]+50])
ylabel('Gridness score')
title('Vel. modulation onto E cells')

subplot(2, 1, 2)
hold('on')
for it in range(len(Iplace)):
    plot([Iplace[it]]*Ntrials, gridnessScores_i[it, :], '.', color=[0.1, 0.1, 0.1])
    avg = np.mean(gridnessScores_i[it, :])
    std = np.std(gridnessScores_i[it, :])
    avg_gridness_i.append(avg)
    std_gridness_i.append(std)
    Iplace_individual += [Iplace[it]]*Ntrials
    errorbar(Iplace[it], avg, std, fmt='ob')
xlim([-25, Iplace[-1]+50])
xlabel('Place cell input amp. (pA)')
ylabel('Gridness score')
title('Vel. modulation onto I cells')

savefig(dirName + '/fig_gridness_Iplace.pdf')


h5file = openFile(dirName + '/fig_gridness_Iplace.h5', mode = "w")

h5file.createArray(h5file.root, 'Iplace', Iplace)
h5file.createArray(h5file.root, 'Iplace_individual', Iplace_individual)
h5file.createArray(h5file.root, 'Iplace_individual_e_0',
        [0]*len(gridnessScores_e_no_place))


h5file.createArray(h5file.root, 'gridnessScores_e', gridnessScores_e.ravel())
h5file.createArray(h5file.root, 'gridnessScores_e_0', gridnessScores_e_no_place)
h5file.createArray(h5file.root, 'gridnessScores_i', gridnessScores_i.ravel())

h5file.createArray(h5file.root, 'avg_gridness_e', avg_gridness_e)
h5file.createArray(h5file.root, 'avg_gridness_e_0',
        np.mean(gridnessScores_e_no_place))
h5file.createArray(h5file.root, 'avg_gridness_i', avg_gridness_i)
h5file.createArray(h5file.root, 'std_gridness_e', std_gridness_e)
h5file.createArray(h5file.root, 'std_gridness_i', std_gridness_i)
h5file.createArray(h5file.root, 'std_gridness_e_0',
        np.std(gridnessScores_e_no_place))


h5file.close()
