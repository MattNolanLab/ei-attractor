#
#   fig_gridness_Iplace.py
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
import pp

from scipy.io           import loadmat
from scipy.io           import savemat
from scipy              import stats
from scipy.stats        import sem
from matplotlib.pyplot  import *
from tables             import *

from analysis.grid_cells import SNSpatialRate2D, cellGridnessScore


jobStart = 2000
loadGridness = True

Iplace_vec = [0, 10, 20, 30, 40]
Ntrials = 10

trialNum = 0
dumpNum = 7

arenaDiam = 180.0     # cm
h = 3.0

# Neuron to extract spikes from
neuronNum = 10
spikeType = 'excitatory'


rcParams['font.size'] = 10


dirName = "output/"
fileNamePrefix = ''
spikesFileNameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"
gridnessFileNameTemp = "{0}/job{1:04}_gridness_scores_excitatory"


def loadAndCalcGridness(spikesFileNameTemp, dirName, fileNamePrefix, jobNum,
        trialNum, dumpNum, spikeType, neuronNum, arenaDiam, h):

    print 'jobNum: ' + str(jobNum)

    spikesFileName = spikesFileNameTemp
    spikesFileName = spikesFileName.format(dirName, fileNamePrefix, jobNum,
                                trialNum, dumpNum)
    try:
        data = scipy.io.loadmat(spikesFileName +  '_output.mat')
    except:
        print "warning: could not open: " + spikesFileName + '_output.mat'
        return (None, None)
    
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
    spikes = numpy.delete(spikes, numpy.nonzero(spikes < 0)[0])
    
    rateMap, xedges, yedges = SNSpatialRate2D(spikes, pos_x, pos_y, rat_dt, arenaDiam, h)
    maxRate = numpy.max(rateMap.flatten())
    G, crossCorr, angles = cellGridnessScore(rateMap, arenaDiam, h, corr_cutRmin)
    return (G, maxRate)



gridnessScores = np.ndarray((len(Iplace_vec)), dtype=object)
maxRates = np.ndarray((len(Iplace_vec)), dtype=object)

jobs = []

# vel. modulation of E cells
if loadGridness == True:
    data = loadmat(gridnessFileNameTemp.format(dirName, jobStart) +  '.mat')
    gridnessScores = data['gridnessScores'].flatten()
    maxRates = data['maxRates'].flatten()
else:
    job_server = pp.Server(ppservers=())

    jobNum = jobStart - 1
    for Iplace_it in range(len(Iplace_vec)):
        Iplace = Iplace_vec[Iplace_it]
        jobs.append([])
        for nTrial in range(Ntrials):
            jobNum += 1
            jobs[Iplace_it].append(job_server.submit(
                loadAndCalcGridness,
                (spikesFileNameTemp, dirName, fileNamePrefix, jobNum,
                    trialNum, dumpNum, spikeType, neuronNum, arenaDiam, h),
                (),
                ("numpy", "scipy.io", "grid_cell_analysis")))

    for Iplace_it in range(len(Iplace_vec)):
        gridnessScores[Iplace_it] = []
        maxRates[Iplace_it] = []
        for nTrial in range(Ntrials):
            G, maxRate = jobs[Iplace_it][nTrial]()
            if G is not None:
                gridnessScores[Iplace_it].append(G)
                maxRates[Iplace_it].append(maxRate)

    gridnessDataOut = {
            "Iplace_vec"        : Iplace_vec,
            "Ntrials"           : Ntrials,
            "gridnessScores"    : gridnessScores,
            "maxRates"          : maxRates
    }
    savemat(gridnessFileNameTemp.format(dirName, jobStart) +  '.mat', gridnessDataOut,
            oned_as='row')


avg_gridness = []
std_gridness = []
avg_rate = []
std_rate = []

figure(figsize=(3, 3))
subplot(2, 1, 1)
hold('on')
for it in range(len(Iplace_vec)):
    G_vec = np.array(gridnessScores[it]).flatten()
    plot([Iplace_vec[it]]*len(G_vec), G_vec, '.', color=[0.1, 0.1, 0.1])
    avg = np.mean(G_vec)
    std = sem(G_vec)
    avg_gridness.append(avg)
    std_gridness.append(std)

errorbar(Iplace_vec, avg_gridness, std_gridness, fmt='o-b')


xlim([-5, Iplace_vec[-1]+5])
ylabel('Gridness score')
#subplots_adjust(left=0.2)

subplot(2, 1, 2)
hold('on')
for it in range(len(Iplace_vec)):
    rate_vec = np.array(maxRates[it]).flatten()
    plot([Iplace_vec[it]]*len(rate_vec), rate_vec, '.', color=[0.1, 0.1, 0.1])
    avg = np.mean(rate_vec)
    std = sem(rate_vec)
    avg_rate.append(avg)
    std_rate.append(std)

errorbar(Iplace_vec, avg_rate, std_rate, fmt='o-b')
ylabel('Max. firing rate (Hz)')
xlabel('Place cell input amplitude (pA)')
xlim([-5, Iplace_vec[-1]+5])
tight_layout()
savefig('fig_gridness_Iplace.eps')

show()

#savefig(dirName + '/fig_gridness_Iplace.pdf')

#h5file = openFile(dirName + '/fig_gridness_Iplace.h5', mode = "w")
#
#h5file.createArray(h5file.root, 'Iplace', Iplace)
