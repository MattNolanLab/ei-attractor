#
#   fig_bump_drift.py
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


rcParams['font.size'] = 10


arenaDiam = 180.0     # cm
h = 3.0
Ne = 34

# Neuron to extract spikes from

def createIplaceLegend(Iplace_vec):
    leg = []
    for Iplace in Iplace_vec:
        leg.append(str(Iplace) + " pA")

    return leg

Iplace_all = range(0, 50, 10)
#Iplace_all = [0]
nTrials = 20
leg = createIplaceLegend(Iplace_all)

jobStart = 1000
trialNum = 0
dumpNum = 0

dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}_dump{4:03}"

bumpPos_std_all = []
bumpPos_var_all = []

jobNum = jobStart - 1
for Iplace in Iplace_all:
    bumpPos_all = []
    
    for trial_it in range(nTrials):
        jobNum += 1
        print 'jobNum: ' + str(jobNum)
    
        fileName = fileNameTemp +  '_output.mat'
        fileName = fileName.format(dirName, fileNamePrefix, jobNum, trialNum, dumpNum)
        try:
            data = loadmat(fileName)
        except:
            print "warning: could not open: " + fileName
            continue
    
        
        bumpPos = np.unwrap(data['bumpPos'][:, 0].ravel(), discont=Ne)
        #bumpPos = data['bumpPos'][:, 0].ravel()
        times   = data['bumpPos_times'].ravel()
        bumpPos_all.append(bumpPos)

    bumpPos_all = np.array(bumpPos_all).T
    bumpPos_std = np.std(bumpPos_all, 1) 
    bumpPos_var = np.var(bumpPos_all, 1)
    bumpPos_std_all.append(bumpPos_std)
    bumpPos_var_all.append(bumpPos_var)

    figure(figsize=(3,3))
    subplot(2, 1, 1)
    plot(times, bumpPos_all)
    ylabel('Bump position\n(neurons)')

    subplot(2, 1, 2)
    plot(times, bumpPos_var)
    ylabel('Variance\n(neurons)')
    xlabel('Time (s)')

    #title('Iplace = ' + str(Iplace) + ' pA')
    tight_layout()
    savefig('bump_drift_' + str(Iplace) + ' pA.eps')

figure(figsize=(4, 3))
plot(times, np.array(bumpPos_var_all).T)
xlabel('Time (s)')
ylabel('Variance (neurons^2)')
legend(leg, loc="best")
tight_layout()
savefig("fig_bump_drifts_Iplace.eps")
    
#show()
