#!/usr/bin/env python
#
#   fig_last_spike.py
#
#   Plot histograms of the last E spike time.
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
import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LinearLocator, MaxNLocator, \
        ScalarFormatter

from parameters  import JobTrialSpace2D
from figures.EI_plotting import plotGridTrial
from plotting.global_defs import globalAxesSettings

import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
plt.rcParams['font.size'] = 11

noise_sigma_vec = [0, 150, 300]
NTrials = 10
dirs = ('EI_param_sweep_{0}pA',    (30, 30))
varList = ['spikeMon_e', 'events', 'times']
nbins = 40


    

def aggregateLastSpike(sp, varList, trialNumList):
    func = lambda x: x[-1]
    return sp.aggregateData(varList, trialNumList, funReduce=func,
            saveData=True).flatten()


################################################################################
plt.figure(figsize=(2.9, 2.9))
ax = plt.gca()                                                                                                               
plt.hold('on')                                                                                                               
globalAxesSettings(ax)                                                                                                       

for noise_sigma in noise_sigma_vec:
    dir = dirs[0].format(int(noise_sigma))
    rootDir = "output/grids/{0}".format(dir)
    shape   = dirs[1]
    sp = JobTrialSpace2D(shape, rootDir)
    lastSpikeT = aggregateLastSpike(sp, varList, range(NTrials))
    filtIdx = np.logical_not(np.isnan(lastSpikeT))
    ax.hist(lastSpikeT[filtIdx]*1e-3, bins=nbins, histtype='step', align='mid')

leg = []                                                                                                                     
for s in noise_sigma_vec:                                                                                                        
    leg.append("{0}".format(int(s)))                                                                                         
ax.legend(leg, loc='best', title='$\sigma$ (pA)', frameon=False,                                                         
        fontsize='small', ncol=2)                                                                                            
                                                                                                                             
ax.set_xlabel("Gridness score")                                                                                              
ax.set_ylabel('Count')
                                                                                                                             
ax.spines['top'].set_visible(False)                                                                                          
ax.spines['right'].set_visible(False)                                                                                        
#ax.xaxis.set_major_locator(MultipleLocator(50))                                                                             
#ax.yaxis.set_major_locator(MultipleLocator(50))                                                                              
#ax.xaxis.set_minor_locator(AutoMinorLocator(2))                                                                              
#ax.yaxis.set_minor_locator(AutoMinorLocator(3))                                                                              
ax.margins(0.05, 0.025)                    


###############################################################################
plt.tight_layout()
dir = dirs[0].format(0)
rootDir = "output/grids/{0}".format(dir)
plt.savefig(rootDir + '/../last_spike_T.pdf')

