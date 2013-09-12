#!/usr/bin/env python
#
#   plot_grids.py
#
#   Grid fields in the EI parameter sweep.
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
import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt

from parameters  import JobTrialSpace2D
from figures.EI_plotting import plotGridTrial
import logging as lg
#lg.basicConfig(level=lg.WARN)
lg.basicConfig(level=lg.INFO)


# Other
plt.rcParams['font.size'] = 11

noise_sigma = 0
dirs = \
    ('EI_param_sweep_{0}pA',    (30, 30))

NTrials = 1
dir = dirs[0].format(int(noise_sigma))
rootDir = "output/grids_50pA_Ivel/{0}".format(dir)
shape   = dirs[1]

sp = JobTrialSpace2D(shape, rootDir)
iterList  = ['g_AMPA_total', 'g_GABA_total']
    
################################################################################
plt.figure(figsize=(2.9, 2.9))
G = plotGridTrial(sp, ['gridnessScore'], iterList,
        trialNumList=range(NTrials),
        r=0,
        c=5,
        xlabel="I (nS)",
        ylabel='E (nS)',
        clBarLabel = "Gridness score",
        clbarNTicks=3)
###############################################################################
plt.tight_layout()
plt.savefig(sp.rootDir +
        '/../analysis_EI_gridness_{0}pA.png'.format(int(noise_sigma)))

