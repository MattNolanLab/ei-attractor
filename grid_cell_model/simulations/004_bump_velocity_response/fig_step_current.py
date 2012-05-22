#
#   fig_step_current.py
#
#   Step current response figures.
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
from scipy.optimize     import leastsq
from matplotlib.pyplot  import *



jobRange = [1000, 1015]
genRange = [0, 9]

leg = (
        '4 nrns')

rcParams['font.size'] = 14


dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_gen{3:04}"

sheetSize = 68  # !!! MUST be set properly



jobN = jobRange[1] - jobRange[0] + 1
genN = genRange[1] - genRange[0] + 1

t_start = 1
dt = 0.1e-3


for job_it in range(jobN):
    jobNum = job_it + jobRange[0]

    bumpPos_all = []

    figure()
    suptitle('Vel. I: ' + str(Ivel) + ' pA')
    subplot2grid((4, 1), (0, 0), rowspan=3)
    hold('on')

    for gen_it in range(genN):
        genNum   = genRange[0] + gen_it

        fileName = fileNameTemp
        fileName = fileName.format(dirName, fileNamePrefix, jobNum,
                                    genNum)
        try:
            data = loadmat(fileName + '_output.mat')
        except:
            print "warning: could not open: " + fileName
            continue


        bumpPos           = data['bumpPos'][:, 0].ravel()
        bump_times        = data['bumpPos_times'].ravel()
        Iext_vel_e_times  = data['stateMon_Iext_vel_e_times'].ravel()
        Iext_vel_e_values = data['stateMon_Iext_vel_e_values'][0, :]
        Ivel              = data['Ivel'][0][0]
        stepStart_t       = data['stepStart_t'][0][0]
        stepEnd_t         = data['stepEnd_t'][0][0]
        bumpPos_all.append(bumpPos)

        plot(bump_times, bumpPos, color='black', alpha=0.15)
        hold('on')
    
    plot(bump_times, np.mean(bumpPos_all, 0), color='black', linewidth=2)
    ylim([-sheetSize/2, sheetSize/2])
    axvline(x=stepStart_t/1e3, linestyle='--', color='grey')
    axvline(x=stepEnd_t/1e3, linestyle='--', color='grey')
    ylabel('Bump position (neurons)')

    subplot2grid((4, 1), (3, 0))
    plot(Iext_vel_e_times[t_start/dt:], Iext_vel_e_values[t_start/dt:]*1e12)
    ylabel('Velocity current (pA)')
    xlabel('time (s)')
    #gcf().subplots_adjust(bottom=0.4)
    tight_layout()

    savefig(fileName + '_step_current.pdf')


