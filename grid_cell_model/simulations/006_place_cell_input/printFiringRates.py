#
#   printFiringRates.py
#
#   Print firing rate plots of the simulated population of neurons for each of
#   the jobs selected.
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
import numpy
import pp
import scipy.io

from scipy.io           import loadmat
from scipy.optimize     import leastsq
from matplotlib.pyplot  import *

import spike_analysis



jobs = np.arange(9171, 9192)
genRange = [0, 9]

rcParams['font.size'] = 14


class FileOptions(object):
    dirName = "output_local/"
    fileNamePrefix = ''
    fileNameTemp = "{0}/{1}job{2:04}_gen{3:04}"


def processData(jobNum, gen_it, fo):
    print jobNum, gen_it
    fileName = fo.fileNameTemp +  '_output.mat'
    fileName = fileName.format(fo.dirName, fo.fileNamePrefix, jobNum, gen_it)
    try:
        data = scipy.io.loadmat(fileName, variable_names=["spikeCell_e",
            "spikeCell_i", "options"])
    except:
        print "warning: could not open: " + fileName
        return False
    
    spikeCell_e = data['spikeCell_e'].flatten()
    spikeCell_i = data['spikeCell_i'].flatten()
    T = int(data['options']['time'])
    
    F_tstart = 0
    F_tend = T*1e-3
    F_dt = 0.05
    F_winLen = 0.25
    Fe, Fe_t = spike_analysis.firingRateSlidingWindow(spikeCell_e, F_tstart, F_tend, F_dt, F_winLen) 
    Fi, Fi_t = spike_analysis.firingRateSlidingWindow(spikeCell_i, F_tstart, F_tend, F_dt, F_winLen) 

    from matplotlib.pyplot import figure, subplot, pcolormesh, xlabel, ylabel, savefig
    figure()
    subplot(2, 1, 1)
    T, FR = numpy.meshgrid(Fe_t, numpy.arange(len(Fe)))
    pcolormesh(T, FR, Fe)
    ylabel('E neuron #')
    subplot(2, 1, 2)
    T, FR = numpy.meshgrid(Fi_t, numpy.arange(len(Fi)))
    pcolormesh(T, FR, Fi)
    xlabel('Time (s)')
    ylabel('I neuron #')

    savefig('{0}/job{1:04}_gen{2:04}_firingRate.png'.format(fo.dirName, jobNum,
        gen_it))

    return True




fo = FileOptions()

parallel = True

if parallel:
    jobsParallel = []
    job_server = pp.Server(ppservers=('jupiter1', 'localhost', 'jupiter3'))
    for jobNum in jobs:
        for gen_it in xrange(genRange[0], genRange[1] + 1):
            jobsParallel.append(job_server.submit(
                    processData,
                    (jobNum, gen_it, fo),
                    (),
                    ("numpy", "scipy.io", "spike_analysis", "matplotlib.pyplot")))

    it = 0
    res = []
    for jobNum in jobs:
        for gen_it in xrange(genRange[0], genRange[1] + 1):
            res.append(jobsParallel[it]())
            it = it+1
                    
else:
    for jobNum in jobs:
        for gen_it in xrange(genRange[0], genRange[1] + 1):
            process(jobNum, gen_it, fo)

    


