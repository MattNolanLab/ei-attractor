#
#   bump_velocity_response.py
#
#   Estimation of bump velocity reponse to velocity current injection.
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



jobRange_all = [
        [100, 119],
        [120, 139],
        [140, 159],
        [160, 179],
        [180, 199],
        [200, 219],
        [220, 249],
        [240, 259]]
fitLine = [
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False]
genRange = [0, 9]

leg = (
        '0 nrns',
        '1 nrns',
        '2 nrns',
        '3 nrns',
        '4 nrns',
        '5 nrns',
        '6 nrns',
        '7 nrns')

rcParams['font.size'] = 14



def fitCircularSlope(bumpPos, times, normFac):
    '''
    Fit a (circular) slope line onto velocity response of the bump and extract
    the slope (velocity), in neurons/s
    '''
    t = np.array(times) - times[0]
    bumpPos_norm = np.unwrap(1.0 * bumpPos / normFac * 2 * np.pi) # normalise to 2*Pi
    func = lambda X: X[0]*t - bumpPos_norm
    x0 = np.array([0.0])  # slope
    x = leastsq(func, x0)
    return x[0][0] / 2. / np.pi * normFac
    

def getLineFit(Y):
    '''
    Fit a line to data
    '''
    X = np.arange(len(Y))
    
    func = lambda P: P[0]*X  - Y
    P0 = np.array([0.0]) # slope
    P = leastsq(func, P0)
    return P[0][0]*X, P[0][0]
    

figure()
subplot(111)
hold('on')

for Ivel_types in range(len(jobRange_all)):
    jobRange = jobRange_all[Ivel_types]
    print "jobRange: " + str(jobRange)

    jobN = jobRange[1] - jobRange[0] + 1
    genN = genRange[1] - genRange[0] + 1
    
    sheetSize = 68  # !!! MUST be set properly
    
    t_start = 3
    
    dirName = "output/"
    fileNamePrefix = ''
    fileNameTemp = "{0}/{1}job{2:04}_gen{3:04}"
    
    
    res         = []
    gen_avg     = []
    gen_stderr  = []
    
    for job_it in range(jobN):
        jobNum = job_it + jobRange[0]
    
        res.append([])
        for gen_it in range(genN):
            genNum   = genRange[0] + gen_it
            fileName = fileNameTemp +  '_output.mat'
            fileName = fileName.format(dirName, fileNamePrefix, jobNum,
                                        genNum)
            try:
                data = loadmat(fileName)
            except:
                print "warning: could not open: " + fileName
                continue
    
    
            bumpPos = data['bumpPos'][:, 0].ravel()
            times   = data['bumpPos_times'].ravel()
            normFac = sheetSize
            res[job_it].append(fitCircularSlope(bumpPos, times, normFac))
    
        gen_avg.append(np.mean(res[job_it]))
        gen_stderr.append(np.std(res[job_it]) / np.sqrt(len(res[job_it])))
    
    
    # TODO: export this in simulations
    Ivel = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190]
    errorbar(Ivel, gen_avg, gen_stderr, fmt='o-')
    xlabel('Velocity current (pA)')
    ylabel('Bump velocity (neurons/s)')

    if fitLine[Ivel_types]:
        line, slope = getLineFit(gen_avg)
        slope = slope/(Ivel[1] - Ivel[0])
        plot(Ivel, line)
        title("Line fit slope: " + str(slope) + ' nrns/s/pA')
    

legend(leg, loc='best')
savefig(dirName + '/bump_velocity_modulation.pdf')

