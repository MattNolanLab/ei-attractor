#
#   fig_bump_velocity_response.py
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
        [5000, 5015]]
fitLine = [
        True]
genRange = [0, 9]

leg = (
        '4 nrns',)

stderr_th = 5   # neurons

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
    Ivel = range(0, 160, 10)
    errorbar(Ivel, gen_avg, gen_stderr, fmt='o-')
    xlabel('Velocity current (pA)')
    ylabel('Bump velocity (neurons/s)')

    print gen_stderr

    if fitLine[Ivel_types]:
        ## Use data points only until the std.error is reasonably small
        #it = 0
        #while it < len(Ivel):
        #    if gen_stderr[it] > stderr_th:
        #        break
        #    it += 1
        it = 10
        line, slope = getLineFit(gen_avg[0:it])
        slope = slope/(Ivel[1] - Ivel[0])
        plot(Ivel[0:it], line)
        title("Line fit slope: " + str(slope) + ' nrns/s/pA')
    

legend(leg, loc='best')
savefig('{0}/job{1:04}_bump_velocity_estimation.pdf'.format(dirName, jobRange[0]))

