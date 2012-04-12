#
#   bump_variance.py
#
#   Bump variance estimation
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
from scipy.io import loadmat
from matplotlib.pyplot import *

import numpy as np

from data_analysis import bumpVariancePos

# Generate bump stability figures based on simulated stationary attractor
# network


jobRange = [0, 49]
trialRange = [0, 24]
jobN = jobRange[1] - jobRange[0] + 1
trialN = trialRange[1] - trialRange[0] + 1

t_start = 3

dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}"

hist_nbins= 50

res = bumpVariancePos(dirName, fileNamePrefix, jobRange, trialRange, t_start)

for job_it in range(jobN):
    jobNum = job_it + jobRange[0]
    for trial_it in range(trialN):
        if res.pos_x_vec[job_it, trial_it] is None:
            continue
        trialNum = trialRange[0] + trial_it
        fileName = fileNameTemp.format(dirName, fileNamePrefix, jobNum,
                                    trialNum)
        pos_x = res.pos_x_vec[job_it, trial_it]
        pos_y = res.pos_y_vec[job_it, trial_it]
        Ne = res.o_vec[job_it, trial_it]['Ne'][0][0][0][0]

        #figure()
        #subplot(121)
        #hist(pos_x/Ne, hist_nbins, range=[-0.5, 0.5], normed=False)
        #xlabel('Position x (normalised)')
        #ylabel('Frequency')
        #subplot(122)
        #hist(pos_y/Ne, hist_nbins, range=[-0.5, 0.5], normed=False)
        #xlabel('Position y (normalised)')
        #savefig(fileName + '_bump_position_hist.pdf')
        #close(gcf())

job_its = np.arange(jobN)

figure()
hold(True)
errorbar(job_its, res.pos_var_mean[:, 0], res.pos_var_std[:, 0], fmt='o')
errorbar(job_its, res.pos_var_mean[:, 1], res.pos_var_std[:, 1], fmt='o')
margins(0.05, 0.05)
legend(('X', 'Y'))
xlabel('Network generation trial')
ylabel('Average bump variance (normalised)')
savefig(dirName + '/bump_position_variance_errbars.pdf')
hold(False)

figure()
subplot(211)
errorbar(job_its, res.pos_var_mean[:, 0], res.pos_var_std[:, 0], fmt='o')
margins(0.05, 0.05)
ylabel('Average bump variance X (normalised)')
subplot(212)
errorbar(job_its, res.pos_var_mean[:, 1], res.pos_var_std[:, 1], fmt='o')
margins(0.05, 0.05)
xlabel('Network generation trial')
ylabel('Average bump variance Y (normalised)')
savefig(dirName + '/bump_position_variance_errbars_sep.pdf')

figure()
hold(True)
plot(job_its, res.pos_var_mean[:, 0], 'o')
plot(job_its, res.pos_var_mean[:, 1], 'o')
margins(0.05, 0.05)
legend(('X', 'Y'))
xlabel('Network generation trial')
ylabel('Average bump variance (normalised)')
savefig(dirName + '/bump_position_variance_lines.pdf')
hold(False)

