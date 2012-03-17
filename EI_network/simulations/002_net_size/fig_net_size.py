from scipy.io import loadmat
from matplotlib.pyplot import *

import numpy as np

from data_analysis import bumpVariancePos

# Generate bump stability figures based on simulated stationary attractor
# network


jobRange = [100, 108]
trialRange = [0, 9]
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

        print('Histogram for ' + fileName)
        figure()
        subplot(121)
        hist(pos_x/Ne, hist_nbins, range=[-0.5, 0.5], normed=False)
        xlabel('Position x (normalised)')
        ylabel('Frequency')
        subplot(122)
        hist(pos_y/Ne, hist_nbins, range=[-0.5, 0.5], normed=False)
        xlabel('Position y (normalised)')
        savefig(fileName + '_bump_position_hist.pdf')
        close(gcf())


figure()
hold(True)
errorbar(np.array(res.pos_N)**2, res.pos_var_mean[:, 0], res.pos_var_std[:, 0], fmt='o-')
errorbar(np.array(res.pos_N)**2, res.pos_var_mean[:, 1], res.pos_var_std[:, 1], fmt='o-')
legend(('X', 'Y'))
xlabel('Network size E and I')
ylabel('Average bump variance (normalised)')
savefig(dirName + '/bump_position_variance_errbars.pdf')
hold(False)

figure()
subplot(211)
errorbar(np.array(res.pos_N)**2, res.pos_var_mean[:, 0], res.pos_var_std[:, 0], fmt='o-')
ylabel('Average bump variance X (normalised)')
subplot(212)
errorbar(np.array(res.pos_N)**2, res.pos_var_mean[:, 1], res.pos_var_std[:, 1], fmt='o-')
xlabel('Network size E and I')
ylabel('Average bump variance Y (normalised)')
savefig(dirName + '/bump_position_variance_errbars_sep.pdf')

figure()
hold(True)
plot(np.array(res.pos_N)**2, res.pos_var_mean[:, 0], 'o-')
plot(np.array(res.pos_N)**2, res.pos_var_mean[:, 1], 'o-')
legend(('X', 'Y'))
xlabel('Network size E and I')
ylabel('Average bump variance (normalised)')
savefig(dirName + '/bump_position_variance_lines.pdf')
hold(False)

