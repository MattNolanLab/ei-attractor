from scipy.io import loadmat
from matplotlib.pyplot import *

import numpy as np

# Generate bump stability figures based on simulated stationary attractor
# network

def circularVariance(x, range):
    ''' Returns circular moments of a vector of circular variable x, defined in
    the 'range' (this will be mapped to a circle)'''
    c = np.exp(1j*2*np.pi*x/range)
    avg = np.mean(c)
    theta_avg = np.angle(avg)
    theta_var = 1 - np.abs(avg)
    return (theta_avg, theta_var)
    

jobRange = [0, 6]
trialRange = [0, 9]

dirName = "output/"
fileNamePrefix = ''
fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}"

jobN = jobRange[1] - jobRange[0]+1
pos_var = []
pos_N = []
pos_var_mean = np.ndarray((jobN, 2))
pos_var_std = np.ndarray((jobN, 2))

hist_nbins= 50

for job_it in range(jobN):
    jobNum = job_it + jobRange[0]
    pos_var.append([])
    for trialNum in range(trialRange[0], trialRange[1]+1):
        try:
            fileName = fileNameTemp.format(dirName, fileNamePrefix, jobNum,
                trialNum)
            print "Processing file: " + fileName + "_output.mat"
            data = loadmat(fileName + "_output.mat")
        except:
            print "Warning: could not open: " + fileName
            continue

        o = data['options']
        #spikes = data['spikeCell_e']
        pos_x = data['bumpPos'][:, 0]
        pos_y = data['bumpPos'][:, 1]
        
        Ne = o['Ne'][0][0][0][0] # Ugly but these bastards cannot export to
                                 # matlab
        mean_x, var_x = circularVariance(pos_x, Ne)
        mean_y, var_y = circularVariance(pos_y, Ne)
        pos_var[job_it].append([var_x, var_y])

        figure
        subplot(121)
        hist(pos_x/Ne, hist_nbins, range=[-0.5, 0.5], normed=False)
        xlabel('Position x (normalised)')
        ylabel('Frequency')
        subplot(122)
        hist(pos_y/Ne, hist_nbins, range=[-0.5, 0.5], normed=False)
        xlabel('Position y (normalised)')
        savefig(fileName + '_bump_position_hist.pdf')
        close(gcf())

    pos_var_mean[job_it, :] = np.mean(pos_var[job_it], 0)
    pos_var_std[job_it, :] = np.std(pos_var[job_it], 0)
    pos_N.append(Ne)

figure()
hold(True)
errorbar(np.array(pos_N)**2, pos_var_mean[:, 0], pos_var_std[:, 0], fmt='o-')
errorbar(np.array(pos_N)**2, pos_var_mean[:, 1], pos_var_std[:, 1], fmt='o-')
legend(('X', 'Y'))
xlabel('Network size E and I')
ylabel('Average bump variance (normalised)')
savefig(dirName + '/bump_position_variance_errbars.pdf')
hold(False)

figure()
subplot(211)
errorbar(np.array(pos_N)**2, pos_var_mean[:, 0], pos_var_std[:, 0], fmt='o-')
ylabel('Average bump variance X (normalised)')
subplot(212)
errorbar(np.array(pos_N)**2, pos_var_mean[:, 1], pos_var_std[:, 1], fmt='o-')
xlabel('Network size E and I')
ylabel('Average bump variance Y (normalised)')
savefig(dirName + '/bump_position_variance_errbars_sep.pdf')

figure()
hold(True)
plot(np.array(pos_N)**2, pos_var_mean[:, 0], 'o-')
plot(np.array(pos_N)**2, pos_var_mean[:, 1], 'o-')
legend(('X', 'Y'))
xlabel('Network size E and I')
ylabel('Average bump variance (normalised)')
savefig(dirName + '/bump_position_variance_lines.pdf')
hold(False)

