from scipy.io import loadmat
from matplotlib.pyplot import *

import numpy as np


class ReturnStruct(object):
    pass

def circularVariance(x, range):
    ''' Returns circular moments of a vector of circular variable x, defined in
    the 'range' (this will be mapped to a circle)'''
    c = np.exp(1j*2*np.pi*x/range)
    avg = np.mean(c)
    theta_avg = np.angle(avg)
    theta_var = 1 - np.abs(avg)
    return (theta_avg, theta_var)
    

def bumpVariancePos(dirName, fileNamePrefix, jobRange, trialRange, t_start):
    '''t_start - in discrete time'''
    ret = ReturnStruct()

    fileNameTemp = "{0}/{1}job{2:04}_trial{3:04}"
    F_dt = 0.2
    it_start = t_start/F_dt
    
    jobN = jobRange[1] - jobRange[0]+1
    trialN = trialRange[1] - trialRange[0] + 1
    pos_var = []
    pos_N = []
    pos_var_mean = np.ndarray((jobN, 2))
    pos_var_std = np.ndarray((jobN, 2))
    pos_x_vec = np.ndarray((jobN, trialN), dtype=object)
    pos_y_vec = np.ndarray((jobN, trialN), dtype=object)
    o_vec = np.ndarray((jobN, trialN), dtype=object)
    
    for job_it in range(jobN):
        jobNum = job_it + jobRange[0]
        pos_var.append([])
        for trial_it in range(trialN):
            trialNum = trialRange[0] + trial_it
            try:
                fileName = fileNameTemp.format(dirName, fileNamePrefix, jobNum,
                    trialNum)
                print "Processing file: " + fileName + "_output.mat"
                data = loadmat(fileName + "_output.mat")
            except:
                print "Warning: could not open: " + fileName
                continue
    
            o = data['options']
            o_vec[job_it, trial_it] = o
            pos_x = data['bumpPos'][:, 0]
            pos_y = data['bumpPos'][:, 1]
            pos_x_vec[job_it, trial_it] = pos_x
            pos_y_vec[job_it, trial_it] = pos_y
            pos_x = pos_x[it_start:]
            pos_y = pos_y[it_start:]
            
            Ne = o['Ne'][0][0][0][0] # Ugly but these bastards cannot export to
                                     # matlab
            mean_x, var_x = circularVariance(pos_x, Ne)
            mean_y, var_y = circularVariance(pos_y, Ne)
            pos_var[job_it].append([var_x, var_y])
    
        pos_var_mean[job_it, :] = np.mean(pos_var[job_it], 0)
        pos_var_std[job_it, :] = np.std(pos_var[job_it], 0)
        pos_N.append(Ne)
    
    ret.pos_x_vec = pos_x_vec
    ret.pos_y_vec = pos_y_vec
    ret.o_vec = o_vec
    ret.pos_var = pos_var
    ret.pos_var_mean = pos_var_mean
    ret.pos_var_std = pos_var_std
    ret.pos_N = pos_N
    return ret
