#
#   other.py
#
#   Other routines that do not fit into generic structure of the package.
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

__all__ = ['bumpVariancePos']


def bumpVariancePos(dirName, fileNamePrefix, jobRange, trialRange, t_start):
    '''t_start - in discrete time'''
    ret = Position2D()

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



def createIgorSpikeRaster(spikes, yvals=None):
    '''
    spikes  row-wise 2d array of spike times
    '''
    if yvals is None:
        yvals = range(len(spikes)) + 1

    raster_x = np.ndarray((0))
    raster_y = np.ndarray((0))
    for it in range(len(spikes)):
        raster_x = np.hstack((raster_x, spikes[it]))
        raster_y = np.hstack((raster_y, np.zeros((len(spikes[it]))) + yvals[it]))

    return (raster_x, raster_y)


def createJobDir(options):
    # Create job directory in options.output_dir/jobXXXX
    outputDir = options.output_dir + '/job{0:04}'.format(options.job_num) +'/'
    ec = system('mkdir ' + outputDir)
    if ec != 0:
        print "Could not create output directory: " + outputDir
        ec = system('ls ' + outputDir)
        if ec == 0:
            print "But it can be listed --> continuing!"
        else:
            print "And it cannot be listed. Check your permissions and rerun!"
            exit(1)
    return outputDir





