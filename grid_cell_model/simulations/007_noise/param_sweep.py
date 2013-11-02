#
#   param_sweep.py
#
#   Parameter sweep (2D): shared procedures.
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
from submitting.factory   import SubmitterFactory
from submitting.arguments import ArgumentCreator
from data_storage         import DataStorage
from otherpkg.log         import log_info


def submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel,
        appName, rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
        extraIterparams={}):
    ac = ArgumentCreator(p, printout=True)

    GArr = np.linspace(startG, endG, Nvals)
    #GArr = [1.0]
    print(GArr)

    g_AMPA_total_arr     = []
    g_GABA_total_arr     = []
    for E_coupling in GArr:
        for I_coupling in GArr:
            g_AMPA_total_arr.append(E_coupling)
            g_GABA_total_arr.append(I_coupling)


    iterparams = {
        'g_AMPA_total'      : np.array(g_AMPA_total_arr),
        'g_GABA_total'      : np.array(g_GABA_total_arr),
        #'g_AMPA_total'      : [1400],
        #'g_GABA_total'      : [2160]
    }
    iterparams.update(extraIterparams)
    ac.insertDict(iterparams, mult=False)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(ac, appName, envType=ENV,
            rtLimit=rtLimit, output_dir=simRootDir, label=simLabel,
            blocking=blocking, timePrefix=timePrefix, numCPU=numCPU)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run)
    submitter.saveIterParams(iterparams, dry_run=dry_run)


###############################################################################

def getBumpCurrentSlope(noise_sigma, threshold=0):
    fileName = 'bump_slope_data/bump_slope_{0}pA.h5'.format(int(noise_sigma))
    ds = DataStorage.open(fileName, 'r')
    slopes = np.abs(ds['lineFitSlope'].flatten())
    ds.close()
    slopes[slopes < threshold] = np.nan
    return slopes

def getSpeedPercentile(p, path, grid_lambda, Nx):
    '''
    Retrieve the file containing animal positions and calculate the bump
    speed value at the specified percentile.

    Parameters
    ----------
    p : float
        The specified percentile.
    path : string
        Path to the file containing rat velocities
    grid_lambda : float
        Grid field spacing (cm)
    Nx : int
        Neural sheet size (neurons). THe bump has to travel this distance (in
        units of neurons) in order to return back to its original position,
        i.e. form a grid field.
    output : float
        The bump speed at the p-th percentile
    '''
    from scipy.io import loadmat
    data  = loadmat(path)
    dt    = float(data['dt'])
    pos_x = data['pos_x'].flatten()
    pos_y = data['pos_y'].flatten()
    vel_x = np.diff(pos_x)/dt
    vel_y = np.diff(pos_y)/dt
    animal_s = np.abs(np.hstack((vel_x, vel_y)))
    bump_s = float(Nx) / grid_lambda * animal_s
    res   = np.percentile(bump_s, p)

    msg = "Loaded velocity data from: {0}".format(path)
    log_info("getAnimalSpeedPercentile", msg)
    msg = "{0:.2f}th percentile: {1:.3f}".format(p, res)
    log_info("getAnimalSpeedPercentile", msg)

    return res


