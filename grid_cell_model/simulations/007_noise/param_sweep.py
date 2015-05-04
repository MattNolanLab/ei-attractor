'''Parameter sweep (2D): shared procedures.'''
from __future__ import absolute_import, print_function

import numpy as np
from grid_cell_model.submitting.factory   import SubmitterFactory
from grid_cell_model.submitting.arguments import ArgumentCreator
from grid_cell_model.submitting.noise.slopes import (DefaultSelector,
                                                     NoThetaSelector)
from grid_cell_model.data_storage import DataStorage
from grid_cell_model.otherpkg.log import log_info


def submitParamSweep(p, startG, endG, Nvals, ENV, simRootDir, simLabel, appName,
                     rtLimit, numCPU, blocking, timePrefix, numRepeat, dry_run,
                     extraIterparams=None, rc=None, **kwargs):
    '''Submit and save metadata for the gE vs gI parameter sweep.'''
    printout = kwargs.pop('printout', True)
    if extraIterparams is None:
        extraIterparams = {}
    ac = ArgumentCreator(p, printout=printout)

    GArr = np.linspace(startG, endG, Nvals)
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
    }
    dimension_labels = ['g_AMPA_total', 'g_GABA_total']
    dimensions = [Nvals, Nvals]
    iterparams.update(extraIterparams)
    ac.insertDict(iterparams, mult=False)

    ###############################################################################
    submitter = SubmitterFactory.getSubmitter(
        ac, appName, envType=ENV, rtLimit=rtLimit, output_dir=simRootDir,
        label=simLabel, blocking=blocking, timePrefix=timePrefix, numCPU=numCPU,
        **kwargs)
    ac.setOption('output_dir', submitter.outputDir())
    startJobNum = 0
    filter = rc[0]*len(GArr) + rc[1] if rc is not None else None
    submitter.submitAll(startJobNum, numRepeat, dry_run=dry_run, filter=filter)
    submitter.saveIterParams(iterparams, dimension_labels, dimensions,
                             dry_run=dry_run)


###############################################################################

def getBumpCurrentSlope(noise_sigma, threshold=0, type=None):
    '''
    Parameters
    ----------
    noise_sigma : int
        Noise level (sigma of the Gaussian)
    threshold : float
        Threshold below which slope values will be replaced with ``NaN``.
    type : string, optional
        If ``None`` the regular bump slope files will be used. If ``no_theta``,
        the bump slope files specific for the simulations wihtout theta
        oscillations will be used.
    '''
    data_root = 'bump_slope_data'
    selector_cls = None

    if type is None:
        selector_cls = DefaultSelector
    elif type == 'no_theta':
        selector_cls = NoThetaSelector
    else:
        raise ValueError('Invalid bump slope type.')

    selector = selector_cls(data_root, threshold)
    return selector.get_slopes(noise_sigma)

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
