#
#   aggregate.py
#
#   Data aggregation, mainly parameter sweeps
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
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
import numpy.ma as ma

from parameters import DataSpace
from otherpkg.log import log_warn


def aggregate2DTrial(sp, varList, trialNumList, fReduce=np.mean,
        ignoreNaNs=False):
    '''
    Aggregate all the data from a 2D ParamSpace, applying fReduce on the trials
    of all the data sets in the parameter space.

    Parameters
    ----------
    sp : ParamSpace
        A parameter space to apply the reduction on.
    varList : list of strings
        A variable list, specifying the location of the variable to reduce.
        This will be prefixed with ['analysis']
    trialNumList : list of ints
        A list specifying exactly which trials are to be processed.
    fReduce : a function f(data, axis)
        A reduction function.
    ignoreNaNs : bool, optional
        If True, mask the NaN values.
    output : a 2D numpy array of the reduced values
    '''
    varList = ['analysis'] + varList
    retVar = sp.aggregateData(varList, trialNumList, funReduce=np.mean,
            saveData=False)
    if (ignoreNaNs):
        nans = np.isnan(retVar)
        retVar = ma.MaskedArray(retVar, mask=nans)
    return fReduce(retVar, 2)


def aggregate2D(sp, varList, funReduce=None):
    '''
    Aggregate all the data from a 2D ParamSpace, applying fReduce on the trials
    of all the data sets in the parameter space, however the data is retrieved
    from the top-level of the data hierarchy, i.e. sp['analysis']. funReduce is
    applied on the data of each item in the ParamSpace (not necessarily
    trials).

    Parameters
    ----------
    sp : ParamSpace
        A parameter space to apply the reduction on.
    varList : list of strings
        A variable list, specifying the location of the variable to reduce.
        This will be prefixed with ['analysis']
    funReduce : a function f(data, axis)
        A reduction function.
    output : a 2D numpy array of the reduced values
    '''
    varList = ['analysis'] + varList
    return sp.aggregateData(varList, funReduce=funReduce,
            trialNumList='all-at-once', saveData=True)



def computeYX(sp, iterList, r=0, c=0, trialNum=0, normalize=True, **kw):
    E, I = sp.getIteratedParameters(iterList)
    if (normalize):
        Ne = DataSpace.getNetParam(sp[r][c][trialNum].data, 'net_Ne')
        Ni = DataSpace.getNetParam(sp[r][c][trialNum].data, 'net_Ni')
    else:
        Ne = 1
        Ni = 1
    return E/Ne, I/Ni



def computeVelYX(sp, iterList, r=0, c=0, trialNum=0, normalize=True, **kw):
    E, I = sp.getIteratedParameters(iterList)
    if (normalize):
        Ne = DataSpace.getNetParam(sp[r][c][trialNum].data['IvelData'][0],
                'net_Ne')
        Ni = DataSpace.getNetParam(sp[r][c][trialNum].data['IvelData'][0],
                'net_Ni')
    else:
        Ne = 1.0
        Ni = 1.0

    return E/Ne, I/Ni



def aggregateType(sp, iterList, types, NTrials, ignoreNaNs=False, **kw):
    '''
    Automatically aggregate data according to the type of the data.
    '''
    type, subType = types
    vars          = ['analysis']
    output_dtype  = 'array'
    funReduce     = None

    if (type == 'gamma'):
        # Gamma oscillation analyses
        if (subType == 'acVal'):
            # Autocorrelation first local maximum
            vars += ['acVal']
        elif (subType == 'freq'):
            # Gamm frequency
            vars += ['freq']
        elif (subType == 'acVec'):
            # All the autocorrelations
            vars += ['acVec']
            output_dtype = 'list'
        else:
            raise ValueError('Unknown gamma subtype: {0}'.format(subType))
        trialNumList  = np.arange(NTrials)

    elif (type == 'bump'):
        trialNumList  = np.arange(NTrials)
        if (subType == 'sigma'):
            vars += ['bump_e', 'sigma']
        else:
            raise ValueError('Unknown bump subtype: {0}'.format(subType))

    elif (type == 'velocity'):
        if (subType == 'slope'):
            vars += ['lineFitSlope']
        elif (subType == 'fitErr'):
            vars += ['lineFitErr']
            funReduce = np.sum
        else:
            raise ValueError('Unknown velocity subtype: {0}'.format(subType))
        trialNumList = 'all-at-once'

    elif (type == 'grids'):
        trialNumList  = np.arange(NTrials)
        if (subType == 'gridnessScore'):
            vars += ['gridnessScore']
            funReduce = None
        else:
            raise ValueError('Unknown grids subtype: {0}'.format(subType))
    elif type == 'FR':
        trialNumList = np.arange(NTrials)
        if subType == 'E':
            vars += ['FR_e', 'avg']
            funReduce = None
        elif subType == 'I_10': # user should be aware of 10 neurons limitation
            vars += ['FR_i', 'all']
            funReduce = None
        else:
            raise ValueError('Unknown FR subtype: {0}'.format(subType))

    else:
        raise ValueError('Unknown aggregation type: {0}'.format(type))


    data = sp.aggregateData(vars, trialNumList, output_dtype=output_dtype,
            loadData=True, saveData=False, funReduce=funReduce)
    if (ignoreNaNs):
        log_warn('aggregateType', 'Ignoring NaNs')
        nans = np.isnan(data)
        data = ma.MaskedArray(data, mask=nans)

    if (type == 'velocity'):
        Y, X = computeVelYX(sp, iterList, normalize=False, **kw)
    else:
        Y, X = computeYX(sp, iterList, normalize=False, **kw)
        data = np.mean(data, axis=2) # TODO: fix the trials, stack them
        # bump sigma is a reciprocal
        if (type == 'bump'):
            if (subType == 'sigma'):
                data = 1./data

    return data, X, Y



