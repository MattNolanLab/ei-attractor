#!/usr/bin/env python
'''
Aggregate grid field data into the reductions file.
'''
import numpy as np
from grid_cell_model.parameters import JobTrialSpace2D
from grid_cell_model.submitting import flagparse

evenSpacingType = 'even-spacing'
detailedNoiseType = 'detailed-noise'
allowedTypes = [evenSpacingType, detailedNoiseType]

# Positions
allowedPositions = ['EI-1_3', 'EI-3_1']
detailedShape = (31, 9)

# Even spacing
evenShape     = (31, 31)

parser = flagparse.FlagParser()
parser.add_argument("type",      type=str, choices=allowedTypes,
        metavar='types', help='Type of the aggregation. Can be one of %s' %
        str(allowedTypes))
parser.add_argument("where",     type=str, help='Root directory')
parser.add_argument("--ns",      type=str, choices=["0pA", "150pA", "300pA"])
parser.add_argument('--ntrials', type=int, default=3)
parser.add_argument('--noLoadData', action='store_true')
parser.add_argument('--position',type=str, choices=allowedPositions)
parser.add_flag('--gridFields')
parser.add_flag('--FR')
args = parser.parse_args()

ns_all = ['0pA', '150pA', '300pA']
trialNumList = range(args.ntrials)
varListBase = ['analysis']
loadData = not args.noLoadData

################################################################################

# determine the iterator
if args.type == evenSpacingType:
    shape = evenShape
    subDirs = ns_all if args.ns is None  else [args.ns]
else:
    shape = detailedShape
    subDirs = allowedPositions if args.position is None else [args.position]

for subDir in subDirs:
    rootDir = '{0}/{1}'.format(args.where, subDir)
    print rootDir, shape

    sp = JobTrialSpace2D(shape, rootDir)

    if args.gridFields or args.all:
        sp.aggregateData(varListBase + ['rateMap_e'], trialNumList,
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='list')
        sp.aggregateData(varListBase + ['rateMap_e_X'], [trialNumList[0]],
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='list')
        sp.aggregateData(varListBase + ['rateMap_e_Y'], [trialNumList[0]],
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='list')
        sp.aggregateData(varListBase + ['corr'], trialNumList, funReduce=None,
                saveData=True, loadData=loadData, output_dtype='list')
        sp.aggregateData(varListBase + ['corr_X'], [trialNumList[0]],
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='list')
        sp.aggregateData(varListBase + ['corr_Y'], [trialNumList[0]],
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='list')
        sp.aggregateData(varListBase + ['gridnessScore'], trialNumList,
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='array')
        sp.aggregateData(['options', 'arenaSize'], [trialNumList[0]],
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='array')

    if args.FR or args.all:
        sp.aggregateData(varListBase + ['FR_e', 'avg'], trialNumList,
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='array')
        sp.aggregateData(varListBase + ['FR_e', 'all'], trialNumList,
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='list')

        # An ugly method to extract only 10 recorded I cells and average their
        # firing rate AND store as 'all' dataset
        def extractFR_i(x):
            N = 10 # 10 neurons recorded, ugly, I know
            return np.mean(x[0:N])
        sp.aggregateData(varListBase + ['FR_i', 'avg'], trialNumList,
                funReduce=None, saveData=True, loadData=loadData,
                output_dtype='array')
        sp.aggregateData(varListBase + ['FR_i', 'all'], trialNumList,
                funReduce=extractFR_i, saveData=True, loadData=False,
                output_dtype='array')

