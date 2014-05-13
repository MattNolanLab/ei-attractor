#!/usr/bin/env python
'''
Aggregate bump data into the reductions file.
'''
import numpy as np
from grid_cell_model.parameters  import JobTrialSpace2D
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
        metavar='type', help='Type of the aggregation. Can be one of %s' %
        str(allowedTypes))
parser.add_argument("where",     type=str, help='Root directory')

parser.add_argument("--ns",      type=str, choices=['0pA', '150pA', '300pA'])
parser.add_argument('--ntrials', type=int, default=5)
parser.add_argument('--noLoadData', action='store_true')
parser.add_argument('--position',type=str, choices=allowedPositions)

parser.add_flag('--bump')
parser.add_flag('--positions')
parser.add_flag('--isBump')
parser.add_flag('--uniformML')
parser.add_flag('--AC')
parser.add_flag('--FR_e')
parser.add_flag('--FR_i')
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

    if args.bump or args.all:
        for suffix in ['', '_full']:
            sp.aggregateData(varListBase + ['bump_e'+suffix, 'sigma'],
                    trialNumList, funReduce=None, loadData=loadData,saveData=True,
                    output_dtype='array')
            sp.aggregateData(varListBase + ['bump_e'+suffix, 'err2'], trialNumList,
                    funReduce=None, loadData=loadData, saveData=True,
                    output_dtype='array')
            sp.aggregateData(varListBase + ['bump_e'+suffix, 'bump_e_rateMap'],
                    trialNumList, funReduce=None, loadData=loadData, saveData=True,
                    output_dtype='list')
            sp.aggregateData(varListBase + ['bump_i'+suffix, 'bump_i_rateMap'],
                    trialNumList, funReduce=None, loadData=loadData, saveData=True,
                    output_dtype='list')

    if args.positions or args.all:
        bumpPosVars = varListBase + ['bump_e', 'positions']
        sp.aggregateData(bumpPosVars + ['A'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['mu_x'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['mu_y'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['sigma'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['err2Sum'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['ln_L'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['lh_precision'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['times'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')

    if args.isBump or args.all:
        bumpPosVars = varListBase + ['bump_e', 'isBump']
        sp.aggregateData(bumpPosVars + ['isBumpFrames'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['fracTotal'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='array')

    if args.uniformML or args.all:
        bumpPosVars = varListBase + ['bump_e', 'uniformML']
        sp.aggregateData(bumpPosVars + ['mu'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['sigma2'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['ln_L'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['err2'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')
        sp.aggregateData(bumpPosVars + ['times'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')

    if args.AC or args.all:
        sp.aggregateData(varListBase + ['acVal'], trialNumList, funReduce=np.mean,
                loadData=loadData, saveData=True, output_dtype='array')
        sp.aggregateData(varListBase + ['acVec'], trialNumList, funReduce=None,
                loadData=loadData, saveData=True, output_dtype='list')
        sp.aggregateData(varListBase + ['freq'], trialNumList, funReduce=np.mean,
                loadData=loadData, saveData=True, output_dtype='array')

    if args.FR_e or args.all:
        FRVars = varListBase + ['FR_e']
        sp.aggregateData(FRVars + ['popSliding'], trialNumList, funReduce=None,
                         loadData=loadData, saveData=True, output_dtype='list')

    if args.FR_i or args.all:
        FRVars = varListBase + ['FR_i']
        sp.aggregateData(FRVars + ['popSliding'], trialNumList, funReduce=None,
                         loadData=loadData, saveData=True, output_dtype='list')
