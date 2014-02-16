#!/usr/bin/env python
'''
Aggregate bump data into the reductions file.
'''
import numpy as np
from parameters  import JobTrialSpace2D

from submitting import flagparse

parser = flagparse.FlagParser()
parser.add_flag('--bump')
parser.add_flag('--positions')
parser.add_flag('--AC')
parser.add_argument('--ntrials', type=int, default=5)
parser.add_argument('--noLoadData', action='store_true')
args = parser.parse_args()


ns_all  = [0, 150, 300]
ns_none = [-100]
dirs = \
    ("output/even_spacing/const_position/{0}pA",     (31, 31), ns_all)
    #("output/even_spacing/gamma_bump/{0}pA",     (31, 31), ns_all)
    #("output/detailed_noise/gamma_bump/EI-3_1",  (31, 9),  ns_none)
    #("output/detailed_noise/gamma_bump/EI-1_3",  (31, 9),  ns_none)

trialNumList = xrange(args.ntrials)
varListBase = ['analysis']
loadData = not args.noLoadData

################################################################################
shape        = dirs[1]
noise_sigmas = dirs[2]
for noise_sigma in noise_sigmas:
    rootDir = dirs[0].format(int(noise_sigma))
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
        sp.aggregateData(bumpPosVars + ['errSum'],
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
        sp.aggregateData(bumpPosVars + ['t'],
                trialNumList, funReduce=None, loadData=loadData,saveData=True,
                output_dtype='list')


    if args.AC or args.all:
        sp.aggregateData(varListBase + ['acVal'], trialNumList, funReduce=np.mean,
                loadData=loadData, saveData=True, output_dtype='array')
        sp.aggregateData(varListBase + ['acVec'], trialNumList, funReduce=None,
                loadData=loadData, saveData=True, output_dtype='list')
        sp.aggregateData(varListBase + ['freq'], trialNumList, funReduce=np.mean,
                loadData=loadData, saveData=True, output_dtype='array')
