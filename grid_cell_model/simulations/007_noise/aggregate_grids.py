#!/usr/bin/env python
'''
Aggregate grid field data into the reductions file.
'''
import numpy as np
from parameters  import JobTrialSpace2D
from submitting import flagparse

parser = flagparse.FlagParser()
parser.add_argument("--where",   type=str, required=True)
parser.add_argument("--ns",      type=int, choices=[0, 150, 300])
parser.add_argument('--ntrials', type=int, default=3)
parser.add_argument('--noLoadData', action='store_true')

parser.add_flag('--gridFields')
parser.add_flag('--FR')
args = parser.parse_args()

ns_all = [0, 150, 300]
shape = (31, 31)
noise_sigmas = ns_all if args.ns is None  else [args.ns]
trialNumList = range(args.ntrials)
varListBase = ['analysis']
loadData = not args.noLoadData

################################################################################
for noise_sigma in noise_sigmas:
    rootDir = '{0}/{1}pA'.format(args.where, noise_sigma)

    if args.gridFields or args.all:
        sp = JobTrialSpace2D(shape, rootDir)
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
        sp.aggregateData(varListBase + ['FR_i', 'all'], trialNumList,
                funReduce=extractFR_i, saveData=True, loadData=False,
                output_dtype='array')

