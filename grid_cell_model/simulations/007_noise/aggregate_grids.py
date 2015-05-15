#!/usr/bin/env python
'''
Aggregate grid field data into the reductions file.
'''
import numpy as np
from grid_cell_model.parameters import JobTrialSpace2D

from aggregate import AggregationParser

parser = AggregationParser()
parser.add_flag('--gridFields')
parser.add_flag('--gridFields_i_cells')
parser.add_flag('--FR')
args = parser.parse_args()

trialNumList = range(args.ntrials)
varListBase = ['analysis']
varListBase_i = ['analysis', 'i_fields']

################################################################################
for subDir in parser.subdirs:
    rootDir = '{0}/{1}'.format(args.where, subDir)
    print rootDir, parser.shape

    sp = JobTrialSpace2D(parser.shape, rootDir)

    if args.gridFields or args.all:
        sp.aggregateData(varListBase + ['rateMap_e'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase + ['rateMap_e_X'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase + ['rateMap_e_Y'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase + ['corr'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase + ['corr_X'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase + ['corr_Y'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase + ['gridnessScore'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='array')
        sp.aggregateData(['options', 'arenaSize'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='array')

    if args.gridFields_i_cells or args.all:
        sp.aggregateData(varListBase_i + ['rateMap_i'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase_i + ['rateMap_i_X'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase_i + ['rateMap_i_Y'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase_i + ['corr_i'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase_i + ['corr_X'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase_i + ['corr_Y'],
                         [trialNumList[0]], funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')
        sp.aggregateData(varListBase_i + ['gridnessScore'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='array')

    if args.FR or args.all:
        sp.aggregateData(varListBase + ['FR_e', 'avg'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='array')
        sp.aggregateData(varListBase + ['FR_e', 'all'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='list')

        # An ugly method to extract only 10 recorded I cells and average their
        # firing rate AND store as 'all' dataset
        def extractFR_i(x):
            N = 10 # 10 neurons recorded, ugly, I know
            return np.mean(x[0:N])
        sp.aggregateData(varListBase + ['FR_i', 'avg'],
                         trialNumList, funReduce=None, saveData=True,
                         loadData=parser.load_data, output_dtype='array')
        sp.aggregateData(varListBase + ['FR_i', 'all'],
                         trialNumList, funReduce=extractFR_i, saveData=True,
                         loadData=False, output_dtype='array')
