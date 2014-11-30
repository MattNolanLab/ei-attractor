#!/usr/bin/env python
'''Aggregate velocity analysis outputs.'''
import numpy as np
from grid_cell_model.parameters import JobTrialSpace2D

from aggregate import AggregationParser

parser = AggregationParser()
parser.add_flag('--line_fits')
args = parser.parse_args()

varListBase  = ['analysis']

################################################################################
for subDir in parser.subdirs:
    rootDir = '{0}/{1}'.format(args.where, subDir)
    print rootDir, parser.shape

    sp = JobTrialSpace2D(parser.shape, rootDir)

    if args.line_fits or args.all:
        sp.aggregateData(varListBase + ['lineFitErr'],
                         funReduce=np.sum, trialNumList='all-at-once',
                         saveData=True, output_dtype='array',
                         loadData=parser.load_data)
        sp.aggregateData(varListBase + ['lineFitSlope'],
                         funReduce=None, trialNumList='all-at-once',
                         saveData=True, output_dtype='array',
                         loadData=parser.load_data)
        sp.aggregateData(varListBase + ['bumpVelAll'],
                         funReduce=None, trialNumList='all-at-once',
                         saveData=True, output_dtype='list',
                         loadData=parser.load_data)
