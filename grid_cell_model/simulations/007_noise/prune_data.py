#!/usr/bin/env python
from parameters  import JobTrialSpace2D
from submitting import flagparse
from visitors import data_manipulation as dm

###############################################################################
parser = flagparse.FlagParser()
parser.add_argument('-r', '--row',  type=int)
parser.add_argument('-c', '--col',  type=int)
parser.add_argument("--where", type=str, required=True)
parser.add_argument("--ns", type=int, choices=[0, 150, 300])
parser.add_argument("--prune_what", type=str, required=True)
parser.add_argument('--type', type=str, choices=['velocity'], required=True)
parser.add_argument('--output_dir', type=str, help='Unused')
parser.add_flag("--ns_all")
parser.add_flag("--repack", help='Whether to repack data after pruning')
o = parser.parse_args()

if not o.ns_all and o.ns is None:
    raise RuntimeError("Must specify either --ns or --ns_all!")
if (o.row is None) ^ (o.col is None):
    raise RuntimeError("Must specify either both --row and --col or none!")



ns_all  = [0, 150, 300]
shape   = (31, 31)

if o.row is not None and o.col is not None:
    dataPoints = [(o.row, o.col)] 
else:
    dataPoints = None

################################################################################
noise_sigmas = ns_all if o.ns_all  else [o.ns]
for noise_sigma in noise_sigmas:
    rootDir = '{0}/{1}pA'.format(o.where, int(noise_sigma))

    sp = JobTrialSpace2D(shape, rootDir, dataPoints=dataPoints)

    visitor = None
    if o.type == 'velocity':
        visitor = dm.VelocityPruningVisitor(what=o.prune_what)
        sp.visit(visitor, trialList='all-at-once')
        if o.repack:
            if o.row is not None:
                sp.repackItem(o.row, o.col)
            else:
                sp.repackAllItems()
    else:
        raise RuntimeError("Whoops! We shouldn't be here!")

