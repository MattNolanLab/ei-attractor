#
#   submit_job.py
#
#   Submit job(s) to the cluster/workstation
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

from default_params import defaultParameters
from common         import *


EDDIE = False  # if eddie, submit on a cluster using qsub

QSUB_PARAMS = "-P inf_ndtc -cwd -l h_rt=06:00:00 -pe memory-2G 2"

net_generations=5

parameters = defaultParameters

ac = ArgumentCreator(parameters)

iterparams = {
        'Iext_e'    : [1, 2, 3, 4],
        'Iext_i'    : [5, 6, 7, 8]}

ac.insertDict(iterparams, mult=False)
print ac
