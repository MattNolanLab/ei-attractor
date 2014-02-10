#!/bin/sh
#
#   cluster_submit.sh
#
#   Submit job to the cluster.
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

#$ -P inf_ndtc
#$ -cwd
#$ -j y

# Initialise environment module
. /etc/profile.d/modules.sh

module load python/2.7.5


BASE=../../
export LOCAL_DIR=/exports/work/inf_ndtc/lsolanka
export PYTHONPATH=$PYTHONPATH:$LOCAL_DIR/GridCells/grid_cell_model
export PYTHONUSERBASE=$LOCAL_DIR/usr/local/

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCAL_DIR/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCAL_DIR/usr/local/lib/nest
export LD_LIBRARY_PATH

trap 'echo catch signal USR2 at `date +"%D %T"`' usr2

# Run the program
python $* 
