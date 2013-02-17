#!/bin/sh
#
#   eddie_submit.sh
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

# Initialise environment module

. /etc/profile.d/modules.sh

# Use python 2.6

module load python/2.6.6
module load openmpi-gcc


BASE=../../
export LOCAL_DIR=/exports/work/inf_ndtc/lsolanka
export PYTHONPATH="$LOCAL_DIR/usr/local/lib/python2.6/site-packages:$LOCAL_DIR/usr/local/lib64/python2.6/site-packages:$BASE"

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$LOCAL_DIR/usr/local/lib:$LOCAL_DIR/usr/local/lib/nest"


# Run the program
mpirun -np $NSLOTS python2.6 $* 
