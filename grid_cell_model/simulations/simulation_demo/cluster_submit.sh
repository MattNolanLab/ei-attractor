#!/bin/sh
#
#   cluster_submit.sh
#
#   Submit job to the cluster.
#
#   This file is called by qsub when running simulations on an SGE cluster.
#
#   While this file is part of the repository for this directory, it is not
#   generic and needs to be adjusted to reflect the settings for the particular
#   SGE environment you are using. It is not going to work for you if you leave
#   the values as they are committed here.
#
#   Please consult the documentation for your SGE environment so that you know
#   what each particular setting means.
#

#$ -P inf_ndtc
#$ -cwd
#$ -j y

# Initialise environment module
. /etc/profile.d/modules.sh

module load python/2.7.5


BASE=../../
export LOCAL_DIR=/exports/work/inf_ndtc/lsolanka

# This sets paths to shared dynamic libraries, i.e. especially the grid cell
# NEST module.
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCAL_DIR/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCAL_DIR/usr/local/lib/nest
export LD_LIBRARY_PATH

# virtualenvwrapper
# This is not strictly necessary but is a convenience wrapper for python
# virtual environments. See https://virtualenvwrapper.readthedocs.org for more
# information.
#
# Provides the ``workon`` command, see below.
export WORKON_HOME=$LOCAL_DIR/Envs
source $LOCAL_DIR/usr/local/bin/virtualenvwrapper.sh

# Do not remove this line!
trap 'echo catch signal USR2 at `date +"%D %T"`' usr2

# Run the program.
# This assumes you have a virtual environment called 'noise' and the project
# has been installed into this virtual environment.
#
# The next step essentially takes the arguments supplied to the script on the
# command line (by qsub) and runs the python simulation.
workon noise
echo "Virtual environment: $VIRTUAL_ENV"
python $* 
