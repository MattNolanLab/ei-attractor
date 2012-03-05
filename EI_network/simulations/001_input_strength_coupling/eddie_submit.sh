#!/bin/sh
########################################
#                                      #
# SGE job script for ECDF Cluster      #
#                                      #
# by ECDF System Team                  #
# ecdf-systems-team@lists.ed.ac.uk     #
#                                      #
########################################

# Initialise environment module

. /etc/profile.d/modules.sh

# Use python 2.6

module load python/2.6.3


BASE=../../src/mult_bump_spiking_net/
export PYTHONPATH="/exports/work/informatics/s0966762/python-modules/lib/python2.6/site-packages:$BASE"


# Run the program
python2.6 simulation.py $* 
