#!/bin/sh
########################################
#                                      #
# SGE MPI job script for ECDF Cluster  #
#                                      #
# by ECDF System Team                  #
# ecdf-systems-team@lists.ed.ac.uk     #
#                                      #
########################################

# Grid Engine options

#$ -N plotStatistics
#$ -cwd
#$ -l h_rt=00:02:00

# Initialise environment module

. /etc/profile.d/modules.sh

# Use python 2.6

module load matlab/4.0-r2008b

F_JOB_ID=$1


# Run the program
matlab -nodisplay -r "d = dir('results/*job${F_JOB_ID}_*.mat'); plotStatistics(['results/' d(end).name], [2048]); exit"
