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

#$ -N fiete_path_integration
#$ -cwd
#$ -l h_rt=04:00:00

# Initialise environment module

. /etc/profile.d/modules.sh

# Use python 2.6

module load python/2.6.3
module load matlab/4.0-r2008b

export PYTHONPATH=/exports/work/informatics/s0966762/python-modules/lib/python2.6/site-packages

F_SHEET_SIZE="$1"
F_TIME="$2"
F_ALPHA="$3"
F_CONN_MULT="$4"
F_JOB_ID=$5


# Run the program
python2.6 fiete_path_integration.py -w -n $F_JOB_ID -s $F_SHEET_SIZE -t $F_TIME --alpha=$F_ALPHA -c $F_CONN_MULT
#matlab -nodisplay -r "d = dir('results/*job${F_JOB_ID}_*.mat'); plotStatistics(['results/' d(end).name], [2048]); exit"
