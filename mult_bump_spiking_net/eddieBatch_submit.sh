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

export PYTHONPATH=/exports/work/informatics/s0966762/python-modules/lib/python2.6/site-packages

F_SHEET_SIZE="$1"
F_TIME="$2"
F_ALPHA="$3"
F_CONN_MULT="$4"
F_JOB_ID=$5
F_INPUT=$6
F_TAUM=$7
F_TAUI=$8
F_LAMBDA_NET=$9
F_THRESHOLD=${10}
F_L=${11}
F_NOISE_SIGMA=${12}


# Run the program
python2.6 fiete_path_integration.py -w -n $F_JOB_ID -s $F_SHEET_SIZE -t $F_TIME \
    --alpha=$F_ALPHA -c $F_CONN_MULT -i $F_INPUT --taum $F_TAUM --taui $F_TAUI  \
    --lambda-net=$F_LAMBDA_NET --threshold $F_THRESHOLD -l $F_L \
    --noise-sigma $F_NOISE_SIGMA
