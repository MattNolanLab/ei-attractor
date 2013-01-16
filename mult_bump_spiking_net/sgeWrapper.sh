#!/bin/sh
###################################
#                                 #
# Matlab SGE job wrapper          #
# Copyright 2007 ECDF System Team #
# University of Edinburgh         #
#                                 #
###################################

#$ -S /bin/sh

ulimit -v 2048000
echo `ulimit -v` virtual memory can be used

echo "Executing: ${MDCE_MATLAB_COMMAND}"
echo "${MDCE_MATLAB_EXE}" ${MDCE_MATLAB_ARGS}
exec "${MDCE_MATLAB_EXE}" ${MDCE_MATLAB_ARGS}
exit 0 
