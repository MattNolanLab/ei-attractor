#!/bin/sh

########################################
# Submit fiete_path_integration.py jobs
# with different network parameters
# into the cluster
########################################

# Constant parameters
F_SHEET_SIZE="64"
F_TIME="10"

# Parameters different for each job
CONN_MULT="20"
#CONN_MULT="
#    0.5000    1.0000    1.5000    2.0000    2.5000    3.0000    3.5000    4.0000    4.5000    5.0000    5.5000    6.0000    6.5000
#    7.0000    7.5000    8.0000    8.5000    9.0000    9.5000   10.0000   10.5000   11.0000   11.5000   12.0000   12.5000   13.0000
#   13.5000   14.0000   14.5000   15.0000   15.5000   16.0000   16.5000   17.0000   17.5000   18.0000   18.5000   19.0000   19.5000
#   20.0000   20.5000   21.0000   21.5000   22.0000   22.5000   23.0000   23.5000   24.0000   24.5000   25.0000   25.5000   26.0000
#   26.5000   27.0000   27.5000   28.0000   28.5000   29.0000   29.5000   30.0000   30.5000   31.0000   31.5000   32.0000   32.5000
#   33.0000   33.5000   34.0000   34.5000   35.0000   35.5000   36.0000   36.5000   37.0000   37.5000   38.0000   38.5000   39.0000
#   39.5000   40.0000   40.5000   41.0000   41.5000   42.0000   42.5000   43.0000   43.5000   44.0000   44.5000   45.0000   45.5000
#   46.0000   46.5000   47.0000   47.5000   48.0000   48.5000   49.0000   49.5000   50.0000"

#F_ALPHA="0.025 0.027 0.03 0.035 0.04 0.045 0.050 0.055 0.060 0.07 0.08 0.09 0.100"
F_ALPHA="0.015"

#F_INPUT="0.3000    0.4000    0.5000    0.6000    0.7000    0.8000"
F_INPUT="0.3"

F_TAUM="10"

F_TAUI="9"

REPEAT=1

F_LAMBDA_NET="14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30"
#F_LAMBDA_NET="20"

job_id=10000
for alpha in $F_ALPHA; do
    for conn_mult in $CONN_MULT; do
        for input in $F_INPUT; do
            for taum in $F_TAUM; do
                for taui in $F_TAUI; do
                    for lambda_net in $F_LAMBDA_NET; do

#####################
repeat=1
while [ $repeat -le $REPEAT ]; do
    echo "conn_mult=$conn_mult"
    echo "alpha=$alpha"
    echo "job_id=$job_id"
    echo "lambda_net=$lambda_net"

    qsub ./eddieBatch_fiete_path_integration.sh $F_SHEET_SIZE $F_TIME $alpha $conn_mult $job_id $input $taum $taui $lambda_net

    echo

    let job_id=$job_id+1
    let repeat=$repeat+1
done
#####################

                    done
                done
            done
        done
    done
done