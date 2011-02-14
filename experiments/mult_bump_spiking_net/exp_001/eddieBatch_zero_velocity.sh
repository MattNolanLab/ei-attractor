#!/bin/sh

########################################
# Submit fiete_path_integration.py jobs
# with different network parameters
# into the cluster
########################################

# different lambda_net parameters

F_SHEET_SIZE="96"

F_TIME="20"

CONN_MULT="20"

F_ALPHA="0.013"

F_INPUT="0.3"

F_TAUM="10"

F_TAUI="10"

REPEAT=10

F_LAMBDA_NET="13"

F_THRESHOLD="-20"

F_L="2"

F_NOISE_SIGMA="0.01"

job_id=1000
for alpha in $F_ALPHA; do
    for conn_mult in $CONN_MULT; do
        for input in $F_INPUT; do
            for taum in $F_TAUM; do
                for taui in $F_TAUI; do
                    for lambda_net in $F_LAMBDA_NET; do
                        for threshold in $F_THRESHOLD; do
                            for l_param in $F_L; do
                                for sheet_size in $F_SHEET_SIZE; do

#####################
repeat=1
while [ $repeat -le $REPEAT ]; do
    echo "sheet_size=$sheet_size"
    echo "conn_mult=$conn_mult"
    echo "alpha=$alpha"
    echo "job_id=$job_id"
    echo "lambda_net=$lambda_net"
    echo "threshold=$threshold"

    qsub ./eddieBatch_fiete_path_integration.sh $sheet_size $F_TIME $alpha $conn_mult $job_id $input $taum $taui $lambda_net $threshold $l_param $F_NOISE_SIGMA

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
        done
    done
done
