#!/bin/sh

########################################
# Submit fiete_path_integration.py jobs
# with different network parameters
# into the cluster
########################################

# Produce recordings of membrane potentials of all single neurons
# This requires changes in the python codes

F_SHEET_SIZE="96"

F_TIME="10"

CONN_MULT="20"

F_ALPHA="0.015"

F_INPUT="0.3"

F_TAUM="10"

F_TAUI="10"

REPEAT=4

F_LAMBDA_NET="20"

F_THRESHOLD="-20"

F_L="2"

job_id=23500
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

    nice python2.6 fiete_path_integration.py --record-sn --record-sn-row -w -n $job_id -s $sheet_size -t $F_TIME \
        --alpha=$alpha -c $conn_mult -i $input --taum $taum --taui $taui  \
        --lambda-net=$lambda_net --threshold $threshold -l $l_param&

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
