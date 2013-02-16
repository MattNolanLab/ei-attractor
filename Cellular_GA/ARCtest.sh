#! /bin/bash

# Author: Lukas Solanka <xsolan00@stud.fit.vutbr.cz>
# Counts Effectivity of the Cellular GA (MPI version must be compiled)

APP_NAME=${APP_NAME-"testOneMax"}

#NUM_PROCS="1 4 9 16 25 36"
NUM_PROCS=${NUM_PROCS-"1 4 9"}  # The first value must be 1

# Do some cleanup&compile
#make clean
USE_MPI=1 make $APP_NAME || exit 1

logfile=mpitest.log

#
# Select GA parameters here (testGA -h)
# To change the chromosome lengthyou must edit testGA.cc
#
xsize=${xsize-"10"}	# Number of population matrix columns
ysize=${ysize-"10"}	# Number of population matrix rows
gen=${gen-"100"}	# Max generations
rlen=${rlen-"1"}	# Random walk length

echo -n > $logfile

for procnum in $NUM_PROCS; 
do
    echo "-----------------------------------------------------------------------"
    echo "Processors: $procnum"
    ./mpi $procnum $APP_NAME -g $gen -xsize $xsize -ysize $ysize -r $rlen 2>> $logfile
    sleep 1
done

cat $logfile | grep "<gastat>" | cut -d" " -f 3,5,7 |
{
    echo "#Processors/cores	Time [s]	Speedup		Effectivity"

    while read procnum gennum time
	do
		if [ "$procnum" -eq 1 ]
		then
			time_one=$time
		fi
		
		speedup=`bc <<< "scale=5; $time_one/$time"`
		eff=`bc <<< "scale=5; $speedup/$procnum"`
		# Print statistics
		echo "$procnum			$time		$speedup		$eff"
	done
}
