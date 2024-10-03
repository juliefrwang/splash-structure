#!/bin/sh
#SBATCH -p horence,normal
#SBATCH --time=5:00:00
#SBATCH --mem=60000      # In MB
#SBATCH --job-name=test     # job name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Variables to be set
RESULTS_DIR=/oak/stanford/groups/horence/juliew/github/splash-structure/test_runs # This folder has to be made ahead and will contain results of all runs.
DATA_HANDLE=juliew  # Specifically, for the current run, result will save to folder named RESULTS_DIR/${DATA_HANDLE}_results/
SPLASH_OUT_FILE=/oak/stanford/groups/horence/juliew/github/splash-structure/test_runs/test.after_correction.scores.tsv # SPLASH significant anchor output file
COMPACTOR_FILE=/oak/stanford/groups/horence/juliew/github/splash-structure/test_runs/test_compactor.tsv # compactor output file

# Load the environment
source /oak/stanford/groups/horence/juliew/envs/structure_run/bin/activate

# Go to the results directory
cd $RESULTS_DIR

# run structure on target
time python3 /oak/stanford/groups/horence/juliew/github/splash-structure/scripts/main_structure_target.py $SPLASH_OUT_FILE $DATA_HANDLE &
pid1=$!
# run structure on compactor 40 mers
time python3 /oak/stanford/groups/horence/juliew/github/splash-structure/scripts/main_structure_compactor_40mers.py $COMPACTOR_FILE $DATA_HANDLE &
pid2=$!


# wait for all programs to finish 
wait $pid1
status1=$?
wait $pid2
status2=$?

# check the status of all python jobs
if [ $status1 -ne 0 ]; then
    echo "Error: Structure on target for $DATA_HANDLE failed with status $status1"
    exit 1
fi

if [ $status2 -ne 0 ]; then
    echo "Error: Structure on compactor_40 for $DATA_HANDLE failed with status $status2"
    exit 1
fi
