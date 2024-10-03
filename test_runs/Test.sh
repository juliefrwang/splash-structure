#!/bin/sh

# Get the current directory (where the user runs the script)
BASE_DIR=$(pwd)

# Define the results directory relative to the user's repo path
DATA_HANDLE="default_test_run"  # Specifically, for the current run, results will save to folder named RESULTS_DIR/${DATA_HANDLE}_results/
RESULTS_DIR="${BASE_DIR}/${DATA_HANDLE}_results" 

# Define paths for files relative to the user's repo path
SPLASH_OUT_FILE="${BASE_DIR}/test_data/test.after_correction.scores.tsv"  # SPLASH significant anchor output file
COMPACTOR_FILE="${BASE_DIR}/test_data/test_compactor.tsv"  # Compactor output file

# Print the paths (for debugging purposes)
echo "Results Directory: ${RESULTS_DIR}"
echo "SPLASH Output File: ${SPLASH_OUT_FILE}"
echo "Compactor File: ${COMPACTOR_FILE}"

# Load the environment
source /oak/stanford/groups/horence/juliew/envs/structure_run/bin/activate

# Go to the results directory
cd $BASE_DIR

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
