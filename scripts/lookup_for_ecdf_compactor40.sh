#!/bin/sh
#SBATCH -p horence,normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --job-name=lookup_query     # job name
#SBATCH --output=./lookup/lookup-query-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=frwang

DATA_HANDLE=$1
LOOKUP_IDX="/oak/stanford/groups/horence/juliew/structure/lookup_table/lookup_list2/lookup"

mkdir lookup
cd lookup
STRUCTURE_OUT_ALL="/oak/stanford/groups/horence/juliew/structure/0906_results_rerun_com40/${DATA_HANDLE}_results/structure_on_compactors_40mers.tsv"
tail -n +2 $STRUCTURE_OUT_ALL | awk '{print ">Seq" NR; print $2}' > all_compactors.fasta

INPUT_FASTA="/oak/stanford/groups/horence/juliew/structure/0906_results_rerun_com40/${DATA_HANDLE}_results/lookup/all_compactors.fasta"
OUTPUT_FILE="/oak/stanford/groups/horence/juliew/structure/0906_results_rerun_com40/${DATA_HANDLE}_results/lookup/out_com40_lookup.txt"
/oak/stanford/groups/horence/juliew/workflows/splash-2.4.2/lookup_table query --stats_fmt with_stats $INPUT_FASTA $LOOKUP_IDX $OUTPUT_FILE &
pid3=$!

wait $pid3
status3=$?

if [ $status3 -ne 0 ]; then
    echo "Error: Lookup table on compactor_40 for $DATA_HANDLE failed with status $status3"
fi
