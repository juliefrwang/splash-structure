#!/bin/sh
#SBATCH -p horence,normal,owners
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=80000
#SBATCH --job-name=rfam     # job name
#SBATCH --output=./rfam/run_rfam-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

# This script is used to annotate the structure result with Rfam database using software infernal.
# Input Argumnet:
# $1: path to file of "structure_on_target.tsv" or "structure_on_compactors_40mers.tsv"
# $2: "extendor" or "compactor_40" or "compactor_20", denoting whether the input structure result is from extendors or compactors

# STEP 1: create fasta file from structure result

# python3 /oak/stanford/groups/horence/juliew/structure/src/prep_query_for_rfam.py $1 $2

# echo "juliew: Fasta file created."

# STEP 2: run infernal
Rfam=/oak/stanford/groups/horence/juliew/Rfam_software/infernal-1.1.4
TOTAL_RESI=$(esl-seqstat ${2}.fasta | grep "Total # residues" | awk '{print $NF}')
TOTAL_RESI=$(echo "$TOTAL_RESI*2/1000000" | bc -l)

echo "juliew: Start rfam annotation."    

cmscan -Z $TOTAL_RESI --notextw --rfam --cut_ga --nohmmonly --cpu 32 \
    --tblout ${2}.fasta_RFAM.tblout --fmt 2 --clanin $Rfam/Rfam.clanin $Rfam/Rfam.cm ${2}.fasta

echo "juliew: rfam annotation finished."  

# STEP 3: parse infernal output
# Extract relevant lines and reformat with tab separators
head -n -10 ${2}.fasta_RFAM.tblout | 
    awk 'BEGIN {OFS="\t"} NR>2 {$1=$1; print}' | 
    awk 'BEGIN {FS=OFS="\t"} {for (i=28; i<=NF; i++) $27 = $27 " " $i; NF=27; print}' > temp.tsv

echo "juliew: temp file created."

# Create the header for the output file
echo -e "#idx\ttarget name\taccession\tquery name\taccession\tclan name\tmd\tmdl from\tmdl to\tseq from\tseq to\tstran\ttrun\tpass\tgc\tbias\tscore\tE-value\tinc\tolp\tanyidx\tafrct1\tafrct2\twinidx\twfrct1\twfrct2\tdescription of target" > ${2}.fasta_RFAM.tsv
echo "juliew: start append"

# Append the contents of temp.tsv to the output file
cat temp.tsv >> ${2}.fasta_RFAM.tsv

# Clean up temporary file
rm temp.tsv
echo "juliew: temp remove."

# STEP 4: merge structure result with rfam annotation
# python3 /oak/stanford/groups/horence/juliew/structure/src/merge_structure_rfam.py $1 ${2}.fasta_RFAM.tsv
