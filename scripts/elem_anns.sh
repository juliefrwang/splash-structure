#!/bin/bash
#
#SBATCH --job-name=nextflowscript
#SBATCH --output=launch_nf.%j.out
#SBATCH --error=launch_nf.%j.err
#SBATCH --time=12:00:00
#SBATCH -p horence,normal
#SBATCH --nodes=1
#SBATCH --mem=16Gb
#SBATCH --requeue

# This script is for running element annotations on the cluster. 
# Two arguments are required: 
# 1. the Full path to extendors/compactors list file 
# 2. the output directory.

ml java/18.0.2

anchor=$1
structure_result_dir=$2
outdir=$3
annotation_samplesheet="/oak/stanford/groups/horence/juliew/structure/annotation_nf/element_annotations_samplesheet_5.csv"

cd "$structure_result_dir"/elem_anns
mkdir "$outdir"_work
cd "$outdir"_work
mkdir "$outdir"

/oak/stanford/groups/horence/juliew/workflows/nextflow run salzmanlab/nomad \
    -r element_annotations \
    -profile singularity,horence_no_quake \
    -c /oak/stanford/groups/horence/juliew/structure/annotation_nf/nextflow_element_ann.config \
    --anchors $anchor \
    --outdir $outdir/ \
    -latest \
    -resume \
    --element_annotations_samplesheet $annotation_samplesheet
