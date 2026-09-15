#!/usr/bin/env bash

#SBATCH --array=1-9
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30GB
#SBATCH --job-name=clair3
#SBATCH --output=logs/clair3-%A.%a.log
#SBATCH --error=logs/clair3-%A.%a.err

module load micromamba
mamba activate clair3


set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <speciesfile> <input_dir> <output_dir> <model_name>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage

speciesfile=${1}
input_dir=${2}
output_dir=${3}
model_name=${4}

# Get species
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$speciesfile")

mkdir -p ${output_dir}/${species}_output

INPUT_DIR=${input_dir}
OUTPUT_DIR=${output_dir}/${species}_output
THREADS=$SLURM_CPUS_PER_TASK
MODEL_PATH=/scratch/amarcion/Long_Read_Mapping/SNPs_clair3/models/${model_name}

run_clair3.sh --bam_fn=${INPUT_DIR}/${species}.minimap2.sorted.bam \
	--ref_fn=${INPUT_DIR}/${species}.Assembly_v2.fa \
	--threads=${THREADS} \
	--platform=hifi \
	--model_path=${MODEL_PATH} \
	--output=${OUTPUT_DIR} \
	--include_all_ctgs

