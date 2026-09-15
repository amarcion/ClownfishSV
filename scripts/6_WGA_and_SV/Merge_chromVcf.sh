#!/usr/bin/env bash

#SBATCH --array=1-23
#SBATCH --partition=cpu
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name=merge
#SBATCH --output=logs/mergevcf-%A.%a.log
#SBATCH --export=None
#SBATCH --error=logs/mergevcf-%A.%a.err

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <filename> <wd>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage

filename=$1
wd=$2
inputPrefix=$3
outputPrefix=$4

conda activate tabix

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$filename")

# Skip empty lines
if [ -z "$chrom" ]; then
    echo "Empty line at index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing chromosome: $chrom"

cd $wd

# Main command
vcf-merge "$inputPrefix"."$chrom".vcf.gz > "$outputPrefix"."$chrom".vcf

