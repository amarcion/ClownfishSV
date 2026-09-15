#!/usr/bin/env bash

#SBATCH --array=1-16
#SBATCH --partition=cpu
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20GB
#SBATCH --job-name=sam2bam
#SBATCH --output=logs/sam2bam-%A.%a.log
#SBATCH --error=logs/sam2bam-%A.%a.err

module load gcc/12.3.0 samtools/1.19.2

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <speciesfile> <wd>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage

speciesfile=${1}
wd=${2}

# Get species
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$speciesfile")

cd ${wd}

# sam 2 bam
samtools view -@ 16 -bS ${species}.minimap2.sam > ${species}.minimap2.bam
samtools sort -@ 16 ${species}.minimap2.bam > ${species}.minimap2.sorted.bam
samtools index -@ 16 ${species}.minimap2.sorted.bam