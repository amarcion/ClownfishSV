#!/usr/bin/env bash

#SBATCH --array=1-16
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50GB
#SBATCH --job-name=minimap
#SBATCH --output=logs/minimap2-%A.%a.log
#SBATCH --error=logs/minimap2-%A.%a.err

module load gcc/12.3.0 minimap2/2.28 samtools/1.19.2

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <speciesfile> <wd>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage

speciesfile=${1}
wd=${2}

# Get species
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$speciesfile")


# Align the reads

minimap2 -ax map-hifi -t 16 Final_Assemblies/${species}.Assembly_v2.fa hifi_reads/${species}.hifi_reads.fastq -o minimap2_Align/${species}.minimap2.sam

