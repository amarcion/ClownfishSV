#!/usr/bin/env bash

#SBATCH --array=1-9
#SBATCH --job-name=FastQC
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=outscripts/FastQC_%j.out
#SBATCH --error=outscripts/FastQC_%j.err
#SBATCH --time=2:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <SpeciesFile>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
SpeciesFile=${2}

##############################

species=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SpeciesFile)

module load gcc fastqc

cd "$workindir"

mkdir -p FastQC_report_"$species"

fastqc -o ./FastQC_report_"$species"/ ../HiFiReads/"$species".hifi_reads.fastq.gz


