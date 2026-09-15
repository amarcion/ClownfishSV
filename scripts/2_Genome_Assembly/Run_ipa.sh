#!/usr/bin/env bash

#SBATCH --job-name=ipa
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=15
#SBATCH --mem=80G
#SBATCH --output=outscripts/ipa_%j.out
#SBATCH --error=outscripts/ipa_%j.err
#SBATCH --time=72:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <FastqcFile>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
FastqcFile=${2}			# Fastqc file

##############################

module load gcc/10.4.0
conda activate ipa

echo "${FastqcFile}"

cd "$workindir"

ipa local --nthreads 15 --njobs 2 --tmp-dir /scratch/amarcion/tmp  -i "${FastqcFile}"
