#!/usr/bin/env bash

#SBATCH --job-name=flye
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=80G
#SBATCH --output=outscripts/flye_%j.out
#SBATCH --error=outscripts/flye_%j.err
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

echo "${FastqcFile}"

cd "$workindir"

~/Flye/bin/flye --pacbio-hifi "${FastqcFile}" -o . --scaffold -t 15

