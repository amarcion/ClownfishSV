#!/usr/bin/env bash

#SBATCH --job-name=Jelly
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --output=outscripts/Jelly_%j.out
#SBATCH --error=outscripts/Jelly_%j.err
#SBATCH --time=48:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <FastqFile> <OutPrefix> <kpar>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
FastqFile=${2}			# Fastqc file
OutPrefix=${3}      # OutputPrefix
kpar=${4}

##############################

module load gcc jellyfish/2.2.7

echo "${OutPrefix}"
echo "${FastqFile}"

cd "$workindir"

jellyfish count -C -m "$kpar" -s 100M -t 12 -o "${OutPrefix}".jf <(zcat "${FastqFile}")
jellyfish histo -t 12 "${OutPrefix}".jf >  "${OutPrefix}".histo


