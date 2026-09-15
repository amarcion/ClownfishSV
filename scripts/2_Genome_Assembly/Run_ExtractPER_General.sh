#!/usr/bin/env bash

#SBATCH --job-name=job
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --output=outscripts/%j.out
#SBATCH --error=outscripts/%j.err
#SBATCH --time=5:00:00

##############################
# Arguments

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <wd> <speciesPrefix> <dvalue> <lvalue> <bvalue> <outprefix> <arg7>" >&2; exit 1; }
[ "$#" -eq 7 ] || usage
mkdir -p outscripts

wd=${1}				# Working directory
speciesPrefix=${2}	# Species Prefix
dvalue=${3}			# Distance between paired ends
lvalue=${4}			# Minimum HiFi read length
bvalue=${5}			# Starting position for read extraction
outprefix=${6}		# Output prefix
scriptPath=${7-/home/amarcion/}

######
SECONDS=0

cd "$wd"

python "$scriptPath"/Extract_PE_from_HiFI.py -w 200 -s 1000 -d "$dvalue" -l "$lvalue" -b "$bvalue" -t Separated "$speciesPrefix".hifi_reads.fastq "$speciesPrefix"."$outprefix"


echo "The command took $SECONDS s."
