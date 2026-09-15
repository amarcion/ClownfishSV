#!/usr/bin/env bash

#SBATCH --job-name=PurgeDub_Step1
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --output=outscripts/PurgeDub_Step1_%j.out
#SBATCH --error=outscripts/PurgeDub_Step1_%j.err
#SBATCH --time=48:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <newDir> <assembly> <FastqReads> <OutNamePrefix>" >&2; exit 1; }
[ "$#" -eq 5 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
newDir=${2}				# New directory for output
assembly=${3}			# Assembly Name
FastqReads=${4}			# Fastq Reads
OutNamePrefix=${5}		# Output prefix

##############################

module load gcc/10.4.0 minimap2/2.14

cd "$workindir"
mkdir -p "$newDir"

minimap2 -t 12 -xasm20 "$assembly" "$FastqReads" | gzip -c - > "$newDir"/"$OutNamePrefix".paf.gz

cd "$newDir"
~/purge_dups/bin/pbcstat "$OutNamePrefix".paf.gz
~/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log




