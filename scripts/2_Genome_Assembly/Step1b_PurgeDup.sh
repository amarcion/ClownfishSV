#!/usr/bin/env bash

#SBATCH --job-name=PurgeDub_Step2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --output=outscripts/PurgeDub_Step2_%j.out
#SBATCH --error=outscripts/PurgeDub_Step2_%j.err
#SBATCH --time=48:00:00


##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <assembly> <SplittedAssembly> <OutNamePrefix>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
assembly=${2}			# Assembly name
SplittedAssembly=${3}	# Name of splitted assembly
OutNamePrefix=${4}		# Outfile prefix
##############################

module load gcc/10.4.0 minimap2/2.14

cd "$workindir"

~/purge_dups/bin/split_fa "$assembly" > "$SplittedAssembly"
minimap2 -t 12 -xasm5 -DP "$SplittedAssembly" "$SplittedAssembly" | gzip -c - > "$OutNamePrefix".self.paf.gz