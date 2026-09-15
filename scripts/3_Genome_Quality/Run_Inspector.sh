#!/usr/bin/env bash

#SBATCH --job-name=Inspector
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --output=outscripts/Inspector_%j.out
#SBATCH --error=outscripts/Inspector_%j.err
#SBATCH --time=12:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <Assembly> <Reads> <OutFolder>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
Assembly=${2}      		# Fasta assembly
Reads=${3}				# HiFi reads
OutFolder=${4}			# output file name index
##############################

# Cluster-specific: point this at your own Inspector install
pathToScript="${INSPECTOR_PATH:-/users/amarcion/work/clownfish/amarcion/Software/Inspector}"

micromamba activate inspector
module load minimap2 samtools
export PATH="$pathToScript":$PATH

cd "$workindir"

"$pathToScript"/inspector.py -c "$Assembly" -r "$Reads" -t 12 -o "$OutFolder"  --datatype hifi

"$pathToScript"/inspector-correct.py -i "$OutFolder" -t 12 --datatype pacbio-hifi 
