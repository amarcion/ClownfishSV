#!/usr/bin/env bash

#SBATCH --job-name=bwa
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --output=outscripts/bwaIndex.%j.out
#SBATCH --error=outscripts/bwaIndex.%j.err
#SBATCH --time=3:00:00

##############################
# Arguments

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <wd> <speciesPrefix> <Assembly>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage
mkdir -p outscripts

wd=${1}
speciesPrefix=${2}
Assembly=${3}

######
SECONDS=0

module load gcc bwa samtools bamtools picard

cd "$wd"

bwa index -p "$speciesPrefix" "$Assembly"


echo "The command took $SECONDS s."


