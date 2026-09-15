#!/usr/bin/env bash

#SBATCH --job-name=Busco
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --output=outscripts/Busco_%j.out
#SBATCH --error=outscripts/Busco_%j.err
#SBATCH --time=10:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <Assembly> <OutPrefix>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
Assembly=${2}      # Fasta assembly
OutPrefix=${3}			# output file name index

##############################

module load gcc/10.4.0 bwa/0.7.17 python/3.9.13 miniconda3/4.10.3
conda activate busco

cd "$workindir"

busco -i "$Assembly" -l actinopterygii_odb10 -o "$OutPrefix" -m genome -c 6 --tar --offline


