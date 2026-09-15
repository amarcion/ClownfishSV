#!/usr/bin/env bash

#SBATCH --job-name=EDTA
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --output=outscripts/EDTA_%j.out
#SBATCH --error=outscripts/EDTA_%j.err
#SBATCH --time=48:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <FastaAssembly>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage
mkdir -p outscripts

workindir=${1}		# Working directory
FastaAssembly=${2}

##############################

# Cluster-specific: point this at your own EDTA install
pathToEDTAscript="${EDTA_PATH:-/users/amarcion/work/clownfish/amarcion/Software/EDTA}"

module load gcc python miniconda3
conda activate edta_envs

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8


cd "$workindir"

perl "$pathToEDTAscript"/EDTA.pl --genome "$FastaAssembly" --overwrite 1 --sensitive 1 --anno 1 --evaluate 0 --threads 12