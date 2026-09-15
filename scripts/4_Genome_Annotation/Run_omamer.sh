#!/usr/bin/env bash

#SBATCH --job-name=job
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=outscripts/%j.out
#SBATCH --error=outscripts/%j.err
#SBATCH --time=1:00:00


##############################

module load micromamba
micromamba activate OMArk

# Arguments

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <inputAA> <OutputFile> <db>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage
mkdir -p outscripts

inputAA=${1}
OutputFile=${2}
db=${3}

echo "$arg"

omamer search --db "$db".h5 --query "$inputAA" --out "$OutputFile"


