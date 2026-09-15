#!/usr/bin/env bash

#SBATCH --job-name=PrepareByChrom
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=outscripts/PrepareByChrom_%j.out
#SBATCH --error=outscripts/PrepareByChrom_%j.err
#SBATCH --time=01:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <Assembly> <Species> <arg4>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
Assembly=${2}			# Assembly file
Species=${3}			# Species 
PythonScript=${4-scripts/2_Genome_Assembly/Create_Assembly_By_Chrom.py}
##############################

conda activate Bio

originalDirectory=$(pwd)
pathToPythonScript="$originalDirectory"/"$PythonScript"
echo "$pathToPythonScript"

cd "$workindir"
mkdir -p "$Species"

python "$pathToPythonScript" "$Assembly" "$Species"/ "$Species"

