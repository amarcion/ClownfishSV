#!/usr/bin/env bash

#SBATCH --job-name=canu
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --output=outscripts/canu_%j.out
#SBATCH --error=outscripts/canu_%j.err
#SBATCH --time=72:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <FastqcFile> <GenomeSize> <Prefix>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
FastqcFile=${2}			# Fastqc file
GenomeSize=${3}			# Expected Genome size
Prefix=${4}			# Predix for output

##############################

module load gcc/10.4.0
conda activate canu

echo "${FastqcFile}"
echo "${GenomeSize}"

cd "$workindir"

canu -p "${Prefix}" -d  genomeSize="${GenomeSize}" -pacbio-hifi "${FastqcFile}"