#!/usr/bin/env bash

#SBATCH --job-name=RepeatMod
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --output=outscripts/RepeatMod_%j.out
#SBATCH --error=outscripts/RepeatMod_%j.err
#SBATCH --time=72:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <FastaAssembly> <DBName>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
FastaAssembly=${2}
DBName=${3}      # Fasta assembly

##############################

conda activate RepMask_envs
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

cd "$workindir"

mkdir -p "$DBName"

# Build database
BuildDatabase -name "$DBName"/"$DBName" "$FastaAssembly"
RepeatModeler -database "$DBName"/"$DBName" -pa 4 -LTRStruct

