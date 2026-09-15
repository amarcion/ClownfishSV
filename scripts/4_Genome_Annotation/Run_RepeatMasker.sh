#!/usr/bin/env bash

#SBATCH --job-name=RepeatMask
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --output=outscripts/RepeatMask_%j.out
#SBATCH --error=outscripts/RepeatMask_%j.err
#SBATCH --time=48:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <FastaAssembly> <Library> <OutDir> <TypeMask>" >&2; exit 1; }
[ "$#" -eq 5 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
FastaAssembly=${2}
Library=${3}      # Fasta assembly
OutDir=${4}
TypeMask=${5}

##############################

conda activate RepMask_envs
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

cd "$workindir"
mkdir -p "$OutDir"

if [[ "$TypeMask" == "SM" ]]; then
    RepeatMasker -lib "$Library" -dir "$OutDir" -xsmall -pa 15 "$FastaAssembly"
else
    RepeatMasker -lib "$Library" -dir "$OutDir" -pa 15 "$FastaAssembly"
fi

