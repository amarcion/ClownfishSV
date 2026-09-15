#!/usr/bin/env bash

#SBATCH --job-name=Scaff
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=outscripts/Scaff.%j.out
#SBATCH --error=outscripts/Scaff.%j.err
#SBATCH --time=10:00:00

##############################
# Arguments

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <wd> <inputBam> <inputAssembly> <ReadPrefix> <OutPrefix>" >&2; exit 1; }
[ "$#" -eq 5 ] || usage
mkdir -p outscripts

wd=${1}
inputBam=${2}
inputAssembly=${3}
ReadPrefix=${4}
OutPrefix=${5}

######
SECONDS=0

module load gcc perl bwa

cd "$wd"

sh ~/Softwares/P_RNA_scaffolder/P_RNA_scaffolder.sh -d ~/Softwares/P_RNA_scaffolder -i "$inputBam" -j "$inputAssembly" -F "$ReadPrefix".w1kb.R1.fastq -R "$ReadPrefix".w1kb.R2.fastq -b no -o "$OutPrefix".out


echo "The command took $SECONDS s."


