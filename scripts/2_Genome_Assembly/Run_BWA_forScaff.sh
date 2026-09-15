#!/usr/bin/env bash

#SBATCH --job-name=bwa
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --output=outscripts/bwaScaff.%j.out
#SBATCH --error=outscripts/bwaScaff.%j.err
#SBATCH --time=15:00:00

##############################
# Arguments

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <wd> <speciesPrefix> <LongReadsPrefix> <outFile>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

wd=${1}
speciesPrefix=${2}
LongReadsPrefix=${3}
outFile=${4}

######
SECONDS=0

module load gcc bwa samtools  bamtools picard

cd "$wd"

bwa mem -M -t 6 "$speciesPrefix".Softmasked "$LongReadsPrefix".w1kb.R1.fastq "$LongReadsPrefix".w1kb.R2.fastq  > "$outFile"

echo "The command took $SECONDS s."


