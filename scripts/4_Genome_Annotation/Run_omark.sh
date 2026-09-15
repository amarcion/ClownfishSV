#!/usr/bin/env bash

#SBATCH --job-name=job
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=outscripts/%j.out
#SBATCH --error=outscripts/%j.err
#SBATCH --time=1:00:00


##############################

module load micromamba
micromamba activate OMArk

# Arguments

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <wd> <species> <db>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage
mkdir -p outscripts

wd=${1}
species=${2}
db=${3}

echo "$arg"

mkdir -p Omark_output/"$species".Annotation_v2."$db"

omark -f Omamer_Out/"$species".Annotation_v2."$db".omamer -d "$db".h5 -o Omark_output/"$species".Annotation_v2."$db" -i Splice_files/"$species".Annotation_v2.splice


