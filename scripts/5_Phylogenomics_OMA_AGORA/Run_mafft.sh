#!/usr/bin/env bash

#SBATCH --job-name=mafft
#SBATCH --partition=cpu
#SBATCH --time=30:00
#SBATCH --output=logs/mafft_%A_%a.out
#SBATCH --error=logs/mafft_%A_%a.err
#SBATCH --array=1-2%100     # 15,000 jobs, max 100 running at once

module load gcc/12.3.0 mafft/7.505 trimal/1.4.1

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <WD> <OG_FILE> <IN_FOLDER> <OUT_FOLDER> <OUT_FOLDER_TRIMAL>" >&2; exit 1; }
[ "$#" -eq 5 ] || usage

WD=${1}
OG_FILE=${2}
IN_FOLDER=${3}
OUT_FOLDER=${4}
OUT_FOLDER_TRIMAL=${5}

cd ${WD}
mkdir -p ${OUT_FOLDER}
mkdir -p ${OUT_FOLDER_TRIMAL}

# Get the OG corresponding to the OG
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${OG_FILE})
BASENAME=$(basename "$FILE" .fa)
OUT_NAME=${OUT_FOLDER}/${BASENAME}.aln.fa
OUT_NAME_TIMAL=${OUT_FOLDER_TRIMAL}/${BASENAME}.Trimaln.fa

mafft ${IN_FOLDER}/$FILE > ${OUT_NAME}

trimal -in ${OUT_NAME} -out ${OUT_NAME_TIMAL} -gt 0.75 -st 0.001

