#!/usr/bin/env bash

#SBATCH --array=1-18
#SBATCH --partition=cpu
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=hetwind
#SBATCH --output=logs/het_windows_selfmap_%A_%a.out
#SBATCH --error=logs/het_windows_selfmap_%A_%a.err

set -euo pipefail

# ============================================================
# PARAMETERS 
# ============================================================

INPUT_PATH=$1
SPECIES_LIST=$2
WINDOW_SIZE=100000


# ============================================================
# GET CHROMOSOME FOR THIS ARRAY JOB
# ============================================================

SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SPECIES_LIST})
echo "[INFO] Processing chromosome: ${SPECIES}"
echo "[INFO] Array task ID: ${SLURM_ARRAY_TASK_ID}"

# ============================================================
# RUN HET by windows
# ============================================================

module load micromamba
mamba activate Bio

python scripts/3_Genome_Quality/het_windows_selfmap.py ${INPUT_PATH}/${SPECIES}.merge_output.vcf.gz ${SPECIES} ${WINDOW_SIZE}

