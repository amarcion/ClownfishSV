#!/usr/bin/env bash

#SBATCH --array=1-16
#SBATCH --partition=cpu
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20GB
#SBATCH --job-name=sniffle2
#SBATCH --output=logs/sniffle2-%A.%a.log
#SBATCH --error=logs/sniffle2-%A.%a.err

module load micromamba
# Cluster-specific: point this at your own Sniffles2 conda/mamba env
SNIFFLE_CONDA_ENV="${SNIFFLE_CONDA_ENV:-/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/einaciom/envs/sniffle}"
mamba activate "$SNIFFLE_CONDA_ENV"

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <speciesfile> <wd>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage

speciesfile=${1}
wd=${2}

# Get species
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$speciesfile")

cd ${wd}

# sam 2 bam
sniffles -i minimap2_Align/${species}.minimap2.sorted.bam -v ${species}.sniffle_out.vcf