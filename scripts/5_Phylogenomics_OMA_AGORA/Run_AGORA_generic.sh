#!/usr/bin/env bash

#SBATCH --job-name=AGORAgeneric
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --output=logs/AGORA.Generic.out
#SBATCH --error=logs/AGORA.Generic.err


set -euo pipefail

module load micromamba
# Cluster-specific: point this at your own AGORA conda/micromamba env
AGORA_CONDA_ENV="${AGORA_CONDA_ENV:-/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/Conda/envs/agora}"
micromamba activate "$AGORA_CONDA_ENV"

cd 7_OMA_and_AncestralGenome

./Agora/src/agora-generic.py data/Species.nwk data/orthologyGroups.%s.list data/genes.%s.list \
-workingDir=4b_AGORA_NewFiles/Generic -nbThreads=20  2> AGORA.Generic.log
