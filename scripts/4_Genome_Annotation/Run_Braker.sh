#!/usr/bin/env bash

#SBATCH --job-name=braker
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=80G
#SBATCH --output=outscripts/Braker_%j.out
#SBATCH --error=outscripts/Braker_%j.err
#SBATCH --time=32:00:00

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <Assembly> <Species>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage
mkdir -p outscripts

workindir=${1}	
Assembly=${2}	
Species=${3}	

module load gcc miniconda3
conda activate braker_envs

# Cluster-specific: point these at your own ProtHint / GeneMark-ES / perl5 installs
PROTHINT_BIN_DIR="${PROTHINT_BIN_DIR:-/work/FAC/FBM/DBC/nsalamin/clownfish/amarcion/Software/ProtHint/bin}"
GENEMARK_DIR="${GENEMARK_DIR:-/work/FAC/FBM/DBC/nsalamin/clownfish/amarcion/Software/gmes_linux_64_4/}"
PERL5LIB_DIR="${PERL5LIB_DIR:-/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/Conda/envs/braker_envs/lib/perl5/5.32/}"
export PATH=$PATH:"$PROTHINT_BIN_DIR":"$GENEMARK_DIR"
export GENEMARK_PATH="$GENEMARK_DIR"
export PERL5LIB="$PERL5LIB_DIR"

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

mkdir -p "$workindir"
cd "$workindir"

braker.pl --genome="$Assembly" --prot_seq=../Eukaryota.fa --species="$Species"  --threads 15


