#!/usr/bin/env bash

#SBATCH --job-name=InterProScan
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G
#SBATCH --output=outscripts/InterProScan_%j.out
#SBATCH --error=outscripts/InterProScan_%j.err
#SBATCH --time=24:00:00

module load gcc miniconda3

conda activate interproscan
# Cluster-specific: point this at your own InterProScan install
INTERPROSCAN_DIR="${INTERPROSCAN_DIR:-/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/amarcion/interproscan-5.66-98.0/}"
export PATH=$PATH:"$INTERPROSCAN_DIR"

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <species> <query> <analyses>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}	
species=${2} 
query=${3}	
analyses=${4}


cd "$workindir"
mkdir -p InterProScan_"$species"
mkdir -p InterProScan_"$species"/tmp

interproscan.sh -appl "$analyses" -i "$query" -cpu 16 -d ./InterProScan_"$species" -goterms -T InterProScan_"$species"/tmp 
