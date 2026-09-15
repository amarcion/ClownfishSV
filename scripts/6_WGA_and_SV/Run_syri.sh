#!/usr/bin/env bash
# Run_syri.sh - whole-genome alignment (MUMmer4 nucmer) + structural-variant
# calling (SyRI) between a reference and one query assembly.
#
# Usage: sbatch Run_syri.sh <workindir> <prefix> <reference.fasta> <query.fasta>
#   workindir        directory to run in (created if missing)
#   prefix           output file prefix (e.g. BIARef_vs_AKA)
#   reference.fasta  reference assembly (SyRI's -r)
#   query.fasta      query assembly (SyRI's -q)
#
# Cluster-specific - edit for your environment:
#   module load gcc mummer4/4.0.0rc1 ; conda activate syri_env
#SBATCH --job-name=syri
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --output=outscripts/syri_%j.out
#SBATCH --error=outscripts/syri_%j.err
#SBATCH --time=12:00:00

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <prefix> <reference.fasta> <query.fasta>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage

workindir=$1
prefix=$2
reference=$3
query=$4

module load gcc mummer4/4.0.0rc1
conda activate syri_env


SYRI_BIN="${SYRI_BIN:-syri}"

mkdir -p "$workindir" outscripts
cd "$workindir"

nucmer --maxmatch -c 500 -b 500 -l 100 -t 12 "$reference" "$query" -p "$prefix"
delta-filter -m -i 90 -l 100 "$prefix".delta > "$prefix".filtered.delta
show-coords -THrd "$prefix".filtered.delta > "$prefix".filtered.coords

"$SYRI_BIN" -c "$prefix".filtered.coords -d "$prefix".filtered.delta -r "$reference" -q "$query"
