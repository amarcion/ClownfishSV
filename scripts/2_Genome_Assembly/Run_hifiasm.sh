#!/usr/bin/env bash
# Run_hifiasm.sh - assemble one species' HiFi reads with Hifiasm, GFA -> FASTA.
#
# Usage: sbatch Run_hifiasm.sh <workindir> <hifi_reads.fastq> <out_prefix>
#   workindir         directory to run in (created if missing)
#   hifi_reads.fastq  path to the HiFi reads, relative to workindir
#   out_prefix        prefix for hifiasm's output files
#
# Cluster-specific - edit for your environment:
#   module load gcc miniconda3 ; conda activate hifiasm
#SBATCH --job-name=hifiasm
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --output=outscripts/hifiasm_%j.out
#SBATCH --error=outscripts/hifiasm_%j.err
#SBATCH --time=72:00:00

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <hifi_reads.fastq> <out_prefix>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage

workindir=$1
fastq_file=$2
out_prefix=$3

module load gcc miniconda3
conda activate hifiasm

mkdir -p "$workindir"
cd "$workindir"
mkdir -p outscripts

echo "reads:  $fastq_file"
echo "prefix: $out_prefix"

hifiasm -o "$out_prefix" -t 15 "$fastq_file"

# GFA -> FASTA: keep only the sequence ("S") lines
awk '/^S/ { print ">"$2; print $3 }' "$out_prefix".bp.p_ctg.gfa > "$out_prefix".bp.p_ctg.fa
