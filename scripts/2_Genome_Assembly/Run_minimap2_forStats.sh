#!/usr/bin/env bash

#SBATCH --job-name=mapping
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --output=outscripts/minimap2_%j.out
#SBATCH --error=outscripts/minimap2_%j.err
#SBATCH --time=48:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <Assembly> <OutindexFileMmi> <FastqReads> <outAlignmentPrefix>" >&2; exit 1; }
[ "$#" -eq 5 ] || usage
mkdir -p outscripts

workindir=${1}				# Working directory
Assembly=${2}      			# Fasta assembly
OutindexFileMmi=${3}		# output file name index
FastqReads=${4}				# Fastq file
outAlignmentPrefix=${5}		# Prefix for alignment file
##############################

module load gcc/10.4.0 minimap2/2.14 samtools

cd "$workindir"

minimap2 -d "$OutindexFileMmi" "$Assembly" -t 6
minimap2 -t 6 -a "$OutindexFileMmi" "$FastqReads" > "$outAlignmentPrefix".sam

samtools view -@ 6 -Sb "$outAlignmentPrefix".sam > "$outAlignmentPrefix".bam
samtools sort -@ 6  "$outAlignmentPrefix".bam > "$outAlignmentPrefix".Sorted.bam
samtools index -@ 6 "$outAlignmentPrefix".Sorted.bam

# Get mapping statistics and Coverage
samtools stats -@ 6 -c 1,100,5 -in "$outAlignmentPrefix".Sorted.bam > "$outAlignmentPrefix".Stats
cat "$outAlignmentPrefix".Stats  |  grep ^COV | cut -f 2- > "$outAlignmentPrefix".Stats.Cov
