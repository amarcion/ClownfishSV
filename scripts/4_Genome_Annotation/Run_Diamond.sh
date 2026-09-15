#!/usr/bin/env bash

#SBATCH --job-name=diamond
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --output=outscripts/Diamond_%j.out
#SBATCH --error=outscripts/Diamond_%j.err
#SBATCH --time=24:00:00

module purge
dcsrsoft use arolle

module load gcc diamond/2.0.15

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <query> <database> <outPrefix>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}	
query=${2}	
database=${3}	
outPrefix=${4}	


cd "$workindir"

diamond blastp -q "$query" -d "$database" -o "$outPrefix".Max1.diamond.blastp -p 6 --outfmt 6 --un "$outPrefix".Max1.unaligned --max-target-seqs 1
diamond blastp -q "$query" -d "$database" -o "$outPrefix".Max5.diamond.blastp -p 6 --outfmt 6 --un "$outPrefix".Max5.unaligned --max-target-seqs 5


