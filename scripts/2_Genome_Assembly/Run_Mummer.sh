#!/usr/bin/env bash

#SBATCH --job-name=mummer
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --output=outscripts/mummerV2_%j.out
#SBATCH --error=outscripts/mummerV2_%j.err
#SBATCH --time=10:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <Reference> <Query> <OutFilePrefix>" >&2; exit 1; }
[ "$#" -eq 4 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
Reference=${2}      	# Reference sequences (Assembly)
Query=${3}				# Query sequences (Assembly)
OutFilePrefix=${4}		# output file name index
##############################

module load gcc mummer4/4.0.0rc1
cd "$workindir"

nucmer -b 500 -g 200 -t 16 -p "$OutFilePrefix" "$Reference" "$Query"

# Filter output
delta-filter -l 100 -q "$OutFilePrefix".delta > "$OutFilePrefix".delta.filter

# Get coordinates in different formats
show-coords -roTlH "$OutFilePrefix".delta.filter > "$OutFilePrefix".Filter.Coord.tsv
show-coords -croTlH -L 5000 "$OutFilePrefix".delta.filter > "$OutFilePrefix".Filter2.Coord.tsv
show-coords -roTlH "$OutFilePrefix".delta > "$OutFilePrefix".Coord.tsv
show-coords -c "$OutFilePrefix".delta.filter > "$OutFilePrefix".delta.filter.coords
show-coords -cT "$OutFilePrefix".delta.filter > "$OutFilePrefix".delta.filter.coords.tsv
