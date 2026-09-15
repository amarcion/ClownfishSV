#!/usr/bin/env bash

#SBATCH --job-name=PurgeDub_Step2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=50G
#SBATCH --output=outscripts/PurgeDub_Step2_%j.out
#SBATCH --error=outscripts/PurgeDub_Step2_%j.err
#SBATCH --time=48:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <workindir> <OutputPrefix> <Assembly>" >&2; exit 1; }
[ "$#" -eq 3 ] || usage
mkdir -p outscripts

workindir=${1}			# Working directory
OutputPrefix=${2}		# Output prefix
Assembly=${3}			# Assembly file name
##############################

module load gcc/10.4.0 minimap2/2.14

cd "$workindir"

~/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov "$OutputPrefix".self.paf.gz > dups.bed 2> purge_dups.log
~/purge_dups/bin/get_seqs -p "$OutputPrefix".DupPurged -e dups.bed "$Assembly"