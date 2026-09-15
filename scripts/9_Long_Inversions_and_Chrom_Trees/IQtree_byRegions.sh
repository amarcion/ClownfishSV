#!/usr/bin/env bash

#SBATCH --array=1-54
#SBATCH --partition=cpu
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30GB
#SBATCH --job-name=iqtree
#SBATCH --output=logs/iqtree-%A.%a.log
#SBATCH --error=logs/iqtree-%A.%a.err

module load gcc/12.3.0  iqtree2/2.2.2.7

set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <filename> <wd>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage

filename=${1}
wd=${2}



fastafile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$filename")

# Skip empty lines
if [ -z "$fastafile" ]; then
    echo "Empty line at index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing file: $fastafile"

cd $wd

# Main command
iqtree2 -s ${fastafile} -m GTR+ASC -B 1000 -alrt 1000 -ntmax 10