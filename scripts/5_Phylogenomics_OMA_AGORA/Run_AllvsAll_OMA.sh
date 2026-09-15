#!/usr/bin/env bash

#SBATCH --array=1-500
#SBATCH --partition=cpu
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --job-name=oma2
#SBATCH --output=logs/oma2-%A.%a.log
#SBATCH --export=None
#SBATCH --error=logs/oma2-%A.%a.err

set -euo pipefail

cd /scratch/amarcion/OMA.2.6.0/

export NR_PROCESSES=500
./bin/oma -s -W 7000
if [[ "$?" == "99" ]] ; then
scontrol requeue \
${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
fi
exit 0
