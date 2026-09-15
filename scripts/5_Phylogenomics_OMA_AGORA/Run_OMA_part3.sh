#!/usr/bin/env bash

#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name=oma3
#SBATCH --output=logs/oma3-%J.log
#SBATCH --export=None
#SBATCH --error=logs/oma2-%J.err

set -euo pipefail

cd /scratch/amarcion/OMA.2.6.0/

./bin/oma
