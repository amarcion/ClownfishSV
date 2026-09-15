#!/usr/bin/env bash

#SBATCH --job-name=IQtree
#SBATCH --partition=cpu
#SBATCH --output=logs/IQtree.out
#SBATCH --error=logs/IQtree.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --time=48:00:00


set -euo pipefail

module load gcc/12.3.0  iqtree2/2.2.2.7

cd Additional/10_1to1OG_Tree/Concat

iqtree2 -s Concat_1to1OG.ClownDams.Trimmed.fasta -st DNA -o DTR -m GTR+G+I -p Partitions.Concat_1to1OG.ClownDams.txt -bb 10000 -alrt 1000 -ntmax 10



