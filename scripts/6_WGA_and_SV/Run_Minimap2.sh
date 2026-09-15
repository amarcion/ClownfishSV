#!/usr/bin/bash --login

#SBATCH --job-name=minimap2
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --output=outscripts/minimap_%j.out
#SBATCH --error=outscripts/minimap_%j.err
#SBATCH --time=48:00:00

##############################
# Arguments
workindir=${1}			# Working directory
Prefix=${2}
Ref=${3}
QueryRef=${4}

module load gcc/12.3.0 minimap2/2.28 micromamba/1.4.2 samtools/1.19.2


cd "$workindir"

minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK "$Ref" "$QueryRef" > "$Prefix".paf

