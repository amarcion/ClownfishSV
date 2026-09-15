#!/usr/bin/env bash

#SBATCH --array=1-50
#SBATCH --job-name=Scaff
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --output=outscripts/ManualChrScaff_%a_%j.out
#SBATCH --error=outscripts/ManualChrScaff_%a_%j.err
#SBATCH --time=12:00:00

##############################
# Argument to pass: 
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <wd> <species> <threads> <ReadsPath> <arg5>" >&2; exit 1; }
[ "$#" -eq 5 ] || usage
mkdir -p outscripts

wd=${1}
species=${2}	
threads=${3}
ReadsPath=${4}
ChromToAnalyzeFile=${5-"$species".ChromToAnalyze.txt}
##############################

module purge
dcsrsoft use arolle
module load gcc python miniconda3

cd "$wd"

echo "$ChromToAnalyzeFile"

originalwd=$(pwd)
ChromToAnalyze=$(sed -n ${SLURM_ARRAY_TASK_ID}p $ChromToAnalyzeFile)


# We prepare the files by joining all scaffolds in a chromosome
cd "$species"/"$ChromToAnalyze"
chromoFastaTemp="$species".temp."$ChromToAnalyze".fa
cat *.fa > "$chromoFastaTemp"
ln -s "$originalwd"/"$ReadsPath"/"$species".hifi_reads.fastq


####
# Scaffolding with long reads
####

conda activate samba_env
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

mkdir -p temp_samba

samba.sh -r "$chromoFastaTemp" -q "$species".hifi_reads.fastq-t "$threads" -d asm -m 2500

# Move temporary files and make soft link to final scaffolds
[ ! -f minimap.err ] || mv minimap.err temp_samba
[ 'scaffold_*' = "$(echo scaffold_*)" ] || mv scaffold_* ./temp_samba 
[ '*.fa.*' = "$(echo *.fa.*)" ] || mv *.fa.* ./temp_samba 
cp ./temp_samba/"$chromoFastaTemp".scaffolds.fa "$species".temp."$ChromToAnalyze".Samba.fa


conda deactivate


####
# Scaffolding with ntLink
####

conda activate Tigmint
mkdir -p temp_ntLink

ntLink scaffold target="$species".temp."$ChromToAnalyze".Samba.fa reads="$species".hifi_reads.fastq t="$threads" k=32 w=100 z=1000 n=2 a=1 conservative=True

# Move temporary files and make soft link to final scaffolds
[ '*.Samba.fa.*' = "$(echo *.Samba.fa.*)" ] || mv *.Samba.fa.* ./temp_ntLink 
cp ./temp_ntLink/"$species".temp."$ChromToAnalyze".Samba.fa.k32.w100.z1000.stitch.abyss-scaffold.fa "$species".temp."$ChromToAnalyze".Samba.NtLink.fa


####
# Statistics
####

conda activate assemblyStats
assembly-stats -t *.fa > Scaffolds_Statistics."$scaffoldID".txt




