#!/usr/bin/env bash

#SBATCH --job-name=AssStats
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=outscripts/AssemblyStats_%a.out
#SBATCH --error=outscripts/AssemblyStats_%a.err
#SBATCH --time=6:00:00

##############################
# Argument to pass: 
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <wd> <species> <purge_status>" >&2; exit 1; }
[ "$#" -ge 2 ] || usage
mkdir -p outscripts

wd=${1}						# Working directory
species=${2}				# Species to analyse
purge_status=${3:-NoPurged}	# Whether to analyse first assemblies (NoPurged, default) or "Purged" ones

##############################

module load gcc miniconda3
conda activate assemblyStats

cd "$wd"

if [ "$purge_status" = "NoPurged" ]; then
#IPA
assembly-stats -t ./IPA/"$species"/final.a_ctg.fasta ./IPA/"$species"/final.p_ctg.fasta > ../2_Assembly_Statistics/AssemblySummary_"$species".IPA.txt
#Hifiasm
assembly-stats -t ./Hifiasm/"$species"/"$species".hifiasm.bp.hap1.p_ctg.fa ./Hifiasm/"$species"/"$species".hifiasm.bp.hap2.p_ctg.fa ./Hifiasm/"$species"/"$species".hifiasm.bp.p_ctg.fa > ../2_Assembly_Statistics/AssemblySummary_"$species".Hifiasm.txt
#HiCanu
assembly-stats -t ./Canu/"$species"/"$species"_canu.contigs.fasta Canu/"$species"/"$species"_canu.unassembled.fasta > ../2_Assembly_Statistics/AssemblySummary_"$species".Canu.txt
#Flye
assembly-stats -t ./Flye/"$species"/assembly.fasta > ../2_Assembly_Statistics/AssemblySummary_"$species".Flye.txt

else

#Hifiasm Purged
assembly-stats -t ./PurgeDup_Hifiasm/"$species"/"$species".Split.DupPurged.purged.fa ./PurgeDup_Hifiasm/"$species"/"$species".Split.DupPurged.hap.fa > ../2_Assembly_Statistics/AssemblySummary_"$species".PurgedAssembly.txt
#Canu Purged
assembly-stats -t ./PurgeDup_Canu/"$species"/"$species".Canu.DupPurded.purged.fa > ../2_Assembly_Statistics/AssemblySummary_"$species".PurgedCanu.Summary.txt


fi
