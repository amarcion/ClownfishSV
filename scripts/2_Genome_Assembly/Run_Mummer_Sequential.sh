#!/usr/bin/env bash

#SBATCH --job-name=MummerSeq
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --output=outscripts/mummerSeq_%j.out
#SBATCH --error=outscripts/mummerSeq_%j.err
#SBATCH --time=32:00:00

##############################
# Arguments
set -euo pipefail

usage() { echo "Usage: sbatch $(basename "$0") <InputFileSequence> <TempDirectory> <out_prefix> <min_coverge_align> <min_Meanidentity> <arg6>" >&2; exit 1; }
[ "$#" -eq 6 ] || usage
mkdir -p outscripts

InputFileSequence=${1}			
TempDirectory=${2}      
out_prefix=${3}
min_coverge_align=${4}
min_Meanidentity=${5}
pythonScript=${6-Sequential_mummer.py}		

##############################

module load gcc mummer4/4.0.0rc1 miniconda3

conda activate Bio

mkdir -p $TempDirectory

python "$pythonScript" $InputFileSequence $TempDirectory $out_prefix $min_coverge_align $min_Meanidentity
#Ex:
#python Sequential_mummer.py CfAKA94.provaScaff.fa ./tmp/ CfAKA94.provaMummer 0.8 80
