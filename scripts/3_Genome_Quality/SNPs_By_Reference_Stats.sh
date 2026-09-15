#!/usr/bin/env bash
#SBATCH --job-name=vcf_stats
#SBATCH --array=1-18
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00

set -euo pipefail

module load micromamba bcftools/1.21 vcftools/0.1.16

# Cluster-specific: point this at your own conda/mamba env with samtools/bcftools/etc.
BIO_CONDA_ENV="${BIO_CONDA_ENV:-/work/FAC/FBM/DBC/nsalamin/software/nsalamin/clownfishes/Conda/envs/Bio/}"
mamba activate "$BIO_CONDA_ENV"


SPECIES=(AKA AKY ALL BIA CLA CRP EPH FRE LAT LAZ MCC OCE OMA POL PRC PRD SAN SEB)
SP=${SPECIES[$SLURM_ARRAY_TASK_ID-1]}

VCF=${SP}_output/merge_output.vcf.gz
BAM=${SP}.minimap2.sorted.bam 
REF=${SP}.Assembly_v2.fa
OUTDIR=${SP}
mkdir -p ${OUTDIR}

# 1. bcftools stats - overall VCF statistics
bcftools stats ${VCF} > ${OUTDIR}/${SP}_bcftools_stats.txt

# 2. vcftools - heterozygosity per site
vcftools --gzvcf ${VCF} \
  --het \
  --out ${OUTDIR}/${SP}

# 3. vcftools - SNP density in 100kb windows
vcftools --gzvcf ${VCF} \
  --SNPdensity 100000 \
  --out ${OUTDIR}/${SP}

# 4. vcftools - allele frequency
vcftools --gzvcf ${VCF} \
  --freq2 \
  --out ${OUTDIR}/${SP}

# 5. vcftools - depth per site
vcftools --gzvcf ${VCF} \
  --site-depth \
  --out ${OUTDIR}/${SP}

# 6. mosdepth - windowed depth from BAM (100kb windows)
mosdepth \
  --threads ${SLURM_CPUS_PER_TASK} \
  --by 100000 \
  --quantize 0:5:30:200: \
  ${OUTDIR}/${SP}_mosdepth \
  ${BAM}