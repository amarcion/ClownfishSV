#!/usr/bin/bash --login

#SBATCH --array=1-18
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30GB
#SBATCH --job-name=breakOr
#SBATCH --output=logs/breakOr-%A.%a.log
#SBATCH --error=logs/breakOr-%A.%a.err

module load gcc/12.3.0 samtools/1.19.2

SPECIES_BREAKS=${1}
MAPPING_DIR=${2}
MAPPING_SUFFIX=${3}
WINDOW=${4}
OUT_FOLDER=${5}

mkdir -p ${OUT_FOLDER}

# Get species
read SPECIES CHR BP1 BP2 <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_BREAKS")"

# Get Break points 1
samtools view -F 4 ${MAPPING_DIR}/${SPECIES}.${MAPPING_SUFFIX} ${CHR}:$((BP1-WINDOW))-$((BP1+WINDOW)) |
	awk '{
            flag=$2
            if(and(flag,16)) strand="REV"
            else strand="FWD"
            print strand
        }' | sort | uniq -c > ${OUT_FOLDER}/${SPECIES}.${CHR}_INV.BreakPoints.txt
    

 # Get Break points 2
samtools view -F 4 ${MAPPING_DIR}/${SPECIES}.${MAPPING_SUFFIX} ${CHR}:$((BP2-WINDOW))-$((BP2+WINDOW)) | \
        awk '{
            flag=$2
            if(and(flag,16)) strand="REV"
            else strand="FWD"
            print strand
        }' | sort | uniq -c >> ${OUT_FOLDER}/${SPECIES}.${CHR}_INV.BreakPoints.txt
        