#!/usr/bin/env bash
# Submit_Jellyfish.sh - submit one Run_JellyFish.sh job per species.
#
# Usage: bash Submit_Jellyfish.sh <species_list.txt>
#   species_list.txt   one species code per line (e.g. metadata/Species.txt)

set -euo pipefail

usage() { echo "Usage: bash $(basename "$0") <species_list.txt>" >&2; exit 1; }
[ "$#" -eq 1 ] || usage
species_list=$1


while IFS= read -r species || [ -n "$species" ]; do
    [ -z "$species" ] && continue
    echo "$species"
    sbatch scripts/1_PacBio_Data_and_Preprocessing/Run_JellyFish.sh \
        1_PacBioData_and_Preprocessing/JellyFish \
        ../HiFiReads/"$species".T2.hifi_reads.fastq.gz \
        "$species"_k21 21
done < "$species_list"
