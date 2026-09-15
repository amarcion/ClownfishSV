#!/usr/bin/env bash

set -euo pipefail

usage() { echo "Usage: bash $(basename "$0") <filename> <arg2>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage

filename=$1

# Read length is the mean HiFi read length for each species
declare -A ReadLength=( [AKA]=14357
[AKY]=18203
[ALL]=13422
[BIA]=8949
[CRP]=13395
[EPH]=17097
[FRE]=13275
[LAT]=16529
[LAZ]=9427
[MCC]=12875
[OMA]=14177
[PRC]=11756
[PRD]=16514
[POL]=17085
[SAN]=10820
[SEB]=15823
)


while IFS= read -r species || [[ -n "$species" ]]; do
    [ -z "$species" ] && continue

    length=${ReadLength[$species]}

    if [ -z "$length" ]; then
        echo "No length found for $species"
        continue
    fi

    echo "Species: $species | Length: $length"

Rscript GenomeScope.R "$species"_k21.histo 21 "$length" "$species"

done<"$filename"
