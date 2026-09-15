#!/usr/bin/env bash

set -euo pipefail

usage() { echo "Usage: bash $(basename "$0") <filename>" >&2; exit 1; }
[ "$#" -eq 1 ] || usage

filename=$1

while IFS= read -r line || [[ -n "$line" ]]; do
    # Skip empty lines
    if [ -z "$line" ]; then
        continue
    fi
    
    cd "$line"/braker
    cp braker.aa ../../"$line".braker.aa
    cp braker.codingseq ../../"$line".braker.codingseq
    cp braker.gtf ../../"$line".braker.gtf

done<"$filename"
