#!/usr/bin/env bash

module load tabix/2013-12-16

set -euo pipefail

usage() {
  echo "Usage: bash $(basename "$0") <species_list> <wd> <outdir> <input_pattern> <svtypes>" >&2
  echo "  species_list   file with one species name per line" >&2
  echo "  wd             working directory containing input VCFs" >&2
  echo "  outdir         output directory name (created under wd)" >&2
  echo "  input_pattern  input VCF filename pattern, use {species} as placeholder" >&2
  echo "                 e.g. 'BIARef_vs_{species}.Splitvcf.SR.vcf.gz'" >&2
  echo "  svtypes        comma-separated SV types to keep (no <>)," >&2
  echo "                 e.g. 'INV,DUP,INVDP,TRANS,INVTR'" >&2
  exit 1
}
[ "$#" -eq 5 ] || usage

species_list=$1
wd=$2
outdirname=$3
input_pattern=$4
svtypes=$5

curr_dir=$(pwd)

cd "$wd"
mkdir -p "$outdirname"

while IFS= read -r species || [[ -n "$species" ]]; do
    [ -z "$species" ] && continue

    echo "Processing: $species"

    infile=${input_pattern//\{species\}/$species}
    if [ ! -f "$infile" ]; then
        echo "WARNING: $infile not found, skipping" >&2
        continue
    fi

    prefix=$(basename "$infile" .vcf.gz)

    zcat "$infile" | \
    awk -v outdir="$outdirname" -v prefix="$prefix" -v svtypes="$svtypes" '
    BEGIN {
        n = split(svtypes, arr, ",")
        for (i = 1; i <= n; i++) allowed["<" arr[i] ">"] = 1
    }

    /^#/ {
        header[++h] = $0
        next
    }

    {
        chr = $1
        svtype = $5
        if (!(svtype in allowed)) next

        gsub(/[<>]/, "", svtype)
        fname = outdir "/" prefix "." chr "." svtype ".vcf"

        if (!(fname in seen)) {
            seen[fname] = 1
            for (i = 1; i <= h; i++) print header[i] > fname
        }

        print >> fname
    }'

    for f in "$outdirname"/*.vcf; do
        [ -e "$f" ] || continue
        echo "Compressing and indexing: $f"
        bgzip "$f"
        tabix -p vcf "$f.gz"
    done

done < "$curr_dir/$species_list"