#!/usr/bin/env bash
#!/bin/sh

set -euo pipefail

usage() { echo "Usage: bash $(basename "$0") <output> <contig>" >&2; exit 1; }
[ "$#" -eq 2 ] || usage

output=$1
contig=$2

pathToSoftware=/scratch/amarcion/Scaffolding_PE/P_RNA_scaffolder


module load gcc perl bwa

sort -k1,1 -k3,3nr $output/reliable.connections > $output/sort.reliable.connection
"$pathToSoftware"/find_end_node $output/sort.reliable.connection $output/start.node
"$pathToSoftware"/select_nodes $output/start.node $output/both.nodes
echo -n "" > $output/intron.txt
perl "$pathToSoftware"/form_path.pl $output/both.nodes $output/intron.txt 100 > $output/both.path
sed 's/->/\n/g' $output/both.path |sed 's/\/r//g' |grep -v "N(" |sort -u > $output/scaffolded.fragment.id
perl "$pathToSoftware"/generate_scaffold.pl $contig $output/both.path P_RNA_scaffold_ > $output/scaffold.fasta & perl $pathToSoftware/generate_unscaffold.pl $contig $output/scaffolded.fragment.id  > $output/unscaffold.fasta
wait
cat $output/scaffold.fasta $output/unscaffold.fasta >$output/P_RNA_scaffold.fasta
