#!/usr/bin/env python3
"""Keep one representative scaffold per chromosome (the ``scfNN.1`` records).

Assembly scaffolds are named ``scf<chrom>.<part>`` (part 1 = longest scaffold of
that chromosome, see Section 2.8). For SyRI we want a single sequence per
chromosome, so this keeps only ``.1`` scaffolds of chromosomes 1..24 and drops the
rest. Chromosomes flagged in the reverse-complement table are flipped so every
assembly has the same orientation.

Usage:
    Extract_longest_chromosomes.py <genome.fasta> <species_id> <out_dir> [reverse_complement_table]

Writes ``<out_dir>/<species_id>.LongestChroms.fa``.

reverse_complement_table (optional): TSV, ``<species_id>\t<scfID,scfID,...>``; the
listed scaffolds of that species are reverse-complemented. 
"""
import argparse
import os

from Bio import SeqIO

MAX_CHROM = 24


def load_reverse_complement_table(path):
    """species_id -> set of scaffold IDs to reverse-complement."""
    table = {}
    if not path:
        return table
    if not os.path.isfile(path):
        print(f"[warn] reverse-complement table not found ({path}); proceeding without it")
        return table
    with open(path) as handle:
        for line in handle:
            species, scaffolds = line.rstrip("\n").split("\t")
            table[species] = set(scaffolds.split(","))
    return table


def split_scaffold_id(scaffold_id):
    """'scf03.1' -> (chromosome_number, part_number) = (3, 1)."""
    dot = scaffold_id.find(".")
    chrom = int(scaffold_id[3:dot])
    part = int(scaffold_id[dot + 1:])
    return chrom, part


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("genome_fasta")
    parser.add_argument("species_id")
    parser.add_argument("out_dir")
    parser.add_argument("reverse_complement_table", nargs="?", default="",
                        help="optional TSV: species_id <TAB> comma-separated scaffold IDs to flip")
    args = parser.parse_args()

    flip_by_species = load_reverse_complement_table(args.reverse_complement_table)
    to_flip = flip_by_species.get(args.species_id, set())

    out_path = os.path.join(args.out_dir, f"{args.species_id}.LongestChroms.fa")
    with open(out_path, "w") as out:
        for record in SeqIO.parse(args.genome_fasta, "fasta"):
            chrom, part = split_scaffold_id(record.id)
            if part != 1 or chrom > MAX_CHROM:
                continue
            seq = record.seq.reverse_complement() if record.id in to_flip else record.seq
            print(f">{record.id}", file=out)
            print(seq, file=out)


if __name__ == "__main__":
    main()
