#!/usr/bin/env python3
"""Prefix the names of scaffolds produced by P_RNA_scaffolder.

After a scaffolding round, records whose ID starts with ``P_RNA_scaffold`` are the
new joins; they are renamed ``<prefix>_<id>`` so the round is traceable. Untouched
scaffolds keep their name.

Usage:
    Rename_ScaffoldedSequences.py <in.fasta> <out.fasta> <prefix>   # e.g. prefix = R1
"""
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_fasta")
    parser.add_argument("out_fasta")
    parser.add_argument("prefix", help="round tag prepended to newly joined scaffolds")
    args = parser.parse_args()

    with open(args.out_fasta, "w") as out:
        for record in SeqIO.parse(args.in_fasta, "fasta"):
            name = f"{args.prefix}_{record.id}" if record.id.startswith("P_RNA_scaffold") else record.id
            print(">" + name, file=out)
            print(record.seq, file=out)


if __name__ == "__main__":
    main()
