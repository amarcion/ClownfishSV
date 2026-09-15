#!/usr/bin/env python3
"""Partition a FASTA into "selected" and "everything else".

Used to pull the mitochondrial contig(s) out of an assembly: the named records go
to one file, the rest to another.

Usage:
    Extract_Sequences_from_fasta.py <in.fasta> <id1,id2,...> <selected.fasta> <remaining.fasta>
"""
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_fasta")
    parser.add_argument("ids_to_select", help="comma-separated record IDs")
    parser.add_argument("selected_fasta", help="output: the named records")
    parser.add_argument("remaining_fasta", help="output: all the other records")
    args = parser.parse_args()

    wanted = set(args.ids_to_select.split(","))

    with open(args.remaining_fasta, "w") as keep, open(args.selected_fasta, "w") as selected:
        for record in SeqIO.parse(args.in_fasta, "fasta"):
            out = selected if record.id in wanted else keep
            print(">" + record.id, file=out)
            print(record.seq, file=out)


if __name__ == "__main__":
    main()
