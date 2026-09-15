#!/usr/bin/env python3
"""Rewrite a FASTA, reverse-complementing a chosen subset of records.

SyRI expects homologous chromosomes to be on the same strand. Some chromosomes of
the individual assemblies are assembled in the opposite orientation; those are 
flipped here before alignment.

Usage:
    Get_reverse_complement.py <in.fasta> <out.fasta> <id1,id2,...>

The third argument is a comma-separated list of record IDs to reverse-complement
(e.g. "CM132556.1,CM132559.1"). Behaviour is identical to the original.
"""
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_fasta")
    parser.add_argument("out_fasta")
    parser.add_argument("ids_to_flip",
                        help="comma-separated record IDs to reverse-complement")
    args = parser.parse_args()

    to_flip = set(args.ids_to_flip.split(","))

    with open(args.out_fasta, "w") as out:
        for record in SeqIO.parse(args.in_fasta, "fasta"):
            seq = record.seq.reverse_complement() if record.id in to_flip else record.seq
            print(f">{record.id}", file=out)
            print(seq, file=out)


if __name__ == "__main__":
    main()
