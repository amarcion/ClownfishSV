#!/usr/bin/env python3
"""Drop scaffolds shorter than a threshold and shorten hifiasm-style names.

Records >= ``min_length`` bp are kept. If a kept record's ID contains ``_`` (the
name a scaffold gets after joining several ``ptgXXXXXX`` contigs), the ``ptg``
pieces are concatenated into a single ``ptg<digits><digits>...`` name; otherwise
the ID is unchanged. 

Usage:
    Remove_Short_Scaffolds.py <in.fasta> <out.fasta> <min_length_bp>
"""
import argparse

from Bio import SeqIO


def shorten_id(record_id):
    if "_" not in record_id:
        return record_id
    new_id = "ptg"
    for piece in record_id.split("_"):
        if piece.startswith("ptg"):
            new_id += piece[5:]          # drop the "ptg" + 2 chars, keep the numeric tail
    print(new_id)
    return new_id


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_fasta")
    parser.add_argument("out_fasta")
    parser.add_argument("min_length", type=int, help="minimum scaffold length to keep (bp)")
    args = parser.parse_args()

    with open(args.out_fasta, "w") as out:
        for record in SeqIO.parse(args.in_fasta, "fasta"):
            if len(record.seq) < args.min_length:
                continue
            print(">" + shorten_id(record.id), file=out)
            print(record.seq, file=out)


if __name__ == "__main__":
    main()
