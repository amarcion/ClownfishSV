#!/usr/bin/env python3
"""Write the length of every sequence in a FASTA (``<id>\\t<length>`` per line)."""
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_fasta")
    parser.add_argument("out_tsv")
    args = parser.parse_args()

    with open(args.out_tsv, "w") as out:
        for record in SeqIO.parse(args.in_fasta, "fasta"):
            print(f"{record.id}\t{len(record.seq)}", file=out)


if __name__ == "__main__":
    main()
