#!/usr/bin/env python3
"""GC content in fixed non-overlapping windows, as a bedGraph-like TSV.

Output columns:  <seq_id>  <start(1-based)>  <end>  <gc_fraction>

Usage:
    Calculate_GC_content.py <in.fasta> <out.bedGraph> <window_bp>

DIFFERENCE FROM THE ORIGINAL: the original script had two bugs affecting only the
final, shorter-than-a-window slice of each sequence - it sliced ``seq[w:window]``
instead of ``seq[w:]`` and reported the *previous* window's GC value. This version
computes the trailing window correctly. Full-window rows are unchanged (including
the original's ``end = start + window - 1`` convention).
"""
import argparse

from Bio import SeqIO


def gc_fraction(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_fasta")
    parser.add_argument("out_bedgraph")
    parser.add_argument("window", type=int, help="window size in bp")
    args = parser.parse_args()

    win = args.window
    with open(args.out_bedgraph, "w") as out:
        for record in SeqIO.parse(args.in_fasta, "fasta"):
            seq = str(record.seq)
            length = len(seq)
            start = 0
            while start < length:
                if start + win < length:
                    window_seq = seq[start:start + win]
                    print(record.id, start + 1, start + win - 1, gc_fraction(window_seq),
                          sep="\t", file=out)
                    start += win
                else:
                    window_seq = seq[start:]                     # fixed: full remaining slice
                    print(record.id, start + 1, length, gc_fraction(window_seq),
                          sep="\t", file=out)
                    break


if __name__ == "__main__":
    main()
