#!/usr/bin/env python3
"""Report the total length of the sequences actually analysed for one species.

Not every scaffold of every assembly enters the SyRI structural-variant analysis
(only full chromosomes / longest scaffolds, chr24 excluded). To turn SV counts and
lengths into per-genome proportions we need the length of exactly that analysed
subset, i.e. the FASTA that was handed to SyRI.

Usage:
    Get_AnalyzedGenomeLength.py <analysed_genome.fasta> <species_id>

Prints one tab-separated line:  <species_id>\t<total_length_bp>
(run it in a loop over species and redirect to a summary table).
"""
import argparse

from Bio import SeqIO


def total_sequence_length(fasta_path):
    """Sum of the lengths of every record in a FASTA file."""
    return sum(len(record.seq) for record in SeqIO.parse(fasta_path, "fasta"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("analysed_genome", help="FASTA of the sequences given to SyRI for this species")
    parser.add_argument("species_id", help="short species code, printed as the first column")
    args = parser.parse_args()

    print(args.species_id, total_sequence_length(args.analysed_genome), sep="\t")


if __name__ == "__main__":
    main()
