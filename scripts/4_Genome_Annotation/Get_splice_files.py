#!/usr/bin/env python3
"""Build an OMArk/OMA ``.splice`` file from a protein FASTA.

Isoform IDs are ``<gene>.<something>``; this groups records by ``<gene>`` (the part
before the first ``.``) and writes one line per gene listing all its isoform IDs
separated by ``;``. Gene order follows first appearance in the FASTA.

Usage:
    Get_splice_files.py <proteins.fasta> <out.splice>
"""
import argparse
from collections import OrderedDict

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("proteins_fasta")
    parser.add_argument("out_splice")
    args = parser.parse_args()

    isoforms_by_gene = OrderedDict()
    for record in SeqIO.parse(args.proteins_fasta, "fasta"):
        gene = record.id[:record.id.find(".")]
        isoforms_by_gene.setdefault(gene, []).append(record.id)

    with open(args.out_splice, "w") as out:
        for isoforms in isoforms_by_gene.values():
            print(";".join(isoforms), file=out)


if __name__ == "__main__":
    main()
