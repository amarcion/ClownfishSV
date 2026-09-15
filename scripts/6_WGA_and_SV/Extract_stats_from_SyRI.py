#!/usr/bin/env python3
"""Flatten a SyRI ``*.syri.summary`` file into one tab-separated row per species.

SyRI writes a fixed-layout summary table (one variant class per line, columns
"count" and "length/bp"). This script pulls the numbers we report in the paper into
a single line so many species can be concatenated into one table.

Usage:
    Extract_stats_from_SyRI.py <sample.syri.summary> <species_id>

Output columns (tab separated), matching the original script exactly:
    species  Synteny_count Synteny_length  Inversion_count Inversion_length
    Transloc_count Transloc_length  Duplication_count Duplication_length
    NotAligned_count NotAligned_length  SNP
    Insertions_count Insertions_length  Deletion_count Deletion_length

NOTE: the row/column indices below are tied to the SyRI v1.6.x summary layout
(the same magic numbers as the original). If SyRI changes its summary format they
must be updated.
"""
import argparse

# (row index in the summary file, column index within that row)
#   column  1 = the "count" field ; column -1 = the last field = "length"
FIELDS = [
    ("Synteny_count",     (2, 1)),
    ("Synteny_length",    (2, -1)),
    ("Inversion_count",   (3, 1)),
    ("Inversion_length",  (3, -1)),
    ("Transloc_count",    (4, 1)),
    ("Transloc_length",   (4, -1)),
    ("Duplication_count", (6, 1)),
    ("Duplication_length",(6, -1)),
    ("NotAligned_count",  (8, 1)),
    ("NotAligned_length", (8, -1)),
    ("SNP",               (13, -1)),
    ("Insertions_count",  (14, 1)),
    ("Insertions_length", (14, -1)),
    ("Deletion_count",    (15, 1)),
    ("Deletion_length",   (15, -2)),  # deletions: length is the second-to-last field
]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("syri_summary", help="a SyRI *.syri.summary file")
    parser.add_argument("species_id", help="short species code, printed as the first column")
    args = parser.parse_args()

    with open(args.syri_summary) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]

    values = [rows[row][col] for _, (row, col) in FIELDS]
    print(args.species_id, *values, sep="\t")


if __name__ == "__main__":
    main()
