#!/usr/bin/env python3
"""Count and sum the length of "long" indels (>= 50 bp) in a SyRI ``*.syri.out`` file.

SyRI reports insertions and deletions among the *short* variants regardless of
size. For the manuscript, indels >= 50 bp are treated as structural variants, so
this script totals them separately from the small ones.

Column layout of ``*.syri.out`` (0-based), for the fields used here:
    1,2   = reference start / end
    6,7   = query start / end
    10    = variant type tag (e.g. "INS", "DEL", "SYNAL", ...)

Length is taken on the side that carries the sequence:
    INS -> query span   (col 7 - col 6)
    DEL -> reference span (col 2 - col 1)

Usage:
    Sum_long_INDELs.py <sample.syri.out> <species_id>

Output (tab separated), identical to the original:
    species  sum_len_INS  count_INS  sum_len_DEL  count_DEL
"""
import argparse

MIN_LEN = 50


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("syri_out", help="a SyRI *.syri.out file")
    parser.add_argument("species_id", help="short species code, printed as the first column")
    args = parser.parse_args()

    sum_len = {"INS": 0, "DEL": 0}
    count = {"INS": 0, "DEL": 0}

    with open(args.syri_out) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            vtype = fields[10]
            if vtype == "INS":
                length = int(fields[7]) - int(fields[6])
            elif vtype == "DEL":
                length = int(fields[2]) - int(fields[1])
            else:
                continue
            if length >= MIN_LEN:
                sum_len[vtype] += length
                count[vtype] += 1

    print(args.species_id, sum_len["INS"], count["INS"],
          sum_len["DEL"], count["DEL"], sep="\t")


if __name__ == "__main__":
    main()
