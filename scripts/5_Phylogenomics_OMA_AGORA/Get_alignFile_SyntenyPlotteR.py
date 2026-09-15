#!/usr/bin/env python3

import sys
from collections import defaultdict

def compute_sign(col4, col9):
    return "+" if int(col4) == int(col9) else "-"


def main():
    if len(sys.argv) != 6:
        print("Usage: script.py <input_file> <output_file> <chrom_output_file> <ref_species> <target_species>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    chrom_output_file = sys.argv[3]
    ref_species = sys.argv[4]
    target_species = sys.argv[5]

    # Target species (col6, length from col7/col8)
    target_chrom_lengths = defaultdict(int)

    # Reference species (col1, length from col3)
    ref_chrom_lengths = defaultdict(int)

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if not line.strip():
                continue

            cols = line.strip().split()

            col1 = cols[0]             # reference chromosome
            col2 = cols[1]
            col3 = int(cols[2])        # reference coordinate
            col4 = cols[3]
            col6 = cols[5]             # target chromosome
            col7 = int(cols[6])        # target start
            col8 = int(cols[7])        # target end
            col9 = cols[8]

            # --- update reference chromosome length (col3) ---
            if col3 > ref_chrom_lengths[col1]:
                ref_chrom_lengths[col1] = col3

            # --- update target chromosome length (max col7/col8) ---
            target_coord = max(col7, col8)
            if target_coord > target_chrom_lengths[col6]:
                target_chrom_lengths[col6] = target_coord

            # compute sign
            new_sign = compute_sign(col4, col9)

            output_line = [
                col1, str(col2), str(col3),
                col6, str(col7), str(col8),
                new_sign,
                ref_species,
                target_species
            ]

            outfile.write("\t".join(output_line) + "\n")

    # --- write chromosome summary file ---
    with open(chrom_output_file, "w") as chrom_out:

        # 1. Target species FIRST
        for chrom in sorted(target_chrom_lengths):
            chrom_out.write(f"{chrom}\t{target_chrom_lengths[chrom]}\t{target_species}\n")

        # 2. Reference species SECOND
        for chrom in sorted(ref_chrom_lengths):
            chrom_out.write(f"{chrom}\t{ref_chrom_lengths[chrom]}\t{ref_species}\n")


if __name__ == "__main__":
    main()