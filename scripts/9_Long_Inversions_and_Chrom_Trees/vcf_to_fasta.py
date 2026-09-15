import sys
import csv
from collections import defaultdict
from intervaltree import IntervalTree

def parse_synteny_file(filepath):
    tree = defaultdict(IntervalTree)
    with open(filepath) as f:
        for line in f:
            if line.strip():
                chrom, start, end = line.strip().split()
                tree[chrom].addi(int(start), int(end))
    return tree

def position_in_synteny(tree, chrom, pos):
    return chrom in tree and tree[chrom].overlaps(pos)

def parse_vcf(vcf_file, chrom, start=None, end=None):
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split('\t')
                species_list = header[9:]
                break

        snps = []
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split('\t')
            vcf_chrom, pos = cols[0], int(cols[1])
            if vcf_chrom != chrom:
                continue
            # If start/end provided, filter by region; otherwise take whole chromosome
            if start is not None and end is not None:
                if not (start <= pos <= end):
                    continue

            ref, alt = cols[3], cols[4]
            genotypes = cols[9:]
            snps.append((pos, ref, alt, genotypes))

    return snps, species_list

def build_alignment(snps, species_list, synteny_trees, chrom):
    alignment = {sp: [] for sp in species_list}
    missing_counts = {sp: 0 for sp in species_list}
    alignment["BIARef"] = []
    # BIARef will never have missing data
    missing_counts["BIARef"] = 0

    for pos, ref, alt, genotypes in snps:

        # Add REF for BIARef always
        alignment["BIARef"].append(ref)

        for sp_idx, sp in enumerate(species_list):
            tree = synteny_trees[sp]
            in_synteny = position_in_synteny(tree, chrom, pos)
            gt = genotypes[sp_idx].split(':')[0]

            if not in_synteny:
                alignment[sp].append('N')
                missing_counts[sp] += 1
            elif gt == '1':
                alignment[sp].append(alt)
            else:
                alignment[sp].append(ref)

    return alignment, missing_counts


def filter_invariant_sites(alignment, species_list):
    """Remove sites where all non-N bases are identical (invariant)."""
    all_species = ["BIARef"] + species_list
    n_sites = len(alignment["BIARef"])

    keep = []
    for i in range(n_sites):
        bases = set()
        for sp in all_species:
            b = alignment[sp][i].upper()
            if b != 'N':
                bases.add(b)
        if len(bases) >= 2:
            keep.append(i)

    n_removed = n_sites - len(keep)
    filtered = {sp: [alignment[sp][i] for i in keep] for sp in all_species}
    return filtered, n_removed


def filter_missing_data(alignment, species_list, max_missing=0.6):
    """Remove sites where more than max_missing fraction of samples have N."""
    # BIARef is the reference and never has missing data, so only count real samples
    n_samples = len(species_list)
    n_sites = len(alignment["BIARef"])
    all_species = ["BIARef"] + species_list

    keep = []
    for i in range(n_sites):
        n_missing = sum(1 for sp in species_list if alignment[sp][i].upper() == 'N')
        if n_missing / n_samples <= max_missing:
            keep.append(i)

    n_removed = n_sites - len(keep)
    filtered = {sp: [alignment[sp][i] for i in keep] for sp in all_species}
    return filtered, n_removed


def write_fasta(alignment, output_file):
    with open(output_file, 'w') as f:
        for sp, seq in alignment.items():
            f.write(f">{sp}\n{''.join(seq)}\n")

def main():
    # USER INPUTS
    # Usage:
    #   Whole chromosome:  script.py <chrom> <outfile> <synteny_dir> <vcf>
    #   Specific region:   script.py <chrom> <start> <end> <outfile> <synteny_dir> <vcf>

    if len(sys.argv) == 5:
        # Whole-chromosome mode: no start/end provided
        chrom                   = sys.argv[1]
        outfile_name            = sys.argv[2]
        path_to_syntenyic_regions = sys.argv[3]
        vcf_file                = sys.argv[4]
        region_start, region_end = None, None
        print(f"No region specified — extracting whole chromosome: {chrom}")

    elif len(sys.argv) == 7:
        # Region mode: start and end provided
        chrom                   = sys.argv[1]
        region_start            = int(sys.argv[2])
        region_end              = int(sys.argv[3])
        outfile_name            = sys.argv[4]
        path_to_syntenyic_regions = sys.argv[5]
        vcf_file                = sys.argv[6]
        print(f"Extracting region: {chrom}:{region_start}-{region_end}")

    else:
        print("Usage:")
        print("  Whole chromosome : script.py <chrom> <outfile> <synteny_dir> <vcf>")
        print("  Specific region  : script.py <chrom> <start> <end> <outfile> <synteny_dir> <vcf>")
        sys.exit(1)

    synteny_files = {
        "AKA": path_to_syntenyic_regions + "/AKA.BIARef_Aligned.txt",
        "AKY": path_to_syntenyic_regions + "/AKY.BIARef_Aligned.txt",
        "ALL": path_to_syntenyic_regions + "/ALL.BIARef_Aligned.txt",
        "BIA": path_to_syntenyic_regions + "/BIA.BIARef_Aligned.txt",
        "CLA": path_to_syntenyic_regions + "/CLA.BIARef_Aligned.txt",
        "CRP": path_to_syntenyic_regions + "/CRP.BIARef_Aligned.txt",
        "DTR": path_to_syntenyic_regions + "/DTR.BIARef_Aligned.txt",
        "EPH": path_to_syntenyic_regions + "/EPH.BIARef_Aligned.txt",
        "FRE": path_to_syntenyic_regions + "/FRE.BIARef_Aligned.txt",
        "LAT": path_to_syntenyic_regions + "/LAT.BIARef_Aligned.txt",
        "LAZ": path_to_syntenyic_regions + "/LAZ.BIARef_Aligned.txt",
        "MCC": path_to_syntenyic_regions + "/MCC.BIARef_Aligned.txt",
        "OMA": path_to_syntenyic_regions + "/OMA.BIARef_Aligned.txt",
        "OCE": path_to_syntenyic_regions + "/OCE.BIARef_Aligned.txt",        
        "POL": path_to_syntenyic_regions + "/POL.BIARef_Aligned.txt",
        "PRC": path_to_syntenyic_regions + "/PRC.BIARef_Aligned.txt",
        "PRD": path_to_syntenyic_regions + "/PRD.BIARef_Aligned.txt",
        "SAN": path_to_syntenyic_regions + "/SAN.BIARef_Aligned.txt",
        "SEB": path_to_syntenyic_regions + "/SEB.BIARef_Aligned.txt",
    }

    # Load synteny data
    synteny_trees = {sp: parse_synteny_file(path) for sp, path in synteny_files.items()}

    # Parse VCF
    snps, species_list = parse_vcf(vcf_file, chrom, region_start, region_end)

    # Limit to species we have synteny info for
    species_list = [sp for sp in species_list if sp in synteny_files]

    # Build alignment
    alignment, missing_counts = build_alignment(snps, species_list, synteny_trees, chrom)

    # Filter invariant sites
    alignment, n_inv_removed = filter_invariant_sites(alignment, species_list)
    after_inv = len(alignment["BIARef"])

    # Filter sites with too much missing data (>60% of samples)
    alignment, n_miss_removed = filter_missing_data(alignment, species_list, max_missing=0.6)
    after_miss = len(alignment["BIARef"])

    # Output FASTA
    write_fasta(alignment, outfile_name)

    # Print summary
    print(f"\nSite filtering summary:")
    print(f"  Total SNPs extracted     : {len(snps)}")
    print(f"  Invariant sites removed  : {n_inv_removed}")
    print(f"  High-missingness removed : {n_miss_removed}  (>60% samples as N)")
    print(f"  Final sites kept         : {after_miss}")
    print("\nMissing data counts (N per species, before site filtering):")
    for sp, count in missing_counts.items():
        print(f"  {sp}: {count} N's")

if __name__ == "__main__":
    main()