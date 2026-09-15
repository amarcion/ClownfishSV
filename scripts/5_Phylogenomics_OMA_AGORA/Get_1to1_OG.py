import sys 

infile_name = sys.argv[1]
outfile_name = sys.argv[1]

def parse_species(gene_id):
    """Extract species code from a gene ID.
    Handles both 'AKA.g123.t1' and 'AMPPE017281'-type identifiers."""
    return gene_id.split('.')[0]


def filter_hogs(filepath, strict_1to1=False):
    """
    Filter HOGs based on species composition and copy-number rules.

    Parameters
    ----------
    filepath : str
        Path to the HOG file (tab-separated, HOGname then gene IDs).
    strict_1to1 : bool, optional
        If True → require exactly 1 copy in *all* species_needed.
        If False → allow at least 15 species_needed with 1 copy (default).

    Returns
    -------
    valid_hogs : list of str
        List of HOG names passing the filtering rules.
    """

    species_needed = [
        "FRE", "EPH", "CLA", "AKA", "SAN", "PRD", "SEB", "POL",
        "OMA", "ALL", "LAT", "CRP", "MCC", "AKY", "LAZ", "PRC",
        "OCE", "BIA"
    ]
    species_at_least_one = ["ACH", "DTR"]
    species_not_necessary_out = ["ORENI", "OREAU"]
    species_not_necessary_in = ["AMPPE", "AMPOC"]

    # threshold changes depending on strictness
    needed_threshold = len(species_needed) if strict_1to1 else 15

    valid_hogs = []

    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            hog = parts[0]
            gene_ids = parts[1:]

            # count genes per species
            species_counts = {}
            for gid in gene_ids:
                sp = parse_species(gid)
                species_counts[sp] = species_counts.get(sp, 0) + 1

            # check needed species
            needed_present = sum(1 for s in species_needed if species_counts.get(s, 0) == 1)
            has_required = any(species_counts.get(s, 0) >= 1 for s in species_at_least_one)
            bad_in_species = any(
                s in species_counts and species_counts[s] != 1 for s in species_not_necessary_in
            )

            # apply thresholds
            if needed_present >= needed_threshold and has_required and not bad_in_species:
                valid_hogs.append(hog)

    print(f"Found {len(valid_hogs)} valid HOGs ({'strict' if strict_1to1 else 'relaxed'} 1-to-1 mode)")
    return valid_hogs



# get 1to1 OG
valid_hogs_strict = filter_hogs(infile_name, strict_1to1=True)

with open(outfile_name, "w") as outfile:
with open(infile_name, "r") as infile:
    for line in infile:
        line_split = line.rstrip().split("\t")
        if line_split[0] in valid_hogs_strict:
            print(line.rstrip(), file=outfile)
                
