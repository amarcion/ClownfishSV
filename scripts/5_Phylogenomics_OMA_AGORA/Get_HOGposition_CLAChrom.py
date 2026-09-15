import sys


gene_positions_file = sys.argv[1]
hogs_file = sys.argv[2]
output_file = sys.argv[3]

# Read CLA gene → chromosome mapping ---
gene_to_chr = {}
with open(gene_positions_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        chr_name = parts[0]
        gene = parts[4]
        gene_to_chr[gene] = chr_name

# Parse HOGs and find CLA gene ---
with open(hogs_file) as f, open(output_file, "w") as out:
    out.write("HOG\tCLA_gene\tCLA_chr\n")
    for line in f:
        parts = line.strip().split("\t")
        if not parts:
            continue
        hog = parts[0]
        genes = parts[1].split(" ")
        
        # find the CLA gene in this HOG
        cla_gene = next((g for g in genes if g.startswith("CLA.")), None)
        
        if cla_gene and cla_gene in gene_to_chr:
            cla_chr = gene_to_chr[cla_gene]
        else:
            cla_chr = "NA"
        
        out.write(f"{hog}\t{cla_gene or 'NA'}\t{cla_chr}\n")

print(f"Results written to {output_file}")