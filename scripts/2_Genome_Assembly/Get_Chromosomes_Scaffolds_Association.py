import pandas as pd
import re
from collections import defaultdict
import sys

file_name = sys.argv[1] 
output_file_name = sys.argv[2]

# Read the file while skipping dynamic header lines
with open(file_name, 'r') as f:
    lines = f.readlines()
    start_idx = next(i for i, line in enumerate(lines) if line.startswith("[S1]")) + 1  # Find data start

# Read the file with proper whitespace handling
df = pd.read_csv(file_name, sep=r'\s+', skiprows=start_idx, header=None,
                 names=["S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "COV_R", "COV_Q", "CHROM", "SCAFFOLD"],
                 engine='python')  # Use 'python' engine for flexible parsing

# Remove non-numeric characters from numeric columns
def clean_numeric(value):
    if isinstance(value, str):
        cleaned_value = re.sub(r'[^0-9]', '', value)  # Remove non-digit characters
        return int(cleaned_value) if cleaned_value else None  # Convert to int or None
    return value

numeric_columns = ["S1", "E1", "S2", "E2", "LEN1", "LEN2"]
df[numeric_columns] = df[numeric_columns].map(clean_numeric)

df.dropna(inplace=True)  # Remove invalid rows
df = df.astype({"S1": int, "E1": int, "S2": int, "E2": int})

# Compute coverage per row
df["COVERAGE"] = df["E2"] - df["S2"]

# Dictionary to store cumulative covered length per (scaffold, chromosome)
cumulative_coverage = defaultdict(lambda: defaultdict(int))

# Aggregate coverage for each (chromosome, scaffold) pair
for _, row in df.iterrows():
    cumulative_coverage[row["SCAFFOLD"]][row["CHROM"]] += row["LEN2"]  # Add up covered lengths

# Determine the major chromosome with the highest cumulative length per scaffold
major_chromosome = {}
covered_lengths = {}
for scaffold, chrom_dict in cumulative_coverage.items():
    if chrom_dict:  # Ensure there is valid data
        major_chr = max(chrom_dict, key=chrom_dict.get)  # Chromosome with highest cumulative coverage
        major_chromosome[scaffold] = major_chr
        covered_lengths[scaffold] = chrom_dict[major_chr]  # Get total covered length

# Convert to DataFrame for better visualization
result_df = pd.DataFrame({
    "Scaffold": list(major_chromosome.keys()),
    "Major_Chromosome": list(major_chromosome.values()),
    "Total_Covered_Length": list(covered_lengths.values())
})
print(result_df)


# Save results to a file
result_df.to_csv(output_file_name,  sep='\t', index=False)
