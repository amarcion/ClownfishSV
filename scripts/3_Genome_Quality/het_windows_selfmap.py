import cyvcf2
import pandas as pd
import sys
from collections import defaultdict

vcf_file = sys.argv[1]      # e.g. AKA.vcf.gz
species = sys.argv[2]        # e.g. AKA
window_size = int(sys.argv[3])  # e.g. 100000

vcf = cyvcf2.VCF(vcf_file)
results = []
windows = defaultdict(lambda: [0, 0, 0])  # [het, total, depth_sum]

for variant in vcf:
    # Skip non-SNPs and missing genotypes
    if not variant.is_snp:
        continue
    chrom = variant.CHROM
    pos = variant.POS
    window = (pos - 1) // window_size

    gt = variant.genotypes[0]  # single individual
    a1, a2 = gt[0], gt[1]
    if a1 < 0 or a2 < 0:      # skip missing
        continue

    windows[(chrom, window)][1] += 1          # total SNPs
    if a1 != a2:
        windows[(chrom, window)][0] += 1      # het SNPs

    # Get depth if available
    try:
        dp = variant.format('DP')[0][0]
        windows[(chrom, window)][2] += dp
    except:
        pass

for (chrom, window), (het, total, depth_sum) in windows.items():
    if total > 0:
        start = window * window_size + 1
        end = (window + 1) * window_size
        results.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "species": species,
            "het_count": het,
            "total_snps": total,
            "het_fraction": het / total,
            "het_per_bp": het / window_size,
            "mean_depth": depth_sum / total if total > 0 else 0
        })

df = pd.DataFrame(results)
df.to_csv(f"{species}_het_windows.csv", index=False)
print(f"{species}: mean het = {df['het_per_bp'].mean():.6f}")