#!/usr/bin/env python3
# summarize_self_mapping_SVs.py

import pysam
import pandas as pd
import sys
from collections import defaultdict

def parse_sv_vcf(vcf_file, species_name):
    """Parse Sniffles2 VCF and extract SV info"""
    
    vcf = pysam.VariantFile(vcf_file)
    
    records = []
    
    for rec in vcf:
        # Skip if not PASS
        if rec.filter.keys() and \
           "PASS" not in rec.filter.keys():
            continue
            
        # Get SV type
        svtype = rec.info.get("SVTYPE", "UNKNOWN")
        
        # Get SV length
        svlen = rec.info.get("SVLEN", None)
        if svlen is None:
            end = rec.info.get("END", rec.pos)
            svlen = abs(end - rec.pos)
        else:
            svlen = abs(svlen)
        
        # Get genotype for this sample
        sample = list(rec.samples.keys())[0]
        gt = rec.samples[sample]["GT"]
        
        # Convert genotype tuple to string
        if gt is None or None in gt:
            gt_str = "./."
        elif gt == (0, 0):
            gt_str = "0/0"
        elif gt in [(0, 1), (1, 0)]:
            gt_str = "0/1"
        elif gt == (1, 1):
            gt_str = "1/1"
        else:
            gt_str = "/".join(map(str, gt))
        
        records.append({
            "species":  species_name,
            "chrom":    rec.chrom,
            "pos":      rec.pos,
            "svtype":   svtype,
            "svlen":    svlen,
            "genotype": gt_str
        })
    
    return pd.DataFrame(records)


def summarize_svs(df):
    """Create summary table per SV type and genotype"""
    
    summary = df.groupby(
        ["species", "svtype", "genotype"]
    ).agg(
        count    = ("svlen", "count"),
        total_len_mb = ("svlen", lambda x: 
                        round(x.sum() / 1e6, 3)),
        mean_len_kb  = ("svlen", lambda x: 
                        round(x.mean() / 1e3, 2)),
        median_len_kb = ("svlen", lambda x: 
                         round(x.median() / 1e3, 2))
    ).reset_index()
    
    return summary


def main():
    import glob
    import os
    
    os.chdir(sys.argv[1])
    
    # Find all VCF files
    # Assumes naming: SPECIES.vcf or SPECIES_sniffles.vcf
    vcf_files = glob.glob("*.vcf") + \
                glob.glob("*.vcf.gz")
    
    all_data = []
    
    for vcf_file in sorted(vcf_files):
        species = os.path.basename(vcf_file)\
                    .replace(".vcf.gz", "")\
                    .replace(".vcf", "")\
                    .replace("_sniffles", "")
        
        print(f"Processing {species}...")
        df = parse_sv_vcf(vcf_file, species)
        all_data.append(df)
    
    # Combine all species
    all_df = pd.concat(all_data, ignore_index=True)
    
    # Save raw data
    all_df.to_csv("all_SVs_raw.csv", index=False)
    
    # Summary per species, svtype, genotype
    summary = summarize_svs(all_df)
    summary.to_csv("SV_summary_by_genotype.csv", 
                   index=False)
    
    # Print formatted summary
    print("\n" + "="*70)
    print("SV SUMMARY BY TYPE AND GENOTYPE")
    print("="*70)
    
    for species in summary["species"].unique():
        print(f"\n{'='*50}")
        print(f"Species: {species}")
        print(f"{'='*50}")
        
        sp_data = summary[
            summary["species"] == species
        ].sort_values(["svtype", "genotype"])
        
        print(sp_data.to_string(index=False))
    
    # Overall summary collapsed across species
    print("\n" + "="*70)
    print("OVERALL SUMMARY (collapsed across species)")
    print("="*70)
    
    overall = all_df.groupby(
        ["svtype", "genotype"]
    ).agg(
        count        = ("svlen", "count"),
        total_len_mb = ("svlen", lambda x: 
                        round(x.sum() / 1e6, 3)),
        mean_len_kb  = ("svlen", lambda x: 
                        round(x.mean() / 1e3, 2))
    ).reset_index()
    
    print(overall.to_string(index=False))
    overall.to_csv("SV_overall_summary.csv", index=False)
    
    # Heterozygosity summary
    print("\n" + "="*70)
    print("HETEROZYGOUS SVs SUMMARY (0/1 genotypes)")
    print("Potential missed heterozygous regions")
    print("="*70)
    
    het_svs = all_df[all_df["genotype"] == "0/1"]
    
    het_summary = het_svs.groupby(
        ["species", "svtype"]
    ).agg(
        het_count    = ("svlen", "count"),
        het_len_mb   = ("svlen", lambda x: 
                        round(x.sum() / 1e6, 3))
    ).reset_index()
    
    print(het_summary.to_string(index=False))
    het_summary.to_csv("heterozygous_SVs_summary.csv",
                       index=False)

if __name__ == "__main__":
    main()