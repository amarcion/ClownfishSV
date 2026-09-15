#!/usr/bin/env python3
"""Split a merged SyRI short-variant (ShV) VCF into SNPs, long indels and short indels.

From one per-chromosome merged ShV VCF this writes, for chromosome ``<chr>``:

    <prefix>.<chr>.SNP.vcf           biallelic SNPs only (multiallelic SNPs dropped)
    <prefix>.<chr>.INDELS.vcf        indels >= 50 bp   ("long" -> treated as SV)
    <prefix>.<chr>.ShortINDELS.vcf   indels < 50 bp
    <prefix>.<chr>.INDELS.Summary.txt   one line per long indel:
        chrom  pos  id  length  comma-separated carrier species

Indel length is taken as the max of two estimates parsed from the INFO field
(``END`` - ``POS`` and ``EndB`` - ``StartB``), exactly as in the original.
A sample "carries" the variant when its GT contains a "1".

Usage:
    Filter_ShV_files.py <merged_ShV.vcf> <chrom_label> <out_prefix>
"""
import argparse

MIN_LONG_INDEL = 50


def _int_between(text, after, before):
    start = text.find(after) + len(after)
    end = text.find(before)
    return int(text[start:end])


def indel_length_from_end(pos, info):
    """abs(END - POS) parsed from INFO (``END=...;EndB=``)."""
    return abs(_int_between(info, "END=", ";EndB=") - int(pos))


def indel_length_from_query(info):
    """EndB - StartB parsed from INFO (``StartB=...;VarType=`` / ``EndB=...;Parent``)."""
    start = _int_between(info, "StartB=", ";VarType=")
    stop = _int_between(info, "EndB=", ";Parent")
    return stop - start


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("merged_shv_vcf")
    parser.add_argument("chrom_label", help="used only to name the output files")
    parser.add_argument("out_prefix")
    args = parser.parse_args()

    header = []
    long_indels, short_indels, summary = [], [], []
    n_biallelic_snp = n_multiallelic_snp = 0

    snp_path = f"{args.out_prefix}.{args.chrom_label}.SNP.vcf"
    with open(args.merged_shv_vcf) as infile, open(snp_path, "w") as out_snp:
        for line in infile:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
                print(line, file=out_snp)
                continue

            fields = line.split("\t")
            var_id, alt, info = fields[2], fields[4], fields[7]

            if var_id.startswith("SNP"):
                if "," in alt:                     # multiallelic -> drop
                    n_multiallelic_snp += 1
                else:
                    n_biallelic_snp += 1
                    print(line, file=out_snp)
                continue

            # --- indel ---
            length = max(indel_length_from_end(fields[1], info),
                         indel_length_from_query(info))
            if length < MIN_LONG_INDEL:
                short_indels.append(line)
                continue

            long_indels.append(line)

            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0
            sample_names = header[-1].split("\t")[9:]
            carriers = [sample_names[i]
                        for i, sample in enumerate(fields[9:])
                        if "1" in sample.split(":")[gt_idx]]
            summary.append("\t".join(fields[0:3]) + f"\t{length}\t" + ",".join(carriers))

    _write(f"{args.out_prefix}.{args.chrom_label}.INDELS.vcf", header + long_indels)
    _write(f"{args.out_prefix}.{args.chrom_label}.ShortINDELS.vcf", header + short_indels)
    _write(f"{args.out_prefix}.{args.chrom_label}.INDELS.Summary.txt", summary)

    print("Biallelic SNPS:", n_biallelic_snp)
    print("SNPs not biallelic:", n_multiallelic_snp)
    print("Small INDELS:", len(short_indels))
    print("LONG INDELS:", len(long_indels))


def _write(path, lines):
    with open(path, "w") as handle:
        for line in lines:
            print(line, file=handle)


if __name__ == "__main__":
    main()
