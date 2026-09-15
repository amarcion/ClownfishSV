#!/usr/bin/env python3
"""Split a reformatted SyRI VCF (from Reformat_SyRI_vcf.py) into four files.

    <prefix>.SR.vcf    structural rearrangements: INV, TRANS, INVTR, DUP, INVDP
    <prefix>.SYN.vcf   syntenic regions (VarType=SR + "SYN"; not real SVs)
    <prefix>.ShV.vcf   short variants: SNP / INS / DEL  (VarType=ShV)
    <prefix>.AL.vcf    all alignment blocks (VarType=.)
    <prefix>.summary   record counts (see keys below)

Classification uses only substrings of the INFO field, exactly as the original:
    "VarType=ShV"  -> short variant (sub-classified by "Parent=SYN/INV/INVTR/TRANS")
    "VarType=.;"   -> alignment block
    "VarType=SR"   -> "SYN" in line -> syntenic, else structural (by "<INV>" etc.)

Usage:
    Split_vcf_SyRI.py -i <reformatted.vcf> -o <out_prefix>
"""
import argparse

COUNT_KEYS = [
    "c_ShV_total", "c_ShV_syn", "c_ShV_INVTR", "c_ShV_IN", "c_ShV_TRANS",
    "c_align_lines", "c_syn",
    "c_SR_total", "c_SR_INVDP", "c_SR_INVTR", "c_SR_DUP", "c_SR_TRANS", "c_SR_INV",
]


def split_vcf_by_variant(input_vcf, out_prefix):
    counts = {key: 0 for key in COUNT_KEYS}

    paths = {
        "ShV": f"{out_prefix}.ShV.vcf",
        "AL": f"{out_prefix}.AL.vcf",
        "SYN": f"{out_prefix}.SYN.vcf",
        "SR": f"{out_prefix}.SR.vcf",
    }
    handles = {name: open(path, "w") for name, path in paths.items()}
    try:
        with open(input_vcf) as infile:
            for raw in infile:
                line = raw.rstrip("\n")

                if line.startswith("#"):                       # header -> every file
                    for handle in handles.values():
                        print(line, file=handle)
                    continue

                if "VarType=ShV" in line:
                    print(line, file=handles["ShV"])
                    counts["c_ShV_total"] += 1
                    if "Parent=SYN" in line:
                        counts["c_ShV_syn"] += 1
                    elif "Parent=INVTR" in line:
                        counts["c_ShV_INVTR"] += 1
                    elif "Parent=INV" in line:
                        counts["c_ShV_IN"] += 1
                    elif "Parent=TRANS" in line:
                        counts["c_ShV_TRANS"] += 1
                    else:
                        print(line)                            # unexpected short-variant parent

                elif "VarType=.;" in line:
                    print(line, file=handles["AL"])
                    counts["c_align_lines"] += 1

                elif "VarType=SR" in line:
                    if "SYN" in line:
                        counts["c_syn"] += 1
                        print(line, file=handles["SYN"])
                    else:
                        counts["c_SR_total"] += 1
                        print(line, file=handles["SR"])
                        for tag, key in (("<INVDP>", "c_SR_INVDP"), ("<INVTR>", "c_SR_INVTR"),
                                         ("<DUP>", "c_SR_DUP"), ("<TRANS>", "c_SR_TRANS"),
                                         ("<INV>", "c_SR_INV")):
                            if tag in line:
                                counts[key] += 1
                                break
                        else:
                            print(line)                        # unexpected SR type
                else:
                    print(line)                                # unexpected record
    finally:
        for handle in handles.values():
            handle.close()
    return counts


def main():
    parser = argparse.ArgumentParser(
        description="Split a reformatted SyRI VCF into syntenic / structural / short-variant / alignment files.")
    parser.add_argument("-i", "--input_vcf", required=True, help="input VCF from Reformat_SyRI_vcf.py")
    parser.add_argument("-o", "--out_prefix", required=True, help="output prefix")
    args = parser.parse_args()

    counts = split_vcf_by_variant(args.input_vcf, args.out_prefix)
    with open(args.out_prefix + ".summary", "w") as summary:
        for key in COUNT_KEYS:
            print(key, counts[key], file=summary)


if __name__ == "__main__":
    main()
