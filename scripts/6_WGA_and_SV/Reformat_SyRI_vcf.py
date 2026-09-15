#!/usr/bin/env python3
"""Convert a SyRI ``*.syri.out`` table into a single-sample VCF.

Adapted from SyRI's own ``writeout``/``getVCF`` routine
(https://github.com/schneebergerlab/syri). Compared with the VCF SyRI writes
itself this version:
  * always emits a FORMAT/GT column (``GT  1``) so the file can be merged across
    species with ``vcf-merge``;
  * tags each record's INFO with ``VarType=SR`` (structural rearrangement) /
    ``ShV`` (short variant) / ``.`` so downstream splitting is trivial;
  * skips records where REF == ALT after IUPAC normalisation.

``*.syri.out`` columns (tab separated):
    achr astart aend  aseq bseq  bchr bstart bend  id parent  vartype dupclass

Usage:
    Reformat_SyRI_vcf.py <in.syri.out> <out.vcf> <work_dir> <sample_name> <reference.fasta>

``in.syri.out`` and ``out.vcf`` are both taken relative to ``work_dir`` (kept from
the original). ``reference.fasta`` is only read for the ``##contig`` lengths.
"""
import argparse
import logging
import os
from datetime import date

import pandas as pd
from Bio import SeqIO

logger = logging.getLogger("Reformat_SyRI_vcf")

SYRI_COLUMNS = ["achr", "astart", "aend", "aseq", "bseq",
                "bchr", "bstart", "bend", "id", "parent", "vartype", "dupclass"]

# IUPAC -> unambiguous, used to detect REF == ALT (same as upstream SyRI)
_IUPAC_TAB = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                           "ACGTNacgtnACCAGACAAAaccagacaaa")

_VCF_HEADER_ALTS = [
    ("SYN", "Syntenic region"), ("INV", "Inversion"), ("TRANS", "Translocation"),
    ("INVTR", "Inverted Translocation"), ("DUP", "Duplication"), ("INVDP", "Inverted Duplication"),
    ("SYNAL", "Syntenic alignment"), ("INVAL", "Inversion alignment"),
    ("TRANSAL", "Translocation alignment"), ("INVTRAL", "Inverted Translocation alignment"),
    ("DUPAL", "Duplication alignment"), ("INVDPAL", "Inverted Duplication alignment"),
    ("HDR", "Highly diverged regions"), ("INS", "Insertion in non-reference genome"),
    ("DEL", "Deletion in non-reference genome"), ("CPG", "Copy gain in non-reference genome"),
    ("CPL", "Copy loss in non-reference genome"), ("SNP", "Single nucleotide polymorphism"),
    ("TDM", "Tandem repeat"), ("NOTAL", "Not Aligned region"),
]
_VCF_HEADER_INFO = [
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position on reference genome\">",
    "##INFO=<ID=ChrB,Number=1,Type=String,Description=\"Chromosome ID on the non-reference genome\">",
    "##INFO=<ID=StartB,Number=1,Type=Integer,Description=\"Start position on non-reference genome\">",
    "##INFO=<ID=EndB,Number=1,Type=Integer,Description=\"End position on non-reference genome\">",
    "##INFO=<ID=Parent,Number=1,Type=String,Description=\"ID of the parent SR\">",
    "##INFO=<ID=VarType,Number=1,Type=String,Description=\"SR for structural arrangements, ShV for short variants, missing otherwise\">",
    "##INFO=<ID=DupType,Number=1,Type=String,Description=\"Copy gain or loss in the non-reference genome\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
]

STRUCTURAL = {"SYN", "INV", "TRANS", "INVTR", "DUP", "INVDP"}
ALIGNMENTS = {"SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDPAL"}
COPY_TANDEM = {"CPG", "CPL", "TDM"}
SHORT = {"SNP", "DEL", "INS", "HDR"}


def load_syri_out(path):
    data = pd.read_table(path, header=None, keep_default_na=False, dtype=object)
    data.columns = SYRI_COLUMNS
    data = data.loc[data["achr"] != "-"].copy()
    for col in ("astart", "aend"):
        data[col] = data[col].astype(int)
    # sort numerically by chromosome when possible, else lexicographically
    try:
        sort_chr = data["achr"].astype(int)
    except ValueError:
        logger.debug("Chromosome values are sorted lexicographically.")
        sort_chr = data["achr"]
    data = data.assign(_sortchr=sort_chr).sort_values(["_sortchr", "astart", "aend"]).drop(columns="_sortchr")
    return data


def info_field(row):
    """Build the INFO string for one syri.out row (row is a namedtuple)."""
    vt = row.vartype
    if vt in STRUCTURAL:
        return ";".join([f"END={row.aend}", f"ChrB={row.bchr}", f"StartB={row.bstart}",
                         f"EndB={row.bend}", "Parent=.", "VarType=SR", f"DupType={row.dupclass}"])
    if vt == "NOTAL":
        return ";".join([f"END={row.aend}", "ChrB=.", "StartB=.", "EndB=.",
                         "Parent=.", "VarType=.", "DupType=."])
    if vt in ALIGNMENTS:
        return ";".join([f"END={row.aend}", f"ChrB={row.bchr}", f"StartB={row.bstart}",
                         f"EndB={row.bend}", f"Parent={row.parent}", "VarType=.", "DupType=."])
    if vt in COPY_TANDEM:
        return ";".join([f"END={row.aend}", f"ChrB={row.bchr}", f"StartB={row.bstart}",
                         f"EndB={row.bend}", f"Parent={row.parent}", "VarType=ShV", "DupType=."])
    # SHORT (SNP/DEL/INS/HDR)
    return ";".join([f"END={row.aend}", f"ChrB={row.bchr}", f"StartB={row.bstart}",
                     f"EndB={row.bend}", f"Parent={row.parent}", "VarType=ShV", "DupType=."])


def write_vcf(data, out_path, sample_name, contig_lengths):
    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.3\n")
        out.write("##fileDate=" + str(date.today()).replace("-", "") + "\n")
        out.write("##source=syri\n")
        for contig, length in contig_lengths.items():
            out.write(f"##contig=<ID={contig},length={length}>\n")
        for alt_id, desc in _VCF_HEADER_ALTS:
            out.write(f'##ALT=<ID={alt_id},Description="{desc}">\n')
        for line in _VCF_HEADER_INFO:
            out.write(line + "\n")
        out.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                             "FILTER", "INFO", "FORMAT", sample_name]) + "\n")

        fmt = "GT\t1"
        for row in data.itertuples(index=False):
            ref, alt = "N", f"<{row.vartype}>"
            record = [row.achr, str(row.astart), row.id, ref, alt, ".", "PASS"]

            if row.vartype in SHORT:
                has_ref, has_alt = row.aseq != "-", row.bseq != "-"
                if has_ref != has_alt:
                    logger.error("Inconsistent annotation: need sequence for both REF and ALT or neither.")
                elif has_ref and has_alt:
                    ref_n, alt_n = row.aseq.translate(_IUPAC_TAB), row.bseq.translate(_IUPAC_TAB)
                    if ref_n.upper() == alt_n.upper():
                        continue                         # REF == ALT -> not a variant
                    record = [row.achr, str(row.astart), row.id, ref_n, alt_n, ".", "PASS"]

            record.append(info_field(row))
            record.append(fmt)
            out.write("\t".join(record) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("syri_out")
    parser.add_argument("out_vcf")
    parser.add_argument("work_dir", help="directory that both syri_out and out_vcf are relative to")
    parser.add_argument("sample_name")
    parser.add_argument("reference_fasta", help="only read for ##contig lengths")
    args = parser.parse_args()

    contig_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(args.reference_fasta, "fasta")}
    data = load_syri_out(os.path.join(args.work_dir, args.syri_out))
    write_vcf(data, os.path.join(args.work_dir, args.out_vcf), args.sample_name, contig_lengths)


if __name__ == "__main__":
    main()
