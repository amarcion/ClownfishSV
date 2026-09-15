#!/usr/bin/env python3
"""
Extract coding genes overlapping inversion breakpoints and genes contained
within an inversion, for a given species/scaffold/interval, and join with
functional annotation (GeneName, GeneDescription, GO terms).
"""
import re
import sys
import csv
import json
import argparse
from collections import defaultdict

ANNOT_DIR = "/Users/amarcion/Documents/PacBio_Sequencing_Assembly/10_Additional/16_Long_Inversions_Final/A_Long_Inv_BreakPoints_perSpecies/Annotations"

GTF_MAP = {
    "AKA": "AKA.Annotation_v2.gtf",
    "AKY": "AKY.Annotation_v2.gtf",
    "ALL": "ALL.Annotation_v2.gtf",
    "BIA": "BIA.Annotation_v2.gtf",
    "CRP": "CRP.Annotation_v2.gtf",
    "EPH": "EPH.Annotation_v2.gtf",
    "FRE": "FRE.Annotation_v2.gtf",
    "LAT": "LAT.Annotation_v2.gtf",
    "LAZ": "LAZ.Annotation_v2.gtf",
    "MCC": "MCC.Annotation_v2.gtf",
    "OMA": "OMA.Annotation_v2.gtf",
    "PRC": "PRC.Annotation_v2.gtf",
    "PRD": "PRD.Annotation_v2.gtf",
    "POL": "POL.Annotation_v2.gtf",
    "SAN": "SAN.Annotation_v2.gtf",
    "SEB": "SEB.Annotation_v2.gtf",
    "CLA": "A.clarkii_Final_Annotation.gtf",
}

FUNC_MAP = {sp: f"{sp}.Annotation_v2.FunctionalAnnotation.txt" for sp in GTF_MAP if sp != "CLA"}
FUNC_MAP["CLA"] = "CLA.Annotation_v2.FunctionalAnnotation.txt"

ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


def load_genes(species, scaffold):
    """Return dict gene_id -> [start, end, strand, seqname]"""
    path = f"{ANNOT_DIR}/{GTF_MAP[species]}"
    genes = {}
    if species == "CLA":
        # No explicit 'gene' feature; aggregate from 'transcript' lines by gene_id
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                f_ = line.rstrip("\n").split("\t")
                if len(f_) < 9:
                    continue
                seqname, source, feature, start, end, score, strand, frame, attr = f_[:9]
                if seqname != scaffold or feature != "transcript":
                    continue
                a = parse_attrs(attr)
                gid = a.get("gene_id")
                if gid is None:
                    continue
                start, end = int(start), int(end)
                if gid not in genes:
                    genes[gid] = [start, end, strand, seqname]
                else:
                    genes[gid][0] = min(genes[gid][0], start)
                    genes[gid][1] = max(genes[gid][1], end)
    else:
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                f_ = line.rstrip("\n").split("\t")
                if len(f_) < 9:
                    continue
                seqname, source, feature, start, end, score, strand, frame, attr = f_[:9]
                if seqname != scaffold or feature != "gene":
                    continue
                # gene line attr is bare gene_id, e.g. "g5795"
                gid = attr.strip().strip(";").strip()
                start, end = int(start), int(end)
                genes[gid] = [start, end, strand, seqname]
    return genes


def load_functional_annotation(species):
    """Return dict gene_id -> best transcript row (prefers IsLongestIsoform == Yes)
    and a dict gene_id -> list of all GO terms across all its isoforms (union)."""
    path = f"{ANNOT_DIR}/{FUNC_MAP[species]}"
    best = {}
    go_union = defaultdict(set)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tid = row["TranscriptID"]
            if species == "CLA":
                gid = tid  # TranscriptID IS the gene id for CLA
            else:
                gid = tid.rsplit(".t", 1)[0] if ".t" in tid else tid
            gos = row.get("GeneOntologies", "")
            if gos and gos != "NA":
                for g in gos.split(","):
                    go_union[gid].add(g.strip())
            is_longest = row.get("IsLongestIsoform", "") == "Yes"
            if gid not in best or is_longest:
                best[gid] = row
    return best, go_union


def classify_genes(genes, start_bp, stop_bp):
    """Return (breakpoint_genes, inside_genes) as dict gene_id -> info"""
    lo, hi = min(start_bp, stop_bp), max(start_bp, stop_bp)
    bp_genes = {}
    inside_genes = {}
    for gid, (gstart, gend, strand, seqname) in genes.items():
        overlaps_start = gstart <= start_bp <= gend
        overlaps_stop = gstart <= stop_bp <= gend
        if overlaps_start or overlaps_stop:
            which = []
            if overlaps_start:
                which.append("Start")
            if overlaps_stop:
                which.append("Stop")
            bp_genes[gid] = {
                "start": gstart, "end": gend, "strand": strand,
                "breakpoint": ";".join(which)
            }
        elif gstart >= lo and gend <= hi:
            inside_genes[gid] = {"start": gstart, "end": gend, "strand": strand}
    return bp_genes, inside_genes


def nearest_flanks(genes, bp_pos):
    """Nearest gene upstream (end < bp) and downstream (start > bp) of a breakpoint position
    that does not directly overlap it."""
    upstream = None  # (distance, gid)
    downstream = None
    for gid, (gstart, gend, strand, seqname) in genes.items():
        if gstart <= bp_pos <= gend:
            continue  # overlapping, handled elsewhere
        if gend < bp_pos:
            d = bp_pos - gend
            if upstream is None or d < upstream[0]:
                upstream = (d, gid)
        elif gstart > bp_pos:
            d = gstart - bp_pos
            if downstream is None or d < downstream[0]:
                downstream = (d, gid)
    return upstream, downstream


def annotate(gene_ids, best, go_union, go_lookup):
    rows = []
    for gid in gene_ids:
        row = best.get(gid)
        gos = sorted(go_union.get(gid, []))
        go_names = []
        for g in gos:
            info = go_lookup.get(g)
            if info:
                ns = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}.get(info["namespace"], "?")
                go_names.append(f"{g}|{ns}|{info['name']}")
            else:
                go_names.append(g)
        if row:
            rows.append({
                "GeneID_local": gid,
                "GeneName": row.get("GeneName", "NA"),
                "GeneDescription": row.get("GeneDescription", "NA"),
                "NCBI_GeneID": row.get("GeneID", "NA"),
                "GO_terms": ";".join(go_names) if go_names else "NA",
            })
        else:
            rows.append({
                "GeneID_local": gid, "GeneName": "NA", "GeneDescription": "NA",
                "NCBI_GeneID": "NA", "GO_terms": "NA"
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("scaffold")
    ap.add_argument("start", type=int)
    ap.add_argument("stop", type=int)
    ap.add_argument("outprefix")
    args = ap.parse_args()

    with open("go_lookup.json") as f:
        go_lookup = json.load(f)

    genes = load_genes(args.species, args.scaffold)
    best, go_union = load_functional_annotation(args.species)
    bp_genes, inside_genes = classify_genes(genes, args.start, args.stop)

    bp_rows = annotate(bp_genes.keys(), best, go_union, go_lookup)
    for r, gid in zip(bp_rows, bp_genes.keys()):
        r.update(bp_genes[gid])

    inside_rows = annotate(inside_genes.keys(), best, go_union, go_lookup)
    for r, gid in zip(inside_rows, inside_genes.keys()):
        r.update(inside_genes[gid])

    fields_bp = ["GeneID_local", "start", "end", "strand", "breakpoint", "GeneName", "GeneDescription", "NCBI_GeneID", "GO_terms"]
    fields_in = ["GeneID_local", "start", "end", "strand", "GeneName", "GeneDescription", "NCBI_GeneID", "GO_terms"]

    with open(f"{args.outprefix}.breakpoint_genes.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields_bp, delimiter="\t")
        w.writeheader()
        for r in sorted(bp_rows, key=lambda x: x["start"]):
            w.writerow({k: r.get(k) for k in fields_bp})

    with open(f"{args.outprefix}.inversion_genes.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields_in, delimiter="\t")
        w.writeheader()
        for r in sorted(inside_rows, key=lambda x: x["start"]):
            w.writerow({k: r.get(k) for k in fields_in})

    # nearest flanking genes for breakpoints with no direct overlap
    with open(f"{args.outprefix}.breakpoint_flanks.tsv", "w", newline="") as f:
        fields_flank = ["breakpoint_name", "position", "side", "distance_bp", "GeneID_local", "GeneName", "GeneDescription", "NCBI_GeneID"]
        w = csv.DictWriter(f, fieldnames=fields_flank, delimiter="\t")
        w.writeheader()
        for bp_name, bp_pos in [("Start", args.start), ("Stop", args.stop)]:
            if bp_pos in (args.start, args.stop) and any(bp_name in v["breakpoint"] for v in bp_genes.values()):
                continue  # directly overlapped, skip flank reporting
            up, down = nearest_flanks(genes, bp_pos)
            for side, val in [("upstream", up), ("downstream", down)]:
                if val is None:
                    continue
                d, gid = val
                info = annotate([gid], best, go_union, go_lookup)[0]
                w.writerow({
                    "breakpoint_name": bp_name, "position": bp_pos, "side": side,
                    "distance_bp": d, "GeneID_local": gid,
                    "GeneName": info["GeneName"], "GeneDescription": info["GeneDescription"],
                    "NCBI_GeneID": info["NCBI_GeneID"],
                })

    print(f"{args.species}\t{args.scaffold}\t{args.start}\t{args.stop}\t"
          f"n_total_genes_on_scaffold={len(genes)}\tn_breakpoint_genes={len(bp_genes)}\tn_inside_genes={len(inside_genes)}")


if __name__ == "__main__":
    main()
