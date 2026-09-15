#!/usr/bin/env python3
"""
Check for repetitive elements at inversion breakpoints, and test for enrichment
of repeat content near breakpoints relative to scaffold-wide background, using
RepeatMasker .out annotations.
"""
import re
import sys
import csv
import random
import argparse
import subprocess
from bisect import bisect_left, bisect_right

REPEAT_DIR = "RepeatMask_scf"

RM_MAP = {sp: f"{sp}.Assembly_v2.fa.out" for sp in
          ["AKA", "AKY", "ALL", "BIA", "CRP", "EPH", "FRE", "LAT", "LAZ",
           "MCC", "OMA", "POL", "PRC", "PRD", "SAN", "SEB"]}

WINDOWS = [5000, 20000, 50000]
N_PERM = 2000
random.seed(42)


def load_repeats(species, scaffold):
    """Extract RepeatMasker records for one scaffold. Returns list of dicts."""
    path = f"{REPEAT_DIR}/{RM_MAP[species]}"
    # fast pre-filter with awk (files are 140+ MB, plain python parse of whole file is slow)
    cmd = ["awk", f'NR>3 && $5=="{scaffold}"', path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    records = []
    scaffold_len = 0
    for line in out.splitlines():
        f = line.split()
        if len(f) < 15:
            continue
        begin, end = int(f[5]), int(f[6])
        left = f[7].strip("()")
        strand = f[8]
        rep_name = f[9]
        rep_class = f[10]
        perc_div = float(f[1])
        records.append({
            "start": begin, "end": end, "strand": strand,
            "name": rep_name, "class": rep_class, "div": perc_div
        })
        try:
            scaffold_len = max(scaffold_len, end + int(left))
        except ValueError:
            pass
    return records, scaffold_len


def merge_intervals(intervals):
    """intervals: list of (start,end) 1-based inclusive. Returns merged sorted list."""
    if not intervals:
        return []
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def coverage_in_window(merged_starts, merged_ends, win_start, win_end):
    """bp of merged intervals overlapping [win_start,win_end], via binary search."""
    i = bisect_right(merged_ends, win_start)
    total = 0
    n = len(merged_starts)
    while i < n and merged_starts[i] <= win_end:
        s = max(merged_starts[i], win_start)
        e = min(merged_ends[i], win_end)
        if e >= s:
            total += (e - s + 1)
        i += 1
    return total


def direct_overlaps(records, pos):
    return [r for r in records if r["start"] <= pos <= r["end"]]


def window_stats(records, merged_starts, merged_ends, scaffold_len, center, half_win):
    win_start = max(1, center - half_win)
    win_end = min(scaffold_len, center + half_win)
    win_len = win_end - win_start + 1
    cov = coverage_in_window(merged_starts, merged_ends, win_start, win_end)
    frac = cov / win_len
    overlapping = [r for r in records if r["end"] >= win_start and r["start"] <= win_end]
    class_bp = {}
    for r in overlapping:
        s, e = max(r["start"], win_start), min(r["end"], win_end)
        cls = r["class"].split("/")[0]
        class_bp[cls] = class_bp.get(cls, 0) + (e - s + 1)
    return {
        "win_start": win_start, "win_end": win_end, "win_len": win_len,
        "repeat_bp": cov, "repeat_frac": frac, "n_elements": len(overlapping),
        "class_bp": class_bp,
    }


def permutation_test(merged_starts, merged_ends, scaffold_len, win_len, observed_frac, n_perm=N_PERM):
    if scaffold_len <= win_len:
        return None
    perm_fracs = []
    for _ in range(n_perm):
        start = random.randint(1, scaffold_len - win_len)
        end = start + win_len - 1
        cov = coverage_in_window(merged_starts, merged_ends, start, end)
        perm_fracs.append(cov / win_len)
    mean_p = sum(perm_fracs) / len(perm_fracs)
    var_p = sum((x - mean_p) ** 2 for x in perm_fracs) / len(perm_fracs)
    sd_p = var_p ** 0.5
    n_ge = sum(1 for x in perm_fracs if x >= observed_frac)
    pval = (n_ge + 1) / (n_perm + 1)
    fold = observed_frac / mean_p if mean_p > 0 else float("inf")
    z = (observed_frac - mean_p) / sd_p if sd_p > 0 else float("inf")
    return {"perm_mean": mean_p, "perm_sd": sd_p, "fold_enrichment": fold, "z": z, "pval": pval}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("scaffold")
    ap.add_argument("start", type=int)
    ap.add_argument("stop", type=int)
    ap.add_argument("outprefix")
    args = ap.parse_args()

    records, scaffold_len = load_repeats(args.species, args.scaffold)
    merged = merge_intervals([(r["start"], r["end"]) for r in records])
    merged_starts = [m[0] for m in merged]
    merged_ends = [m[1] for m in merged]
    total_repeat_bp = sum(e - s + 1 for s, e in merged)
    scaffold_frac = total_repeat_bp / scaffold_len if scaffold_len else 0

    # 1) direct overlap at exact breakpoint coordinate
    with open(f"{args.outprefix}.breakpoint_direct_overlap.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["breakpoint_name", "position", "repeat_name", "repeat_class", "rep_start", "rep_end", "perc_div", "strand"])
        for bp_name, pos in [("Start", args.start), ("Stop", args.stop)]:
            hits = direct_overlaps(records, pos)
            if not hits:
                w.writerow([bp_name, pos, "NONE", "", "", "", "", ""])
            for h in hits:
                w.writerow([bp_name, pos, h["name"], h["class"], h["start"], h["end"], h["div"], h["strand"]])

    # 2) window enrichment
    with open(f"{args.outprefix}.window_enrichment.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["breakpoint_name", "window_bp", "win_start", "win_end", "n_elements",
                    "repeat_bp", "repeat_frac", "scaffold_frac", "perm_mean", "perm_sd",
                    "fold_enrichment", "z", "pval", "top_classes(bp)"])
        for bp_name, pos in [("Start", args.start), ("Stop", args.stop)]:
            for half in WINDOWS:
                st = window_stats(records, merged_starts, merged_ends, scaffold_len, pos, half)
                perm = permutation_test(merged_starts, merged_ends, scaffold_len, st["win_len"], st["repeat_frac"])
                top_classes = ";".join(f"{k}:{v}" for k, v in sorted(st["class_bp"].items(), key=lambda x: -x[1])[:5])
                w.writerow([bp_name, half * 2, st["win_start"], st["win_end"], st["n_elements"],
                            st["repeat_bp"], f"{st['repeat_frac']:.4f}", f"{scaffold_frac:.4f}",
                            f"{perm['perm_mean']:.4f}" if perm else "NA",
                            f"{perm['perm_sd']:.4f}" if perm else "NA",
                            f"{perm['fold_enrichment']:.3f}" if perm else "NA",
                            f"{perm['z']:.2f}" if perm else "NA",
                            f"{perm['pval']:.4f}" if perm else "NA",
                            top_classes])

    print(f"{args.species}\t{args.scaffold}\tscaffold_len={scaffold_len}\tn_repeat_elements={len(records)}\t"
          f"scaffold_repeat_frac={scaffold_frac:.4f}")


if __name__ == "__main__":
    main()
