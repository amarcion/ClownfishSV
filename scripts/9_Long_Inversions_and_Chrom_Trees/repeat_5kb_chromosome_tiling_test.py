#!/usr/bin/env python3
"""
Permutation test: does a 5 kb window centered on an inversion breakpoint show
enrichment/depletion of repetitive-element content (bp) relative to the empirical
distribution of repeat content in equal-sized (5 kb), non-overlapping windows
tiled across the rest of the same chromosome/scaffold?
"""
import re
import csv
import argparse
import subprocess
from bisect import bisect_right

REPEAT_DIR = "RepeatMask_scf"

RM_MAP = {sp: f"{sp}.Assembly_v2.fa.out" for sp in
          ["AKA", "AKY", "ALL", "BIA", "CRP", "EPH", "FRE", "LAT", "LAZ",
           "MCC", "OMA", "POL", "PRC", "PRD", "SAN", "SEB"]}

WINDOW = 5000  # default total window size (bp), centered on breakpoint; overridden by --window


def load_repeats(species, scaffold):
    path = f"{REPEAT_DIR}/{RM_MAP[species]}"
    cmd = ["awk", f'NR>3 && $5=="{scaffold}"', path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    intervals = []
    scaffold_len = 0
    for line in out.splitlines():
        f = line.split()
        if len(f) < 15:
            continue
        begin, end = int(f[5]), int(f[6])
        left = f[7].strip("()")
        intervals.append((begin, end))
        try:
            scaffold_len = max(scaffold_len, end + int(left))
        except ValueError:
            pass
    return intervals, scaffold_len


def merge_intervals(intervals):
    if not intervals:
        return [], []
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    starts = [m[0] for m in merged]
    ends = [m[1] for m in merged]
    return starts, ends


def coverage(starts, ends, win_start, win_end):
    i = bisect_right(ends, win_start)
    total = 0
    n = len(starts)
    while i < n and starts[i] <= win_end:
        s = max(starts[i], win_start)
        e = min(ends[i], win_end)
        if e >= s:
            total += (e - s + 1)
        i += 1
    return total


def tile_chromosome(starts, ends, scaffold_len, window=WINDOW):
    """Non-overlapping windows of `window` bp tiled from position 1; partial trailing window dropped."""
    contents = []
    pos = 1
    while pos + window - 1 <= scaffold_len:
        win_end = pos + window - 1
        contents.append(coverage(starts, ends, pos, win_end))
        pos += window
    return contents


def breakpoint_window(center, window=WINDOW):
    half = window // 2
    win_start = center - half
    win_end = win_start + window - 1
    return win_start, win_end


def main():
    global WINDOW
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("scaffold")
    ap.add_argument("start", type=int)
    ap.add_argument("stop", type=int)
    ap.add_argument("outprefix")
    ap.add_argument("--window", type=int, default=5000, help="total window size (bp), centered on breakpoint")
    args = ap.parse_args()
    WINDOW = args.window

    intervals, scaffold_len = load_repeats(args.species, args.scaffold)
    starts, ends = merge_intervals(intervals)

    background = tile_chromosome(starts, ends, scaffold_len, window=WINDOW)
    n_bg = len(background)
    bg_mean = sum(background) / n_bg
    bg_var = sum((x - bg_mean) ** 2 for x in background) / n_bg
    bg_sd = bg_var ** 0.5

    with open(f"{args.outprefix}.{WINDOW//1000}kb_tiling_enrichment.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["breakpoint_name", "position", "win_start", "win_end", "observed_bp",
                    "observed_frac", "n_background_windows", "background_mean_bp",
                    "background_mean_frac", "background_sd_bp", "fold_enrichment", "z",
                    "p_enrichment_one_sided", "p_depletion_one_sided", "p_two_sided"])
        for bp_name, pos in [("Start", args.start), ("Stop", args.stop)]:
            win_start, win_end = breakpoint_window(pos, window=WINDOW)
            win_start_c, win_end_c = max(1, win_start), min(scaffold_len, win_end)
            observed = coverage(starts, ends, win_start_c, win_end_c)
            n_ge = sum(1 for x in background if x >= observed)
            n_le = sum(1 for x in background if x <= observed)
            p_enr = (n_ge + 1) / (n_bg + 1)
            p_dep = (n_le + 1) / (n_bg + 1)
            p_two = min(1.0, 2 * min(p_enr, p_dep))
            fold = observed / bg_mean if bg_mean > 0 else float("inf")
            z = (observed - bg_mean) / bg_sd if bg_sd > 0 else float("inf")
            w.writerow([bp_name, pos, win_start, win_end, observed,
                        f"{observed/WINDOW:.4f}", n_bg, f"{bg_mean:.1f}",
                        f"{bg_mean/WINDOW:.4f}", f"{bg_sd:.1f}", f"{fold:.3f}", f"{z:.2f}",
                        f"{p_enr:.4f}", f"{p_dep:.4f}", f"{p_two:.4f}"])

    print(f"{args.species}\t{args.scaffold}\tscaffold_len={scaffold_len}\t"
          f"n_background_{WINDOW//1000}kb_windows={n_bg}\tbackground_mean_frac={bg_mean/WINDOW:.4f}")


if __name__ == "__main__":
    main()
