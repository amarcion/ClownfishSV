#!/usr/bin/env python3
"""
Re-run the original random-permutation enrichment test (repeat_breakpoint_analysis.py)
but retain the raw 2,000 permuted values per species/breakpoint/window size, save them
to disk, and plot the null distribution with the observed breakpoint value marked.
"""
import csv
import random
import pickle
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from repeat_breakpoint_analysis import (
    load_repeats, merge_intervals, coverage_in_window, WINDOWS, N_PERM
)

random.seed(42)

CHR18_SPECIES = {
    "AKA": (7601449, 25905007), "AKY": (7738486, 26215418), "ALL": (6932024, 25441217),
    "CRP": (7694539, 26503477), "LAT": (6639757, 25125272), "MCC": (1, 17707047),
    "OMA": (7098827, 25662019), "POL": (1892946, 20248721), "SAN": (2216801, 20651451),
    "SEB": (7241588, 25361503),
}
CHR9_SPECIES = {"POL": (40332099, 40645956)}


def run_species(species, scaffold, start, stop):
    records, scaffold_len = load_repeats(species, scaffold)
    merged = merge_intervals([(r["start"], r["end"]) for r in records])
    merged_starts = [m[0] for m in merged]
    merged_ends = [m[1] for m in merged]

    out = {}  # (bp_name, half_win) -> {"observed":..., "perm":[...]}
    for bp_name, pos in [("Start", start), ("Stop", stop)]:
        for half in WINDOWS:
            win_start = max(1, pos - half)
            win_end = min(scaffold_len, pos + half)
            win_len = win_end - win_start + 1
            observed_bp = coverage_in_window(merged_starts, merged_ends, win_start, win_end)
            observed_frac = observed_bp / win_len

            perm_fracs = []
            for _ in range(N_PERM):
                s = random.randint(1, scaffold_len - win_len)
                e = s + win_len - 1
                cov = coverage_in_window(merged_starts, merged_ends, s, e)
                perm_fracs.append(cov / win_len)

            out[(bp_name, half)] = {"observed": observed_frac, "perm": perm_fracs}
    return out


def save_raw_values(all_results, outpath):
    with open(outpath, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species", "breakpoint", "half_window_bp", "perm_index", "perm_frac"])
        for sp, res in all_results.items():
            for (bp_name, half), d in res.items():
                for i, val in enumerate(d["perm"]):
                    w.writerow([sp, bp_name, half, i, f"{val:.6f}"])
    print(f"wrote raw permutation values -> {outpath}")


def plot_grid(all_results, species_list, half_win, outpath, title):
    n_rows = len(species_list)
    fig, axes = plt.subplots(n_rows, 2, figsize=(9, 2.1 * n_rows), sharex=False)
    if n_rows == 1:
        axes = axes.reshape(1, 2)
    for row, sp in enumerate(species_list):
        for col, bp_name in enumerate(["Start", "Stop"]):
            ax = axes[row, col]
            d = all_results[sp][(bp_name, half_win)]
            ax.hist(d["perm"], bins=40, color="#7fa8c9", edgecolor="none")
            ax.axvline(d["observed"], color="crimson", linewidth=2)
            ax.set_title(f"{sp} — {bp_name}", fontsize=9)
            ax.tick_params(labelsize=7)
            if row == n_rows - 1:
                ax.set_xlabel("repeat fraction", fontsize=8)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"wrote plot -> {outpath}")


if __name__ == "__main__":
    print("Running chr18 species...")
    chr18_results = {}
    for sp, (start, stop) in CHR18_SPECIES.items():
        chr18_results[sp] = run_species(sp, "scf18.1", start, stop)
        print(f"  done {sp}")

    save_raw_values(chr18_results, "chr18_permutation_raw/chr18_permutation_values.tsv")
    with open("chr18_permutation_raw/chr18_permutation_results.pkl", "wb") as f:
        pickle.dump(chr18_results, f)

    species18 = list(CHR18_SPECIES.keys())
    for half, label in zip(WINDOWS, ["5kb_half", "20kb_half", "50kb_half"]):
        plot_grid(chr18_results, species18, half,
                   f"figures/chr18_permutation_distributions_{label}.png",
                   f"Chr18: repeat-content null distribution (window = center +/- {half} bp) — observed value in red")

    print("\nRunning chr9 (POL)...")
    chr9_results = {}
    for sp, (start, stop) in CHR9_SPECIES.items():
        chr9_results[sp] = run_species(sp, "scf09.1", start, stop)
        print(f"  done {sp}")

    save_raw_values(chr9_results, "chr9_permutation_raw/chr9_permutation_values.tsv")
    with open("chr9_permutation_raw/chr9_permutation_results.pkl", "wb") as f:
        pickle.dump(chr9_results, f)

    for half, label in zip(WINDOWS, ["5kb_half", "20kb_half", "50kb_half"]):
        plot_grid(chr9_results, list(CHR9_SPECIES.keys()), half,
                   f"figures/chr9_permutation_distributions_{label}.png",
                   f"Chr9 (POL): repeat-content null distribution (window = center +/- {half} bp) — observed value in red")
