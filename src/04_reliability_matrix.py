#!/usr/bin/env python3
"""
04_reliability_matrix.py — Inter-run correlation matrix (Figure 1)

Heatmap of Pearson correlations of per-session pct_top between all six runs,
for fish with complete data across all runs.

Inputs:  data/features/features_by_minute.csv
Outputs: output/figures/04_reliability_matrix.png
"""

import warnings
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings("ignore")


ROOT     = Path(__file__).resolve().parents[1]
FEAT_DIR = ROOT / "data" / "features"
OUT_DIR  = ROOT / "output" / "figures"


def run_label(run: int) -> str:
    day = (run - 1) // 2 + 1
    tod = "AM" if run % 2 == 1 else "PM"
    return f"Day{day}-{tod}"


def plot_correlation_matrix(df_run: pd.DataFrame, out_path: Path):
    fish_full = df_run.groupby("fish_id")["run"].nunique()
    fish_full = fish_full[fish_full == 6].index

    if len(fish_full) < 4:
        print(f"  Skipping correlation matrix (only {len(fish_full)} fish with all 6 runs)")
        return

    wide = (
        df_run[df_run["fish_id"].isin(fish_full)]
        .pivot(index="fish_id", columns="run", values="pct_top_full")
        .dropna()
    )
    corr = wide.corr()

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap="RdBu_r")
    plt.colorbar(im, ax=ax, label="Pearson r")
    labels = [run_label(r) for r in corr.columns]
    ax.set_xticks(range(6)); ax.set_yticks(range(6))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(6):
        for j in range(6):
            ax.text(j, i, f"{corr.values[i,j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title(f"Inter-run correlations (N={len(fish_full)} fish)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()
    print(f"  Saved {out_path.name}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading feature data…")
    df_min = pd.read_csv(FEAT_DIR / "features_by_minute.csv")
    df_min["pct_top"] = 1.0 - df_min["pct_bottom"]

    # Per-fish × per-run pct_top
    df_run = (
        df_min.groupby(["fish_id", "run"])["pct_top"].mean()
        .rename("pct_top_full").reset_index()
    )

    plot_correlation_matrix(df_run, OUT_DIR / "04_reliability_matrix.png")


if __name__ == "__main__":
    main()
