#!/usr/bin/env python3
"""
12_run1_curve.py — within-session pct_top trajectory, Run 1 vs Runs 2-6

Makes the within-session novelty response visible: pooled across fish, pct_top
by minute for the first (novel) exposure vs the familiar-tank re-exposures.
Run 1 shows the classic dive-then-explore rise; Runs 2-6 are flat. This is the
figure that pre-empts the misreading that an "essentially flat" AMONG-run mean
implies the assay elicited no novelty response (the among-run mean is flat
because Run 1's 20-min average ~= later sessions'; the within-session curve is
not).

Inputs:  data/features/features_by_minute.csv
Outputs: output/figures/run1_minute_curve.png
"""

import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FEAT = ROOT / "data" / "features"
FIG  = ROOT / "output" / "figures"


def curve(df):
    """Per-minute mean pct_top and 95% CI (across sessions)."""
    g = df.groupby("minute")["pct_top"]
    m  = g.mean()
    se = g.std() / np.sqrt(g.count())
    return m.index.values, m.values, (1.96 * se).values


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(FEAT / "features_by_minute.csv")
    m = m[m["pct_bottom"].notna()].copy()
    m["pct_top"] = 1.0 - m["pct_bottom"]

    r1  = m[m["run"] == 1]
    r26 = m[m["run"].isin([2, 3, 4, 5, 6])]

    fig, ax = plt.subplots(figsize=(7, 4.6))
    for sub, label, color in [(r1, "Run 1 (Day 1 AM, novel)", "#d62728"),
                              (r26, "Runs 2–6 (familiar)", "#1f77b4")]:
        x, y, ci = curve(sub)
        ax.plot(x, y, "-o", ms=3, color=color, label=label)
        ax.fill_between(x, y - ci, y + ci, color=color, alpha=0.15)

    ax.axhline(0.5, color="grey", lw=0.7, ls=":")
    ax.set_xlabel("Minute into session")
    ax.set_ylabel("pct_top (mean ± 95% CI across sessions)")
    ax.set_title("Within-session pct_top: novel vs. familiar exposures")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    out = FIG / "run1_minute_curve.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"Saved → {out}")


if __name__ == "__main__":
    main()
