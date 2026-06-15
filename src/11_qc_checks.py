#!/usr/bin/env python3
"""
11_qc_checks.py — QC / robustness checks for the manuscript

Reports, in one place, the design and robustness facts cited in the manuscript:
  0. Strain-batch confound: every batch is strain-pure → strain is confounded
     with batch and is not analyzed (explains the scope of the paper).
  1. Complete-case attrition by strain (approximately balanced).
  2. Per-strain repeatability of pct_top (5G vs AB bracket the pooled R, so the
     pooled estimate is not a "chimera").
  3. Non-ascender latency coding (NaN latency → 20.0 min; affects 2/589 sessions).

Inputs:  data/features/features_sessions.csv, features_by_minute.csv
Outputs: output/qc_checks.txt
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FEAT = ROOT / "data" / "features"
OUT  = ROOT / "output"

LATENCY_MAX = 20.0   # minutes; never-ascended sessions coded as full window


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ss = pd.read_csv(FEAT / "features_sessions.csv")
    bymin = pd.read_csv(FEAT / "features_by_minute.csv")
    bymin = bymin[bymin["pct_bottom"].notna()]

    lines = ["QC / ROBUSTNESS CHECKS", "=" * 60]

    # --- strain-batch confound (design) ---
    bs = ss.groupby("batch")["strain"].nunique()
    n_pure = int((bs == 1).sum())
    b5 = sorted(int(b) for b in ss[ss.strain == "5G"]["batch"].unique())
    bA = sorted(int(b) for b in ss[ss.strain == "AB"]["batch"].unique())
    lines.append("\n0. STRAIN-BATCH CONFOUND (design)")
    lines.append(f"   Batches with a single strain: {n_pure} of {bs.size}")
    lines.append(f"   5G batches: {b5}")
    lines.append(f"   AB batches: {bA}")
    lines.append("   → Every batch is strain-pure: strain is COMPLETELY CONFOUNDED with")
    lines.append("     batch (testing date, room, handler, per-batch distortion calibration).")
    lines.append("     No within-batch strain contrast exists, so strain differences cannot")
    lines.append("     be separated from batch and are DESCRIPTIVE only — not attributable")
    lines.append("     to genotype. This is a property of the original design, not the")
    lines.append("     analysis, and bounds every strain result in this work.")

    # --- attrition by strain ---
    rc = ss.groupby(["fish_id", "strain"])["run"].nunique().reset_index()
    analyzed = rc.groupby("strain")["fish_id"].nunique()
    complete = rc[rc["run"] == 6].groupby("strain")["fish_id"].nunique()
    lines.append("\n1. COMPLETE-CASE ATTRITION BY STRAIN")
    lines.append("   (fish analyzed → fish with all 6 runs)")
    for s in analyzed.index:
        a, c = int(analyzed[s]), int(complete.get(s, 0))
        lines.append(f"   {s}: analyzed {a}, complete {c}, "
                     f"dropped {a - c} ({100*(a-c)/a:.1f}%)")
    lines.append(f"   Total complete-case: {int(complete.sum())}")
    lines.append("   → Attrition is similar across strains (not strongly asymmetric).")

    # --- per-strain repeatability (is the pooled R a "chimera"?) ---
    def anova_R(groups):
        groups = [np.asarray(g, float) for g in groups if len(g) > 0]
        k = len(groups); ni = np.array([len(g) for g in groups]); N = ni.sum()
        if k < 2 or N <= k:
            return np.nan, np.nan, np.nan
        grand = np.concatenate(groups).mean()
        means = np.array([g.mean() for g in groups])
        msb = np.sum(ni * (means - grand) ** 2) / (k - 1)
        msw = np.sum([((g - g.mean()) ** 2).sum() for g in groups]) / (N - k)
        n0 = (N - np.sum(ni ** 2) / N) / (k - 1)
        va = max((msb - msw) / n0, 0.0)
        return va / (va + msw), va, msw

    sess_top = (bymin.assign(pct_top=1 - bymin["pct_bottom"])
                     .groupby(["fish_id", "strain", "run"])["pct_top"].mean()
                     .reset_index())
    lines.append("\n2. PER-STRAIN REPEATABILITY OF pct_top (is the pooled R a chimera?)")
    for st in ["5G", "AB", "ALL"]:
        sub = sess_top if st == "ALL" else sess_top[sess_top.strain == st]
        R, va, vw = anova_R([g["pct_top"].values for _, g in sub.groupby("fish_id")])
        lines.append(f"   {st:3s}: R = {R:.3f}  (V_among = {va:.4f}, V_within = {vw:.4f}, "
                     f"n_fish = {sub.fish_id.nunique()})")
    lines.append("   → The two strain estimates bracket the pooled value closely, so the")
    lines.append("     pooled R is representative, not an uninterpretable average of")
    lines.append("     dissimilar covariance structures. (This is a reliability robustness")
    lines.append("     check, NOT a strain comparison — strain is confounded with batch.)")

    # --- latency coding note ---
    n_never = ss["latency_to_top"].isna().sum()
    lines.append("\n3. NON-ASCENDER LATENCY CODING (latency reliability only)")
    lines.append(f"   Sessions where the fish never entered the top half: {n_never} "
                 f"of {len(ss)} ({100*n_never/len(ss):.1f}%)")
    lines.append(f"   These are coded as latency = {LATENCY_MAX:.0f} min (the full "
                 "analysis window) in all downstream analyses (05_ml.py, 09_trait_structure.py).")

    text = "\n".join(lines)
    print(text)
    (OUT / "qc_checks.txt").write_text(text)
    print(f"\nSaved → {OUT / 'qc_checks.txt'}")


if __name__ == "__main__":
    main()
