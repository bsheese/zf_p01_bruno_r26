#!/usr/bin/env python3
"""
09_measure_reliability.py — short-term reliability of the five novel-tank measures

Manuscript Table 2: single-session repeatability R (one-way ANOVA, all sessions)
and six-session composite reliability α (Cronbach, complete cases) for pct_top,
latency, transitions, x_sd, and velocity.

Also writes the per-session trait table used by the (exploratory) supplemental
analyses in `supplemental/`.

Inputs:  data/features/features_by_minute.csv, features_sessions.csv,
         tracking_interp.csv
Outputs: data/features/trait_sessions.csv   — per fish×run, 5 measures
         output/measure_reliability.txt
         output/measure_reliability.csv      — Table 2
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FEAT = ROOT / "data" / "features"
OUT  = ROOT / "output"

LATENCY_MAX = 20.0   # NaN latency → full analysis window (matches 05_ml/ supplemental)
MEASURES = ["pct_top", "latency", "transitions", "x_sd", "velocity"]


def cronbach_alpha(wide: pd.DataFrame) -> float:
    k = wide.shape[1]
    return (k / (k - 1)) * (1 - wide.var(ddof=1).sum() / wide.sum(axis=1).var(ddof=1))


def anova_repeatability(groups):
    """One-way ANOVA single-measure repeatability (Lessells & Boag 1987)."""
    groups = [np.asarray(g, float) for g in groups if len(g) > 0]
    k  = len(groups)
    ni = np.array([len(g) for g in groups])
    N  = ni.sum()
    if k < 2 or N <= k:
        return np.nan
    grand = np.concatenate(groups).mean()
    means = np.array([g.mean() for g in groups])
    msb = np.sum(ni * (means - grand) ** 2) / (k - 1)
    msw = np.sum([((g - g.mean()) ** 2).sum() for g in groups]) / (N - k)
    n0  = (N - np.sum(ni ** 2) / N) / (k - 1)
    va  = max((msb - msw) / n0, 0.0)
    return va / (va + msw) if (va + msw) > 0 else np.nan


def build_trait_sessions() -> pd.DataFrame:
    """One row per fish × run with the five session-level measures."""
    m = pd.read_csv(FEAT / "features_by_minute.csv")
    m = m[m["pct_bottom"].notna()]
    sess = (m.groupby(["fish_id", "batch", "strain", "run"])
              .agg(pct_top=("pct_bottom", lambda s: 1.0 - s.mean()),
                   velocity=("mean_velocity", "mean"))
              .reset_index())

    ss = pd.read_csv(FEAT / "features_sessions.csv")
    ss["latency_to_top"] = ss["latency_to_top"].fillna(LATENCY_MAX)
    sess = sess.merge(ss[["fish_id", "run", "latency_to_top", "n_transitions"]],
                      on=["fish_id", "run"], how="left")

    trk = pd.read_csv(FEAT / "tracking_interp.csv",
                      usecols=["fish_id", "run", "x_interp", "detected"])
    det = trk[trk["detected"] == 1].dropna(subset=["x_interp"])
    xsd = det.groupby(["fish_id", "run"])["x_interp"].std().reset_index(name="x_sd")
    sess = sess.merge(xsd, on=["fish_id", "run"], how="left")

    return sess.rename(columns={"latency_to_top": "latency", "n_transitions": "transitions"})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sess = build_trait_sessions()
    sess.to_csv(FEAT / "trait_sessions.csv", index=False)

    rc = sess.groupby("fish_id")["run"].nunique()
    cc = sess[sess["fish_id"].isin(rc[rc == 6].index)]   # complete-case for α

    lines = ["MEASURE RELIABILITY (manuscript Table 2)",
             "=" * 50,
             f"{len(sess)} sessions; {sess.fish_id.nunique()} fish; "
             f"{cc.fish_id.nunique()} complete-case (all 6 runs)\n",
             "  R = single-session repeatability (one-way ANOVA, all sessions)",
             "  alpha = 6-session composite reliability (Cronbach, complete case)\n",
             f"  {'measure':12s} {'R':>7s} {'alpha':>7s}"]
    rows = []
    for meas in MEASURES:
        R = anova_repeatability([g[meas].dropna().to_numpy()
                                 for _, g in sess.groupby("fish_id")])
        wide = cc.pivot(index="fish_id", columns="run", values=meas).dropna()
        a = cronbach_alpha(wide)
        lines.append(f"  {meas:12s} {R:7.3f} {a:7.3f}")
        rows.append({"measure": meas, "R_single": round(R, 4),
                     "alpha_6session": round(a, 4)})

    text = "\n".join(lines)
    print(text)
    (OUT / "measure_reliability.txt").write_text(text)
    pd.DataFrame(rows).to_csv(OUT / "measure_reliability.csv", index=False)
    print(f"\nSaved → {FEAT / 'trait_sessions.csv'}")
    print(f"Saved → {OUT / 'measure_reliability.txt'}")
    print(f"Saved → {OUT / 'measure_reliability.csv'}")


if __name__ == "__main__":
    main()
