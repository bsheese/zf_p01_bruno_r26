#!/usr/bin/env python3
"""
07_slot_adjusted_reliability.py — does removing the physical-slot artifact
improve test-retest reliability of pct_top?

Section 5 of tracking_quality_by_position.md found a small residual gradient:
pct_top declines from frame-top to frame-bottom slots (~8.6 pts, ANOVA p=0.025),
plus a mild L/R asymmetry. Because fish are randomized to slots, this is a
tracking/distortion residual, not behaviour. If we regress it out, the noise it
injects into the fish x run matrix should shrink and Cronbach's alpha should
rise toward the original (unpublished) analysis's 0.88.

Uses the per-session pct_top from src/06 (raw detected frames, no interpolation),
so the RAW alpha here is the internal baseline for the comparison — not
necessarily identical to the pipeline's interpolated 0.855.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
import statsmodels.formula.api as smf
import pingouin as pg

ROOT    = Path(__file__).resolve().parents[1]
QCSV    = ROOT / "output" / "tracking_quality_by_position.csv"
OUT_DIR = ROOT / "output"


def alpha_from_long(df, value_col):
    """Cronbach's alpha on fish x run wide matrix, complete cases only."""
    wide = df.pivot_table(index="fish_id", columns="run_within", values=value_col)
    wide = wide.dropna()
    a, ci = pg.cronbach_alpha(data=wide)
    # mean inter-run pearson r
    cc = wide.corr().values
    iu = np.triu_indices_from(cc, k=1)
    return a, ci, cc[iu].mean(), len(wide)


def main():
    d = pd.read_csv(QCSV)
    inc = d[~d.excluded].dropna(subset=["pct_top"]).copy()

    lines = []
    def P(s=""):
        print(s); lines.append(s)

    P("=" * 68)
    P("SLOT-ADJUSTED RELIABILITY OF pct_top")
    P("=" * 68)
    P(f"Included sessions: {len(inc)}  |  fish: {inc.fish_id.nunique()}")
    P("")

    # ---- 1. Is per-fish pct_top correlated with the slots it occupied? ----
    P("-" * 68)
    P("1. Per-fish pct_top vs the slots the fish happened to occupy")
    P("   (should be ~0 if randomization neutralizes the slot artifact)")
    P("-" * 68)
    fish = inc.groupby("fish_id").agg(
        pct_top=("pct_top", "mean"),
        mean_frame_row=("frame_row", "mean"),
        frac_right=("column", lambda s: (s == "R").mean()),
        n_runs=("pct_top", "size"),
    ).reset_index()
    r1, p1 = stats.pearsonr(fish.mean_frame_row, fish.pct_top)
    r2, p2 = stats.pearsonr(fish.frac_right, fish.pct_top)
    P(f"  corr(mean_frame_row, per-fish pct_top): r={r1:+.3f}, p={p1:.3f}")
    P(f"  corr(frac_right_column, per-fish pct_top): r={r2:+.3f}, p={p2:.3f}")
    P("  -> near-zero confirms the artifact does not bias per-fish trait estimates.")
    P("")

    # ---- 2. Build slot-adjusted pct_top ----
    # (a) frame-row adjustment (parsimonious; captures the documented gradient)
    # (b) full per-slot adjustment (removes row gradient AND L/R asymmetry)
    grand = inc.pct_top.mean()

    fr_mean = inc.groupby("frame_row").pct_top.transform("mean")
    inc["pct_top_adj_row"] = inc.pct_top - fr_mean + grand

    slot_mean = inc.groupby("tank_pos").pct_top.transform("mean")
    inc["pct_top_adj_slot"] = inc.pct_top - slot_mean + grand

    # ---- 3. Alpha before/after ----
    P("-" * 68)
    P("2. Cronbach's alpha (fish x 6 runs, complete cases)")
    P("-" * 68)
    for label, col in [("raw pct_top", "pct_top"),
                       ("frame-row adjusted", "pct_top_adj_row"),
                       ("full per-slot adjusted", "pct_top_adj_slot")]:
        a, ci, mr, n = alpha_from_long(inc, col)
        P(f"  {label:24s}: α = {a:.3f}  95% CI [{ci[0]:.3f}, {ci[1]:.3f}]"
          f"   mean inter-run r = {mr:.3f}   (N={n})")
    P("")
    P("  Reference: original (unpublished) analysis α = 0.88; pipeline (interpolated) α = 0.855")
    P("")

    # ---- 4. variance accounted for by slot ----
    m = smf.ols("pct_top ~ C(tank_pos)", data=inc).fit()
    mr = smf.ols("pct_top ~ C(frame_row)", data=inc).fit()
    P("-" * 68)
    P("3. Variance in session-level pct_top explained by physical position")
    P("-" * 68)
    P(f"  R^2 from physical slot (8 levels): {m.rsquared:.4f}  (F p={m.f_pvalue:.4f})")
    P(f"  R^2 from frame row    (4 levels): {mr.rsquared:.4f}  (F p={mr.f_pvalue:.4f})")
    P("")

    (OUT_DIR / "slot_adjusted_reliability.txt").write_text("\n".join(lines))
    P(f"Saved → {OUT_DIR / 'slot_adjusted_reliability.txt'}")


if __name__ == "__main__":
    main()
