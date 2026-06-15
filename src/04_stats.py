#!/usr/bin/env python3
"""
04_stats.py — Reliability and behavioral pattern analyses

Inputs:  p01bruno_r26/data/features/features_by_minute.csv
         p01bruno_r26/data/features/features_summary.csv
         p01bruno_r26/data/features/features_sessions.csv

Outputs: p01bruno_r26/output/stats_main.txt       — all text results
         p01bruno_r26/output/stats_correlations.csv — inter-run r matrix
         p01bruno_r26/output/stats_icc.csv          — ICC table by comparison type
         p01bruno_r26/output/stats_habituation.csv  — minute-level mixed model coefs

Session structure
-----------------
  6 runs per fish = 3 days × 2 sessions (AM ≈ 11:00, PM ≈ 14:30)
  run 1 = day 1 AM, run 2 = day 1 PM, run 3 = day 2 AM, ...

Analyses
--------
  1. Cronbach's alpha — overall test-retest reliability of pct_top (6 runs)
  2. Inter-run Pearson correlations — all 15 pairings
  3. Within-day (AM vs. PM) reliability — r for each of 3 day pairs, mean r
  4. Day-to-day reliability — ICC across 3 days, AM-only and PM-only
  5. Strain differences — 5G (inbred) vs. AB (wild type): t-test on overall
     pct_top, and mixed model with strain × run interaction
  6. Within-session habituation — linear mixed model: pct_top ~ minute
     (does pct_top increase over the 20-minute session, as expected?)
     Tested overall and separately by strain.
  7. Summary statistics — per run, per strain, per day, per time-of-day
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path

import pingouin as pg
from scipy import stats
from statsmodels.formula.api import mixedlm

ROOT    = Path(__file__).resolve().parents[1]
FEAT    = ROOT / "data" / "features"
OUT_DIR = ROOT / "output"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_map(run: int) -> dict:
    """Map run 1-6 → day (1-3) and time_of_day (AM/PM)."""
    return {
        "day":         (run - 1) // 2 + 1,
        "time_of_day": "AM" if run % 2 == 1 else "PM",
    }


def cronbach_alpha(df_wide: pd.DataFrame) -> float:
    """Cronbach's alpha from wide DataFrame (rows=subjects, cols=items)."""
    k = df_wide.shape[1]
    return (k / (k - 1)) * (1 - df_wide.var(ddof=1).sum() / df_wide.sum(axis=1).var(ddof=1))


# Physical stand layout: tank_pos (1-8) is a fixed slot. Frame row 1 = top of
# camera frame, 4 = bottom; columns L (slots 1-4) and R (slots 5-8). Recovered
# from the tank polygons in P1-BatchSettingsCombined-v5.csv; see
# tracking_quality_by_position.md (the gradient is captured fully by frame row).
def frame_row_of(tank_pos) -> int:
    return ((int(tank_pos) - 1) % 4) + 1


PCT_TOP_COLS = ["pct_top_full", "pct_top_first10", "pct_top_last10"]


def add_slot_adjustment(df: pd.DataFrame, group_col: str = "frame_row") -> pd.DataFrame:
    """
    Add slot-adjusted pct_top columns ("<col>_adj").

    Fish are randomized across physical slots, so any slot effect on pct_top is a
    tracking/distortion residual rather than behaviour (tracking_quality_by_position.md
    §4-6). The adjustment removes each slot group's mean deviation from the grand
    mean, leaving per-fish trait estimates unbiased while stripping the residual
    noise — which lifts test-retest reliability (~half the 0.855→0.88 gap).

    Slot-group means are estimated from ALL available sessions in `df` for
    stability, then applied row-wise.  group_col="frame_row" (4 levels) is the
    primary adjustment; "tank_pos" (8 levels) is available and gives an identical
    alpha but is less parsimonious.
    """
    df = df.copy()
    df["frame_row"] = df["tank_pos"].map(frame_row_of)
    for col in PCT_TOP_COLS:
        grand    = df[col].mean()
        grp_mean = df.groupby(group_col)[col].transform("mean")
        df[col + "_adj"] = df[col] - grp_mean + grand
    return df


def build_run_table(df_min: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse minute-level data to one row per fish × run.
    Returns pct_top_full, pct_top_first10, pct_top_last10, mean_velocity, tank_pos.

    tank_pos is the physical stand slot, constant within a fish×run; carried
    through so reliability can be computed on slot-adjusted pct_top.
    """
    rows = []
    for (fish_id, batch, run, day, tod, strain, tank_pos), grp in df_min.groupby(
        ["fish_id", "batch", "run", "day", "time_of_day", "strain", "tank_pos"]
    ):
        rows.append({
            "fish_id":         fish_id,
            "batch":           batch,
            "strain":          strain,
            "run":             run,
            "day":             day,
            "time_of_day":     tod,
            "tank_pos":        tank_pos,
            "pct_top_full":    1.0 - grp["pct_bottom"].mean(),
            "pct_top_first10": 1.0 - grp[grp["minute"] < 10]["pct_bottom"].mean(),
            "pct_top_last10":  1.0 - grp[grp["minute"] >= 10]["pct_bottom"].mean(),
            "mean_velocity":   grp["mean_velocity"].mean(),
        })
    return pd.DataFrame(rows)


def section(title: str) -> str:
    return "\n" + "=" * 60 + "\n" + title + "\n" + "=" * 60


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------

def analysis_cronbach(df: pd.DataFrame, fish_full: pd.Index, lines: list):
    lines.append(section("1. CRONBACH'S ALPHA — pct_top across 6 runs"))
    lines.append(f"Fish with all 6 runs: {len(fish_full)} / {df['fish_id'].nunique()}")

    def _alpha(col):
        wide = (
            df[df["fish_id"].isin(fish_full)]
            .pivot(index="fish_id", columns="run", values=col)
            .dropna()
        )
        a, ci = pg.cronbach_alpha(data=wide)
        return a, ci, len(wide)

    lines.append("\n  Raw pct_top:")
    for label, col in [
        ("full 20 min",   "pct_top_full"),
        ("first 10 min",  "pct_top_first10"),
        ("last 10 min",   "pct_top_last10"),
    ]:
        a, ci, n = _alpha(col)
        lines.append(f"    {label} (N={n}): α = {a:.3f}, 95% CI [{ci[0]:.3f}, {ci[1]:.3f}]")

    lines.append("\n  Slot-adjusted pct_top (physical-position residual removed; see")
    lines.append("  tracking_quality_by_position.md). Fish are randomized across slots,")
    lines.append("  so this removes tracking noise without biasing per-fish trait estimates:")
    for label, col in [
        ("full 20 min",   "pct_top_full_adj"),
        ("first 10 min",  "pct_top_first10_adj"),
        ("last 10 min",   "pct_top_last10_adj"),
    ]:
        a, ci, n = _alpha(col)
        lines.append(f"    {label} (N={n}): α = {a:.3f}, 95% CI [{ci[0]:.3f}, {ci[1]:.3f}]")

    lines.append("\n  Reference: original (unpublished) analysis α = .88 (full 20 min)")


def analysis_interrun_corr(df: pd.DataFrame, fish_full: pd.Index, lines: list) -> pd.DataFrame:
    lines.append(section("2. INTER-RUN PEARSON CORRELATIONS — pct_top"))

    sub = df[df["fish_id"].isin(fish_full)][["fish_id", "run", "pct_top_full"]]
    corr_rows = []
    for r1 in range(1, 7):
        for r2 in range(r1 + 1, 7):
            a = sub[sub["run"] == r1].set_index("fish_id")["pct_top_full"]
            b = sub[sub["run"] == r2].set_index("fish_id")["pct_top_full"]
            ab = pd.concat([a, b], axis=1, keys=["r1", "r2"]).dropna()
            if len(ab) > 3:
                res = pg.corr(ab["r1"], ab["r2"])
                r = float(res["r"].iloc[0])
                p = float(res["p_val"].iloc[0])
                m1, m2 = run_map(r1), run_map(r2)
                label = f"R{r1}(D{m1['day']}{m1['time_of_day']}) vs R{r2}(D{m2['day']}{m2['time_of_day']})"
                corr_rows.append({"run1": r1, "run2": r2, "label": label,
                                   "r": round(r, 3), "p": round(p, 4), "n": len(ab)})
                lines.append(f"  {label}: r = {r:.3f}, p = {p:.4f}, n = {len(ab)}")

    corr_df = pd.DataFrame(corr_rows)
    if len(corr_df):
        lines.append(f"\n  Mean inter-run r = {corr_df['r'].mean():.3f}")
    return corr_df


def analysis_within_day_reliability(df: pd.DataFrame, fish_full: pd.Index, lines: list):
    lines.append(section("3. WITHIN-DAY (AM vs. PM) RELIABILITY"))
    lines.append("  Pearson r between AM and PM sessions within each day.")

    sub = df[df["fish_id"].isin(fish_full)][["fish_id", "run", "pct_top_full"]]
    rs = []
    for day, (am_run, pm_run) in enumerate([(1, 2), (3, 4), (5, 6)], start=1):
        am = sub[sub["run"] == am_run].set_index("fish_id")["pct_top_full"]
        pm = sub[sub["run"] == pm_run].set_index("fish_id")["pct_top_full"]
        ab = pd.concat([am, pm], axis=1, keys=["am", "pm"]).dropna()
        if len(ab) > 3:
            res = pg.corr(ab["am"], ab["pm"])
            r = float(res["r"].iloc[0])
            p = float(res["p_val"].iloc[0])
            rs.append(r)
            lines.append(f"  Day {day} (runs {am_run} AM vs {pm_run} PM): "
                         f"r = {r:.3f}, p = {p:.4f}, n = {len(ab)}")

    if rs:
        lines.append(f"\n  Mean within-day r = {np.mean(rs):.3f}")


def analysis_day_to_day_reliability(df: pd.DataFrame, fish_full: pd.Index,
                                    lines: list) -> pd.DataFrame:
    lines.append(section("4. DAY-TO-DAY RELIABILITY"))
    lines.append("  ICC(2,1) across the 3 days, computed separately for AM and PM sessions.")

    sub = df[df["fish_id"].isin(fish_full)].copy()
    icc_rows = []

    for tod, runs in [("AM", [1, 3, 5]), ("PM", [2, 4, 6])]:
        tod_df = sub[sub["run"].isin(runs)][["fish_id", "day", "pct_top_full"]].dropna()
        if tod_df["fish_id"].nunique() < 5:
            continue
        icc_res = pg.intraclass_corr(
            data=tod_df, targets="fish_id", raters="day", ratings="pct_top_full"
        )
        # ICC(A,1) = two-way random, absolute agreement, single measures
        row = icc_res[icc_res["Type"] == "ICC(A,1)"].iloc[0]
        ci = row["CI95"]
        lines.append(f"\n  {tod} sessions (runs {runs}):")
        lines.append(f"    ICC(A,1) = {row['ICC']:.3f}, "
                     f"95% CI [{ci[0]:.3f}, {ci[1]:.3f}], "
                     f"F({row['df1']:.0f},{row['df2']:.0f}) = {row['F']:.2f}, "
                     f"p = {row['pval']:.4f}")
        icc_rows.append({
            "comparison": f"day-to-day_{tod}",
            "ICC":    round(row["ICC"], 3),
            "CI_lo":  round(ci[0], 3),
            "CI_hi":  round(ci[1], 3),
            "F":      round(row["F"], 3),
            "df1":    row["df1"],
            "df2":    row["df2"],
            "p":      round(row["pval"], 4),
        })

    return pd.DataFrame(icc_rows)


def analysis_strain(df: pd.DataFrame, lines: list):
    lines.append(section("5. STRAIN DIFFERENCES — 5G (inbred) vs. AB (wild type)"))

    # Per-fish mean pct_top across all available runs
    fish_means = df.groupby(["fish_id", "strain"])["pct_top_full"].mean().reset_index()
    g5 = fish_means[fish_means["strain"] == "5G"]["pct_top_full"].dropna()
    ab = fish_means[fish_means["strain"] == "AB"]["pct_top_full"].dropna()

    lines.append(f"\n  5G: M = {g5.mean():.3f}, SD = {g5.std():.3f}, N = {len(g5)} fish")
    lines.append(f"  AB: M = {ab.mean():.3f}, SD = {ab.std():.3f}, N = {len(ab)} fish")

    if len(g5) >= 3 and len(ab) >= 3:
        t, p = stats.ttest_ind(g5, ab, equal_var=False)
        d = (g5.mean() - ab.mean()) / np.sqrt((g5.std()**2 + ab.std()**2) / 2)
        lines.append(f"\n  Welch t-test: t = {t:.3f}, p = {p:.4f}, Cohen's d = {d:.3f}")

    # Mixed model: pct_top_full ~ strain + run (with fish as random effect)
    # run coded as continuous to capture linear trend across sessions
    lm_df = df[["fish_id", "strain", "run", "pct_top_full"]].dropna().copy()
    lm_df["run_c"] = lm_df["run"] - lm_df["run"].mean()  # center run
    try:
        md = mixedlm("pct_top_full ~ C(strain) + run_c + C(strain):run_c",
                     lm_df, groups=lm_df["fish_id"])
        mdf = md.fit(reml=True, method="lbfgs")
        lines.append("\n  Mixed model: pct_top_full ~ strain * run (fish random intercept)")
        lines.append("  Fixed effects:")
        for name, coef, p in zip(mdf.params.index, mdf.params, mdf.pvalues):
            if name != "Group Var":
                lines.append(f"    {name:40s}  b = {coef:7.4f}, p = {p:.4f}")
    except Exception as e:
        lines.append(f"\n  Mixed model error: {e}")


def analysis_habituation(df_min: pd.DataFrame, lines: list) -> pd.DataFrame:
    lines.append(section("6. WITHIN-SESSION HABITUATION — pct_top over 20 minutes"))
    lines.append("  Does pct_top increase across the session (classic novel-tank pattern)?")
    lines.append("  Model: pct_top ~ minute + (1 + minute | fish_id), per run and overall.")

    df_min = df_min.copy()
    df_min["pct_top"] = 1.0 - df_min["pct_bottom"]
    df_min = df_min[df_min["pct_top"].notna()].copy()
    df_min["minute_c"] = df_min["minute"] - df_min["minute"].mean()

    coef_rows = []

    # Overall model (all runs pooled)
    try:
        md = mixedlm("pct_top ~ minute_c", df_min, groups=df_min["fish_id"])
        mdf = md.fit(reml=True, method="lbfgs")
        slope = mdf.params["minute_c"]
        p_slope = mdf.pvalues["minute_c"]
        se = mdf.bse["minute_c"]
        lines.append(f"\n  Overall (all runs):")
        lines.append(f"    minute slope = {slope:.5f} per min (SE={se:.5f}), p = {p_slope:.4f}")
        lines.append(f"    Intercept (at mean minute) = {mdf.params['Intercept']:.3f}")
        lines.append(f"    → Over 20 min: expected pct_top change = {slope * 19:.3f}")
        coef_rows.append({"group": "all", "slope": round(slope, 6),
                          "se": round(se, 6), "p": round(p_slope, 4)})
    except Exception as e:
        lines.append(f"\n  Overall model error: {e}")

    # By strain
    for strain in ["5G", "AB"]:
        sub = df_min[df_min["strain"] == strain]
        if sub["fish_id"].nunique() < 5:
            continue
        try:
            md = mixedlm("pct_top ~ minute_c", sub, groups=sub["fish_id"])
            mdf = md.fit(reml=True, method="lbfgs")
            slope = mdf.params["minute_c"]
            p_slope = mdf.pvalues["minute_c"]
            se = mdf.bse["minute_c"]
            lines.append(f"\n  Strain {strain}:")
            lines.append(f"    minute slope = {slope:.5f} per min (SE={se:.5f}), p = {p_slope:.4f}")
            lines.append(f"    Intercept = {mdf.params['Intercept']:.3f}")
            coef_rows.append({"group": strain, "slope": round(slope, 6),
                              "se": round(se, 6), "p": round(p_slope, 4)})
        except Exception as e:
            lines.append(f"\n  Strain {strain} model error: {e}")

    # Strain interaction (does slope differ by strain?)
    try:
        md = mixedlm("pct_top ~ minute_c * C(strain)", df_min, groups=df_min["fish_id"])
        mdf = md.fit(reml=True, method="lbfgs")
        lines.append("\n  Strain × minute interaction:")
        for name, coef, p in zip(mdf.params.index, mdf.params, mdf.pvalues):
            if "strain" in name.lower() or "minute" in name.lower():
                lines.append(f"    {name:45s}  b = {coef:7.5f}, p = {p:.4f}")
    except Exception as e:
        lines.append(f"\n  Strain interaction model error: {e}")

    return pd.DataFrame(coef_rows)


def analysis_extra_metrics(sess: pd.DataFrame, fish_full: pd.Index, lines: list):
    """Reliability and strain analyses for latency_to_top and n_transitions."""
    lines.append(section("8. LATENCY TO TOP AND TRANSITIONS — reliability and strain"))
    lines.append("  latency_to_top : minutes to first detected frame in top half (NaN = never entered)")
    lines.append("  n_transitions  : top↔bottom zone crossings in detected frames")

    for metric, label, invert in [
        ("latency_to_top", "Latency to top (min)", False),
        ("n_transitions",  "N transitions",        False),
    ]:
        lines.append(f"\n  --- {label} ---")

        sub = sess[sess["fish_id"].isin(fish_full)][["fish_id","run","strain",metric]].copy()

        # descriptives per run
        lines.append("  Per run (M ± SD, N non-missing):")
        for run in range(1, 7):
            m = run_map(run)
            v = sub[sub["run"] == run][metric].dropna()
            lines.append(f"    Run {run} (D{m['day']}{m['time_of_day']}): "
                         f"M = {v.mean():.2f}, SD = {v.std():.2f}, N = {len(v)}")

        # strain comparison (per-fish mean across runs)
        fish_means = sub.groupby(["fish_id","strain"])[metric].mean().reset_index()
        g5 = fish_means[fish_means["strain"]=="5G"][metric].dropna()
        ab = fish_means[fish_means["strain"]=="AB"][metric].dropna()
        lines.append(f"\n  Strain means (across all runs):")
        lines.append(f"    5G: M = {g5.mean():.2f}, SD = {g5.std():.2f}, N = {len(g5)}")
        lines.append(f"    AB: M = {ab.mean():.2f}, SD = {ab.std():.2f}, N = {len(ab)}")
        if len(g5) >= 3 and len(ab) >= 3:
            t, p = stats.ttest_ind(g5, ab, equal_var=False)
            d = (g5.mean() - ab.mean()) / np.sqrt((g5.std()**2 + ab.std()**2) / 2)
            lines.append(f"    Welch t = {t:.3f}, p = {p:.4f}, d = {d:.3f}")

        # Cronbach's alpha (fish × run wide)
        wide = sub.pivot(index="fish_id", columns="run", values=metric).dropna()
        if len(wide) >= 5 and wide.shape[1] >= 2:
            try:
                alpha, ci = pg.cronbach_alpha(data=wide)
                lines.append(f"\n  Cronbach's α (all 6 runs, N={len(wide)}): "
                             f"{alpha:.3f}, 95% CI [{ci[0]:.3f}, {ci[1]:.3f}]")
            except Exception as e:
                lines.append(f"\n  Cronbach's α error: {e}")

        # Inter-run correlations
        lines.append("  Inter-run correlations (r):")
        rs = []
        for r1 in range(1, 7):
            for r2 in range(r1+1, 7):
                a = sub[sub["run"]==r1].set_index("fish_id")[metric]
                b = sub[sub["run"]==r2].set_index("fish_id")[metric]
                ab_df = pd.concat([a, b], axis=1, keys=["r1","r2"]).dropna()
                if len(ab_df) > 3:
                    res = pg.corr(ab_df["r1"], ab_df["r2"])
                    r = float(res["r"].iloc[0])
                    p = float(res["p_val"].iloc[0])
                    rs.append(r)
                    m1, m2 = run_map(r1), run_map(r2)
                    lines.append(f"    R{r1}(D{m1['day']}{m1['time_of_day']}) vs "
                                 f"R{r2}(D{m2['day']}{m2['time_of_day']}): "
                                 f"r = {r:.3f}, p = {p:.4f}")
        if rs:
            lines.append(f"  Mean inter-run r = {np.mean(rs):.3f}")

        # ICC day-to-day (AM and PM separately)
        lines.append("  Day-to-day ICC(A,1):")
        for tod, runs in [("AM", [1,3,5]), ("PM", [2,4,6])]:
            tod_df = (sess[sess["fish_id"].isin(fish_full) & sess["run"].isin(runs)]
                      [["fish_id","day",metric]].dropna())
            # keep only fish with complete data across all 3 days
            tod_df = tod_df.groupby("fish_id").filter(lambda x: len(x) == 3)
            if tod_df["fish_id"].nunique() < 5:
                continue
            try:
                icc_res = pg.intraclass_corr(
                    data=tod_df, targets="fish_id", raters="day", ratings=metric)
                row = icc_res[icc_res["Type"] == "ICC(A,1)"].iloc[0]
                ci = row["CI95"]
                n = tod_df["fish_id"].nunique()
                lines.append(f"    {tod} (N={n}): ICC(A,1) = {row['ICC']:.3f}, "
                             f"95% CI [{ci[0]:.3f}, {ci[1]:.3f}], "
                             f"p = {row['pval']:.4f}")
            except Exception as e:
                lines.append(f"    {tod}: ICC error: {e}")


def analysis_summary(df: pd.DataFrame, lines: list):
    lines.append(section("7. SUMMARY STATISTICS"))

    lines.append("\n  Per run (M ± SD, N fish):")
    for run, grp in df.groupby("run"):
        m = run_map(run)
        v = grp["pct_top_full"].dropna()
        lines.append(f"    Run {run} (Day {m['day']} {m['time_of_day']}): "
                     f"M = {v.mean():.3f}, SD = {v.std():.3f}, N = {len(v)}")

    lines.append("\n  Per strain (mean across all runs, N fish):")
    fish_means = df.groupby(["fish_id", "strain"])["pct_top_full"].mean().reset_index()
    for strain, grp in fish_means.groupby("strain"):
        v = grp["pct_top_full"].dropna()
        lines.append(f"    {strain}: M = {v.mean():.3f}, SD = {v.std():.3f}, N = {len(v)}")

    lines.append("\n  AM vs. PM (mean across all days, N obs):")
    for tod, grp in df.groupby("time_of_day"):
        v = grp["pct_top_full"].dropna()
        lines.append(f"    {tod}: M = {v.mean():.3f}, SD = {v.std():.3f}, N = {len(v)}")

    lines.append("\n  Per day (mean across AM+PM, N obs):")
    for day, grp in df.groupby("day"):
        v = grp["pct_top_full"].dropna()
        lines.append(f"    Day {day}: M = {v.mean():.3f}, SD = {v.std():.3f}, N = {len(v)}")

    lines.append("\n  First 10 min vs. last 10 min (overall):")
    for label, col in [("First 10 min", "pct_top_first10"), ("Last 10 min", "pct_top_last10")]:
        v = df[col].dropna()
        lines.append(f"    {label}: M = {v.mean():.3f}, SD = {v.std():.3f}, N = {len(v)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df_min  = pd.read_csv(FEAT / "features_by_minute.csv")
    df_min  = df_min[df_min["pct_bottom"].notna()].copy()
    df_sess = pd.read_csv(FEAT / "features_sessions.csv")

    print(f"Loaded {len(df_min)} minute-rows, "
          f"{df_min['fish_id'].nunique()} fish, "
          f"{df_min.groupby(['fish_id','run']).ngroups} fish×run sessions")

    df = build_run_table(df_min)
    df = add_slot_adjustment(df, group_col="frame_row")
    df.to_csv(OUT_DIR / "run_table_slot_adjusted.csv", index=False)

    # Fish with all 6 runs present (needed for alpha / correlation analyses)
    run_counts = df.groupby("fish_id")["run"].nunique()
    fish_full  = run_counts[run_counts == 6].index
    print(f"  Fish with all 6 runs: {len(fish_full)} / {df['fish_id'].nunique()}")

    lines = []

    analysis_cronbach(df, fish_full, lines)
    corr_df   = analysis_interrun_corr(df, fish_full, lines)
    analysis_within_day_reliability(df, fish_full, lines)
    icc_df    = analysis_day_to_day_reliability(df, fish_full, lines)
    analysis_strain(df, lines)
    hab_df    = analysis_habituation(df_min, lines)
    analysis_extra_metrics(df_sess, fish_full, lines)
    analysis_summary(df, lines)

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    out_text = "\n".join(lines)
    print("\n" + out_text)

    (OUT_DIR / "stats_main.txt").write_text(out_text)
    print(f"\nSaved → {OUT_DIR / 'stats_main.txt'}")

    if len(corr_df):
        corr_df.to_csv(OUT_DIR / "stats_correlations.csv", index=False)
        print(f"Saved → {OUT_DIR / 'stats_correlations.csv'}")

    if len(icc_df):
        icc_df.to_csv(OUT_DIR / "stats_icc.csv", index=False)
        print(f"Saved → {OUT_DIR / 'stats_icc.csv'}")

    if len(hab_df):
        hab_df.to_csv(OUT_DIR / "stats_habituation.csv", index=False)
        print(f"Saved → {OUT_DIR / 'stats_habituation.csv'}")


if __name__ == "__main__":
    main()
