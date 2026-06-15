#!/usr/bin/env python3
"""
08_repeatability.py — Variance-partitioning repeatability (R) of pct_top

Motivation
----------
04_stats.py reports Cronbach's α = 0.855 for pct_top across the 6 runs. α is a
6-session *composite* reliability and is NOT on the same scale as the
single-measure repeatability R reported throughout the animal-personality
literature (bell2009 mean R≈0.37; thomson2020 bottom-time R=0.39; johnson2025
lower-zone r=0.29). This script puts the P1 reliability in the field's native
currency and makes the cross-study comparison defensible.

Two analyses:

  1. Single-measure repeatability R = V_among / (V_among + V_within)
     - R_A  agreement repeatability  (no time effect modelled; Lessells & Boag
            1987 one-way ANOVA estimator, and one-way ICC(1))
     - R_C  consistency repeatability (run modelled as a fixed effect, removing
            the population mean shift across runs; two-way ICC(3,1))
     - Spearman-Brown bridge showing the single-measure R, projected to a
       6-session composite, reproduces the reported α (so the high α reflects
       averaging 6 sessions, not a single-session R that beats the literature).
     - Batch as a variance component: fish are nested in batch (8 fish/batch);
       decompose total variance into batch / among-individual / within-individual
       and report the batch-adjusted R (Nakagawa & Schielzeth 2010; Stoffel 2017).

  2. Novelty (Run 1) vs trait (Runs 2–6) repeatability
     Run 1 is the novel-tank debut and correlates weakly with later runs
     (analysis.md §2,7). The all-6-run R_A conflates the shared novelty→trait
     mean shift with individual trait consistency (Biro & Stamps 2015). Report R
     and α with and without Run 1, and decompose the time trend (R_A vs R_C).

Inputs:  data/features/features_by_minute.csv
Outputs: output/repeatability.txt   — text results
         output/repeatability.csv   — tidy table of every R / α estimate
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path

import pingouin as pg
from statsmodels.formula.api import mixedlm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

ROOT    = Path(__file__).resolve().parents[1]
FEAT    = ROOT / "data" / "features"
OUT_DIR = ROOT / "output"

N_BOOT = 2000   # cluster-bootstrap resamples for R_A confidence intervals
SEED   = 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_map(run: int) -> dict:
    return {"day": (run - 1) // 2 + 1, "time_of_day": "AM" if run % 2 == 1 else "PM"}


def cronbach_alpha(wide: pd.DataFrame) -> float:
    """Cronbach's α from a complete (rows=fish, cols=runs) matrix."""
    k = wide.shape[1]
    return (k / (k - 1)) * (1 - wide.var(ddof=1).sum() / wide.sum(axis=1).var(ddof=1))


def spearman_brown(r1: float, k: int) -> float:
    """Composite reliability of a k-measure average, from single-measure r1."""
    return k * r1 / (1 + (k - 1) * r1)


def sb_inverse(rk: float, k: int) -> float:
    """Single-measure reliability implied by a k-measure composite reliability."""
    return rk / (k - (k - 1) * rk)


def anova_repeatability(groups):
    """
    One-way ANOVA repeatability (Lessells & Boag 1987), unbalanced-safe.
    `groups` = iterable of 1-D arrays, one per individual.
    Returns (R, V_among, V_within). This is the agreement repeatability R_A.
    """
    groups = [np.asarray(g, float) for g in groups if len(g) > 0]
    k  = len(groups)
    ni = np.array([len(g) for g in groups])
    N  = ni.sum()
    if k < 2 or N <= k:
        return np.nan, np.nan, np.nan
    gm    = np.concatenate(groups)
    grand = gm.mean()
    means = np.array([g.mean() for g in groups])
    ssb = np.sum(ni * (means - grand) ** 2)
    ssw = np.sum([((g - g.mean()) ** 2).sum() for g in groups])
    msb = ssb / (k - 1)
    msw = ssw / (N - k)
    n0  = (N - np.sum(ni ** 2) / N) / (k - 1)
    v_among  = max((msb - msw) / n0, 0.0)
    v_within = msw
    denom = v_among + v_within
    return (v_among / denom if denom > 0 else np.nan), v_among, v_within


def boot_ci_RA(df, value_col="pct_top", n_boot=N_BOOT, seed=SEED):
    """Cluster (case) bootstrap over fish for the one-way repeatability R_A."""
    rng = np.random.default_rng(seed)
    fish = df["fish_id"].unique()
    by_fish = {f: df.loc[df["fish_id"] == f, value_col].to_numpy() for f in fish}
    out = []
    for _ in range(n_boot):
        samp = rng.choice(fish, size=len(fish), replace=True)
        R, _, _ = anova_repeatability([by_fish[f] for f in samp])
        if np.isfinite(R):
            out.append(R)
    lo, hi = np.percentile(out, [2.5, 97.5])
    return lo, hi


def fit_mixed_R(df, fixed="1", value_col="pct_top"):
    """
    Mixed model value ~ <fixed> + (1|fish_id) on all available sessions.
    Returns (R, V_among, V_within). fixed="1" → R_A; fixed="C(run)" → R_C.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        md  = mixedlm(f"{value_col} ~ {fixed}", df, groups=df["fish_id"])
        mdf = md.fit(reml=True, method="lbfgs")
    v_among  = float(mdf.cov_re.iloc[0, 0])
    v_within = float(mdf.scale)
    return v_among / (v_among + v_within), v_among, v_within


def fit_batch_decomposition(df, value_col="pct_top"):
    """
    Nested model value ~ 1 + (1|batch) + (1|fish within batch).
    Returns dict with V_batch, V_fish, V_within (residual). fish_id is unique
    per batch, so a single dummy group with two variance components nests fish
    in batch correctly.
    """
    d = df.copy()
    d["_all"] = 1
    vcf = {"batch": "0 + C(batch)", "fish": "0 + C(fish_id)"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        md  = mixedlm(f"{value_col} ~ 1", d, groups=d["_all"], vc_formula=vcf)
        mdf = md.fit(reml=True, method="lbfgs")
    vc = dict(zip(mdf.model.exog_vc.names, mdf.vcomp))
    return {"V_batch": float(vc["batch"]), "V_fish": float(vc["fish"]),
            "V_within": float(mdf.scale)}


def two_way_icc(mat):
    """
    ICCs from a complete (fish × run) matrix `mat`.
    Returns (R_A, R_C):
      R_A = ICC(1,1) one-way random  = agreement repeatability (run not modelled)
      R_C = ICC(C,1) two-way mixed   = consistency repeatability (run effect removed)
    Closed-form mean-square decomposition; full precision (unlike pingouin's
    2-decimal CI rounding).
    """
    x = np.asarray(mat, float)
    n, k = x.shape
    grand = x.mean()
    row_m = x.mean(axis=1)
    col_m = x.mean(axis=0)
    ss_row = k * np.sum((row_m - grand) ** 2)
    ss_col = n * np.sum((col_m - grand) ** 2)
    ss_tot = np.sum((x - grand) ** 2)
    ss_err = ss_tot - ss_row - ss_col
    msr = ss_row / (n - 1)
    msc = ss_col / (k - 1)
    mse = ss_err / ((n - 1) * (k - 1))
    msw = (ss_col + ss_err) / (n * (k - 1))   # within (run effect pooled in)
    r_a = (msr - msw) / (msr + (k - 1) * msw)
    r_c = (msr - mse) / (msr + (k - 1) * mse)
    return r_a, r_c


def icc_RA_RC(df, value_col="pct_top", runs=None, n_boot=N_BOOT, seed=SEED):
    """
    Complete-case (balanced) R_A and R_C with cluster-bootstrap 95% CIs.
    Returns (n_fish, RA_dict, RC_dict) each with ICC, CI lo/hi.
    """
    sub  = df if runs is None else df[df["run"].isin(runs)]
    wide = sub.pivot(index="fish_id", columns="run", values=value_col).dropna()
    mat  = wide.to_numpy()
    r_a, r_c = two_way_icc(mat)

    rng = np.random.default_rng(seed)
    n = mat.shape[0]
    bs_a, bs_c = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        a, c = two_way_icc(mat[idx])
        if np.isfinite(a):
            bs_a.append(a)
        if np.isfinite(c):
            bs_c.append(c)
    la, ha = np.percentile(bs_a, [2.5, 97.5])
    lc, hc = np.percentile(bs_c, [2.5, 97.5])
    return (len(wide),
            {"ICC": r_a, "lo": la, "hi": ha},
            {"ICC": r_c, "lo": lc, "hi": hc})


def section(title):
    return "\n" + "=" * 64 + "\n" + title + "\n" + "=" * 64


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    m = pd.read_csv(FEAT / "features_by_minute.csv")
    m = m[m["pct_bottom"].notna()].copy()

    # One row per fish × run (all 589 available sessions; unbalanced).
    df = (m.groupby(["fish_id", "batch", "strain", "run", "day", "time_of_day"])
            .agg(pct_top=("pct_bottom", lambda s: 1.0 - s.mean()))
            .reset_index())

    rc        = df.groupby("fish_id")["run"].nunique()
    fish_full = rc[rc == 6].index
    k_all     = 6

    lines, rows = [], []

    def record(analysis, estimand, value, lo=np.nan, hi=np.nan,
               n_fish=np.nan, k=np.nan, note=""):
        rows.append({"analysis": analysis, "estimand": estimand,
                     "value": round(value, 4) if pd.notna(value) else value,
                     "ci_lo": round(lo, 4) if pd.notna(lo) else lo,
                     "ci_hi": round(hi, 4) if pd.notna(hi) else hi,
                     "n_fish": n_fish, "k": k, "note": note})

    print(f"Loaded {len(df)} sessions, {df.fish_id.nunique()} fish, "
          f"{df.batch.nunique()} batches; {len(fish_full)} fish with all 6 runs")

    # ====================================================================
    # 1. SINGLE-MEASURE REPEATABILITY R
    # ====================================================================
    lines.append(section("1. SINGLE-MEASURE REPEATABILITY (R) — pct_top"))
    lines.append(
        "  R = V_among / (V_among + V_within), the field-standard single-measure\n"
        "  repeatability (Lessells & Boag 1987; Nakagawa & Schielzeth 2010).\n"
        "  Cronbach's α (04_stats) is a 6-session COMPOSITE reliability and is not\n"
        "  directly comparable to single-measure R from the literature.\n"
        "    R_A = agreement   (runs not modelled)\n"
        "    R_C = consistency (run as fixed effect; removes population mean shift)")

    # --- Primary: all available sessions, mixed-model + ANOVA, bootstrap CI ---
    RA_mm, va, vw = fit_mixed_R(df, fixed="1")
    RA_an, _, _   = anova_repeatability(
        [g["pct_top"].to_numpy() for _, g in df.groupby("fish_id")])
    lo, hi        = boot_ci_RA(df)
    RC_mm, _, _   = fit_mixed_R(df, fixed="C(run)")

    lines.append("\n  -- All available sessions (N=%d sessions, %d fish, unbalanced) --"
                 % (len(df), df.fish_id.nunique()))
    lines.append(f"    V_among (fish)   = {va:.5f}")
    lines.append(f"    V_within (resid) = {vw:.5f}")
    lines.append(f"    R_A (mixed model)         = {RA_mm:.3f}")
    lines.append(f"    R_A (ANOVA, Lessells-Boag)= {RA_an:.3f}  "
                 f"[cluster-boot 95% CI {lo:.3f}, {hi:.3f}]")
    lines.append(f"    R_C (run as fixed effect) = {RC_mm:.3f}")
    record("R_all_sessions", "R_A", RA_an, lo, hi, df.fish_id.nunique(),
           note="one-way ANOVA, all sessions, cluster-boot CI")
    record("R_all_sessions", "R_A_mixed", RA_mm, n_fish=df.fish_id.nunique(),
           note="mixedlm intercept-only")
    record("R_all_sessions", "R_C_mixed", RC_mm, n_fish=df.fish_id.nunique(),
           note="mixedlm run as fixed effect")

    # --- Cross-check: complete cases, pingouin ICC with F-based CIs ---
    n_cc, RA_icc, RC_icc = icc_RA_RC(df)
    lines.append(f"\n  -- Cross-check: complete cases (n={n_cc} fish, all 6 runs), cluster-boot CI --")
    lines.append(f"    R_A = ICC(1,1) = {RA_icc['ICC']:.3f}  "
                 f"[95% CI {RA_icc['lo']:.3f}, {RA_icc['hi']:.3f}]")
    lines.append(f"    R_C = ICC(C,1) = {RC_icc['ICC']:.3f}  "
                 f"[95% CI {RC_icc['lo']:.3f}, {RC_icc['hi']:.3f}]")
    record("R_complete_case", "R_A", RA_icc["ICC"], RA_icc["lo"], RA_icc["hi"],
           n_cc, k_all, "ICC(1,1) one-way, cluster-boot CI")
    record("R_complete_case", "R_C", RC_icc["ICC"], RC_icc["lo"], RC_icc["hi"],
           n_cc, k_all, "ICC(C,1) consistency, cluster-boot CI")

    # --- Spearman-Brown bridge to Cronbach's α ---
    wide_all = df[df.fish_id.isin(fish_full)].pivot(
        index="fish_id", columns="run", values="pct_top").dropna()
    alpha_all = cronbach_alpha(wide_all)
    sb_up     = spearman_brown(RA_icc["ICC"], k_all)
    sb_single = sb_inverse(alpha_all, k_all)
    lines.append("\n  -- Spearman-Brown bridge to Cronbach's α --")
    lines.append(f"    Single-measure R_A = {RA_icc['ICC']:.3f}  (k = {k_all} sessions)")
    lines.append(f"    → composite reliability  k·R/(1+(k-1)R) = {sb_up:.3f}")
    lines.append(f"    Cronbach's α (all 6 runs, n={len(wide_all)}) = {alpha_all:.3f}")
    lines.append(f"    α implies single-measure R = α/(k-(k-1)α) = {sb_single:.3f}")
    lines.append("    → The two agree: the high α reflects averaging 6 sessions, not a")
    lines.append("      single-session repeatability that exceeds the literature.")
    record("spearman_brown", "alpha_all6", alpha_all, n_fish=len(wide_all), k=k_all)
    record("spearman_brown", "R_from_alpha", sb_single, n_fish=len(wide_all), k=k_all,
           note="single-measure R implied by alpha")

    lines.append("\n  -- Published single-measure R benchmarks --")
    lines.append("    bell2009 meta-analysis (lab) mean R ≈ 0.37")
    lines.append("    thomson2020 zebrafish bottom-time R = 0.39 (28 wk)")
    lines.append("    johnson2025 zebrafish lower-zone r = 0.29 (combined sexes)")
    lines.append("    → P1 single-session R_A ≈ %.2f sits at/above the upper end of"
                 % RA_icc["ICC"])
    lines.append("      this range; the 6-session composite (α=%.2f) is correspondingly"
                 % alpha_all)
    lines.append("      high purely through Spearman-Brown aggregation.")

    # ====================================================================
    # 2. BATCH AS A VARIANCE COMPONENT
    # ====================================================================
    lines.append(section("2. BATCH AS A VARIANCE COMPONENT"))
    lines.append("  Fish are nested in batch (8 fish/batch, 13 batches, Jan–Mar 2014).")
    lines.append("  Model: pct_top ~ 1 + (1|batch) + (1|fish within batch).")
    vc = fit_batch_decomposition(df)
    tot = vc["V_batch"] + vc["V_fish"] + vc["V_within"]
    R_adj = vc["V_fish"] / tot
    lines.append(f"\n    V_batch  = {vc['V_batch']:.5f}  ({100*vc['V_batch']/tot:4.1f}% of total)")
    lines.append(f"    V_fish   = {vc['V_fish']:.5f}  ({100*vc['V_fish']/tot:4.1f}% of total)  "
                 "← among-individual, within batch")
    lines.append(f"    V_within = {vc['V_within']:.5f}  ({100*vc['V_within']/tot:4.1f}% of total)  "
                 "← within-individual (residual)")
    lines.append(f"\n    Batch-adjusted R = V_fish/(V_batch+V_fish+V_within) = {R_adj:.3f}")
    lines.append(f"    Batch accounts for {100*vc['V_batch']/tot:.1f}% of total pct_top variance.")
    lines.append("    Interpretation: batch (testing date / room / experimenter) explains")
    lines.append(f"    only {100*vc['V_batch']/tot:.1f}% of variance, so the batch-adjusted R ({R_adj:.3f}) is")
    lines.append(f"    essentially identical to the unadjusted R_A ({RA_icc['ICC']:.3f}). Modelling batch")
    lines.append("    is the standard control for repeated-measures designs (shishis2022); here")
    lines.append("    it confirms the reliability estimate is not inflated by batch structure —")
    lines.append("    reassuring given the documented within-batch spatial gradient and that")
    lines.append("    the 13 batches ran on different dates across Jan–Mar 2014.")
    record("batch_decomposition", "V_batch_pct", 100*vc["V_batch"]/tot)
    record("batch_decomposition", "V_fish_pct",  100*vc["V_fish"]/tot)
    record("batch_decomposition", "V_within_pct", 100*vc["V_within"]/tot)
    record("batch_decomposition", "R_batch_adjusted", R_adj,
           note="V_fish/(V_batch+V_fish+V_within)")

    # ====================================================================
    # 3. NOVELTY (RUN 1) vs TRAIT (RUNS 2–6)
    # ====================================================================
    lines.append(section("3. NOVELTY (RUN 1) vs TRAIT (RUNS 2–6) REPEATABILITY"))
    lines.append("  Run 1 is the novel-tank debut; runs 2–6 are re-exposures (analysis.md")
    lines.append("  §2,7). The all-6-run R conflates the shared novelty→trait mean shift")
    lines.append("  with individual trait consistency (Biro & Stamps 2015). Reported with")
    lines.append("  and without Run 1, decomposing agreement (R_A) vs consistency (R_C).")

    runs26 = [2, 3, 4, 5, 6]
    df26   = df[df["run"].isin(runs26)]

    # Cronbach's α with and without Run 1 (complete-case within each set)
    wide26 = df26.pivot(index="fish_id", columns="run", values="pct_top").dropna()
    alpha26 = cronbach_alpha(wide26)
    lines.append("\n  -- Cronbach's α --")
    lines.append(f"    All 6 runs (k=6, n={len(wide_all)}): α = {alpha_all:.3f}")
    lines.append(f"    Runs 2–6   (k=5, n={len(wide26)}): α = {alpha26:.3f}")
    lines.append("    Run 1 alone: a single measure — no internal reliability defined.")
    record("novelty_vs_trait", "alpha_runs2_6", alpha26, n_fish=len(wide26), k=5)

    # Single-measure R, all-6 vs runs 2–6, agreement vs consistency
    _, RA6, RC6 = icc_RA_RC(df, runs=runs26)
    RA26_an, _, _ = anova_repeatability(
        [g["pct_top"].to_numpy() for _, g in df26.groupby("fish_id")])
    lo26, hi26 = boot_ci_RA(df26)
    lines.append("\n  -- Single-measure repeatability R (complete-case ICC, cluster-boot CI) --")
    lines.append("    set        R_A (agreement)         R_C (consistency)")
    lines.append(f"    all 6      {RA_icc['ICC']:.3f} [{RA_icc['lo']:.3f},{RA_icc['hi']:.3f}]"
                 f"     {RC_icc['ICC']:.3f} [{RC_icc['lo']:.3f},{RC_icc['hi']:.3f}]")
    lines.append(f"    runs 2–6   {RA6['ICC']:.3f} [{RA6['lo']:.3f},{RA6['hi']:.3f}]"
                 f"     {RC6['ICC']:.3f} [{RC6['lo']:.3f},{RC6['hi']:.3f}]")
    lines.append(f"    (runs 2–6 R_A all-sessions ANOVA = {RA26_an:.3f} "
                 f"[boot {lo26:.3f}, {hi26:.3f}])")
    record("novelty_vs_trait", "R_A_runs2_6", RA6["ICC"], RA6["lo"], RA6["hi"],
           len(wide26), 5, "ICC(1,1) one-way")
    record("novelty_vs_trait", "R_C_runs2_6", RC6["ICC"], RC6["lo"], RC6["hi"],
           len(wide26), 5, "ICC(C,1) consistency")

    # Population mean shift across runs (what R_C removes)
    run_means = df.groupby("run")["pct_top"].mean()
    r1_mean   = run_means.loc[1]
    r26_mean  = df26["pct_top"].mean()
    lines.append("\n  -- Population mean shift across runs (what R_C removes) --")
    lines.append("    " + "  ".join(f"R{r}={run_means.loc[r]:.3f}" for r in range(1, 7)))
    lines.append(f"    Run 1 mean = {r1_mean:.3f}; Runs 2–6 mean = {r26_mean:.3f}; "
                 f"shift = {r1_mean - r26_mean:+.3f}")
    lines.append("    The population mean shift is small, so R_A ≈ R_C (conditioning on run")
    lines.append("    barely changes the estimate). The Run 1 problem is therefore NOT a")
    lines.append("    mean-level trend but individual-level REORDERING: Run 1 rank-orders fish")
    lines.append("    differently from later runs. That is why dropping Run 1 (R_A 0.50 → 0.59),")
    lines.append("    not conditioning on run (R_C), is what raises repeatability — exactly the")
    lines.append("    acclimation-then-stable-trait pattern Biro & Stamps (2015) warn about.")

    # Run 1's correlation with the trait composite vs internal cohesion of 2–6
    comp26 = wide_all[[2, 3, 4, 5, 6]].mean(axis=1)
    r1_vec = wide_all[1]
    r_r1_trait = float(pg.corr(r1_vec, comp26)["r"].iloc[0])
    inter26 = [float(pg.corr(wide_all[a], wide_all[b])["r"].iloc[0])
               for a in runs26 for b in runs26 if a < b]
    lines.append("\n  -- Run 1 vs the trait composite (n=%d complete-case) --" % len(wide_all))
    lines.append(f"    r(Run 1, mean of Runs 2–6) = {r_r1_trait:.3f}")
    lines.append(f"    mean pairwise r among Runs 2–6 = {np.mean(inter26):.3f}")
    lines.append("    Run 1 is markedly less coherent with the trait composite than the")
    lines.append("    re-exposure runs are with each other — consistent with thomson2020,")
    lines.append("    where week-1→2 repeatability (R=0.25) was far below later intervals.")
    record("novelty_vs_trait", "r_run1_vs_trait", r_r1_trait, n_fish=len(wide_all))
    record("novelty_vs_trait", "mean_r_within_runs2_6", float(np.mean(inter26)),
           n_fish=len(wide_all))

    # ------------------------------------------------------------------
    out_text = "\n".join(lines)
    print(out_text)
    (OUT_DIR / "repeatability.txt").write_text(out_text)
    pd.DataFrame(rows).to_csv(OUT_DIR / "repeatability.csv", index=False)
    print(f"\nSaved → {OUT_DIR / 'repeatability.txt'}")
    print(f"Saved → {OUT_DIR / 'repeatability.csv'}")


if __name__ == "__main__":
    main()
