#!/usr/bin/env python3
"""
trait_structure.py — EXPLORATORY multivariate structure (NOT in the manuscript)

Among-individual correlations among the five measures, PCA dimensionality, and
the 5G/AB difference on the principal axes. Demoted to supplemental because:
  - strain is completely confounded with batch (every batch is strain-pure), so
    the strain/PC results are descriptive of batches, not genotype (see
    src/11_qc_checks.py §0); and
  - the bivariate-mixed-model among-individual correlation does not separate
    within- from among-individual covariance, so it is an approximation only.
Retained as a record of the exploratory analysis. The manuscript uses only the
univariate reliabilities (src/09_measure_reliability.py, Table 2).

Inputs:  data/features/trait_sessions.csv  (from src/09_measure_reliability.py)
Outputs: output/supplemental/trait_structure.txt
         output/supplemental/trait_correlations_{raw,among}.csv
         output/supplemental/trait_pca.csv
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path

from scipy import stats
from sklearn.decomposition import PCA
from statsmodels.formula.api import mixedlm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

ROOT = Path(__file__).resolve().parents[1]
FEAT = ROOT / "data" / "features"
OUT  = ROOT / "output" / "supplemental"
TRAITS = ["pct_top", "latency", "transitions", "x_sd", "velocity"]


def cronbach_alpha(wide):
    k = wide.shape[1]
    return (k / (k - 1)) * (1 - wide.var(ddof=1).sum() / wide.sum(axis=1).var(ddof=1))


def fisher_ci(r, n, alpha=0.05):
    if not np.isfinite(r) or n < 4 or abs(r) >= 1:
        return np.nan, np.nan
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    return np.tanh(z - zc * se), np.tanh(z + zc * se)


def among_corr_mm(sess, ta, tb):
    d = sess[["fish_id", ta, tb]].dropna().copy()
    for t in (ta, tb):
        sd = d[t].std()
        d[t] = (d[t] - d[t].mean()) / (sd if sd > 0 else 1.0)
    long = d.melt(id_vars="fish_id", value_vars=[ta, tb], var_name="trait", value_name="val")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        md = mixedlm("val ~ 0 + C(trait)", long, groups=long["fish_id"],
                     re_formula="0 + C(trait)")
        mdf = md.fit(reml=True, method="lbfgs")
    c = mdf.cov_re.values
    denom = np.sqrt(c[0, 0] * c[1, 1])
    return c[0, 1] / denom if denom > 0 else np.nan


def section(t):
    return "\n" + "=" * 66 + "\n" + t + "\n" + "=" * 66


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sess = pd.read_csv(FEAT / "trait_sessions.csv")
    rc = sess.groupby("fish_id")["run"].nunique()
    cc = sess[sess["fish_id"].isin(rc[rc == 6].index)]
    fm = cc.groupby("fish_id")[TRAITS].mean()
    n_fish = len(fm)
    rel = {t: cronbach_alpha(cc.pivot(index="fish_id", columns="run", values=t).dropna())
           for t in TRAITS}

    lines = ["EXPLORATORY multivariate structure — NOT in the manuscript",
             "(strain confounded with batch; among-individual correlation approximate)"]
    raw_mat = fm.corr()

    # raw per-fish-mean correlations
    lines.append(section("RAW PER-FISH-MEAN CORRELATIONS (n=%d)" % n_fish))
    raw_rows = []
    for i, a in enumerate(TRAITS):
        for b in TRAITS[i + 1:]:
            r = raw_mat.loc[a, b]
            lo, hi = fisher_ci(r, n_fish)
            p = stats.pearsonr(fm[a], fm[b])[1]
            lines.append(f"  {a:12s} × {b:12s}  r = {r:+.3f}  [{lo:+.3f}, {hi:+.3f}]  p = {p:.4f}")
            raw_rows.append({"trait_a": a, "trait_b": b, "r_mean": round(r, 4),
                             "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                             "p": round(p, 4), "n_fish": n_fish})

    # among-individual (repeatability-corrected) correlations
    lines.append(section("AMONG-INDIVIDUAL CORRELATIONS (disattenuated / mixed-model)"))
    among_rows = []
    for i, a in enumerate(TRAITS):
        for b in TRAITS[i + 1:]:
            r = raw_mat.loc[a, b]
            f = np.sqrt(rel[a] * rel[b])
            r_dis = r / f
            lo, hi = fisher_ci(r, n_fish)
            r_mm = among_corr_mm(sess, a, b)
            lines.append(f"  {a:12s} × {b:12s}  r_dis={r_dis:+.3f} [{lo/f:+.3f},{hi/f:+.3f}]  r_mm={r_mm:+.3f}")
            among_rows.append({"trait_a": a, "trait_b": b, "r_disattenuated": round(r_dis, 4),
                               "ci_lo": round(lo / f, 4), "ci_hi": round(hi / f, 4),
                               "r_mixed_model": round(r_mm, 4)})

    # PCA + strain overlay
    lines.append(section("PCA on z-scored per-fish means; 5G/AB on the axes"))
    Z = (fm - fm.mean()) / fm.std()
    pca = PCA().fit(Z.values)
    ev, evr = pca.explained_variance_, pca.explained_variance_ratio_
    n_keep = max(2, int((ev > 1).sum()))
    load = pca.components_.copy()
    for j in range(load.shape[0]):
        if load[j, np.argmax(np.abs(load[j]))] < 0:
            load[j] *= -1
    for j in range(len(ev)):
        lines.append(f"  PC{j+1}: eigenvalue={ev[j]:.3f}, %var={100*evr[j]:.1f}")
    lines.append("  Loadings (oriented): " + ", ".join(
        f"PC{j+1}[" + ",".join(f"{t}={load[j,i]:+.2f}" for i, t in enumerate(TRAITS)) + "]"
        for j in range(n_keep)))
    pca_rows = [{"trait": t, **{f"PC{j+1}": round(load[j, i], 4) for j in range(len(ev))}}
                for i, t in enumerate(TRAITS)]

    scores = pca.transform(Z.values)
    strain = cc.groupby("fish_id")["strain"].first().reindex(fm.index).values
    lines.append("  5G/AB difference on each component (Welch t; DESCRIPTIVE — batch-confounded):")
    for j in range(n_keep):
        s = scores[:, j] * (1 if load[j].sum() >= 0 else 1)  # sign already oriented in `load`
        if pca.components_[j, np.argmax(np.abs(pca.components_[j]))] < 0:
            s = -s
        g5, ab = s[strain == "5G"], s[strain == "AB"]
        t, p = stats.ttest_ind(g5, ab, equal_var=False)
        d = (g5.mean() - ab.mean()) / np.sqrt((g5.std()**2 + ab.std()**2) / 2)
        lines.append(f"    PC{j+1}: t={t:.2f}, p={p:.4f}, d={d:+.2f}")

    text = "\n".join(lines)
    print(text)
    (OUT / "trait_structure.txt").write_text(text)
    pd.DataFrame(raw_rows).to_csv(OUT / "trait_correlations_raw.csv", index=False)
    pd.DataFrame(among_rows).to_csv(OUT / "trait_correlations_among.csv", index=False)
    pd.DataFrame(pca_rows).to_csv(OUT / "trait_pca.csv", index=False)
    print(f"\nSaved → {OUT}/trait_structure.txt (+ correlation/pca CSVs)")


if __name__ == "__main__":
    main()
