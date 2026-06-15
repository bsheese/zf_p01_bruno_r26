#!/usr/bin/env python3
"""
05_ml.py — Strain classification from behavioral features

Binary classification: 5G (inbred, label=0) vs AB (wild type, label=1)
Unit of analysis: individual fish (not session)
Evaluation: Leave-one-out cross-validation (LOOCV)

Feature sets
------------
  A  Aggregate  (N≈103): pct_top_mean, pct_top_sd, latency_mean,
                          transitions_mean, x_sd_mean
  B  Trajectory (N=82):  pct_top_r1 … pct_top_r6
  C  Theory-driven (N≈95): delta_r1r2, novelty_effect,
                             x_sd_mean, transitions_mean, latency_mean

Classifiers
-----------
  Logistic Regression (L2, C=1, scaled)
  Linear Discriminant Analysis (scaled)
  Random Forest (200 trees, unscaled)

Outputs
-------
  output/ml_results.txt
  figures/06a_ml_roc.png        — ROC curves per feature set
  figures/06b_ml_features.png   — feature importance / coefficients
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import (roc_auc_score, accuracy_score,
                             balanced_accuracy_score, roc_curve)
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT    = Path(__file__).resolve().parents[1]
FEAT    = ROOT / "data" / "features"
OUT_DIR = ROOT / "output" / "supplemental"   # exploratory — NOT in the manuscript
FIG_DIR = OUT_DIR / "figures"

LATENCY_MAX = 20.0   # impute "never entered top" as full session length


# ---------------------------------------------------------------------------
# Feature engineering
# ---------------------------------------------------------------------------

def build_feature_matrix() -> pd.DataFrame:
    """Return one row per fish with all derived features + strain label."""

    # --- pct_top per run (from minute-level data) ---
    df_min = pd.read_csv(FEAT / "features_by_minute.csv")
    df_min["pct_top"] = 1.0 - df_min["pct_bottom"]
    run_means = (df_min.groupby(["fish_id", "strain", "run"])["pct_top"]
                       .mean().reset_index())
    wide_pct = run_means.pivot(index="fish_id", columns="run", values="pct_top")
    wide_pct.columns = [f"pct_top_r{c}" for c in wide_pct.columns]

    # --- latency + transitions per run ---
    df_sess = pd.read_csv(FEAT / "features_sessions.csv")
    df_sess["latency_to_top"] = df_sess["latency_to_top"].fillna(LATENCY_MAX)
    lat_wide  = df_sess.pivot(index="fish_id", columns="run",
                              values="latency_to_top")
    lat_wide.columns  = [f"latency_r{c}"     for c in lat_wide.columns]
    trans_wide = df_sess.pivot(index="fish_id", columns="run",
                               values="n_transitions")
    trans_wide.columns = [f"transitions_r{c}" for c in trans_wide.columns]

    # --- within-session x SD per run ---
    df_trk = pd.read_csv(FEAT / "tracking_interp.csv",
                         usecols=["fish_id", "run", "x_interp", "detected"])
    det = df_trk[df_trk["detected"] == 1].dropna(subset=["x_interp"])
    x_sd = (det.groupby(["fish_id", "run"])["x_interp"]
               .std().reset_index(name="x_sd"))
    x_sd_wide = x_sd.pivot(index="fish_id", columns="run", values="x_sd")
    x_sd_wide.columns = [f"x_sd_r{c}" for c in x_sd_wide.columns]

    # --- strain lookup ---
    strain_map = (df_sess.groupby("fish_id")["strain"].first())

    # --- merge everything ---
    df = pd.DataFrame(index=strain_map.index)
    df["strain"] = strain_map
    df = df.join(wide_pct).join(lat_wide).join(trans_wide).join(x_sd_wide)

    # --- aggregate features (mean / sd across runs) ---
    pct_cols   = [f"pct_top_r{r}"     for r in range(1, 7)]
    lat_cols   = [f"latency_r{r}"     for r in range(1, 7)]
    trans_cols = [f"transitions_r{r}" for r in range(1, 7)]
    x_sd_cols  = [f"x_sd_r{r}"        for r in range(1, 7)]

    df["pct_top_mean"]     = df[pct_cols].mean(axis=1)
    df["pct_top_sd"]       = df[pct_cols].std(axis=1)
    df["latency_mean"]     = df[lat_cols].mean(axis=1)
    df["transitions_mean"] = df[trans_cols].mean(axis=1)
    df["x_sd_mean"]        = df[x_sd_cols].mean(axis=1)

    # --- derived / theory-driven features ---
    # delta: how much pct_top DROPS from run 1 to run 2 (5G drops, AB doesn't)
    df["delta_r1r2"]      = df["pct_top_r1"] - df["pct_top_r2"]
    # novelty effect: run 1 vs mean of re-exposure runs
    re_cols = [f"pct_top_r{r}" for r in range(2, 7)]
    df["novelty_effect"]  = df["pct_top_r1"] - df[re_cols].mean(axis=1)
    # cross-session slope (linear trend across runs 1–6)
    run_nums = np.arange(1, 7)
    def traj_slope(row):
        vals = row[[f"pct_top_r{r}" for r in range(1, 7)]].values.astype(float)
        mask = ~np.isnan(vals)
        if mask.sum() < 3:
            return np.nan
        return np.polyfit(run_nums[mask], vals[mask], 1)[0]
    df["traj_slope"] = df.apply(traj_slope, axis=1)

    df["label"] = (df["strain"] == "AB").astype(int)   # 0=5G, 1=AB
    return df.reset_index()


# ---------------------------------------------------------------------------
# LOOCV evaluation
# ---------------------------------------------------------------------------

def loocv_eval(X: np.ndarray, y: np.ndarray,
               clf_pipeline) -> dict:
    """
    Run LOOCV and return dict with probs, preds, and scalar metrics.
    """
    loo   = LeaveOneOut()
    probs = np.full(len(y), np.nan)
    preds = np.full(len(y), np.nan)

    for train_idx, test_idx in loo.split(X):
        clf_pipeline.fit(X[train_idx], y[train_idx])
        probs[test_idx] = clf_pipeline.predict_proba(X[test_idx])[:, 1]
        preds[test_idx] = clf_pipeline.predict(X[test_idx])

    auc  = roc_auc_score(y, probs)
    acc  = accuracy_score(y, preds)
    bacc = balanced_accuracy_score(y, preds)
    fpr, tpr, _ = roc_curve(y, probs)
    return dict(probs=probs, preds=preds, fpr=fpr, tpr=tpr,
                auc=auc, acc=acc, bacc=bacc)


def make_pipelines() -> dict:
    lr  = Pipeline([("imp", SimpleImputer(strategy="mean")),
                    ("sc",  StandardScaler()),
                    ("clf", LogisticRegression(C=1.0, max_iter=1000,
                                               solver="liblinear",
                                               random_state=42))])
    lda = Pipeline([("imp", SimpleImputer(strategy="mean")),
                    ("sc",  StandardScaler()),
                    ("clf", LinearDiscriminantAnalysis())])
    rf  = Pipeline([("imp", SimpleImputer(strategy="mean")),
                    ("clf", RandomForestClassifier(n_estimators=200,
                                                   max_features="sqrt",
                                                   random_state=42))])
    return {"LR": lr, "LDA": lda, "RF": rf}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = build_feature_matrix()
    print(f"Feature matrix: {len(df)} fish  "
          f"(5G={( df['strain']=='5G').sum()}, AB={(df['strain']=='AB').sum()})")

    # ------------------------------------------------------------------
    # Define feature sets
    # ------------------------------------------------------------------
    pct_run_cols   = [f"pct_top_r{r}"     for r in range(1, 7)]
    lat_run_cols   = [f"latency_r{r}"     for r in range(1, 7)]
    trans_run_cols = [f"transitions_r{r}" for r in range(1, 7)]

    FEATURE_SETS = {
        "A — Aggregate": {
            "cols": ["pct_top_mean", "pct_top_sd",
                     "latency_mean", "transitions_mean", "x_sd_mean"],
            "label": "Aggregate\n(means across runs)",
        },
        "B — pct_top trajectory": {
            "cols": pct_run_cols,
            "label": "pct_top trajectory\n(per-run values r1–r6)",
        },
        "C — Theory-driven": {
            "cols": ["delta_r1r2", "novelty_effect",
                     "x_sd_mean", "transitions_mean", "latency_mean"],
            "label": "Theory-driven\n(delta, novelty, activity)",
        },
        "D — Full trajectory": {
            "cols": pct_run_cols + trans_run_cols + lat_run_cols,
            "label": "Full trajectory\n(pct_top + transitions + latency r1–r6)",
        },
    }

    pipelines = make_pipelines()
    lines     = []
    all_results = {}   # (fset_name, clf_name) → result dict

    lines.append("=" * 65)
    lines.append("STRAIN CLASSIFICATION — 5G vs AB")
    lines.append("Evaluation: Leave-one-out cross-validation (LOOCV)")
    lines.append("=" * 65)

    for fset_name, fset in FEATURE_SETS.items():
        lines.append(f"\n{'─'*65}")
        lines.append(f"Feature set {fset_name}")
        lines.append(f"  Features: {fset['cols']}")

        cols = fset["cols"]
        sub  = df.dropna(subset=cols + ["label"])
        X    = sub[cols].values.astype(float)
        y    = sub["label"].values

        lines.append(f"  N fish: {len(sub)}  "
                     f"(5G={int((y==0).sum())}, AB={int((y==1).sum())})")
        lines.append(f"  {'Classifier':<8}  {'AUC':>6}  {'Acc':>6}  {'BalAcc':>8}")

        for clf_name, pipe in pipelines.items():
            res = loocv_eval(X, y, pipe)
            all_results[(fset_name, clf_name)] = {**res, "fset": fset_name,
                                                   "clf": clf_name, "n": len(sub)}
            lines.append(f"  {clf_name:<8}  {res['auc']:>6.3f}  "
                         f"{res['acc']:>6.3f}  {res['bacc']:>8.3f}")

    lines.append("\n" + "=" * 65)
    lines.append("BEST MODEL PER FEATURE SET (by AUC, RF)")
    lines.append("=" * 65)
    for fset_name in FEATURE_SETS:
        rf_res = all_results[(fset_name, "RF")]
        lines.append(f"  {fset_name[:30]:<30}  AUC={rf_res['auc']:.3f}  "
                     f"Acc={rf_res['acc']:.3f}  BalAcc={rf_res['bacc']:.3f}")

    # ------------------------------------------------------------------
    # Feature importance: fit on full data
    # ------------------------------------------------------------------
    lines.append("\n" + "=" * 65)
    lines.append("FEATURE IMPORTANCE — Theory-driven set C (RF, full fit)")
    lines.append("=" * 65)
    fset_c = FEATURE_SETS["C — Theory-driven"]
    sub_c  = df.dropna(subset=fset_c["cols"] + ["label"])
    X_c    = sub_c[fset_c["cols"]].values.astype(float)
    y_c    = sub_c["label"].values
    imp_pipe = make_pipelines()["RF"]
    imp_pipe.fit(X_c, y_c)
    importances = imp_pipe.named_steps["clf"].feature_importances_
    for feat, imp in sorted(zip(fset_c["cols"], importances),
                            key=lambda x: -x[1]):
        lines.append(f"  {feat:<25}  {imp:.4f}")

    lines.append("\n" + "=" * 65)
    lines.append("FEATURE IMPORTANCE — pct_top trajectory set B (RF, full fit)")
    lines.append("=" * 65)
    sub_b  = df.dropna(subset=pct_run_cols + ["label"])
    X_b    = sub_b[pct_run_cols].values.astype(float)
    y_b    = sub_b["label"].values
    imp_b  = make_pipelines()["RF"]
    imp_b.fit(X_b, y_b)
    for feat, imp in sorted(zip(pct_run_cols,
                                imp_b.named_steps["clf"].feature_importances_),
                            key=lambda x: -x[1]):
        lines.append(f"  {feat:<25}  {imp:.4f}")

    lines.append("\n" + "=" * 65)
    lines.append("LR COEFFICIENTS — Theory-driven set C (standardised)")
    lines.append("=" * 65)
    lr_pipe = make_pipelines()["LR"]
    lr_pipe.fit(X_c, y_c)
    coefs = lr_pipe.named_steps["clf"].coef_[0]
    for feat, coef in sorted(zip(fset_c["cols"], coefs), key=lambda x: -abs(x[1])):
        lines.append(f"  {feat:<25}  {coef:+.4f}  "
                     f"({'→ AB' if coef > 0 else '→ 5G'})")

    out_txt = "\n".join(lines)
    print(out_txt)
    (OUT_DIR / "ml_results.txt").write_text(out_txt)
    print(f"\nSaved → {OUT_DIR / 'ml_results.txt'}")

    # ------------------------------------------------------------------
    # Figure 1 — ROC curves
    # ------------------------------------------------------------------
    CLF_STYLE = {"LR":  ("-",  "#2171b5"),
                 "LDA": ("--", "#e6550d"),
                 "RF":  (":",  "#31a354")}

    fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), sharey=True)
    for ax, (fset_name, fset) in zip(axes, FEATURE_SETS.items()):
        ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.4)
        for clf_name, style in CLF_STYLE.items():
            res = all_results[(fset_name, clf_name)]
            ax.plot(res["fpr"], res["tpr"],
                    linestyle=style[0], color=style[1], lw=1.8,
                    label=f"{clf_name}  AUC={res['auc']:.3f}")
        ax.set_title(fset["label"], fontsize=9)
        ax.set_xlabel("False positive rate", fontsize=9)
        ax.legend(fontsize=7.5, loc="lower right")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_aspect("equal")
    axes[0].set_ylabel("True positive rate", fontsize=9)
    fig.suptitle("ROC curves — strain classification (5G vs AB), LOOCV",
                 fontsize=11)
    plt.tight_layout()
    roc_path = FIG_DIR / "06a_ml_roc.png"
    plt.savefig(roc_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved → {roc_path}")

    # ------------------------------------------------------------------
    # Figure 2 — Feature importance + LR coefficients
    # ------------------------------------------------------------------
    fig2, axes2 = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: RF importances — theory-driven set C
    ax = axes2[0]
    feats_c = fset_c["cols"]
    imp_c   = imp_pipe.named_steps["clf"].feature_importances_
    order   = np.argsort(imp_c)
    ax.barh([feats_c[i] for i in order], imp_c[order], color="#2171b5")
    ax.set_title("RF importance — Set C\n(theory-driven)", fontsize=10)
    ax.set_xlabel("Mean decrease in impurity", fontsize=9)

    # Panel B: LR coefficients — theory-driven set C
    ax = axes2[1]
    order2 = np.argsort(coefs)
    bar_colors = ["#e6550d" if coefs[i] > 0 else "#2171b5" for i in order2]
    ax.barh([feats_c[i] for i in order2], coefs[order2], color=bar_colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title("LR coefficients — Set C\n(orange→AB, blue→5G)", fontsize=10)
    ax.set_xlabel("Standardised coefficient", fontsize=9)

    # Panel C: RF importances — pct_top trajectory set B
    ax = axes2[2]
    imp_b_vals = imp_b.named_steps["clf"].feature_importances_
    x_pos = np.arange(1, 7)
    ax.bar(x_pos, imp_b_vals, color="#31a354", alpha=0.85)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"r{r}" for r in range(1, 7)], fontsize=9)
    ax.set_title("RF importance — Set B\n(pct_top trajectory)", fontsize=10)
    ax.set_xlabel("Run", fontsize=9)
    ax.set_ylabel("Mean decrease in impurity", fontsize=9)

    plt.tight_layout()
    feat_path = FIG_DIR / "06b_ml_features.png"
    plt.savefig(feat_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved → {feat_path}")


if __name__ == "__main__":
    main()
