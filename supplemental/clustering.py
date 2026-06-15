#!/usr/bin/env python3
"""
10_clustering.py — Unsupervised behavioral phenotype clustering

Motivation
----------
The supervised classifier (05_ml.py, analysis.md §13) assumes strain is the
partition and asks how well behavior predicts it. This script asks the
complementary, unsupervised question: do natural behavioral clusters exist, and
do they reproduce the 5G/AB boundary or cut across it? rajput2022 found four
exploratory phenotypes (bold, shy, wall-huggers, active explorers) that cross
strain and sex; johnson2025 found stable high/medium/low anxiety groups.

Section 15 established two dissociable dimensions (PC1 activity, PC2 vertical
anxiety) and showed the strain difference lives on PC2 (d=0.87) while PC1 does
NOT separate strains — so clusters may well cut across strain.

Design
------
- Cluster on all FIVE z-scored traits (pct_top, latency, transitions, x_sd,
  velocity), not just the 2 PCs — keep the 30% of variance in PC3–5. PC1/PC2
  (refit identically to §15) are used only for visualization.
- Complete-case fish (all 6 runs, n=82), matching §15 so the PC axes coincide.
- k = 2, 3, 4 with k-means (n_init=25) and Ward hierarchical; pick k by silhouette.
- Stability: bootstrap-resample fish, refit k-means, score the bootstrap
  partition against the full-data partition by adjusted Rand index (ARI).
- Clusters vs strain: contingency table, chi-square, and ARI(cluster, strain).
  Low ARI ⇒ clusters do not reproduce strain.
- Cluster profiles: mean trait z-scores per cluster (phenotype description).

Inputs:  data/features/trait_sessions.csv  (from 09_trait_structure.py)
Outputs: output/clustering.txt
         output/cluster_assignments.csv
         output/figures/07a_clusters_pca.png
"""

import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score, adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
FEAT = ROOT / "data" / "features"
OUT  = ROOT / "output" / "supplemental"   # exploratory — NOT in the manuscript
FIG  = OUT / "figures"

TRAITS = ["pct_top", "latency", "transitions", "x_sd", "velocity"]
KS     = [2, 3, 4]
SEED   = 0
N_BOOT = 500


def section(t):
    return "\n" + "=" * 66 + "\n" + t + "\n" + "=" * 66


def kmeans(k):
    return KMeans(n_clusters=k, n_init=25, random_state=SEED)


def bootstrap_stability(X, k, n_boot=N_BOOT, seed=SEED):
    """Mean ARI between the full-data k-means partition and partitions fit on
    bootstrap resamples (predicted back onto all original points)."""
    rng  = np.random.default_rng(seed)
    base = kmeans(k).fit_predict(X)
    n    = len(X)
    aris = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        km  = kmeans(k).fit(X[idx])
        aris.append(adjusted_rand_score(base, km.predict(X)))
    return float(np.mean(aris)), float(np.std(aris))


def main():
    FIG.mkdir(parents=True, exist_ok=True)

    sess = pd.read_csv(FEAT / "trait_sessions.csv")
    rc   = sess.groupby("fish_id")["run"].nunique()
    cc   = sess[sess["fish_id"].isin(rc[rc == 6].index)]
    fm   = cc.groupby("fish_id")[TRAITS].mean()
    strain = cc.groupby("fish_id")["strain"].first().reindex(fm.index)

    X = StandardScaler().fit_transform(fm.values)          # cluster space (5-D)
    pca = PCA().fit(X)
    pcs = pca.transform(X)[:, :2]                          # viz space (PC1, PC2)
    # orient PCs identically to §15 (largest-magnitude loading positive)
    for j in range(2):
        if pca.components_[j, np.argmax(np.abs(pca.components_[j]))] < 0:
            pcs[:, j] *= -1

    lines = []
    print(f"Clustering {len(fm)} complete-case fish on {len(TRAITS)} z-scored traits")
    lines.append("Unsupervised clustering of %d complete-case fish (all 6 runs)" % len(fm))
    lines.append("Features: 5 z-scored per-fish trait means %s" % TRAITS)
    lines.append("PC1/PC2 (from §15) used for visualization only.")

    # ====================================================================
    # 1. CHOOSING k — silhouette + stability
    # ====================================================================
    lines.append(section("1. CHOOSING k — silhouette and bootstrap stability"))
    lines.append("  Silhouette: higher = better separation. Stability: mean ARI between")
    lines.append("  full-data and %d bootstrap k-means partitions (1=identical)." % N_BOOT)
    lines.append("\n    k   KMeans sil   Ward sil   stability ARI (mean±sd)")
    labels = {}
    for k in KS:
        km_lab   = kmeans(k).fit_predict(X)
        ward_lab = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(X)
        sil_km   = silhouette_score(X, km_lab)
        sil_wd   = silhouette_score(X, ward_lab)
        stab, sd = bootstrap_stability(X, k)
        labels[k] = km_lab
        lines.append(f"    {k}   {sil_km:9.3f}   {sil_wd:8.3f}   {stab:.3f} ± {sd:.3f}")
    best_k = max(KS, key=lambda k: silhouette_score(X, labels[k]))
    lines.append(f"\n  Best k by KMeans silhouette: k = {best_k}")

    # ====================================================================
    # 2. CLUSTERS vs STRAIN
    # ====================================================================
    lines.append(section("2. CLUSTERS vs STRAIN — do clusters reproduce 5G/AB?"))
    for k in KS:
        lab = labels[k]
        ct  = pd.crosstab(pd.Series(lab, name="cluster"),
                          pd.Series(strain.values, name="strain"))
        chi2, p, _, _ = stats.chi2_contingency(ct)
        ari = adjusted_rand_score(strain.values, lab)
        lines.append(f"\n  k = {k}:")
        for cl in sorted(ct.index):
            row = ct.loc[cl]
            n5, na = int(row.get("5G", 0)), int(row.get("AB", 0))
            tot = n5 + na
            lines.append(f"    cluster {cl}: n={tot:2d}  5G={n5:2d} AB={na:2d}  "
                         f"({100*na/tot:.0f}% AB)")
        lines.append(f"    chi-square p = {p:.4f};  ARI(cluster, strain) = {ari:+.3f}")
    lines.append("\n  ARI near 0 ⇒ clusters are largely independent of strain (cut across it);")
    lines.append("  ARI near 1 ⇒ clusters reproduce the strain boundary.")

    # ====================================================================
    # 3. CLUSTER PROFILES (best k)
    # ====================================================================
    lines.append(section(f"3. CLUSTER PROFILES (k = {best_k}) — mean trait z-scores"))
    lab = labels[best_k]
    Zdf = pd.DataFrame(X, columns=TRAITS, index=fm.index)
    Zdf["cluster"] = lab
    Zdf["PC1"] = pcs[:, 0]; Zdf["PC2"] = pcs[:, 1]
    prof = Zdf.groupby("cluster")[TRAITS + ["PC1", "PC2"]].mean()
    hdr = "    cluster  n   " + " ".join(f"{t[:8]:>8}" for t in TRAITS) + "   PC1    PC2"
    lines.append(hdr)
    for cl in sorted(prof.index):
        n = int((lab == cl).sum())
        vals = " ".join(f"{prof.loc[cl, t]:+8.2f}" for t in TRAITS)
        lines.append(f"    {cl:^7d} {n:2d}   {vals}   "
                     f"{prof.loc[cl,'PC1']:+5.2f}  {prof.loc[cl,'PC2']:+5.2f}")
    lines.append("\n  (z-scores: + above the sample mean, − below. pct_top↑/latency↓ = bolder;")
    lines.append("   transitions/x_sd/velocity↑ = more active.)")

    # ====================================================================
    # Figure: PC1/PC2 colored by cluster (best k) and by strain
    # ====================================================================
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True, sharey=True)
    palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    for cl in sorted(np.unique(lab)):
        m = lab == cl
        ax[0].scatter(pcs[m, 0], pcs[m, 1], s=42, alpha=0.8,
                      color=palette[cl % len(palette)], label=f"cluster {cl}")
    ax[0].set_title(f"k-means clusters (k={best_k})")
    ax[0].legend(frameon=False, fontsize=9)
    for st, c in [("5G", "#9467bd"), ("AB", "#8c564b")]:
        m = strain.values == st
        ax[1].scatter(pcs[m, 0], pcs[m, 1], s=42, alpha=0.8, color=c, label=st)
    ax[1].set_title("Same points, colored by strain")
    ax[1].legend(frameon=False, fontsize=9)
    for a in ax:
        a.axhline(0, color="grey", lw=0.6, ls=":"); a.axvline(0, color="grey", lw=0.6, ls=":")
        a.set_xlabel("PC1 (general activity)")
    ax[0].set_ylabel("PC2 (vertical anxiety: pct_top↑ / latency↓)")
    fig.suptitle("Behavioral phenotype clusters vs. strain", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "07a_clusters_pca.png", dpi=140)
    plt.close(fig)

    # ------------------------------------------------------------------
    out_text = "\n".join(lines)
    print(out_text)
    (OUT / "clustering.txt").write_text(out_text)

    assign = pd.DataFrame({"fish_id": fm.index, "strain": strain.values,
                           "PC1": pcs[:, 0], "PC2": pcs[:, 1]})
    for k in KS:
        assign[f"cluster_k{k}"] = labels[k]
    assign.to_csv(OUT / "cluster_assignments.csv", index=False)

    print(f"\nSaved → {OUT / 'clustering.txt'}")
    print(f"Saved → {OUT / 'cluster_assignments.csv'}")
    print(f"Saved → {FIG / '07a_clusters_pca.png'}")


if __name__ == "__main__":
    main()
