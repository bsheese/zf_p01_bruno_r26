#!/usr/bin/env python3
"""
05_plots.py — Behavioral visualization plots

Produces the following figures in p01bruno_r26/output/figures/:

  01_trajectory_grid.png    — x/y trajectory for each tank, one session (batch/run)
  02_minute_by_minute.png   — pct_top per minute, mean ± SE across fish
  03_run_means.png          — mean pct_top per run (Day × AM/PM)
  04_reliability_matrix.png — inter-run correlation matrix heatmap
  05_cronbach_by_batch.png  — Cronbach's alpha per batch (when ≥ 2 batches present)

Usage
-----
  python 05_plots.py                  # uses all available data
  python 05_plots.py --batch 1 --run 1  # trajectory plot for one session only
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from pathlib import Path

warnings.filterwarnings("ignore")


ROOT     = Path(__file__).resolve().parents[1]
FEAT_DIR = ROOT / "data" / "features"
SETTINGS = ROOT / "data" / "source" / "P1-BatchSettingsCombined-v5.csv"
OUT_DIR  = ROOT / "output" / "figures"

FPS          = 29.97
SAMPLE_EVERY = 5


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_polygons(row_index: int) -> list[np.ndarray]:
    cfg = pd.read_csv(SETTINGS, header=None)
    row = cfg.iloc[row_index]
    return [
        np.array([[float(row[2+t*8+i*2]), float(row[2+t*8+i*2+1])] for i in range(4)])
        for t in range(8)
    ]


def run_label(run: int) -> str:
    day = (run - 1) // 2 + 1
    tod = "AM" if run % 2 == 1 else "PM"
    return f"Day{day}-{tod}"


def sem(x):
    x = x.dropna()
    return x.std() / np.sqrt(len(x)) if len(x) > 1 else 0.0


# ---------------------------------------------------------------------------
# Plot 1 — trajectory grid for one session
# ---------------------------------------------------------------------------

def plot_trajectory_grid(df_interp: pd.DataFrame, batch: int, run: int,
                          row_index: int, out_path: Path):
    polys = load_polygons(row_index)
    sub   = df_interp[(df_interp["batch"] == batch) & (df_interp["run"] == run)].copy()

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flat

    colors = plt.cm.plasma(np.linspace(0, 1, 8))

    for t_idx in range(8):
        ax = axes[t_idx]
        tank_data = sub[sub["tank_pos"] == t_idx + 1].sort_values("frame")

        # draw tank outline
        poly_xy = polys[t_idx]
        patch = MplPolygon(poly_xy, closed=True, fill=False, edgecolor="gray", lw=1.5)
        ax.add_patch(patch)

        # draw midline (y_mid)
        if "y_mid" in tank_data.columns and len(tank_data):
            ymid = tank_data["y_mid"].iloc[0]
            xmin = poly_xy[:, 0].min(); xmax = poly_xy[:, 0].max()
            ax.axhline(ymid, xmin=0, xmax=1, color="gray", ls="--", lw=0.8, alpha=0.5)

        if len(tank_data) == 0:
            ax.set_title(f"Tank {t_idx+1} (excluded)", fontsize=8)
            ax.axis("off")
            continue

        xi = tank_data["x_interp"].values
        yi = tank_data["y_interp"].values
        valid = ~(np.isnan(xi) | np.isnan(yi))

        # colour trajectory by time
        time_norm = np.linspace(0, 1, valid.sum())
        scatter_c = plt.cm.viridis(time_norm)

        ax.scatter(xi[valid], yi[valid], c=scatter_c, s=1.5, lw=0)
        ax.set_title(f"Tank {t_idx+1} (fish {tank_data['fish_id'].iloc[0]})", fontsize=8)
        ax.set_xlim(poly_xy[:, 0].min() - 5, poly_xy[:, 0].max() + 5)
        ax.set_ylim(poly_xy[:, 1].max() + 5, poly_xy[:, 1].min() - 5)   # y flipped
        ax.set_aspect("equal")
        ax.axis("off")

    batch_label = f"Batch {batch}, Run {run} ({run_label(run)})"
    fig.suptitle(f"Fish trajectories — {batch_label}", fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=120)
    plt.close()
    print(f"  Saved {out_path.name}")


# ---------------------------------------------------------------------------
# Plot 2 — minute-by-minute pct_top mean ± SE, overall + by strain
# ---------------------------------------------------------------------------

STRAIN_COLOR = {"5G": "#2171b5", "AB": "#e6550d"}
STRAIN_LABEL = {"5G": "5G (inbred)", "AB": "AB (wild type)"}


def _minute_curve(ax, df, color, label, zorder=2):
    """Draw mean ± SE ribbon for one group onto ax."""
    grp   = df.groupby("minute")["pct_top"]
    mins  = sorted(df["minute"].unique())
    means = np.array([grp.get_group(m).mean() for m in mins])
    sems  = np.array([sem(grp.get_group(m)) for m in mins])
    n     = int(df["fish_id"].nunique())
    ax.plot(mins, means, color=color, lw=2, label=f"{label} (n={n})", zorder=zorder)
    ax.fill_between(mins, means - sems, means + sems, color=color, alpha=0.18, zorder=zorder-1)


def plot_minute_by_minute(df_min: pd.DataFrame, out_path: Path):
    """Overall grand-mean curve (all runs, all fish)."""
    grp = df_min.groupby("minute")["pct_top"]
    mins  = sorted(df_min["minute"].unique())
    means = [grp.get_group(m).mean() for m in mins]
    sems  = [sem(grp.get_group(m)) for m in mins]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.errorbar(mins, means, yerr=sems, fmt="-o", ms=4, lw=1.5, capsize=3)
    ax.axhline(0.5, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("Minute")
    ax.set_ylabel("% time in top of tank")
    ax.set_title(f"Minute-by-minute vertical position (N={df_min['fish_id'].nunique()} fish)")
    ax.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()
    print(f"  Saved {out_path.name}")


def plot_minute_by_strain(df_min: pd.DataFrame, out_path: Path):
    """
    Grand-average minute-by-minute curve, 5G vs AB, all 6 runs pooled.
    Each minute point = mean ± SE across fish × runs for that strain.
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    for strain in ["5G", "AB"]:
        sub = df_min[df_min["strain"] == strain]
        _minute_curve(ax, sub, STRAIN_COLOR[strain], STRAIN_LABEL[strain])

    # overall (thin gray behind)
    _minute_curve(ax, df_min, "gray", "Overall", zorder=1)

    ax.axhline(0.5, color="black", ls="--", lw=0.7, alpha=0.4)
    ax.set_xlabel("Minute in session")
    ax.set_ylabel("Proportion of time in top half")
    ax.set_title("Within-session habituation by strain (all 6 runs pooled, mean ± SE)")
    ax.set_ylim(0.25, 0.75)
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved {out_path.name}")


def plot_minute_by_run_and_strain(df_min: pd.DataFrame, out_path: Path):
    """
    2×3 grid: one panel per run (Day 1 AM … Day 3 PM), each showing 5G vs AB.
    Lets you see whether the within-session curve changes across test days.
    """
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=True, sharex=True)
    axes = axes.flat

    for ax, run in zip(axes, range(1, 7)):
        day = (run - 1) // 2 + 1
        tod = "AM" if run % 2 == 1 else "PM"
        run_df = df_min[df_min["run"] == run]

        for strain in ["5G", "AB"]:
            sub = run_df[run_df["strain"] == strain]
            if sub["fish_id"].nunique() < 3:
                continue
            _minute_curve(ax, sub, STRAIN_COLOR[strain], STRAIN_LABEL[strain])

        ax.axhline(0.5, color="black", ls="--", lw=0.6, alpha=0.35)
        ax.set_title(f"Run {run}  Day {day} {tod}", fontsize=9)
        ax.set_ylim(0.2, 0.85)

        if run in (1, 4):
            ax.set_ylabel("Prop. time in top half", fontsize=8)
        if run in (4, 5, 6):
            ax.set_xlabel("Minute", fontsize=8)

    # shared legend on first panel only
    axes[0].legend(fontsize=7, loc="upper left")

    fig.suptitle("Minute-by-minute top-half proportion by strain and run (mean ± SE)",
                 fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved {out_path.name}")


# ---------------------------------------------------------------------------
# Plot 2d — Run 1 detailed: spaghetti + first-vs-last scatter + slope dist
# ---------------------------------------------------------------------------

def plot_run1_detail(df_min: pd.DataFrame, out_path: Path):
    """
    Three-panel deep-dive on Run 1 (the true novel-tank exposure).

    Panel A — spaghetti plot: individual fish minute curves (semi-transparent)
               with bold mean ± SE per strain.
    Panel B — first-5-min vs last-5-min per fish: one dot per fish, colour =
               strain; diagonal = no change; above diagonal = habituated upward.
    Panel C — distribution of per-fish linear slopes (minute → pct_top) by
               strain, shown as violin + individual points.
    """
    from scipy import stats as spstats

    r1 = df_min[df_min["run"] == 1].copy()
    r1["pct_top"] = 1.0 - r1["pct_bottom"]
    r1 = r1[r1["pct_top"].notna()]

    strains = ["5G", "AB"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # ------------------------------------------------------------------
    # Panel A — spaghetti
    # ------------------------------------------------------------------
    ax = axes[0]
    for strain in strains:
        sub = r1[r1["strain"] == strain]
        color = STRAIN_COLOR[strain]
        # individual fish (thin, transparent)
        for fid, fgrp in sub.groupby("fish_id"):
            fgrp = fgrp.sort_values("minute")
            ax.plot(fgrp["minute"], fgrp["pct_top"],
                    color=color, alpha=0.12, lw=0.8, zorder=1)
        # mean ± SE (bold)
        _minute_curve(ax, sub, color, STRAIN_LABEL[strain], zorder=3)

    ax.axhline(0.5, color="black", ls="--", lw=0.7, alpha=0.35)
    ax.set_xlabel("Minute in session")
    ax.set_ylabel("Proportion of time in top half")
    ax.set_title("Run 1: individual fish + mean ± SE")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)

    # ------------------------------------------------------------------
    # Panel B — first 5 vs last 5 minutes scatter
    # ------------------------------------------------------------------
    ax = axes[1]
    slope_rows = []
    for strain in strains:
        sub = r1[r1["strain"] == strain]
        for fid, fgrp in sub.groupby("fish_id"):
            fgrp = fgrp.sort_values("minute")
            first5 = fgrp[fgrp["minute"] < 5]["pct_top"].mean()
            last5  = fgrp[fgrp["minute"] >= 15]["pct_top"].mean()
            mins   = fgrp["minute"].values.astype(float)
            pts    = fgrp["pct_top"].values
            if len(mins) >= 5:
                slope, *_ = spstats.linregress(mins, pts)
            else:
                slope = np.nan
            slope_rows.append({"fish_id": fid, "strain": strain,
                                "first5": first5, "last5": last5, "slope": slope})
        color = STRAIN_COLOR[strain]
        fish_df = pd.DataFrame(slope_rows)
        s_df = fish_df[fish_df["strain"] == strain]
        ax.scatter(s_df["first5"], s_df["last5"],
                   color=color, alpha=0.65, s=28, label=STRAIN_LABEL[strain], zorder=3)

    lim = [0, 1]
    ax.plot(lim, lim, "k--", lw=0.8, alpha=0.4, zorder=1)   # diagonal = no change
    ax.set_xlabel("First 5 min — prop. time in top half")
    ax.set_ylabel("Last 5 min — prop. time in top half")
    ax.set_title("Run 1: first vs. last 5 minutes (per fish)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.legend(fontsize=8)

    slope_df = pd.DataFrame(slope_rows)

    # ------------------------------------------------------------------
    # Panel C — per-fish slope distribution
    # ------------------------------------------------------------------
    ax = axes[2]
    positions = [1, 2]
    for pos, strain in zip(positions, strains):
        slopes = slope_df[slope_df["strain"] == strain]["slope"].dropna()
        vp = ax.violinplot(slopes, positions=[pos], widths=0.6,
                           showmedians=True, showextrema=False)
        for pc in vp["bodies"]:
            pc.set_facecolor(STRAIN_COLOR[strain])
            pc.set_alpha(0.55)
        vp["cmedians"].set_color(STRAIN_COLOR[strain])
        vp["cmedians"].set_linewidth(2)
        # jitter
        jx = np.random.default_rng(42).uniform(pos - 0.15, pos + 0.15, len(slopes))
        ax.scatter(jx, slopes, color=STRAIN_COLOR[strain], alpha=0.45, s=18, zorder=3)

        m, s, n = slopes.mean(), slopes.std(), len(slopes)
        ax.text(pos, slopes.max() + 0.002, f"M={m:.4f}\nn={n}", ha="center",
                fontsize=7.5, color=STRAIN_COLOR[strain])

    ax.axhline(0, color="black", ls="--", lw=0.8, alpha=0.4)
    ax.set_xticks(positions)
    ax.set_xticklabels([STRAIN_LABEL[s] for s in strains], fontsize=8)
    ax.set_ylabel("Per-fish slope (prop./min)")
    ax.set_title("Run 1: habituation slope by strain")

    # t-test on slopes
    s5g = slope_df[slope_df["strain"] == "5G"]["slope"].dropna()
    sab = slope_df[slope_df["strain"] == "AB"]["slope"].dropna()
    if len(s5g) >= 3 and len(sab) >= 3:
        t, p = spstats.ttest_ind(s5g, sab, equal_var=False)
        ax.set_xlabel(f"Welch t = {t:.2f}, p = {p:.3f}", fontsize=8)

    fig.suptitle("Run 1 — novel-tank habituation (first exposure)", fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved {out_path.name}")


# ---------------------------------------------------------------------------
# Plot 3 — run means bar chart
# ---------------------------------------------------------------------------

def plot_run_means(df_run: pd.DataFrame, out_path: Path):
    stats = df_run.groupby("run").agg(
        mean=("pct_top_full", "mean"),
        se=("pct_top_full", lambda x: sem(x)),
        n=("pct_top_full", "count"),
    ).reset_index()

    labels = [run_label(r) for r in stats["run"]]
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.bar(x, stats["mean"], yerr=stats["se"], capsize=4,
                  color=plt.cm.Blues(np.linspace(0.4, 0.8, len(x))), edgecolor="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("% time in top of tank")
    ax.set_title(f"Mean % top per run (N≤{stats['n'].max()} fish per run)")
    ax.set_ylim(0, 1)
    ax.axhline(0.5, color="gray", ls="--", lw=0.8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()
    print(f"  Saved {out_path.name}")


# ---------------------------------------------------------------------------
# Plot 4 — inter-run correlation heatmap
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(batch_traj=None, run_traj=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading feature data…")
    df_min   = pd.read_csv(FEAT_DIR / "features_by_minute.csv")
    df_interp = pd.read_csv(FEAT_DIR / "tracking_interp.csv")

    df_min["pct_top"] = 1.0 - df_min["pct_bottom"]

    # Build per-fish × per-run table
    run_rows = []
    for (fish_id, batch, run, day, tod, strain), grp in df_min.groupby(
        ["fish_id", "batch", "run", "day", "time_of_day", "strain"]
    ):
        run_rows.append({
            "fish_id":      fish_id,
            "batch":        batch,
            "run":          run,
            "day":          day,
            "time_of_day":  tod,
            "strain":       strain,
            "pct_top_full": grp["pct_top"].mean(),
        })
    df_run = pd.DataFrame(run_rows)

    n_fish = df_min["fish_id"].nunique()
    print(f"  {n_fish} unique fish, {df_min['batch'].nunique()} batches")

    # Plot 1 — trajectory
    if batch_traj is None:
        batch_traj = int(df_interp["batch"].min())
        run_traj   = 1
    row_idx = (batch_traj - 1) * 6 + (run_traj - 1)
    plot_trajectory_grid(df_interp, batch_traj, run_traj, row_idx,
                         OUT_DIR / "01_trajectory_grid.png")

    # Plot 2 — minute-by-minute (overall)
    plot_minute_by_minute(df_min, OUT_DIR / "02_minute_by_minute.png")

    # Plot 2b — minute-by-minute by strain (grand average)
    plot_minute_by_strain(df_min, OUT_DIR / "02b_minute_by_strain.png")

    # Plot 2c — minute-by-minute by strain × run (small multiples)
    plot_minute_by_run_and_strain(df_min, OUT_DIR / "02c_minute_by_run_strain.png")

    # Plot 2d — Run 1 detail: spaghetti + first/last scatter + slope violin
    plot_run1_detail(df_min, OUT_DIR / "02d_run1_detail.png")

    # Plot 3 — run means
    plot_run_means(df_run, OUT_DIR / "03_run_means.png")

    # Plot 4 — correlation matrix
    plot_correlation_matrix(df_run, OUT_DIR / "04_reliability_matrix.png")

    print("\nAll plots saved to", OUT_DIR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=int)
    parser.add_argument("--run",   type=int)
    args = parser.parse_args()
    main(batch_traj=args.batch, run_traj=args.run)
