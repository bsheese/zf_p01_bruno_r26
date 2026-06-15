#!/usr/bin/env python3
"""
06_tracking_quality.py — Mathematica tracking quality as a function of
physical position on the 8-tank stand.

Question: does the Mathematica tracking quality depend on physical position on
the stand? The frame-top tanks are a concern because that region of the video is
prone to overexposure, which could in principle bias detection. This script asks
whether the Mathematica output shows any positional weakness, or whether it
handled all stand positions comparably well.

Everything is keyed to the physical tank position (1-8), which is a FIXED slot
on the stand (only fish assignment is randomized across sessions). The physical
layout (which slot is frame-top vs frame-bottom, left vs right column) is
recovered from the per-session tank polygons in P1-BatchSettingsCombined-v5.csv.

Works from the SOURCE file so that the 35 QC-excluded sessions are included in
the positional breakdown (the merged file has already dropped them).

Outputs:
  output/tracking_quality_by_position.csv   per-session metrics + layout
  output/tracking_quality_summary.txt       printed tables
  figures/06a_quality_by_position.png        detection rate / extremes by slot
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT      = Path(__file__).resolve().parents[1]
SOURCE    = ROOT / "data" / "source" / "tracking_triallevel_distortion_corrected.csv"
SETTINGS  = ROOT / "data" / "source" / "P1-BatchSettingsCombined-v5.csv"
INCLUDE_CSV = ROOT / "data" / "source" / "tracking_aggregate_uncorrected.csv"
OUT_DIR   = ROOT / "output"
FIG_DIR   = ROOT / "figures"

CUTOFF_FRAME = 1805
MAX_FRAME    = 37805
SAMPLE_EVERY = 5
EXPECTED_SAMPLES = (MAX_FRAME - CUTOFF_FRAME) // SAMPLE_EVERY  # ~7200
Y_MID = 61.0


def physical_layout() -> pd.DataFrame:
    """Mean polygon centroid per physical tank slot (1-8), across all sessions.

    The stand is fixed, so averaging the per-session polygons gives a stable
    picture of where each slot sits in the camera frame. Returns a DataFrame
    indexed by tank_pos with centroid_x/centroid_y, frame_row (1=top..4=bottom),
    and column (L/R).
    """
    cfg = pd.read_csv(SETTINGS, header=None)
    cx = {t: [] for t in range(1, 9)}
    cy = {t: [] for t in range(1, 9)}
    for _, row in cfg.iterrows():
        for t in range(8):
            base = 2 + t * 8
            xs = [float(row[base + i * 2]) for i in range(4)]
            ys = [float(row[base + i * 2 + 1]) for i in range(4)]
            cx[t + 1].append(np.mean(xs))
            cy[t + 1].append(np.mean(ys))
    lay = pd.DataFrame({
        "tank_pos": list(range(1, 9)),
        "centroid_x": [np.mean(cx[t]) for t in range(1, 9)],
        "centroid_y": [np.mean(cy[t]) for t in range(1, 9)],
    })
    # column: split on overall median x
    xmid = lay["centroid_x"].median()
    lay["column"] = np.where(lay["centroid_x"] < xmid, "L", "R")
    # frame_row: rank by centroid_y within each column (1 = top of frame)
    lay["frame_row"] = (
        lay.groupby("column")["centroid_y"].rank(method="dense").astype(int)
    )
    # half: frame-top half (rows 1-2) vs frame-bottom half (rows 3-4)
    lay["frame_half"] = np.where(lay["frame_row"] <= 2, "top_half", "bottom_half")
    return lay


def load_exclude_set() -> set:
    df = pd.read_csv(INCLUDE_CSV)
    df = df[df["include"] == 0].copy()
    df["batch"]      = df["session"].astype(int)
    df["run_within"] = df["tank"].astype(int) // 10
    df["fish_id"]    = df["fish"].astype(int)
    return set(zip(df["fish_id"], df["batch"], df["run_within"]))


def session_metrics() -> pd.DataFrame:
    print(f"Reading {SOURCE.name} …")
    df = pd.read_csv(SOURCE, usecols=["run", "dayafternoon", "tank", "fishid", "framenumber", "x", "y"])
    df = df[(df["framenumber"] >= CUTOFF_FRAME) & (df["framenumber"] < MAX_FRAME)].copy()

    df["batch"]    = df["run"].astype(int)
    df["run_within"] = (df["dayafternoon"].astype(int) // 10)
    df["tank_pos"] = (df["dayafternoon"].astype(int) % 10)
    df["fish_id"]  = df["fishid"].astype(int)
    df["detected"] = df["x"].notna()
    df["in_top"]   = df["y"] >= Y_MID

    g = df.groupby(["batch", "run_within", "tank_pos", "fish_id"])
    rows = []
    for (batch, run_w, pos, fid), sub in g:
        det = sub["detected"]
        n_samp = len(sub)
        n_det = int(det.sum())
        det_rate = n_det / EXPECTED_SAMPLES
        xy = sub.loc[det, ["x", "y"]]
        if n_det > 0:
            x_sd = xy["x"].std()
            y_sd = xy["y"].std()
            pct_top = float(sub.loc[det, "in_top"].mean())
        else:
            x_sd = y_sd = pct_top = np.nan
        rows.append(dict(
            batch=batch, run_within=run_w, tank_pos=pos, fish_id=fid,
            n_samp=n_samp, n_det=n_det, det_rate=det_rate,
            x_sd=x_sd, y_sd=y_sd, pct_top=pct_top,
        ))
    out = pd.DataFrame(rows)

    excl = load_exclude_set()
    out["excluded"] = [
        (r.fish_id, r.batch, r.run_within) in excl for r in out.itertuples()
    ]
    # static-artifact flag: well-detected but spatially pinned (tiny SD)
    out["static_flag"] = (out["det_rate"] > 0.30) & (out["y_sd"] < 5)
    # extreme behavioural value (all-top or all-bottom) — suspicious if also low n_det
    out["extreme_pct"] = (out["pct_top"] <= 0.02) | (out["pct_top"] >= 0.98)
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    lay = physical_layout()
    sess = session_metrics()
    sess = sess.merge(lay, on="tank_pos", how="left")
    sess.to_csv(OUT_DIR / "tracking_quality_by_position.csv", index=False)

    lines = []
    def P(s=""):
        print(s); lines.append(s)

    P("=" * 70)
    P("PHYSICAL STAND LAYOUT (recovered from tank polygons)")
    P("=" * 70)
    P(lay.round(1).to_string(index=False))
    P("")
    P("  frame_row 1 = top of camera frame (the slots most prone to")
    P("  overexposure). frame_row 4 = bottom of frame.")
    P("")

    inc = sess[~sess["excluded"]]

    P("=" * 70)
    P("1. DETECTION RATE BY PHYSICAL SLOT  (all 624 sessions)")
    P("=" * 70)
    tab = (sess.groupby(["tank_pos", "column", "frame_row"])
              .agg(n_sess=("det_rate", "size"),
                   det_rate_mean=("det_rate", "mean"),
                   det_rate_med=("det_rate", "median"),
                   n_excluded=("excluded", "sum"))
              .reset_index()
              .sort_values(["column", "frame_row"]))
    P(tab.round(3).to_string(index=False))
    P("")

    P("=" * 70)
    P("2. DETECTION RATE BY FRAME ROW / HALF  (included sessions only)")
    P("=" * 70)
    for col, lab in [("frame_row", "frame row (1=top..4=bottom)"), ("frame_half", "frame half")]:
        t = (inc.groupby(col)
                .agg(n=("det_rate", "size"),
                     det_rate_mean=("det_rate", "mean"),
                     det_rate_med=("det_rate", "median"),
                     y_sd_mean=("y_sd", "mean"))
                .reset_index())
        P(f"  by {lab}:")
        P(t.round(3).to_string(index=False))
        P("")

    P("=" * 70)
    P("3. WHERE DO THE QC-EXCLUDED SESSIONS FALL?")
    P("=" * 70)
    ex = sess[sess["excluded"]]
    P(f"  {len(ex)} excluded sessions total")
    t = (ex.groupby(["column", "frame_row"]).size()
            .rename("n_excluded").reset_index())
    P(t.to_string(index=False))
    P("")
    P("  excluded-session rate by frame half:")
    t = (sess.groupby("frame_half")["excluded"].agg(["sum", "size", "mean"])
            .rename(columns={"sum": "n_excl", "size": "n_sess", "mean": "excl_rate"}))
    P(t.round(3).to_string())
    P("")

    P("=" * 70)
    P("4. SUSPECT-TRACKING SIGNATURES BY SLOT  (all sessions)")
    P("=" * 70)
    P("  static_flag = det_rate>0.30 AND y_sd<5  (pinned like an artifact)")
    P("  extreme_pct = pct_top <=.02 or >=.98")
    t = (sess.groupby(["tank_pos", "frame_row"])
            .agg(n=("static_flag", "size"),
                 n_static=("static_flag", "sum"),
                 n_extreme=("extreme_pct", "sum"),
                 y_sd_mean=("y_sd", "mean"))
            .reset_index().sort_values("frame_row"))
    P(t.round(2).to_string(index=False))
    P("")
    P(f"  TOTAL static-artifact-like sessions: {int(sess['static_flag'].sum())}")
    P(f"    of which already QC-excluded:      {int((sess['static_flag'] & sess['excluded']).sum())}")
    P("")

    P("=" * 70)
    P("5. INCLUDED-SESSION QUALITY: is it uniform across slots?")
    P("=" * 70)
    t = (inc.groupby("tank_pos")
            .agg(n=("det_rate", "size"),
                 det_rate_mean=("det_rate", "mean"),
                 det_rate_min=("det_rate", "min"),
                 y_sd_mean=("y_sd", "mean"),
                 n_static=("static_flag", "sum"))
            .reset_index())
    P(t.round(3).to_string(index=False))
    P("")
    cv = inc.groupby("tank_pos")["det_rate"].mean()
    P(f"  Across-slot detection-rate range (included): "
      f"{cv.min():.3f} – {cv.max():.3f}  (spread {cv.max()-cv.min():.3f})")
    P("")

    (OUT_DIR / "tracking_quality_summary.txt").write_text("\n".join(lines))

    # ---- figure ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    order = lay.sort_values(["column", "frame_row"])["tank_pos"].tolist()

    m = sess.groupby("tank_pos")["det_rate"].mean().reindex(order)
    mi = inc.groupby("tank_pos")["det_rate"].mean().reindex(order)
    ax = axes[0]
    xs = np.arange(len(order))
    ax.bar(xs - 0.2, m.values, 0.4, label="all sessions", color="#c0392b")
    ax.bar(xs + 0.2, mi.values, 0.4, label="included only", color="#2980b9")
    ax.set_xticks(xs); ax.set_xticklabels([f"T{t}" for t in order])
    ax.set_ylabel("mean detection rate"); ax.set_title("Detection rate by slot")
    ax.axhline(inc["det_rate"].mean(), ls="--", c="k", lw=0.8)
    ax.legend(fontsize=8)

    ax = axes[1]
    exr = sess.groupby("tank_pos")["excluded"].mean().reindex(order)
    ax.bar(xs, exr.values, color="#8e44ad")
    ax.set_xticks(xs); ax.set_xticklabels([f"T{t}" for t in order])
    ax.set_ylabel("fraction QC-excluded"); ax.set_title("Exclusion rate by slot")

    ax = axes[2]
    for col, c in [("L", "#16a085"), ("R", "#d35400")]:
        d = inc[inc["column"] == col]
        gg = d.groupby("frame_row")["det_rate"].mean()
        ax.plot(gg.index, gg.values, "o-", color=c, label=f"col {col}")
    ax.set_xlabel("frame row (1=top → 4=bottom)")
    ax.set_ylabel("mean detection rate (included)")
    ax.set_title("Detection rate vs frame row")
    ax.set_xticks([1, 2, 3, 4]); ax.legend(fontsize=8)

    fig.suptitle("Mathematica tracking quality by physical stand position", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "06a_quality_by_position.png", dpi=130)
    P(f"Saved figure → {FIG_DIR / '06a_quality_by_position.png'}")


if __name__ == "__main__":
    main()
