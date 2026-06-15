#!/usr/bin/env python3
"""
03_features.py — Interpolate trajectories and compute behavioral features

Takes p01bruno_r26/data/merged/tracking_merged.csv and produces:

  tracking_interp.csv   — frame-level data with interpolated positions
  features_by_minute.csv — per-fish × per-minute aggregate features
  features_summary.csv  — per-fish summary (whole-recording means)

Features computed
-----------------
  at frame level (tracking_interp.csv):
    x_interp, y_interp   — linearly interpolated positions (NaN when impossible)
    detected             — 1 if original tracking fired, 0 if interpolated
    at_bottom            — 1 if y_interp > y_mid (bottom half of tank)
    velocity             — Euclidean distance (px/sample) between consecutive frames
                           (NaN if either position is missing)

  per-minute (features_by_minute.csv):
    pct_bottom       — fraction of interpolated frames in bottom half
    mean_velocity    — mean velocity across detected frames in that minute
    n_detected       — raw count of detection events (non-zero original)
    n_frames         — total sampled frames in that minute

  per-fish summary (features_summary.csv):
    All per-minute features averaged over all 20 minutes (minute 0 = first
    full minute after the 60 s cutoff)

Notes on interpolation
----------------------
  Carry-forward is NOT used.  Instead, only frames that are "bracketed" by
  detections on both sides are linearly interpolated.  Isolated runs with no
  surrounding detections remain NaN.

  Velocity is computed only between consecutive detected-or-interpolated frames
  that are both non-NaN.

  The known Mathematica velocity sign bug (averaging signed Δx, Δy → ≈ 0) is
  fixed: velocity = √(Δx² + Δy²).
"""

import numpy as np
import pandas as pd
from pathlib import Path


ROOT    = Path(__file__).resolve().parents[1]
MERGED  = ROOT / "data" / "merged" / "tracking_merged.csv"
OUT_DIR = ROOT / "data" / "features"

SAMPLE_EVERY = 5       # frames between tracked samples
FPS          = 29.97
SAMPLE_FPS   = FPS / SAMPLE_EVERY   # ~5.994 samples per second
MIN_PER_BIN  = 1.0     # aggregate over 1-minute windows


def interpolate_series(x_arr: np.ndarray) -> np.ndarray:
    """
    Linearly interpolate a series where NaN means "not detected."
    Returns array with NaN where interpolation is not possible (leading/trailing gaps).
    """
    x = x_arr.astype(float).copy()

    not_nan = np.where(~np.isnan(x))[0]
    if len(not_nan) < 2:
        return x   # can't interpolate with fewer than 2 anchors

    # Only interpolate interior gaps; leave leading/trailing NaN as is
    first, last = not_nan[0], not_nan[-1]
    # Create interpolated version using pandas (handles non-uniform spacing)
    s = pd.Series(x)
    s_interp = s.interpolate(method="index", limit_area="inside")
    x[first:last+1] = s_interp.values[first:last+1]
    return x


def compute_velocity(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Euclidean distance between consecutive valid positions.
    Result[i] = distance from sample i-1 to sample i.
    NaN when either position is missing.
    """
    vel = np.full(len(x), np.nan)
    for i in range(1, len(x)):
        if not (np.isnan(x[i]) or np.isnan(y[i]) or
                np.isnan(x[i-1]) or np.isnan(y[i-1])):
            vel[i] = np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
    return vel


def process_track(sub: pd.DataFrame) -> pd.DataFrame:
    """
    Interpolate and compute velocity for one fish × one recording (run).
    sub must be sorted by frame.
    """
    sub = sub.sort_values("frame").copy()

    xi = interpolate_series(sub["x"].values)
    yi = interpolate_series(sub["y"].values)

    sub["x_interp"] = xi
    sub["y_interp"] = yi
    sub["detected"] = sub["x"].notna().astype(int)
    # y = 0 at tank floor, increases upward; bottom = y < y_mid
    sub["at_bottom"] = (sub["y_interp"] < sub["y_mid"]).astype(float)
    sub.loc[sub["y_interp"].isna(), "at_bottom"] = np.nan

    sub["velocity"] = compute_velocity(xi, yi)

    return sub


def compute_session_metrics(interp: pd.DataFrame, y_mid: float) -> dict:
    """
    Session-level metrics computed from detected frames only.

    latency_to_top : minutes from analysis window start to first detected
                     frame with y_interp >= y_mid.  NaN if fish never
                     entered the top half during the session.
    n_transitions  : number of top↔bottom zone crossings counted across
                     consecutive detected frames.  Undercount relative to
                     true crossings because undetected periods are invisible;
                     consistent across fish/sessions.
    """
    det = interp[interp["detected"] == 1].dropna(subset=["y_interp", "time_min"])

    # latency
    top = det[det["y_interp"] >= y_mid]
    latency_to_top = float(top["time_min"].min()) if len(top) > 0 else np.nan

    # transitions
    if len(det) < 2:
        n_transitions = 0
    else:
        in_top = (det["y_interp"].values >= y_mid).astype(int)
        n_transitions = int(np.abs(np.diff(in_top)).sum())

    return {"latency_to_top": latency_to_top, "n_transitions": n_transitions}


def aggregate_by_minute(sub: pd.DataFrame) -> pd.DataFrame:
    """Aggregate frame-level features into 1-minute bins."""
    sub = sub.copy()
    sub["minute"] = (sub["time_min"] // MIN_PER_BIN).astype(int)

    rows = []
    for minute, grp in sub.groupby("minute"):
        rows.append({
            "minute":        minute,
            "pct_bottom":    grp["at_bottom"].mean(),   # NaN-safe mean
            "mean_velocity": grp["velocity"].mean(),
            "n_detected":    grp["detected"].sum(),
            "n_frames":      len(grp),
        })
    return pd.DataFrame(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Reading {MERGED} …")
    df = pd.read_csv(MERGED)
    print(f"  {len(df)} rows, {df['fish_id'].nunique()} unique fish, "
          f"{df['row_idx'].nunique()} recordings")

    interp_chunks   = []
    min_chunks      = []
    session_chunks  = []

    groups = list(df.groupby(["row_idx", "tank_pos"]))
    total  = len(groups)
    print(f"  Processing {total} fish×recording groups …")

    for i, ((row_idx, tank_pos), grp) in enumerate(groups):
        if i % 100 == 0:
            print(f"    {i}/{total}")

        meta_cols = ["row_idx","batch","run","day","time_of_day","fish","fish_id","strain","y_mid","tank_pos"]
        interp    = process_track(grp)
        interp_chunks.append(interp)

        y_mid_val = float(grp["y_mid"].iloc[0])

        # per-minute aggregation
        by_min = aggregate_by_minute(interp)
        for col in meta_cols:
            by_min[col] = grp[col].iloc[0]
        min_chunks.append(by_min)

        # session-level metrics
        sess = compute_session_metrics(interp, y_mid_val)
        for col in meta_cols:
            sess[col] = grp[col].iloc[0]
        session_chunks.append(sess)

    # -----------------------------------------------------------------
    # Save frame-level interpolated tracks
    # -----------------------------------------------------------------
    df_interp = pd.concat(interp_chunks, ignore_index=True)
    out_interp = OUT_DIR / "tracking_interp.csv"
    df_interp.to_csv(out_interp, index=False)
    print(f"\nSaved frame-level interpolated data → {out_interp}")

    # -----------------------------------------------------------------
    # Save minute-by-minute features
    # -----------------------------------------------------------------
    df_min = pd.concat(min_chunks, ignore_index=True)
    col_order = [
        "fish_id","fish","batch","run","day","time_of_day","strain","minute",
        "pct_bottom","mean_velocity","n_detected","n_frames","y_mid","row_idx","tank_pos",
    ]
    df_min = df_min[[c for c in col_order if c in df_min.columns]]
    out_min = OUT_DIR / "features_by_minute.csv"
    df_min.to_csv(out_min, index=False)
    print(f"Saved minute-by-minute features   → {out_min}")

    # -----------------------------------------------------------------
    # Save per-fish summary (mean across all minutes)
    # -----------------------------------------------------------------
    summary = (
        df_min.groupby(["fish_id","fish","batch","strain"])
        .agg(  # day/time_of_day vary within fish; not aggregated at summary level
            pct_bottom_mean   = ("pct_bottom",    "mean"),
            pct_bottom_sd     = ("pct_bottom",    "std"),
            velocity_mean     = ("mean_velocity", "mean"),
            velocity_sd       = ("mean_velocity", "std"),
            n_detected_total  = ("n_detected",    "sum"),
            n_frames_total    = ("n_frames",       "sum"),
        )
        .reset_index()
    )
    out_sum = OUT_DIR / "features_summary.csv"
    summary.to_csv(out_sum, index=False)
    print(f"Saved per-fish summary            → {out_sum}")

    # -----------------------------------------------------------------
    # Save session-level metrics (latency to top, n_transitions)
    # -----------------------------------------------------------------
    sess_col_order = [
        "fish_id","fish","batch","run","day","time_of_day","strain",
        "latency_to_top","n_transitions","y_mid","row_idx","tank_pos",
    ]
    df_sess = pd.DataFrame(session_chunks)
    df_sess = df_sess[[c for c in sess_col_order if c in df_sess.columns]]
    out_sess = OUT_DIR / "features_sessions.csv"
    df_sess.to_csv(out_sess, index=False)
    print(f"Saved session-level metrics       → {out_sess}")

    # Quick stats
    print(f"\nQuick summary:")
    print(f"  Mean % bottom (all fish): {100*df_min['pct_bottom'].mean():.1f}%")
    print(f"  Mean velocity: {df_min['mean_velocity'].mean():.2f} px/sample")
    strain_stats = df_min.groupby("strain")["pct_bottom"].mean()
    print(f"  % bottom by strain: {strain_stats.round(3).to_dict()}")


if __name__ == "__main__":
    main()
