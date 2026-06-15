#!/usr/bin/env python3
"""
01_ingest.py — Ingest distortion-corrected Mathematica tracking data

Reads the distortion-corrected frame-by-frame tracking output from the
original Mathematica pipeline and writes tracking_merged.csv, ready for
03_features.py.  This is step 1 of the analysis pipeline; the reanalysis
starts from the Mathematica distortion-corrected tracking output.

Input
-----
  p01bruno_r26/data/source/tracking_triallevel_distortion_corrected.csv

Output
------
  p01bruno_r26/data/merged/tracking_merged.csv

Source-file column encoding
---------------------------
  run          batch number (1–13)
  dayafternoon run_within_batch × 10 + physical_tank_position
               e.g. 13 → run_within=1, physical_tank=3
  tank         fish number within batch (= fishid % 10)
  fishid       fish ID = batch × 10 + fish_within_batch
  framenumber  video frame; 0 = sentinel row; 1805 = first real sample
               (60 s cutoff at 29.97 fps); step = 5 frames (~6 fps)
  x, y         distortion-corrected scaled coordinates
               y = 0 at tank floor, increases upward; NaN = no detection
  deltax/y     signed frame-to-frame displacement (not used downstream)

Top/bottom convention
---------------------
  y_mid = 61 for every row (fixed threshold after distortion correction)
  bottom = y < 61   (closer to floor)
  top    = y >= 61  (closer to water surface)

  This is INVERTED relative to raw-video pixel convention, where
  at_bottom = y > y_mid.  03_features.py is updated accordingly.
"""

import pandas as pd
import numpy as np
from pathlib import Path

ROOT      = Path(__file__).resolve().parents[1]
SOURCE    = ROOT / "data" / "source" / "tracking_triallevel_distortion_corrected.csv"
INCLUDE_CSV = ROOT / "data" / "source" / "tracking_aggregate_uncorrected.csv"
DECODE_XL = ROOT / "data" / "source" / "P1-Tank Randomization Decoding v2.xlsx"
OUT_DIR   = ROOT / "data" / "merged"

FPS          = 29.97
CUTOFF_FRAME = 1805    # first 60 s excluded (experimenter visible)
MAX_FRAME    = 37805   # 20 min post-cutoff (matches original SPSS: framenumber < 37805)
Y_MID        = 61.0    # fixed top/bottom threshold in distortion-corrected coords


def load_exclude_set() -> set:
    """
    Return a set of (fish_id, batch, run_within_batch) tuples flagged include=0
    in the original v5 quality-control CSV.

    The v5 'tank' column encodes run_within_batch*10 + tank_pos; we only need
    (fish_id, batch, run_within_batch) to match against the merged data.

    35 sessions are excluded: 7 missing-video NaNs, 6 empty-tank (fish 138),
    4 sick-fish (fish 116), and 18 sessions with catastrophically low detection
    counts (all below p10 of included sessions; most under 5% of the median).
    Rationale and full session list: p01bruno_r26/manuscript/methods.md
    """
    df = pd.read_csv(INCLUDE_CSV)
    df = df[df["include"] == 0].copy()
    df["batch"]       = df["session"].astype(int)
    df["run_within"]  = df["tank"].astype(int) // 10
    df["fish_id"]     = df["fish"].astype(int)
    return set(zip(df["fish_id"], df["batch"], df["run_within"]))


def load_strain_map() -> dict:
    """Return {batch: strain_label} from the tank randomization workbook."""
    df = pd.read_excel(DECODE_XL, sheet_name="Sheet1")
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={"Batch": "batch", "File name": "file_stem"})
    df["batch"] = df["batch"].astype(int)

    def _strain(stem: str) -> str:
        s = str(stem).upper()
        if "5GWT" in s:
            return "5G"
        if "1GAB" in s or "AB" in s:
            return "AB"
        return "unknown"

    strain_map: dict = {}
    for _, row in df.iterrows():
        b = int(row["batch"])
        if b not in strain_map:
            strain_map[b] = _strain(row["file_stem"])
    return strain_map


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Reading {SOURCE.name} …")
    df = pd.read_csv(SOURCE)
    print(f"  {len(df):,} rows in source file")

    # Drop sentinel rows and enforce analysis window (matches SPSS: 1805 ≤ frame < 37805)
    df = df[df["framenumber"].notna() & (df["framenumber"] > 0) & (df["framenumber"] < MAX_FRAME)].copy()
    print(f"  {len(df):,} rows within analysis window (frames 1805–{MAX_FRAME})")

    # Decode batch/run early so we can apply the include filter
    df["batch"]     = df["run"].astype(int)
    df["run"]       = (df["dayafternoon"].astype(int) // 10)
    df["fish_id"]   = df["fishid"].astype(int)

    # Apply original quality-control exclusions from v5 CSV (include=0).
    # See p01bruno_r26/manuscript/methods.md for the full breakdown and justification.
    exclude_set = load_exclude_set()
    before = len(df)
    key = list(zip(df["fish_id"], df["batch"], df["run"]))
    df = df[[k not in exclude_set for k in key]].copy()
    n_excluded_rows = before - len(df)
    n_excluded_sessions = len(exclude_set)
    print(f"  Excluded {n_excluded_sessions} flagged sessions → dropped {n_excluded_rows:,} rows")

    # -----------------------------------------------------------------------
    # Remaining column decodes (batch, run, fish_id already set above)
    # -----------------------------------------------------------------------
    df["tank_pos"]   = (df["dayafternoon"].astype(int) % 10)    # physical tank position 1-8
    df["fish"]       = df["tank"].astype(int)                   # fish within batch 1-8
    df["frame"]      = df["framenumber"].astype(int)

    # Session structure: 6 runs = 3 days × 2 sessions (AM ≈ 11am, PM ≈ 2:30pm)
    df["day"]        = (df["run"] - 1) // 2 + 1                # 1, 1, 2, 2, 3, 3
    df["time_of_day"] = df["run"].apply(lambda r: "AM" if r % 2 == 1 else "PM")

    # row_idx for compatibility with settings CSV indexing: (batch-1)*6 + (run-1)
    df["row_idx"] = (df["batch"] - 1) * 6 + (df["run"] - 1)

    # Time relative to start of analysis window
    df["time_s"]   = (df["frame"] - CUTOFF_FRAME) / FPS
    df["time_min"] = df["time_s"] / 60.0

    # Fixed top/bottom midpoint (same for all tanks after distortion correction)
    df["y_mid"] = Y_MID

    # Strain lookup
    strain_map = load_strain_map()
    df["strain"] = df["batch"].map(strain_map).fillna("unknown")

    # -----------------------------------------------------------------------
    # Select and order output columns
    # -----------------------------------------------------------------------
    out = df[[
        "row_idx", "batch", "run", "day", "time_of_day", "fish", "fish_id", "strain",
        "frame", "time_s", "time_min", "tank_pos", "x", "y", "y_mid",
    ]].sort_values(["batch", "run", "tank_pos", "frame"]).reset_index(drop=True)

    out_path = OUT_DIR / "tracking_merged.csv"
    out.to_csv(out_path, index=False)

    n_fish    = out["fish_id"].nunique()
    n_batches = out["batch"].nunique()
    n_records = out.groupby(["batch", "run"]).ngroups
    det_rate  = out["x"].notna().mean()

    print(f"\nSaved → {out_path}")
    print(f"  {n_batches} batches, {n_fish} unique fish, {n_records} recordings")
    print(f"  {len(out):,} rows total, detection rate: {100 * det_rate:.1f}%")
    print(f"  Batches present: {sorted(out['batch'].unique().tolist())}")
    print(f"  Runs present:    {sorted(out['run'].unique().tolist())}")


if __name__ == "__main__":
    main()
