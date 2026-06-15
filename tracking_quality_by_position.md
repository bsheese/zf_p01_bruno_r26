# Mathematica Tracking Quality by Stand Position

**Date:** 2026-06-10
**Script:** `src/06_tracking_quality.py` → `output/tracking_quality_by_position.csv`, `output/tracking_quality_summary.txt`, `figures/06a_quality_by_position.png`
**Question:** Does the original Mathematica tracking quality depend on physical position on the 8-tank stand? In particular, do the frame-top slots — where the video is most prone to overexposure — show degraded tracking?

## TL;DR

**The Mathematica tracking is sound across all stand positions.** The failure modes one might fear in the overexposure-prone frame-top slots — tanks rendered untrackable, or detections pinned to static lighting artifacts — **do not appear in the Mathematica output.** Detection is healthy at every slot, the static-artifact failure mode is essentially absent (1 of 624 sessions), and — the decisive test — `pct_top` is *not* depressed in the frame-top slots the overexposure hypothesis predicted would fail. There is one genuine but minor residual (a small frame-row gradient in `pct_top`) that randomization largely washes out.

## Method

Everything is keyed to **physical `tank_pos` (1–8)**, a fixed slot on the stand — only fish assignment is randomized across the 6 runs. The physical layout was recovered from the per-session tank polygons in `P1-BatchSettingsCombined-v5.csv` (mean centroid per slot):

| Slot | Column | Frame row | Position in camera frame |
|------|--------|-----------|--------------------------|
| 1 | L | 1 | top-left (overexposure suspect) |
| 2 | L | 2 | upper-mid-left |
| 3 | L | 3 | lower-mid-left |
| 4 | L | 4 | bottom-left |
| 5 | R | 1 | top-right (overexposure suspect) |
| 6 | R | 2 | upper-mid-right |
| 7 | R | 3 | lower-mid-right |
| 8 | R | 4 | bottom-right |

Frame row 1 = top of the camera frame — the region most prone to overexposure (in the worst tanks, up to ~68% of top-half pixels saturated), which is where detection would be expected to fail if position mattered. Analysis runs on the **source** file (all 624 sessions, including the 35 QC-excluded), so positional weakness can't hide inside the exclusions.

## Findings

### 1. No catastrophic positional failure

Mean detection rate by slot (included sessions): **0.55–0.73**, healthy everywhere. The frame-top slots predicted to fail are *among the better ones* — T5 (top-right) is 0.71, the second-highest of all eight. T1 (top-left) is the weakest at 0.60 but nowhere near failure.

Detection rate by frame row (included): 0.65 / 0.72 / 0.69 / 0.61 (rows 1→4). No row drops out.

### 2. The static-artifact failure mode is essentially absent

A failing tracker can pin detections to a fixed lighting artifact — a static bright spot — yielding many detections clustered at one location with near-zero positional variance. Flagging that signature (det_rate > 0.30 **and** y_sd < 5) across all 624 Mathematica sessions returns **1 session** — and it's borderline, not even one the original QC excluded. Per-slot `y_sd` is a healthy 23–26 units everywhere: fish are moving, not pinned. The frame-top slots' y-distributions are well spread (T1 mean y=66, SD=26; T5 mean y=70, SD=28, with 45% of detections high in the tank) — the opposite of an artifact.

### 3. The "top-row tanks fail" claim was batch-1-specific

`manuscript/methods.md` noted that in batch-1/run-1 the top-row tanks (1,2,5,6) were all excluded — read at the time as a stand-position weakness. Across the full dataset it isn't one: the 35 QC exclusions split **evenly** between frame halves (18 top-half vs 17 bottom-half; exclusion rate 5.8% vs 5.4%). The batch-1 event was a one-session lighting problem, not a systematic property of those slots.

### 4. The overexposure-bias prediction is refuted (the decisive test)

Because fish are randomized to slots, **true `pct_top` must be independent of physical slot** — any slot effect is a tracking artifact. Residual overexposure would specifically kill *top-half* (water-surface) detections in the frame-top slots, depressing their `pct_top`. Observed:

| Frame row | `pct_top` |
|-----------|-----------|
| 1 (top of frame) | **0.563** |
| 2 | 0.542 |
| 3 | 0.493 |
| 4 (bottom of frame) | 0.477 |

`pct_top` is **highest** in the frame-top slots — the exact opposite of the overexposure prediction. The per-session image-adjustment tuning in the Mathematica pipeline evidently handled the saturation in those slots. The slot-level correlation between detection rate and `pct_top` is ~0 (r=0.17, p=0.70), i.e. no detection-driven bias across positions.

### 5. One genuine residual: a small frame-row gradient in `pct_top`

There *is* a statistically detectable slot effect (one-way ANOVA `pct_top ~ slot`: F=2.31, p=0.025): `pct_top` declines monotonically from frame-top to frame-bottom (0.563 → 0.477, a ~8.6-point spread). Since fish are randomized, this is a tracking/distortion-correction residual, not behavior — most likely the fixed `y_mid = 61` threshold not perfectly matching the true physical midpoint across frame rows after the trapezoidal perspective correction. A mild left/right detection asymmetry exists too (L=0.62 vs R=0.71).

**Why it doesn't undermine the results:**
- The gradient (0.086) is ~48% of the between-fish SD (0.177) — present but smaller than the signal.
- Randomization is strong: each fish visits **3.4 of 4 frame rows** across its 6 runs, and 95/103 fish occupy both columns. So the gradient **averages out of per-fish means** rather than biasing them — it adds noise, not bias.
- That added noise is a plausible *partial* contributor to the reliability gap (α = 0.855 vs the original analysis's 0.88): a slot-randomization residual attenuates test-retest correlation slightly.

### 6. Removing the slot artifact recovers ~half the reliability gap

Robustness check (`src/07_slot_adjusted_reliability.py` → `output/slot_adjusted_reliability.txt`). Using session-level `pct_top` from raw detected frames (no interpolation), the slot artifact was regressed out two ways — by frame row (4 levels) and by full physical slot (8 levels) — and Cronbach's alpha recomputed on the fish × 6-run matrix (complete cases, N=82):

| `pct_top` version | Cronbach's α | mean inter-run r |
|-------------------|--------------|------------------|
| raw | 0.852 | 0.497 |
| frame-row adjusted | **0.866** | 0.525 |
| full per-slot adjusted | **0.866** | 0.525 |

Removing the position residual lifts α from 0.852 to **0.866** and mean inter-run r from 0.497 to 0.525 — closing roughly **half** the gap to the original analysis's 0.88. Physical position explains only ~2.4–2.7% of session-level `pct_top` variance (slot F p=0.025; frame-row F p=0.003), but because that variance is pure noise with respect to the fish's trait, removing it still measurably tightens reliability. The remaining gap to 0.88 is most plausibly the interpolation difference (the original analysis's SPSS carry-forward vs. our linear interpolation).

**Now folded into the main pipeline.** `03_features.py` carries `tank_pos` through to the per-minute features; `04_stats.py` derives `frame_row`, builds slot-adjusted `pct_top` columns (`add_slot_adjustment`, frame-row groups), writes `output/run_table_slot_adjusted.csv`, and reports adjusted alpha beside the raw value in Section 1 of `stats_main.txt`. On the pipeline's **interpolated** `pct_top`, slot adjustment lifts full-window α from **0.855 → 0.874** (95% CI [0.826, 0.912], which includes 0.88) and first-10-min α to 0.879 — i.e. interpolation + slot adjustment together converge with the original analysis's 0.88. Raw `pct_top` remains the primary measure for all downstream analyses (strain, habituation, summary are unchanged); the adjusted value is reported as a reliability robustness result.

Two confirmations that the artifact is benign at the trait level:
- Per-fish mean `pct_top` is **uncorrelated** with the slots a fish occupied: vs. mean frame row r=+0.084 (p=0.40); vs. fraction-right-column r=+0.033 (p=0.74). Randomization neutralizes the gradient for per-fish estimates.
- Frame-row and full-slot adjustment give **identical** α (0.866), i.e. essentially all the recoverable artifact is the frame-row gradient; the L/R asymmetry adds nothing once the row gradient is removed.

## Implications

1. **The decision to base the reanalysis on the Mathematica tracking output is well-justified.** The tracking is sound and roughly uniform across all stand positions; it does not carry the overexposure/static-artifact pathologies the frame-top slots might have been expected to show.
2. **The residual frame-row gradient is real but removable, and reclaims ~half the reliability gap** (Section 6). Slot-adjusting `pct_top` lifts α from 0.852 to 0.866; per-fish trait estimates are unbiased by slot. This strengthens the replication and points to interpolation as the likely source of the remaining gap to 0.88.
3. **The frame-differencing "freezing is invisible" limitation is unaffected by this** — it is a property of the algorithm, not of position, and applies equally everywhere.

## Caveats

- Detection *rate* and spatial spread are proxies for tracking quality, not ground-truth accuracy (no hand-labeled positions exist). The strongest evidence here is indirect-but-pointed: the randomization-based `pct_top` test, which would expose systematic mis-tracking regardless of detection counts.
- The 1 static-flagged session and the low-end included sessions (det_rate as low as 0.03 at T5) are individual QC items, not positional patterns.
