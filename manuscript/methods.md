# Methods

## Strain is not analyzed (confounded with batch)

The dataset contains two strains (5G inbred, AB wild type), but each of the 13 testing batches was strain-pure (5G: batches 1, 3, 5, 7, 9, 12; AB: 2, 4, 6, 8, 10, 11, 13; verified by `src/11_qc_checks.py`). There is no within-batch strain variation, so strain is completely confounded with batch and with everything that varied between batches — testing date (January–March 2014), room and water conditions, handler, and per-batch distortion-correction calibration. A "5G vs AB" contrast would be identical to "the batches that ran 5G vs those that ran AB" and could not be attributed to genotype. We therefore do not analyze strain. The batch variance component reported in the reliability model below is a separate, nuisance term and does not bear on this point.

## Fish Identity Across Sessions

### Physical tank randomization

Each batch consisted of 8 fish tested simultaneously in an 8-tank apparatus over 3 consecutive days (2 sessions per day, AM and PM), yielding 6 sessions per fish. Within each session, fish were assigned to physical tank positions 1–8. **The assignment was re-randomized for each session** so that no fish occupied the same physical tank position across all 6 runs. This counterbalances any tank-specific effects (position in the camera frame, lighting gradients, minor water quality differences) on pct_top estimates.

The tank-to-fish assignment for every session of every batch is recorded in the tank randomization decode workbook:

**`p01bruno_r26/data/source/tank_randomization_decode.csv`** (derived from `P1-Tank Randomization Decoding v2.xlsx`)

This workbook lists, for each session (identified by video filename and date), which physical tank position (1–8) contained which fish (identified by fish_within_batch number 1–8). It is the authoritative source for linking tracking data — which records position by physical tank slot — to individual fish identities across sessions.

In the pipeline, `src/01_ingest.py` reads fish identity directly from the `fishid` column of the distortion-corrected source file, which was populated by Bruno's Mathematica pipeline using this same decode workbook. The workbook is not re-read at runtime except to extract strain information; fish identity is already encoded in the source data.

### Known exception: batch 5, runs 1 and 2

Batch 5 is the only batch where runs 1 and 2 (Day 1 AM and Day 1 PM) used **identical** physical tank positions for all 8 fish. Every other batch re-randomized between the AM and PM sessions on Day 1. This is confirmed in both the decode workbook (configurations F and G for batch 5 runs 1–2 list the same tank-to-fish mapping) and the source tracking data. The reason for this exception is not documented. See `manuscript/analysis.md` Section 9 for a discussion of the fish 56/57 anomaly in batch 5 that may or may not be related.

---

## Session Exclusions

### Source of exclusion flags

The pipeline inherits all session-level quality-control decisions from the original Mathematica analysis. The `include` column in `data/source/tracking_aggregate_uncorrected.csv` (originally `Z-P1-Data-Merged Raw No Interpolation-v5 aggregate outcomes.csv` from Attempt 02) marks 35 sessions as `include=0`. All 35 are dropped before any feature extraction. See `src/01_ingest.py → load_exclude_set()`.

### Breakdown of the 35 excluded sessions

| Category | N sessions | Details |
|---|---|---|
| Missing video / NaN data | 7 | Batches 11–13; video files not recovered; NaN in source |
| Empty tank (no fish) | 6 | fish_id 138, all 6 runs — batch 13 |
| Sick fish | 4 | fish_id 116, 4 of 6 runs — batch 11 |
| Tracking failure | 18 | All remaining `include=0` sessions |

**fish_id 116** (batch 11, fish #6): Was sick. Excluded from 4 of 6 runs. Shows pct_top=1.000 in multiple runs — artifact of illness, not exploration behavior.

**fish_id 138** (batch 13, fish #8): Tank was empty. No fish in the water. All 6 sessions excluded; most have 0–3 detections or NaN.

### Why the 18 "tracking failure" sessions are not reinstated

All 18 sessions flagged `include=0` that have any tracking data were examined against the distribution of detection counts in included sessions:

- Median included session: **5,422 detections** over the 20-minute window
- 10th percentile of included sessions: **2,481 detections**

Every one of the 18 sessions falls below that 10th percentile, with most far below it:

| fish_id | batch | run | n_det | note |
|---------|-------|-----|-------|------|
| 13 | 1 | 1 | 89 | < 2% of median |
| 17 | 1 | 1 | 365 | 7% of median |
| 15 | 1 | 1 | 450 | 8% of median |
| 18 | 1 | 1 | 59 | < 2% of median |
| 16 | 1 | 3 | 302 | 6% of median |
| 16 | 1 | 4 | 50 | < 2% of median |
| 18 | 1 | 4 | 109 | 2% of median |
| 37 | 3 | 1 | 123 | 2% of median |
| 38 | 3 | 5 | 109 | 2% of median |
| 35 | 3 | 5 | 367 | 7% of median |
| 37 | 3 | 6 | 227 | 4% of median |
| 71 | 7 | 2 | 196 | 4% of median |
| 84 | 8 | 1 | 132 | 2% of median |
| 86 | 8 | 1 | 172 | 3% of median |
| 83 | 8 | 4 | 73 | 1% of median |
| 88 | 8 | 5 | 295 | 5% of median |
| 86 | 8 | 5 | 349 | 6% of median |
| 113 | 11 | 4 | 115 | 2% of median |
| 126 | 12 | 1 | 142 | 3% of median |
| 123 | 12 | 4 | 1 | clearly bad |

Detection counts this low mean the tracker was essentially non-functional for those sessions. Any pct_top values derived from 50–450 frames out of an expected ~5,000 are noise, not behavior. Several sessions' pct_top values superficially fall within the normal range (0.3–0.7), but that is not grounds for reinstatement — a nearly-failed tracker will sometimes return plausible-looking summary statistics by chance.

**The batch 1 / run 1 spatial pattern** is additional evidence of equipment/physical failure rather than fish behavior: tank positions 1, 2, 5, 6 (the top row of each column in that setup) were all excluded with detection counts of 59–450, while tank positions 3, 4, 7, 8 from the same recording session had normal detection counts. This is consistent with a camera angle or lighting problem affecting the top-row tanks in that specific session.

### Summary

The original researcher's `include=0` QC flags are conservative and well-calibrated. No excluded sessions are reinstated in this reanalysis. The exclusion criterion is tracking quality (insufficient detections), not fish behavior, which means the exclusions are not expected to introduce bias in pct_top estimates.

---

## Behavioral tracking and the freezing limitation

This analysis begins from the distortion-corrected coordinate output of the original Mathematica pipeline; the raw video files were not available for re-tracking (several batches' videos were never recovered). The reanalysis therefore tests the reproducibility of the *analysis* applied to that coordinate stream, not of the upstream tracking, and it inherits the tracking algorithm's biases. This bounds every interpretation below and is stated here, rather than only in the Discussion, because it changes what the primary measure represents.

The tracker is **frame-differencing**: it registers a fish only when the fish moves between consecutive sampled frames, so a completely still (frozen) fish produces no detection (NaN). Freezing — prolonged immobility, typically near the bottom — is a primary anxiety response, so how these gaps are handled determines what the measure represents. Missing positions are filled by **linear interpolation between the detections that bracket each gap**, applied to *interior* gaps only and with **no maximum gap length**; leading and trailing gaps (before the first or after the last detection of a session) have no flanking anchor and are left undefined. **pct_top is then computed over every frame with an interpolated position** (detected plus interior-interpolated frames); the detected-frame-only measures — latency, transition count, velocity — are computed without interpolation.

This makes the freezing blind spot narrower than "freezing is invisible" would imply, but not absent, and the distinction matters for what each measure represents:

- A *within-session* freeze is bracketed by the detections immediately before and after it and is interpolated across — a fish that moves at minute 2, freezes on the bottom, and moves again at minute 18 has those 16 minutes filled by a straight line between the two anchors (there is no gap-length cutoff). Because a bottom freeze is bracketed by bottom detections, the interpolated track stays in the bottom zone and the freeze is **correctly counted as bottom-dwelling** — the only quantity pct_top needs is the top/bottom classification, which is robust whenever both anchors share a zone. Long mid-water gaps, where the straight-line assumption is weakest, are uncommon because vertical movement itself generates detections.
- The residual bias in pct_top is therefore confined to **leading/trailing freezes** — most relevant, a fish that dives and freezes on the bottom at the very start of the window, before its first detection — which are excluded and bias that session's pct_top slightly *upward* (understating early anxiety).
- The **detected-only measures** (latency, transitions, velocity) do not interpolate and so genuinely understate behavior during immobility; they should be read as movement-conditional.

In sum, pct_top is a reasonably faithful index of vertical position that is robust to interior freezing, with a small residual upward bias only for unbracketed start/end freezes; the activity measures are movement-conditional. None of this is recoverable to the absolute time-at-location, because freezing cannot be distinguished from genuine tracking dropout, but for the reliability and individual-differences claims that are this paper's focus a bias that is stable within an individual preserves the rank-ordering of fish.

## Behavioral measures and missing-data coding

Five session-level measures are analyzed: top-half occupancy (pct_top), latency to first top entry, top↔bottom transition count, within-session lateral-position SD (x_sd), and mean swimming velocity. pct_top and velocity are computed over interpolated frames (above); a freeze gap thus enters velocity as slow, constant interpolated motion. x_sd, latency, and transitions are computed over detected frames only. All inherit the movement-conditional limitation above to differing degrees.

**Latency to first top entry** is coded in minutes from the start of the analysis window. Sessions in which the fish never entered the top half are coded as **20.0 minutes** (the full analysis window length); they are not dropped. This affects only **2 of 589 sessions (0.3%)**, so all measures except latency itself are insensitive to the choice. Latency's own (already low) reliability is somewhat sensitive to it, because 20 min is an extreme value relative to the typical sub-minute latency; this is noted where latency reliability is reported.

The complete-case reliability estimates (the six-session composites and the single-measure ICCs) use the 82 fish with usable data in all six runs; the single-measure ANOVA repeatabilities use all 589 sessions. Of the fish analyzed, 21 lacked complete six-run data (11 of 48 5G, 10 of 55 AB), an attrition rate that does not differ materially between the two strains.
