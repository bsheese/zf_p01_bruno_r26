# Zenodo record — copy/paste source

Draft metadata for the Zenodo deposit of the raw tracking data. Fields marked
`‹…›` still need to be filled in before publishing.

---

## Title

Distortion-corrected frame-by-frame tracking coordinates from a zebrafish novel-tank behavioral study (Project P1)

## Description

This dataset contains the per-frame, distortion-corrected positional tracking output from a zebrafish (*Danio rerio*) novel-tank diving study. Individually housed fish were recorded in a novel tank across repeated sessions; an automated tracker localized each fish in every video frame, and the resulting pixel coordinates were scaled and distortion-corrected in custom Wolfram Mathematica software. This file is that coordinate stream — the input to all downstream behavioral analysis, not a summary of it.

The file holds 4,850,062 frame-level rows spanning 13 batches, 104 fish identifiers (8 fish per batch), and 624 recording sessions. Each fish was recorded over six sessions (three days × two sessions per day, morning and afternoon). The primary behavioral measure derived downstream is vertical tank position (time spent in the upper vs. lower half of the tank), a standard index of anxiety-like behavior in the novel-tank test.

The original videos are not included and several batches' source videos were not recovered, so this coordinate stream — not raw video — is the most upstream openly available form of the data. Analyses built on it test the behavioral signal in the tracked coordinates; they inherit the tracking algorithm's properties and cannot re-derive tracking.

Provenance and status. These data were collected and processed by a single laboratory. Findings from this study were presented as a conference poster only; the work has not been submitted to or published in a peer-reviewed journal. This record provides the coordinate data for reproducibility and reuse; it is not a publication and implies no peer review.

## Files

| File | Size | Rows | Format |
|------|------|------|--------|
| `tracking_triallevel_distortion_corrected.csv` | 160,703,963 bytes (≈153 MB) | 4,850,062 | CSV, UTF-8 |

SHA-256: `9f5c88414f86f93079fe100b166dc56ebe5c287dcefc839e4faddb1fb2747a4c`

## Column dictionary

| Column | Description |
|--------|-------------|
| `run` | Batch number (1–13). |
| `dayafternoon` | Encodes `run_within_batch × 10 + physical_tank_position` (e.g. `13` → session 1, tank 3). Sessions 1/3/5 are morning (~11:00), 2/4/6 afternoon (~14:30). |
| `tank` | Fish number within batch (1–8). |
| `fishid` | Global fish identifier = `batch × 10 + fish_within_batch`. |
| `framenumber` | Video frame index. `0` = sentinel row (no data); `1805` = first analyzed frame (first 60 s excluded because the experimenter is in frame); sampling step = 5 frames (~6 Hz at 29.97 fps). |
| `x`, `y` | Distortion-corrected, scaled position. `y = 0` at the tank floor, increasing upward. Empty/`NaN` = no detection in that frame. |
| `deltax`, `deltay` | Signed frame-to-frame displacement. |

## Conventions and notes

- Vertical orientation: `y = 0` is the tank floor and increases toward the water surface — inverted relative to raw-video pixel coordinates. Downstream, the top/bottom threshold is a fixed `y = 61` after distortion correction.
- Missing data: the tracker is frame-differencing, so a `NaN` means the fish did not move between samples; freezing and tracking failure cannot be distinguished from the coordinates alone. Overall detection rate is ≈64%.
- Frame range: real samples run from frame 1805 to 54600; typical analysis windows use frames 1805–37805 (~minutes 1–21).
- No quality-control exclusions are applied here. The 624 sessions include all recordings; the associated analysis applied curated QC flags that removed 35 sessions (empty tank, sick fish, missing video, and catastrophically low detection). This file is the unfiltered coordinate stream.
- Sex was not recorded for individual fish and cannot be recovered.
- Strain composition is batch-pure (each batch is a single strain: an inbred line and a wild-type line), so strain is fully confounded with batch; strain comparisons from these data are descriptive of batches, not genotype.

## Zenodo metadata fields

- Upload type: Dataset
- License: CC BY 4.0
- Language: English
- Keywords: zebrafish; Danio rerio; novel tank test; anxiety-like behavior; behavioral tracking; animal personality; repeatability; ethology; video tracking
- Authors / creators: ‹fill in — list the study investigators and ORCIDs›
- Funding / grants: ‹fill in if applicable›
- Related identifiers: ‹link the GitHub analysis repository once published (relation: "is supplemented by"); and the conference poster if you have a citable record›
