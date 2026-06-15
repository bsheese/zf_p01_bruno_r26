# Repository Contents

This repository presents the following: 
* Real data from an unpublished 2013 study of Zebrafish behavior in the novel-tank paradigm. 103 fish were assayed on six separate occassions for 20 minutes each. So, quite a bit of fish data: 103 fish X 6 sessions X 20 minutes at a sampling rate of 6 frames per second = ~ 4.5 million x, y, pairs. Do with it as you please. 
* A 2026 re-analysis of the 2013 study driven by A.I. which converged with the original human analysis. All the normal A.I. cautions apply here. I'd encourage you to re-derive anything interesting from scratch to verify it. That said, the A.I. has done a lot of work here to translate the raw data into data sets you can work with. You can get a good sense of the data by following the code that does this initial processing. If you find errors, please post use Github Issues or send me a note.


## The Original Study: Developing an Inexpensive High-Throughput Tracking System for Zebrafish 

This repository describes and presents data from a 2013 study that of Zebrafish behavior in the novel-tank paradigm. The primary goal of the study was to develop an inexpensive, high-throughput 
behavior tracking system. The secondary goal was to look at the short-term reliability of responses to the novel-tank paradigm. 

The tracking system involved eight tanks, each with a single fish, stacked on risers and 
recorded with a single camera in an enclosure. Following data collection, video recordings were 
processed using a frame differencing with corrections for optical distortions. 

We recorded for twenty minutes following the introduction of the fish to the novel tanks. For each fish, we repeated this assessment once in the morning and once in the afternoon for three days. In total, we assessed 103 fish on six different occassions. At the time, this was a lot of fish. 

A description of the method and some preliminary findings were presented in a poster session in 2014 (a version of the poster is available in the 'poster' sub-directory): 
- Sheese, B. E., & Deharak, B. (2014, August). Inexpensive high-throughput computer-aided tracking of zebrafish in the novel tank paradigm. Poster presentation at the annual conference of the American Psychological Association. Washington, DC.

The system worked reasonably well, but there were a few significant concerns. First, the frame differencing approach was necessary due to light reflections which varied by tank position. However, frame differencing failed to detect immobile fish and created artifacts in outcomes related to distance. Second, the camera's height, set at the mid-point of the risers, meant the camera's view into the tanks varied slightly in relation to tank's heights on the risers. In particular, it was harder to detect fish in the top row when fish were at the bottom of the tank due to the angle of the camera. Finally, it took about a minute for a trained experimenter to load all eight tanks and close the enclosure. Fish could not be introduced to the novel tanks simultaneously. We systematized how this was done (fish were poured from beakers) and randomized tank position on each assessment, but it still makes interpretation of minute by minute data problematic, since there is a time lag between when the first two fish were introduced and the when the last two fish were introduced.    

These concerns, and others, led to the development of a entirely different system (ATLeS) for susbsequent studies where each fish was individually enclosed. The AtLeS systems is described here: [https://github.com/liffiton/ATLeS] (https://github.com/liffiton/ATLeS)


## 2026 A.I. Re-analysis of the 2014 Data 

> [!IMPORTANT]
> **Status.** This work was presented previously as a conference poster and has
> **not** been peer-reviewed or published in a journal. It is a *within-laboratory*
> A.I. driven reanalysis (a re-implementation of our own earlier human driven analysis), **not** an
> independent replication. Results come from a single cohort at one facility and
> are shared for transparency and reuse, not as settled findings.

This repo does contains the tracking data derived from the original code written in Mathematica. 
That data is sound. Feel free to do whatever with you like with it, but do see the caveats mentioned above. 

Nearly everything else in this repo is A.I. generated and should be treated as such. In particular, the Python code and write-up of that data is a re-analysis conducted by A.I. (largely Claude Opus 4.8, Claude Sonnet 4.6, and Gemini Pro 3.1, with quite a bit of human orchestration). I was interested in seeing if A.I. would replicate the findings from 2013 that were originally conducted using some long-forgotten version of SPSS.  


## 2026 A.I. Driven Re-analysis 
A Python reanalysis of a six-session adult zebrafish (*Danio rerio*) novel-tank
dataset, asking **how reliable the test's behavioral measures are over a short,
within-days interval — and how that reliability should be reported.** The pipeline
re-implements our laboratory's earlier Mathematica analysis from the same
distortion-corrected tracking output and extends it.

---

## What it asks, and what it finds

Most zebrafish repeatability estimates come from intervals of weeks to months;
short-term (multiple sessions within days) reliability is largely unexamined for
adults. Using 103 fish tested twice daily over three days, this reanalysis finds:

- **Top-half occupancy (`pct_top`) is reliably repeatable over the short term.**
  Single-session repeatability is *R* ≈ 0.50; the six-session composite
  (Cronbach's α = 0.855) reproduces our earlier analysis (α ≈ 0.88) and is simply
  the Spearman–Brown step-up of that single-session *R*.
- **Report reliability as a single-measure *R*, not a composite α or a raw
  test–retest correlation.** A composite α looks far above the field's typical
  *R* ≈ 0.37 only because it averages six sessions; on the same scale the estimate
  is unremarkable. This is the main methodological takeaway.
- **The first exposure is different.** It ranks individuals differently from later
  sessions (repeatability rises to *R* ≈ 0.59 when it is excluded) — treat the
  first trial as acclimation, not trait measurement.
- **Reliability is not unique to `pct_top`.** Lateral movement range and swimming
  velocity are equally repeatable; transitions less so; latency least.

"Reliability" here means **short-term repeatability of among-individual
differences** — necessary but not sufficient for a lifelong personality trait, and
over so dense a schedule it may partly reflect a persistent post-handling state.

## Limitations to read first

- **Legacy tracking / freezing.** The original tracker is frame-differencing and
  registers a fish only when it moves, so freezing is not directly observed.
  Interpolation recovers within-session bottom freezes, but immobility is only
  partially observed. We have since adopted freezing-aware tracking; re-examining
  these questions with it is the natural next step.
- **Strain is confounded with batch.** The dataset has two strains (5G, AB) but
  every testing batch was strain-pure, so strain cannot be separated from batch
  (date, handling, per-batch calibration). **No strain comparison is made in the
  paper.** Strain/multivariate/clustering analyses are exploratory only and kept
  in [`supplemental/`](supplemental/README.md).
- **Scope.** Short-term, single cohort, single facility — not generalizable
  without replication.

## Repository layout

```
src/            Manuscript (reliability) pipeline — run in order, 01 → 12
supplemental/   Exploratory strain / multivariate / clustering work (NOT in the
                paper; strain confounded with batch). See supplemental/README.md
manuscript/     The write-up (pandoc → PDF): intro/methods/results/discussion +
                metadata.yaml + Makefile; analysis.md = full lab-notes record
data/
  source/       Original distortion-corrected tracking + decode tables (read-only)
  merged/, features/   Derived data, regenerated by the pipeline
output/         Results tables, text, and the two manuscript figures
  supplemental/ Outputs of the exploratory scripts
lit/, vault/    Reference PDFs (not tracked) and short notes on each
```

## Reproduce

Requires Python 3.10+ and the scientific stack:

```bash
pip install -r requirements.txt
```

Run the manuscript pipeline from the repo root, in order:

```bash
python src/01_ingest.py                      # source → data/merged/
python src/03_features.py                    # → data/features/
python src/04_stats.py                       # inter-run / within-day / day-to-day reliability
python src/05_plots.py                       # → output/figures/04_reliability_matrix.png (Fig 1)
python src/06_tracking_quality.py
python src/07_slot_adjusted_reliability.py
python src/08_repeatability.py               # single-measure R, batch variance, Run1-vs-familiar
python src/09_measure_reliability.py         # per-measure reliability (Table 2)
python src/11_qc_checks.py                   # confound, attrition, per-strain R, latency coding
python src/12_run1_curve.py                  # → output/figures/run1_minute_curve.png (Fig 2)
```

Build the manuscript PDF (requires `pandoc` ≥ 3, a TeX install with `xelatex`, and
the DejaVu Serif font):

```bash
cd manuscript && make        # → manuscript/manuscript.pdf
```

The exploratory analyses are independent and optional:

```bash
python supplemental/trait_structure.py
python supplemental/clustering.py
python supplemental/ml.py
```

## Data availability

The pipeline runs from `data/source/` (frame-by-frame distortion-corrected
positions, QC exclusion flags, and the batch→strain decode tables). All other data
under `data/` is derived and regenerated by the scripts.

The raw frame-by-frame tracking CSV (~153 MB, 4.85 M rows) is too large for Git and
is archived on Zenodo under CC BY 4.0:

> Distortion-corrected frame-by-frame tracking coordinates from a zebrafish
> novel-tank behavioral study (Project P01, 2013–2014). Zenodo.
> DOI: [10.5281/zenodo.20707067](https://doi.org/10.5281/zenodo.20707067)
> (all versions: [10.5281/zenodo.20707066](https://doi.org/10.5281/zenodo.20707066))

Fetch it (downloads into `data/source/` and verifies the SHA-256 checksum), then run
the pipeline:

```bash
python data/source/get_raw_data.py   # → data/source/tracking_triallevel_distortion_corrected.csv
python src/01_ingest.py              # → data/merged/
python src/03_features.py            # → data/features/
```

## Citing

This is unpublished, poster-stage work. If you use the code or data, please cite
the repository and the archived dataset:

> Sheese, B. *Short-term reliability of the zebrafish novel-tank test.* GitHub
> repository, https://github.com/bsheese/zf_p01_bruno_r26
>
> Dataset: Distortion-corrected frame-by-frame tracking coordinates from a
> zebrafish novel-tank behavioral study (Project P01, 2013–2014). Zenodo.
> DOI: [10.5281/zenodo.20707066](https://doi.org/10.5281/zenodo.20707066)

Please do not cite it as a peer-reviewed article. The underlying P1 study was
presented as a conference poster and has not been published.

## License

Dual-licensed (see [`LICENSE`](LICENSE)): **source code** under the MIT License;
**data, text, and figures** under CC-BY-4.0. Third-party reference PDFs under
`lit/` (if present) remain under their publishers' copyright.

## Contact

Brad Sheese — Illinois Wesleyan University — <bsheese@iwu.edu>
