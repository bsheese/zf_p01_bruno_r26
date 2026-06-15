# Analysis Notes — Zebrafish P1 Novel Tank Reanalysis

**Status:** Full lab-notes record (exploratory). The manuscript is reliability-only.
**Last updated:** 2026-06-09
**Pipeline:** `p01bruno_r26/src/` (manuscript) and `supplemental/` (exploratory)

> **Scope map — what backs the manuscript vs. what is exploratory.**
> The manuscript (`manuscript/`, reliability-only) draws on **Sections 1–4, 8, and 14**
> (test-retest reliability, inter-run correlations, within-day and day-to-day reliability,
> single-measure R) plus the per-measure reliability now in `src/09_measure_reliability.py`.
> **Sections 5–7, 9–13, 15–16** (strain differences, within-session habituation by strain,
> spatial distributions, the fish 56/57 and Run-2 investigations, latency/transition strain
> contrasts, ML classification, multivariate structure, clustering) are **exploratory and NOT
> in the manuscript**: strain is completely confounded with batch (every batch is strain-pure),
> so no strain result is attributable to genotype. Those analyses now live in `supplemental/`.
> Wording in the older sections below predates this scoping — e.g. "trait anxiety" should be
> read as "short-term repeatability of among-individual differences," and "independent
> reanalysis" as a within-laboratory re-implementation. See `CLAUDE.md` for the fixed facts.

---

## Overview

This document records findings from the independent Python reanalysis of the P1 zebrafish novel-tank dataset, an **unpublished** study. The goal is to test for **convergence** with the original (unpublished) analysis's reliability results (Cronbach's α ≈ 0.88) starting from the Mathematica distortion-corrected tracking output, and then extend the analysis to cover within-session behavioral patterns, morning-vs-afternoon reliability, day-to-day reliability, and strain differences. The question throughout is whether an independently built pipeline converges on the same conclusions as the original analysis.

The primary behavioral measure throughout is **pct_top**: the proportion of frames (within the analysis window) where the fish is in the upper half of the tank (y ≥ 61 in distortion-corrected coordinates). Higher pct_top = more time exploring the top = less anxiety-like behavior.

### Tracking method limitation: freezing is invisible

The Mathematica tracker uses frame-differencing — it detects a fish only when it moves between consecutive sampled frames. A fish that is completely still produces no detection (NaN). This means **freezing behavior cannot be examined directly**: a stationary fish and a fish that left the tank are indistinguishable in the data. Note (clarified): pct_top is **not** computed over detected frames only — `03_features.py` linearly interpolates interior gaps (bracketed by detections on both sides, no length cap) and computes pct_bottom/pct_top over detected + interior-interpolated frames; only leading/trailing unbracketed gaps are excluded. A within-session bottom freeze is bracketed by bottom detections and is therefore correctly counted as bottom-dwelling. The detected-only metrics are latency, transitions, and velocity. See manuscript/methods.md "Behavioral tracking and the freezing limitation."

This has downstream consequences for other metrics. Distance travelled (velocity) is similarly affected — sessions with more freezing will show artificially low total distance and potentially biased mean velocity, since the stationary periods are simply absent from the record rather than contributing a zero. Any metric derived from the tracking data shares this blind spot. Comparisons of velocity or activity between fish or sessions should be interpreted with this in mind: a fish with low detection count could be inactive (freezing) or could reflect a genuine tracking failure, and the data alone cannot distinguish them.

This is directly relevant to the fish 56/57 anomaly (Section 9) and to fish 41 Run 2, both of which involve sessions with atypically low detection counts. It is also relevant to any future velocity analysis.

---

## Dataset and Session Structure

- **103 fish** across 13 batches (8 fish per batch)
- **6 runs per fish**, yielding up to 618 fish × run sessions; 35 excluded by original QC → 589 sessions analyzed
- **Strains:** 5G (inbred, n=48 fish), AB (wild type, n=55 fish)
- **Analysis window:** frames 1805–37805 (minutes 1–21 of video; first 60 seconds excluded because experimenter is visible in frame)

The 6 runs follow a **3-day × 2-sessions-per-day** design, confirmed from the file names in the tank randomization workbook:

| Run | Date structure | Approx. time |
|-----|---------------|-------------|
| 1   | Day 1 AM      | ~11:00      |
| 2   | Day 1 PM      | ~14:30      |
| 3   | Day 2 AM      | ~11:00      |
| 4   | Day 2 PM      | ~14:30      |
| 5   | Day 3 AM      | ~11:00      |
| 6   | Day 3 PM      | ~14:30      |

Batches ran on different calendar dates across Jan–Mar 2014, but the within-batch structure is always 3 consecutive days. This means "day" refers to day within the testing period for a given batch, not a shared calendar date across batches.

---

## Section 1: Cronbach's Alpha (Test-Retest Reliability)

**Question:** How reliably does pct_top rank-order fish across the 6 runs?

82 of 103 fish had complete data across all 6 runs. Results:

| Window        | α     | 95% CI          |
|---------------|-------|-----------------|
| Full 20 min   | 0.855 | [0.801, 0.899]  |
| First 10 min  | 0.865 | [0.813, 0.906]  |
| Last 10 min   | 0.799 | [0.723, 0.859]  |

The original (unpublished) analysis reported α = .88 for the full 20-minute window. We get 0.855, which is close — the two independent analyses converge. The gap is likely due to small differences in how missing frames are handled (the original used SPSS carry-forward interpolation; our pipeline uses linear interpolation bracketed between actual detections). Both are well into the "good" reliability range.

**The first 10 minutes are more reliable than the last 10 minutes.** This is counterintuitive if you expect the novelty response to make early behavior noisier. One interpretation: early in the session, fish are showing their stable anxiety-like trait (bottom-dwelling) more consistently. By the second half, some fish have habituated, others haven't, introducing more individual variation in the direction of change rather than the level.

---

## Section 2: Inter-Run Correlations

All 15 pairwise Pearson correlations among the 6 runs, computed on fish with complete data (n=82):

| Pair                  | r     | p        |
|-----------------------|-------|----------|
| R1(D1AM) vs R2(D1PM)  | 0.468 | < .001   |
| R1(D1AM) vs R3(D2AM)  | 0.402 | < .001   |
| R1(D1AM) vs R4(D2PM)  | 0.239 | .031     |
| R1(D1AM) vs R5(D3AM)  | 0.351 | .001     |
| R1(D1AM) vs R6(D3PM)  | 0.220 | .047     |
| R2(D1PM) vs R3(D2AM)  | 0.549 | < .001   |
| R2(D1PM) vs R4(D2PM)  | 0.504 | < .001   |
| R2(D1PM) vs R5(D3AM)  | 0.479 | < .001   |
| R2(D1PM) vs R6(D3PM)  | 0.464 | < .001   |
| R3(D2AM) vs R4(D2PM)  | 0.652 | < .001   |
| R3(D2AM) vs R5(D3AM)  | 0.612 | < .001   |
| R3(D2AM) vs R6(D3PM)  | 0.621 | < .001   |
| R4(D2PM) vs R5(D3AM)  | 0.622 | < .001   |
| R4(D2PM) vs R6(D3PM)  | 0.658 | < .001   |
| R5(D3AM) vs R6(D3PM)  | 0.679 | < .001   |

**Mean inter-run r = 0.501**

The standout feature of this table is **Run 1's isolation.** Its correlations with all other runs (0.220–0.468) are substantially weaker than any correlation among runs 2–6 (0.464–0.679). Run 1 is the first-ever exposure to the novel tank, and it is measuring something behaviorally different — or at least noisier — than subsequent re-exposures.

Beyond that, there is a modest tendency for closer-in-time runs to correlate more strongly (R3 vs R4 = 0.652; R3 vs R6 = 0.621 — similar), but the pattern isn't clean enough to claim strong temporal decay. The main temporal story is the Run 1 outlier, not a continuous decay across runs 2–6.

---

## Section 3: Within-Day (AM vs. PM) Reliability

**Question:** How consistent is a fish's pct_top between the morning and afternoon sessions on the same day?

Within-day Pearson r, one per day:

| Day | AM run | PM run | r     | p        | n  |
|-----|--------|--------|-------|----------|----|
| 1   | Run 1  | Run 2  | 0.468 | < .001   | 82 |
| 2   | Run 3  | Run 4  | 0.652 | < .001   | 82 |
| 3   | Run 5  | Run 6  | 0.679 | < .001   | 82 |

**Mean within-day r = 0.600**

Day 1 is noticeably weaker than Days 2 and 3. The obvious reason is that the Day 1 AM session (Run 1) is the novel-tank debut — fish are responding to a genuinely new environment, and their behavior in that session predicts their afternoon behavior less well than later AM-to-PM pairs do. By Days 2 and 3, the tank is familiar and within-day stability improves.

The within-day r of ~0.65–0.68 on Days 2 and 3 is reasonable; 3.5 hours is not a long interval, but a lot can happen physiologically (feeding, social interaction back in home tank, circadian variation). The fact that it's not higher is worth noting.

---

## Section 4: Day-to-Day Reliability

**Question:** How consistent is pct_top across the three testing days, separately for AM and PM sessions?

ICC(A,1) — two-way random, absolute agreement, single measures — computed separately for AM sessions (runs 1, 3, 5) and PM sessions (runs 2, 4, 6):

| Session type | Runs    | ICC(A,1) | 95% CI        | F(81,162) | p        |
|--------------|---------|----------|---------------|-----------|----------|
| AM           | 1, 3, 5 | 0.450    | [0.320, 0.580] | 3.44     | < .001   |
| PM           | 2, 4, 6 | 0.543    | [0.420, 0.660] | 4.54     | < .001   |

Both ICCs are significant, confirming that fish reliably differ from each other across days. But the values themselves (0.45 and 0.54) are only moderate. This is consistent with the inter-run correlation table: the individual run-to-run correlations are in the 0.4–0.7 range, so an ICC aggregating three such correlations should land in that neighborhood.

**PM reliability is higher than AM reliability.** This is almost certainly the Run 1 effect again. The AM ICC pools Run 1 with Runs 3 and 5; since Run 1 is a novelty session that correlates weakly with everything, it drags the AM ICC down. The PM ICC uses only re-exposure sessions (Runs 2, 4, 6) and is correspondingly higher.

**Practical implication:** If you only had one day to measure, a PM session would give you a more reliable rank-ordering of fish than an AM session — though most of that advantage is just the Run 1 effect, not a principled morning-vs.-afternoon difference. If the study were redesigned, you might consider dropping Run 1 from reliability calculations since it's measuring a distinct construct (novelty response vs. trait anxiety).

---

## Section 5: Strain Differences (5G inbred vs. AB wild type)

**Question:** Do 5G and AB fish differ in their overall pct_top, and do they show different patterns across sessions?

### Overall mean pct_top

Per-fish mean pct_top (averaged across all available runs):

| Strain | N fish | M     | SD    |
|--------|--------|-------|-------|
| 5G     | 48     | 0.474 | 0.174 |
| AB     | 55     | 0.526 | 0.208 |

Welch t-test: t = -1.381, p = 0.170, Cohen's d = -0.271

The difference is **not statistically significant**. The effect size is small-to-medium (d = 0.27) in the direction of AB spending more time in the top half, but with this N it is underpowered to detect differences of this magnitude reliably. Both strains are spending roughly equal time in the top and bottom half on average (both near 0.50). This is not the pattern you'd see in a strongly anxious strain — those would be spending 70–80% of the time at the bottom.

### Strain × session trend

A linear mixed model (pct_top_full ~ strain × run_centered, fish as random intercept) found:

- No main effect of strain (b = 0.050, p = 0.188)
- No main effect of run (b = -0.011, p = 0.065; borderline, 5G showing slight decline)
- **Significant strain × run interaction (b = 0.022, p = 0.009)**

The interaction means AB's pct_top trends upward across the 6 sessions while 5G's does not (or trends slightly downward). By the end of 3 days, AB fish are spending somewhat more time in the top half than at the start, while 5G fish are where they started or slightly lower. 

This could reflect:
- AB (outbred) being more behaviorally flexible, adjusting to repeated exposures
- 5G (inbred) showing less phenotypic plasticity across re-exposures
- Or simply that 5G show a Run 1 novelty boost (more exploration initially) that decays, while AB don't show that boost

The within-session habituation analyses below bear on this.

---

## Section 6: Within-Session Habituation (the "classic" novel-tank pattern)

**Question:** Do fish show the expected within-session pattern — starting at the bottom and gradually moving toward the top as they habituate to the tank?

### All runs pooled

Linear mixed model: pct_top ~ minute (centered) + (1 | fish_id), all 6 runs pooled:

- Overall slope = 0.00074 per minute (SE = 0.00039), **p = 0.060** (borderline)
- Over the full 20-minute window, the expected total rise = 0.014 (1.4 percentage points)
- First 10 min mean = 0.502; Last 10 min mean = 0.504 — essentially flat

**Pooled across all 6 runs, the classic pattern is barely detectable.** The mean curve sits right around 0.50 throughout the session with no visible slope.

### By strain (all runs pooled)

| Strain | Slope (per min) | SE      | p     |
|--------|-----------------|---------|-------|
| 5G     | 0.00139         | 0.00057 | 0.015 |
| AB     | 0.00021         | 0.00054 | 0.701 |

5G shows significant habituation; AB does not — at least when all 6 runs are pooled. The strain × minute interaction is in this direction (p = 0.133, borderline).

### The critical finding: all of this is driven by Run 1

When we examined the minute-by-minute plots broken out by individual run (figure `02c_minute_by_run_strain.png`), a completely different picture emerged:

- **Run 1 (Day 1 AM)** — the genuinely novel exposure — shows a clear upward trajectory for both strains, starting around 35–40% top and climbing to ~50–55% by minute 20.
- **Runs 2–6** — all re-exposures — are essentially flat throughout. The within-session habituation pattern has already been extinguished by the time the fish enter the tank a second time.

When you pool all 6 runs, you are averaging one session with a meaningful positive slope against five sessions with zero slope. The pooled result is inevitably near-zero. The statistically significant 5G slope in the pooled model is probably reflecting that 5G shows a larger Run 1 effect than AB (more habituation on the first exposure), not an ongoing within-session process.

---

## Section 7: Run 1 in Detail

Given that Run 1 appears to be qualitatively different from Runs 2–6, we examined it separately with three analyses: individual fish trajectories (spaghetti), first-5-minutes vs. last-5-minutes scatter, and per-fish linear slope distributions. See figure `02d_run1_detail.png`.

### Spaghetti plot observations

The individual fish curves reveal enormous heterogeneity. Some fish enter the tank and immediately go to the top and stay there for 20 minutes; others never leave the bottom. The mean lines (5G starting ~0.35, AB starting ~0.50) are summarizing a very wide distribution, not a tight behavioral cluster.

5G (blue) shows a more pronounced starting-low pattern — the mean begins clearly below 0.5 and rises through 0.5 over the session. AB (orange) starts closer to 0.5 and shows less movement. This visual difference is consistent with the mixed model results.

### First 5 vs. last 5 minutes

Most fish fall above the diagonal (last-5 pct_top > first-5 pct_top), confirming the group-level upward shift is real and not driven by a few outliers. But the scatter is wide — plenty of fish show the opposite pattern, and the correlation between first-5 and last-5 within a single session is not tight. The within-session habituation is a population-level tendency, not a reliable individual-level trait.

### Per-fish slope distributions and strain comparison

Both strains have positive mean slopes in Run 1:
- 5G: M ≈ 0.0032 per minute
- AB: M ≈ 0.0021 per minute

Welch t = -0.36, **p = 0.717** — no significant strain difference in Run 1 slopes.

This is important. The strain × minute interaction we found in the pooled model is **not** because 5G habituates faster than AB within Run 1. Both strains show similar within-session slopes in the first exposure. Whatever is driving the pooled interaction must be happening in the re-exposure sessions (Runs 2–6), where 5G and AB may diverge in ways we haven't yet examined closely.

One hypothesis: 5G fish retain some within-session habituation drift in later runs (perhaps reflecting ongoing familiarity-building across days) while AB fish reach a stable asymptote faster. This would explain the pooled result without requiring a Run 1 difference.

---

## Section 8: Summary of Key Findings

1. **Reliability is good.** Cronbach's α = 0.855, converging with the original analysis's 0.88. The measure is stable enough to use as a trait index.

2. **Run 1 is a distinct session.** Its correlations with all subsequent runs are weaker (r = 0.22–0.47) than correlations among Runs 2–6 (r = 0.46–0.68). Run 1 is measuring a novelty response; Runs 2–6 are measuring something closer to baseline trait anxiety in a familiar environment. This distinction matters for interpreting what the measure means.

3. **Day-to-day reliability is moderate.** ICC(A,1) ≈ 0.45–0.54 across the three days. The PM-only estimate (0.54) is more meaningful than the AM estimate (0.45) because the AM ICC is inflated downward by Run 1.

4. **Within-day reliability improves after Day 1.** AM-to-PM r goes from 0.47 (Day 1) to 0.65–0.68 (Days 2–3), consistent with the first session being behaviorally distinct.

5. **No significant strain difference in overall pct_top.** 5G M=0.47, AB M=0.53, p=0.17, d=0.27. Both strains spend roughly half their time in each half of the tank on average. Neither looks strongly anxious by this measure.

6. **Strain × session interaction is significant (p=0.009).** AB trends upward across sessions; 5G does not. The mechanism is unclear — the effect does not appear in Run 1 slopes, so it must emerge from re-exposure dynamics.

7. **The classic within-session habituation pattern exists in Run 1 only.** Fish start at the bottom and gradually move toward the top in their first exposure. By Run 2 (same afternoon), this pattern is gone. Pooled across all 6 runs, the habituation slope is near-zero and borderline significant (p=0.06).

8. **Within-session habituation in Run 1 is a group tendency, not a reliable individual trait.** Individual fish vary enormously. The population mean rises, but individual fish trajectories span the full range. A fish's Run-1 slope does not appear to be a stable individual characteristic.

---

---

## Section 12: Lateral (x) Movement Variability by Strain

**Question:** The animated position plots suggested 5G fish have wider x error bars than AB. Is that real, and what does it reflect?

The x SD visible in the animation reflects two separable quantities: (1) **within-fish** variability — how much a single fish moves left-right during a session, and (2) **between-fish** variability — whether different fish tend to occupy different lateral positions. These were computed separately.

### Within-fish x SD (lateral movement range per session)

Mean within-fish x SD, averaged across fish, by strain and run:

| Strain | Run 1 | Run 2 | Run 3 | Run 4 | Run 5 | Run 6 | Overall |
|--------|-------|-------|-------|-------|-------|-------|---------|
| 5G     | 58.5  | 52.6  | 55.1  | 51.9  | 53.4  | 50.5  | **53.6 px** |
| AB     | 44.8  | 44.3  | 45.6  | 46.3  | 46.6  | 46.3  | **45.6 px** |

5G fish move laterally about **18% more** than AB fish within a session. The difference is present in every run without exception and does not diminish across days — it is not a novelty effect or a Run 1 artifact. The gap is largest in Run 1 (13.7 px) but remains substantial in Runs 2–6 (6–9 px).

### Between-fish x SD (lateral position spread across fish)

SD of per-fish mean x positions by strain and run:

| Strain | Run 1 | Run 2 | Run 3 | Run 4 | Run 5 | Run 6 | Overall |
|--------|-------|-------|-------|-------|-------|-------|---------|
| 5G     | 16.7  | 22.0  | 21.2  | 19.7  | 17.9  | 24.1  | 20.6 px |
| AB     | 23.0  | 20.0  | 18.3  | 17.5  | 16.9  | 15.2  | 18.7 px |

Between-fish x SD is nearly identical between strains (5G = 20.6, AB = 18.7) and actually reverses direction in Run 1 (AB is wider). Different fish of both strains cluster around similar lateral positions on average. The strain difference visible in the animation is entirely within-fish, not between-fish.

### Interpretation

The x SD difference is a within-fish locomotor phenomenon: **5G fish individually traverse more of the tank's width during a session than AB fish do.** This is not about 5G fish preferring different lateral positions — both strains spread similarly across the tank when you look at where their average positions fall. It is about how much each fish moves back and forth.

This finding connects with the transitions result from Section 11. 5G fish make more top-bottom zone crossings *and* show more lateral movement within sessions. Together these point toward a general locomotor pattern difference: 5G fish exhibit more extensive within-session movement in both axes, while AB fish are calmer movers. AB fish are not simply bold (high pct_top) in a static sense — they appear to reach the top half and remain there with less lateral traversal, which is a qualitatively different behavioral profile from 5G's more wide-ranging movement.

This raises a question about how to interpret strain differences in this dataset. The conventional novel-tank framing treats pct_top as the anxiety measure and implicitly assumes that more movement is better (less freezing, more exploration). But 5G fish show more movement in both the x and y axes (more lateral traversal, more zone crossings) while still ending up with lower pct_top. That combination is not easy to interpret as simply "more anxious" or "less anxious" — 5G fish may be more active but less directional in their exploration, while AB fish move less but more purposefully toward the top half.

Whether this reflects a genuine strain difference in exploratory style (erratic wide-ranging vs. directed goal-oriented) or some other factor (body size, swim speed, visual acuity differences between inbred and outbred lines) cannot be determined from position data alone.

---

---

## Section 13: Machine Learning — Strain Classification

**Question:** Can strain (5G vs AB) be predicted from behavioral features, and which features carry the signal?

This is not a deployment exercise — strain is known from genotype. The goal is to use classification as a feature importance and discriminability analysis: which behavioral measures best separate the two strains, and how much of the strain difference is captured by the conventional pct_top measure vs. the newer metrics identified in Sections 11–12?

### Setup

- **Task:** Binary classification, 5G (label=0) vs AB (label=1), N=103 fish
- **Unit of analysis:** Individual fish (not session). Using sessions as independent observations would inflate N and violate independence
- **Evaluation:** Leave-one-out cross-validation (LOOCV) throughout. With N=103, LOOCV uses 102 fish to train and 1 to test, repeated for every fish. This is conservative and appropriate for small samples
- **Classifiers:** Logistic regression (L2, C=1, standardised features), Linear Discriminant Analysis (LDA, standardised), Random Forest (200 trees). LR and LDA give interpretable coefficients; RF captures nonlinear interactions
- **Latency imputation:** Sessions where a fish never entered the top half are imputed as 20.0 min (full session length) before computing per-fish means. This treats "never entered" as maximum latency rather than as missing data

Four feature sets were tested, ranging from simple aggregates to full session-by-session trajectories. See `src/05_ml.py`.

### Feature sets

**Set A — Aggregate (N=103):** Per-fish means across available runs: `pct_top_mean`, `pct_top_sd`, `latency_mean`, `transitions_mean`, `x_sd_mean`. Uses the full sample because no complete-case restriction is needed.

**Set B — pct_top trajectory (N=82):** The six per-run pct_top values `[r1, r2, r3, r4, r5, r6]` as a feature vector. Requires complete data across all 6 runs. Tests whether the cross-session shape of the pct_top profile is discriminating.

**Set C — Theory-driven (N=94):** Five features derived from prior analysis: `delta_r1r2` (pct_top Run 1 minus Run 2, capturing the 5G-specific drop), `novelty_effect` (Run 1 pct_top minus mean of Runs 2–6), `x_sd_mean`, `transitions_mean`, `latency_mean`. Directly encodes the behavioral patterns identified in Sections 11–12.

**Set D — Full trajectory (N=82):** All 18 per-run values for pct_top, transitions, and latency. Tests whether the full multi-metric session-by-session record adds information beyond the derived features.

### Classification performance (LOOCV)

| Feature set | Classifier | AUC | Accuracy | Balanced acc. | N |
|------------|-----------|-----|----------|--------------|---|
| A — Aggregate | LR | 0.798 | 0.748 | 0.748 | 103 |
| A — Aggregate | LDA | 0.798 | 0.728 | 0.728 | 103 |
| A — Aggregate | RF | 0.801 | 0.757 | 0.754 | 103 |
| B — pct_top trajectory | LR | 0.632 | 0.622 | 0.615 | 82 |
| B — pct_top trajectory | LDA | 0.627 | 0.622 | 0.615 | 82 |
| B — pct_top trajectory | RF | 0.536 | 0.561 | 0.552 | 82 |
| C — Theory-driven | **LR** | **0.837** | **0.777** | **0.772** | 94 |
| C — Theory-driven | LDA | 0.832 | 0.777 | 0.772 | 94 |
| C — Theory-driven | RF | 0.809 | 0.766 | 0.757 | 94 |
| D — Full trajectory | LR | 0.790 | 0.732 | 0.722 | 82 |
| D — Full trajectory | LDA | 0.738 | 0.720 | 0.711 | 82 |
| D — Full trajectory | RF | 0.747 | 0.695 | 0.681 | 82 |

**Best overall:** Set C with LR, AUC=0.837, 77.7% accuracy. All classifiers perform similarly on Set C, suggesting the feature space is well-structured and linearly separable to a reasonable degree.

### Key performance findings

**Set C (theory-driven) outperforms everything,** including the full trajectory set D. Packing behaviorally-motivated derived features into 5 variables beats using 18 raw per-run values. This is partly a degrees-of-freedom issue (fewer features relative to N reduces overfitting) but also reflects that the derived features compress the relevant signal more efficiently than raw trajectories.

**Set B (pct_top trajectory alone) is surprisingly weak** (AUC=0.63 for LR, RF barely above chance at 0.54). The six per-run pct_top values, presented as a raw feature vector, are not very discriminating. This seems counterintuitive given the visible strain × run interaction, but the LOOCV result shows that the raw trajectory is too noisy at the individual fish level for a classifier to reliably extract the population-level pattern with N=82.

**Set D (full trajectory) degrades performance relative to Set A despite having more information.** Adding 12 extra features (per-run transitions and latency) when N=82 hurts generalisation in LOOCV. The information those features carry is better packaged as scalar summaries (Set A) or derived scores (Set C).

**Set A (aggregate means) is nearly as good as Set C** (AUC=0.801 vs 0.837) and uses the full N=103. For practical use, Set A is the simpler and more robust choice; Set C adds ~0.04 AUC by explicitly encoding the cross-session plasticity signal.

### Feature importance

**Theory-driven set C — RF importance and LR coefficients (all features point toward 5G):**

| Feature | RF importance | LR coef (standardised) | Direction |
|---------|--------------|----------------------|-----------|
| x_sd_mean | 0.332 | −1.328 | → 5G |
| latency_mean | 0.238 | −1.252 | → 5G |
| delta_r1r2 | 0.167 | −0.622 | → 5G |
| novelty_effect | 0.151 | +0.039 | (≈ 0) |
| transitions_mean | 0.112 | −0.139 | → 5G |

**`x_sd_mean` is the single strongest predictor in both RF and LR.** A fish with wide lateral movement range is very likely to be 5G. This is the measure identified in Section 12 — the within-session x traversal — which would never appear in a standard novel tank analysis focused on pct_top.

**`latency_mean` is second.** 5G fish take longer on average to first enter the top half. This is a different dimension of the anxiety-related behavior than pct_top measures.

**`delta_r1r2`** (Run 1 minus Run 2 pct_top) contributes meaningfully — the sharp same-day pct_top drop that is characteristic of 5G adds discriminating information beyond the mean level.

**`novelty_effect`** (Run 1 vs mean of Runs 2–6) has near-zero LR coefficient despite moderate RF importance. RF may be capturing a nonlinear interaction with other features; the linear models find it largely redundant once `delta_r1r2` and `x_sd_mean` are included.

**`transitions_mean`** is the weakest contributor. Despite the significant strain difference in Section 11 (5G makes more crossings, p=0.023, d=0.52), the classifier finds it adds little beyond what `x_sd_mean` already captures. This makes sense: both measures reflect overall lateral/vertical movement activity, and `x_sd_mean` carries that signal more cleanly.

**pct_top trajectory (Set B) — run importance:**

RF importance is nearly flat across runs r1–r6, with a slight edge at r1 (0.207) and r2 (0.193). No single run is decisively more informative than others. The classifier cannot reliably exploit the trajectory shape (5G drops, AB stays flat) from the raw values at individual fish resolution — the within-fish noise is too large relative to the between-strain signal in any single run.

### Interpretation

The classification results reframe the strain story from Sections 5–12 in a single coherent picture.

The conventional measure, pct_top, encodes real strain information — Set A using `pct_top_mean` reaches AUC=0.80. But pct_top mean is not the strongest individual discriminator. Two measures that would not appear in a standard novel tank analysis — lateral movement range (`x_sd_mean`) and latency to first top entry — carry at least as much strain signal as the traditional vertical position measure.

The best-performing feature set (Set C, AUC=0.837) works by combining three distinct behavioral dimensions:
1. **Activity style** (`x_sd_mean`, `transitions_mean`): how widely and actively the fish moves
2. **Initial anxiety** (`latency_mean`): how long before the fish first ventures to the top
3. **Behavioral plasticity** (`delta_r1r2`, `novelty_effect`): how much behavior changes between the first and subsequent exposures

5G fish score higher on all three of these dimensions — they move more laterally, take longer to enter the top, and show a sharper behavioral shift between sessions. AB fish score lower — calmer lateral movement, faster initial exploration, more stable re-exposure behavior. The two strains are not simply "more anxious" vs "less anxious" but appear to differ in a broader locomotor and plasticity profile.

### Caveats

- N=103 is small for machine learning. LOOCV is the appropriate evaluation strategy but results should be interpreted with caution — the confidence intervals on AUC=0.837 are wide
- All classifiers use default or lightly-tuned hyperparameters. No inner cross-validation for hyperparameter selection was performed
- The imputation of latency NaN as 20.0 min is a reasonable but arbitrary choice; sensitivity to this choice has not been tested
- The features are correlated (x_sd and transitions in particular), so individual feature importances should not be over-interpreted. The RF importance values are Gini-based and subject to bias for correlated features

---

## Open Questions and Next Steps

- **Why does the strain × run interaction exist if Run 1 slopes don't differ by strain?** Need to look at the per-run minute curves for 5G vs AB separately across Runs 2–6 to see where the divergence begins.

- **Is the first-vs-later-session distinction worth formalizing?** Could compute alpha separately for Run 1 alone vs. Runs 2–6, and report ICCs with and without Run 1. This might be worth including in a methods note or supplementary.

- **Does pct_top in Run 1 correlate differently with pct_top in Run 6 than Run 2?** The inter-run correlation table shows run 1 correlates weakly with all others — but we haven't tested whether Run 1 is more predictive of any particular later run.

- **Are there individual fish who consistently show strong within-session habituation across all 6 runs?** The current analysis averages across fish. Some fish might reliably habituate on every exposure; others might not. Worth computing per-fish mean slope across runs.

- **Velocity.** We compute mean velocity but haven't analyzed it yet. Velocity may add independent information — a fish could be in the top half but barely moving (thigmotaxis up there) vs. actively swimming. Could parallel the pct_top reliability analyses.

- **The Run 2 noise.** → *Investigated. See Section 9 below.*

---

## Section 9: Run 2 Noise Investigation

The per-run small-multiples plot (02c) showed elevated variability in Run 2 (Day 1 PM) for the AB strain, with what appeared to be a late-session spike. Investigation found no systematic methodology issue — the elevated SD in Run 2 is modest (0.260 vs. 0.229–0.255 for other runs) and is driven by three distinct fish-level issues rather than anything structural about the Run 2 session itself. Removing all four affected fish drops the Run 2 SD from 0.260 to 0.256 — effectively unchanged.

### Issue 1: Batch 5 fish 56 and 57 — anomalous behavioral reversal, origin unknown (5G strain)

These two fish show a large, permanent behavioral change between Day 1 (runs 1–2) and Day 2+ (runs 3–6), and the changes are roughly opposite in direction:

| Fish | Run 1 | Run 2 | Run 3 | Run 4 | Run 5 | Run 6 |
|------|-------|-------|-------|-------|-------|-------|
| 56   | 0.237 | 0.163 | 0.769 | 0.793 | 0.693 | 0.844 |
| 57   | 0.800 | 0.520 | 0.010 | 0.010 | 0.070 | 0.050 |

The source data and decode workbook are internally consistent — fish identity assignments agree and were not corrupted by data processing. Both fish have adequate detection counts in all sessions, so this is not a tracking quality issue.

**Why a physical fish swap is one hypothesis:** The reversal is simultaneous (both fish change at the same boundary, Day 1 → Day 2), permanent (4 consecutive sessions), and roughly mirror-image, which would be consistent with a one-time mixup when re-placing fish for Run 3. If a swap occurred, the corrected profiles would be: fish 56 consistently low (0.24, 0.16, 0.01, 0.01, 0.07, 0.05) and fish 57 consistently high (0.80, 0.52, 0.77, 0.79, 0.69, 0.84) — both coherent stable-trait profiles.

**Why it's not clearly a swap:** If it were a simple physical swap, fish 57's runs 3–6 values (recorded as fish 57, but hypothetically fish 56's body) should look like fish 56's Day 1 behavior (~0.20). Instead they're near-zero (0.01–0.07). That's lower than expected — the numbers don't perfectly close. An equally plausible alternative is that fish 57 developed a health or behavioral problem starting Day 2 (stress, injury, illness) that caused it to essentially stop exploring, while fish 56 independently changed in the opposite direction for an unrelated reason. Detection counts for fish 57 in runs 3–6 are in the hundreds to low thousands, so the fish is still physically moving — this isn't the same pattern as fish 116 (immobile sick fish). But a mobile fish can still spend all its time at the bottom.

**The origin cannot be determined from the data alone.** The relevant videos were inspected:
- `ZFish-P1-02-13-14 1430 Main-5GWT` (Run 2, Day 1 PM) — last session before the change; fish 57 in physical tank 2, fish 56 in physical tank 7
- `ZFish-P1-02-14-14 1100 Main-5GWT` (Run 3, Day 2 AM) — first session after the change; fish 57 in physical tank 7, fish 56 in physical tank 8

Video inspection was not conclusive. The footage did not provide clear enough visual evidence to distinguish between a physical fish swap and a genuine behavioral change (e.g., body markings were not distinctive enough to confirm or rule out identity). The cause of the reversal remains unknown.

**Additional anomaly in batch 5:** Batch 5 is the only batch where runs 1 and 2 have identical physical tank positions for every fish. All other batches reassign fish to different physical tanks between AM and PM sessions on Day 1. This is confirmed in both the decode workbook (configurations F and G show the same tank-to-fish mapping) and the source tracking data. Whether this was intentional design for batch 5, a recording error, or related to the fish 56/57 issue is unknown.

**Current status:** Fish 56 and 57 are included in all analyses as-is. Their anomalous profiles inflate within-fish variance. The appropriate response (exclusion, correction, or retention) depends on what the lab records show.

### Issue 2: Fish 41, Run 2 — low-detection session with extreme pct_top (AB strain, batch 4)

Fish 41 is a consistent bottom-dweller: pct_top = 0.052, 0.028, 0.021 in runs 1, 3, and 5 respectively. In Run 2 it has only 1,004 detections (vs. 6,555–7,601 in all other runs) and pct_top = 0.833.

This combination — drastically reduced detection count and a complete behavioral reversal from a strong pattern — is a yellow flag for tracking failure rather than genuine behavior. A likely scenario is that the fish was nearly stationary somewhere in the upper half of the tank for most of the session, with the frame-differencing tracker only firing on rare movements. The result would be a low detection count concentrated in the top half regardless of where the fish actually spent most of its time.

Fish 41 Run 2 was not excluded by the original QC because 1,004 detections clears the implicit detection threshold used (exclusions were generally below ~500). However, 1,004 is approximately 14% of this fish's normal detection rate, which is qualitatively very different from a threshold-based pass.

**Recommendation:** Flag fish 41 Run 2 for sensitivity analysis — report results with and without this session to assess its impact.

### Issue 3: Fish 116, Runs 1–2 — sick fish included for Day 1 (AB strain, batch 11)

Fish 116 is the confirmed sick fish. However, the original QC only excluded Runs 3–6 (where detection counts drop to 14–135), not Runs 1–2. Run 1 has 4,654 detections and pct_top = 0.958; Run 2 has 3,692 detections and pct_top = 1.000. The fish appears to have been active on Day 1 — the sickness manifested as near-immobility starting from Run 3.

Including these two sessions is defensible since the fish was behaviorally active. However, pct_top = 1.000 in Run 2 is the maximum possible value and contributes an extreme data point to the AB strain mean. Removing fish 116 from Run 2 changes the last-5-minute mean pct_top from 0.548 to 0.540 — a small effect.

### Summary

| Issue | Fish | Runs affected | Strain | Recommended action |
|-------|------|---------------|--------|--------------------|
| Anomalous reversal (swap or health event) | 56, 57 | 3–6 | 5G | Check lab notebooks + batch 5 videos; sensitivity analysis |
| Low-detection anomaly | 41 | Run 2 only | AB | Flag; sensitivity analysis |
| Sick fish included | 116 | Runs 1–2 | AB | Low impact; note in methods |
| Batch 5 same-tank Day 1 | All batch 5 | Runs 1 and 2 | 5G | Flag as unusual; investigate source |

The Run 2 late-session spike visible in the per-run plot is real but mild (last-5-min mean = 0.548 for AB, dropping to 0.522 after removing fish 116 and 41). It is not a Run 2-specific methodology problem.

---

## Section 10: Spatial Distribution of Movement (x and y Positions)

**Question:** Do the raw position distributions show the expected unimodal x / bimodal y pattern observed in previous novel-tank work?

Analysis based on 2,822,173 detected frames (movement-only; frame-differencing tracker; all fish × all runs pooled). See figures `03a_xy_distributions.png` and `03b_xy_distributions_by_strain.png`.

### Overall distributions

**Horizontal (x):** Approximately bell-shaped and centered near the tank midpoint (mean = 142.2 px; tank width ≈ 272 px). The distribution is slightly flat-topped rather than sharply peaked, consistent with fish sampling the full width of the tank rather than preferring a lateral position. No strong asymmetry. This replicates the unimodal x pattern reported in previous novel-tank studies.

**Vertical (y):** Clearly bimodal. Two distinct modes:
- Bottom mode ≈ y = 39 px (lower third of tank)
- Top mode ≈ y = 86 px (upper two-thirds of tank)
- Valley near y = 52–58 px (below the tank midpoint of 61)

The bimodal shape reflects a behavioral choice pattern: when fish are moving (the only state the tracker sees), they are predominantly either hugging the floor or hovering in the upper half, with relatively little time in the transitional mid-water zone. This is the classic thigmotaxis signature. The top mode is taller and broader than the bottom mode, consistent with the overall mean pct_top being slightly above 0.50.

### Strain differences

Both strains produce the same structural patterns — unimodal x, bimodal y — but differ in the relative height of the two y modes:

| Strain | x mean | y mean | pct_top (detected frames) |
|--------|--------|--------|--------------------------|
| 5G     | 142.8  | 60.6   | 0.514                    |
| AB     | 141.7  | 62.4   | 0.551                    |

**x:** The two strain distributions are nearly identical. No lateral preference difference by strain.

**y:** Both bimodal, but the balance between modes differs:
- **5G:** bottom mode taller than top mode — when moving, 5G fish are slightly more often in the lower half
- **AB:** top mode taller — when moving, AB fish are slightly more often in the upper half

This frame-level pct_top difference (0.514 vs. 0.551) is a larger version of the same direction as the per-fish summary comparison (0.474 vs. 0.526 from Section 5), and is consistent with the strain × run interaction trend where AB progressively spends more time at the top across sessions.

### Interpretation and caveats

The bimodal y distribution confirms that fish are making a behavioral choice about vertical position, not simply distributing uniformly through the water column. The valley near y ≈ 52–58 px is well below the tank midpoint (y = 61), meaning fish avoid the middle of the water column more than a simple bottom/top dichotomy would suggest — the functional boundary between "bottom-dwelling" and "top-exploring" behavior may be at a lower physical position than the geometric midpoint.

The tracker's movement-only detection means these distributions represent where fish are *while moving*. If fish preferentially freeze at the bottom (as expected during the anxiety phase), the bottom mode is underrepresented relative to true time-at-location. The bimodal pattern is almost certainly stronger in real-time position than what the tracker captures — the bottom mode would be even taller if freezing were included. This reinforces the interpretation that pct_top is a conservative measure of exploration tendency.

---

### By run: the strain crossover

Splitting the y distributions by run reveals the mechanism underlying the strain × run interaction found in Section 5. See figure `03c_y_distribution_by_run_strain.png`.

pct_top computed from detected frames only:

| Run | Day/Time | 5G    | AB    |
|-----|----------|-------|-------|
| 1   | Day 1 AM | 0.571 | 0.542 |
| 2   | Day 1 PM | 0.502 | 0.568 |
| 3   | Day 2 AM | 0.504 | 0.517 |
| 4   | Day 2 PM | 0.507 | 0.546 |
| 5   | Day 3 AM | 0.528 | 0.563 |
| 6   | Day 3 PM | 0.472 | 0.576 |

**The bimodal shape is present in every run for both strains** — it is not a pooling artifact. What changes across runs is the *relative height* of the two modes.

**Run 1 is reversed.** In the first novel-tank exposure, 5G has higher pct_top than AB (0.571 vs 0.542). 5G's top mode is clearly dominant in Run 1; its distribution is the more top-heavy of the two. This is the opposite of every subsequent run.

**The crossover happens between Run 1 and Run 2.** By Run 2 (same afternoon, same batch of fish), 5G drops to 0.502 — essentially at the midline — while AB climbs to 0.568. From Run 2 onward, AB consistently exceeds 5G and the gap is maintained or widens. By Run 6, 5G has drifted below the midline (0.472), meaning 5G fish spend more of their detected movement time in the bottom half by the end of the study.

**5G shows progressive bottom-loading across runs.** Comparing 5G's distribution shape from Run 1 to Run 6, the bottom mode grows and the top mode shrinks. By Run 6 the bottom mode is clearly the taller of the two. The transition is not gradual — the big drop happens between Run 1 and Run 2, and then there is a slower continued drift downward.

**AB is structurally stable.** Its distribution shape changes little across runs. The top mode stays dominant, the bottom mode stays smaller throughout. AB's pct_top rises modestly (0.542 → 0.576) without any reversal or crossover.

**Interpretation.** The pattern is consistent with 5G and AB responding differently to repeated exposure. In Run 1, both strains are encountering a genuinely novel environment; 5G fish are if anything slightly more exploratory at first contact (higher initial pct_top). But 5G's behavior shifts substantially between Run 1 and all subsequent runs — the second exposure extinguishes the novel-tank response and 5G settles into a predominantly bottom-dwelling pattern that persists and gradually deepens across the three days. AB, by contrast, maintains its exploratory tendency across re-exposures with only minor drift.

This is not the picture suggested by the pooled strain comparison (Section 5, p = 0.17, d = 0.27). Pooling across runs conceals the crossover entirely. The strain difference is not simply that one strain is more anxious; rather, the strains respond differently to the transition from novel to familiar environment. 5G may show a stronger initial exploration burst that is then suppressed by familiarity, while AB's behavior is more consistent across exposures regardless of novelty.

Whether this reflects a genuine difference in behavioral flexibility, habituation rate, or some other process is an open question. The run-by-run distribution data provide the clearest view of this effect so far.

---

## Section 11: Latency to Top and Transition Count

Two additional session-level metrics computed from frame-level detected positions:

- **Latency to top**: time (minutes) from the start of the analysis window to the first detected frame in the upper half of the tank (y ≥ 61). NaN if the fish never entered the top during the session. Computed on detected frames; reflects when the fish first moved into the top half, not necessarily when it first crossed the midline.
- **N transitions**: number of top↔bottom zone crossings counted across consecutive detected frames. An undercount of true crossings (gaps between detections are invisible) but consistent across fish and sessions. Captures the frequency of vertical shuttling independently of where the fish spends its time.

### Reliability

| Metric | Cronbach's α | 95% CI | Mean inter-run r | Day-to-day ICC AM | Day-to-day ICC PM |
|--------|-------------|--------|-----------------|-------------------|-------------------|
| pct_top *(reference)* | 0.855 | [0.801, 0.899] | 0.501 | 0.450 | 0.543 |
| Latency to top | 0.616 | [0.470, 0.732] | 0.250 | 0.307 | 0.236 |
| N transitions | 0.764 | [0.676, 0.835] | 0.363 | 0.275 | 0.424 |

Both new metrics are less reliable than pct_top. Transitions hold up reasonably well (α = 0.764), suggesting the overall level of vertical shuttling activity is a moderately stable individual trait. Latency is the weakest (α = 0.616), with a mean inter-run r of only 0.250 — it is a meaningful measure of a state (initial anxiety) but not a consistent trait across repeated exposures.

The day-to-day ICC pattern for latency is inverted relative to pct_top: PM sessions show *lower* ICC (0.236) than AM sessions (0.307), whereas for pct_top PM was higher. This likely reflects a different relationship between the Run 1 effect and latency (discussed below).

### Strain differences

| Metric | 5G mean | AB mean | Welch t | p | Cohen's d |
|--------|---------|---------|---------|---|-----------|
| Latency to top (min) | 1.15 | 0.52 | 2.578 | 0.012 | 0.582 |
| N transitions | 113.8 | 94.1 | 2.322 | 0.023 | 0.524 |

Both metrics show significant strain differences with medium effect sizes — comparable in magnitude to pct_top (d = 0.271, p = 0.170 there). Notably, the strain direction for transitions is **opposite** to pct_top:

- For pct_top, AB > 5G (AB spends more time in the top)
- For latency, 5G > AB (5G takes longer to first enter the top)
- For transitions, **5G > AB** (5G makes more zone crossings)

These three findings together sketch a behaviorally coherent picture: 5G fish are slower to first reach the top, spend less total time in the top, but cross the midline more often when they do move vertically. AB fish get to the top faster, stay there more, and shuttle less. In other words, 5G fish show more erratic vertical movement — repeatedly visiting both zones — while AB fish show more committed top-dwelling once they overcome initial anxiety. The pct_top difference alone would suggest AB is simply bolder, but the transition data complicates that: high transitions in 5G indicates active vertical movement, not just bottom-hugging.

### Per-run patterns

**Latency per run:**

| Run | M (min) | SD  |
|-----|---------|-----|
| 1 (Day 1 AM) | 0.42 | 1.05 |
| 2 (Day 1 PM) | 1.51 | 3.13 |
| 3 (Day 2 AM) | 0.89 | 1.73 |
| 4 (Day 2 PM) | 0.92 | 1.97 |
| 5 (Day 3 AM) | 0.51 | 1.36 |
| 6 (Day 3 PM) | 0.56 | 1.30 |

The most striking feature is **Run 2**: mean latency of 1.51 minutes, nearly triple Run 1 (0.42 min) and roughly double all other runs. This is counterintuitive — fish are returning to the same environment they were in just 3.5 hours earlier, yet they take longer to first enter the top than on the genuinely novel first exposure. The SD in Run 2 (3.13) is also the largest of any run, driven by a subset of fish that were very slow or never entered the top at all on Day 1 PM. This is worth flagging as an open question (see below). By Days 2 and 3, latency is back to near-Run-1 levels.

**Transitions per run:**

| Run | M | SD |
|-----|---|----|
| 1 (Day 1 AM) | 120.9 | 69.8 |
| 2 (Day 1 PM) | 86.4 | 47.2 |
| 3 (Day 2 AM) | 107.9 | 62.7 |
| 4 (Day 2 PM) | 94.5 | 57.1 |
| 5 (Day 3 AM) | 109.5 | 45.6 |
| 6 (Day 3 PM) | 98.6 | 50.6 |

Run 1 shows the most transitions (M = 120.9), consistent with the novelty-driven exploration visible in the pct_top habituation curves. Transitions drop sharply in Run 2 (86.4) and recover partially in Runs 3 and 5 (the AM sessions), with lower values in Runs 4 and 6 (PM sessions). The AM > PM pattern for transitions is consistent across days, suggesting fish are more vertically active in the morning sessions — possibly a circadian or arousal effect, or simply that morning sessions always follow a longer overnight rest period.

### Visual patterns (figure 04a_latency_transitions_by_run_strain.png)

Plotting mean ± SE with individual fish as jittered background points reveals structure not fully visible in the summary numbers.

**Transitions:** The 5G > AB separation is present on every run without exception — this is not a pooled average obscuring run-specific crossings, it is a consistent within-run difference throughout the study. The AM > PM alternation is also clearly visible: both strains show more transitions in the morning sessions (D1AM, D2AM, D3AM) than in the adjacent afternoon sessions. This pattern was noted in the per-run table but the plot makes it unambiguous. Individual variability is large (points range from near-zero to ~450 crossings per session), which accounts for the overlapping error bars despite the consistent directional difference.

**Latency:** The Run 2 spike dominates the plot visually. Both strains show elevated latency on D1PM, but the effect is asymmetric: 5G climbs to a mean of ~2.5 minutes, pulled by a visible cluster of individual fish at 5–18 minutes; AB barely moves from its near-zero baseline. This means the Run 2 latency anomaly is predominantly a **5G phenomenon**. The numerical mean for AB in Run 2 (not shown separately by strain in the per-run table in the stats output) is much closer to its overall baseline than the combined mean of 1.51 minutes suggests. The individual point cloud also clarifies that most fish have near-zero latency in most runs — the run-level means are being driven by a right-skewed tail of slow-to-enter fish, not a general population-level shift. This skew explains both the high SDs and the lower reliability of latency as a trait measure.

### Summary and open questions

The three metrics (pct_top, latency, transitions) are not redundant. They capture different aspects of novel-tank behavior and show different reliability and different strain profiles. pct_top remains the most reliable and interpretable overall measure. Latency captures initial anxiety state but is too noisy for reliable individual trait measurement across sessions. Transitions add information about behavioral style (erratic vs. committed exploration) that pct_top alone misses.

The AM > PM pattern in transitions is a new finding not visible in the pct_top data. It suggests a time-of-day effect on vertical activity that is independent of where in the tank fish spend their time: fish cross the midline more often in morning sessions regardless of strain or day. Whether this reflects circadian regulation, arousal state, or some aspect of home-tank management (feeding schedule, tank disturbance) between sessions is unknown.

The Run 2 latency anomaly is predominantly a 5G effect. AB fish return to the tank in the afternoon with near-normal latency; 5G fish take considerably longer to first enter the top half. Combined with the fact that 5G pct_top drops sharply between Run 1 and Run 2 (Section 10), this suggests that 5G fish show a pronounced same-day re-exposure effect: the first afternoon return to the tank is more anxiety-provoking for them than for AB, even though their initial novel-tank response in Run 1 was actually more exploratory (higher Run 1 pct_top than AB). Whether this reflects faster extinction of the novelty response followed by stronger context-conditioning, or some other mechanism, is an open question that the current data cannot resolve.

---

---

## Section 14: Single-Measure Repeatability (R), Batch Variance, and Novelty vs. Trait

**Pipeline:** `src/08_repeatability.py` → `output/repeatability.txt`, `output/repeatability.csv`. Addresses analysis_plan.md items #2, #4, #6, #8.

**Motivation.** All reliability reported so far is Cronbach's α (Section 1, α = 0.855). α is a *6-session composite* reliability and is not on the same scale as the **single-measure repeatability R** = V_among / (V_among + V_within) reported throughout the animal-personality literature. To place P1 in that literature and to disentangle the Run 1 novelty effect from trait consistency, we computed single-measure R (agreement R_A and consistency R_C), partitioned batch variance, and recomputed reliability with and without Run 1.

### Single-measure repeatability

| Estimate | Method | R | 95% CI |
|----------|--------|-----|--------|
| R_A (all 6 runs, complete case n=82) | one-way ICC(1,1), cluster-boot CI | 0.498 | [0.378, 0.593] |
| R_C (all 6 runs, complete case n=82) | two-way ICC(C,1), cluster-boot CI | 0.496 | [0.381, 0.594] |
| R_A (all 589 sessions, unbalanced) | one-way ANOVA (Lessells & Boag 1987) | 0.511 | [0.409, 0.598] |
| R_A (all 589 sessions, unbalanced) | mixed model (1 \| fish) | 0.517 | — |

**Single-session pct_top repeatability is R ≈ 0.50.** The estimate is stable across complete-case ICC, full-sample ANOVA, and mixed-model approaches. This is the number that should be compared to the literature — **not** α.

**Spearman-Brown bridge.** A single-measure R_A = 0.498 projected to a 6-session composite gives k·R/(1+(k−1)R) = **0.856**, reproducing the reported α = 0.855 almost exactly. Conversely, α = 0.855 implies a single-measure R of 0.496. The high α is therefore not in tension with the literature: it is the arithmetic consequence of averaging 6 sessions, not an unusually reliable single measurement.

**Benchmark placement.** Published single-measure R: bell2009 meta-analysis lab mean R ≈ 0.37; thomson2020 zebrafish bottom-time R = 0.39; johnson2025 lower-zone r = 0.29. P1's single-session R ≈ 0.50 sits at or slightly above the upper end of this range — pct_top is a comparatively (but not extraordinarily) reliable single measure, and the apparent gap between "0.88" and the field disappears once the composite-vs-single distinction is made.

### Batch as a variance component

Fish are nested in batch (8 fish/batch, 13 batches, Jan–Mar 2014). A nested model `pct_top ~ 1 + (1|batch) + (1|fish within batch)` partitions:

| Component | Variance | % of total |
|-----------|----------|-----------|
| Batch | 0.00118 | 1.9% |
| Among-individual (within batch) | 0.03033 | 49.8% |
| Within-individual (residual) | 0.02934 | 48.2% |

**Batch explains only 1.9% of pct_top variance.** The batch-adjusted R (0.498) is identical to the unadjusted R_A. Despite the 13 batches running on different dates and the documented within-batch spatial gradient (tracking_quality_by_position.md), batch-level environmental variance is negligible — the reliability estimate is not inflated by batch structure, and pooling fish across batches is justified.

### Novelty (Run 1) vs. trait (Runs 2–6)

| Set | Cronbach's α | R_A | R_C |
|-----|-------------|-----|-----|
| All 6 runs | 0.855 | 0.498 | 0.496 |
| Runs 2–6 only | 0.876 | 0.586 [0.478, 0.670] | 0.586 [0.484, 0.673] |

**Dropping Run 1 raises single-measure R from 0.50 to 0.59** (and α from 0.855 to 0.876). The mechanism is important and slightly counterintuitive:

- The **population mean shift** across runs is tiny (Run 1 mean = 0.518; Runs 2–6 mean = 0.498; shift = +0.020). Consequently **R_A ≈ R_C** — conditioning on run as a fixed effect (the Biro & Stamps 2015 "consistency" repeatability that removes population-level time trends) barely changes the estimate.
- Therefore the Run 1 problem is **not a mean-level trend** but **individual-level reordering**: Run 1 rank-orders fish differently from the later runs. r(Run 1, mean of Runs 2–6) = 0.411, versus a mean pairwise r of 0.584 among Runs 2–6 themselves.

This is why *dropping* Run 1 — not *conditioning on* run — is what recovers repeatability. The pattern matches thomson2020, where week-1→2 repeatability (R = 0.25) was far below later intervals and was interpreted as the fish needing one tank experience before stable trait behavior emerges. The P1 Run 1 effect is the same documented acclimation phenomenon, not a measurement anomaly: Run 1 measures the novelty response; Runs 2–6 measure trait anxiety in a familiar environment, and the trait is meaningfully more repeatable (R ≈ 0.59) than the all-runs estimate suggests.

**Reporting recommendation.** Report the single-measure R (≈ 0.50 all runs, ≈ 0.59 trait-only) alongside α, with the Spearman-Brown bridge, so the convergence with the original α ≈ 0.88 is stated in a form comparable to the broader literature. Treat Runs 2–6 as the trait-anxiety reliability estimate.

---

---

## Section 15: Multivariate Structure of Individual Differences

**Pipeline:** `src/09_trait_structure.py` → `output/trait_structure.txt`, `trait_correlations_{raw,among}.csv`, `trait_pca.csv`; per-session trait table cached to `data/features/trait_sessions.csv`. Addresses analysis_plan.md items #1 (velocity/x_sd reliability) and #3 (behavioral syndrome structure).

**Question.** Sections 11–13 showed pct_top, latency, transitions, and x_sd are non-redundant and carry strain signal in different directions. How many behavioral dimensions are actually present, and which metrics move together? A behavioral-syndrome claim concerns the **among-individual** correlation (the correlation of true trait values), not the raw correlation of per-fish means, which is attenuated by within-individual noise. Both are reported; the among-individual structure is estimated two independent ways (disattenuation by reliability; bivariate mixed model), which agree.

Five traits per fish × run: pct_top, latency (NaN→20 min), n_transitions, x_sd (within-session lateral SD), and mean velocity (speed *when moving* — the frame-differencing tracker does not see freezing, so velocity is not distance travelled).

### Univariate repeatability of all five traits

| Trait | R (single-session) | α (6-session composite) |
|-------|-------------------|------------------------|
| pct_top | 0.511 | 0.855 |
| x_sd | 0.521 | 0.860 |
| velocity | 0.514 | 0.828 |
| transitions | 0.376 | 0.764 |
| latency | 0.327 | 0.533 |

**x_sd and velocity are as repeatable as pct_top** (single-session R ≈ 0.51–0.52), establishing them as genuine individual traits rather than noise. Velocity's R = 0.51 is consistent with baker2018 swim-speed R = 0.40–0.59 and thomson2020 distance R = 0.46 — though it must be benchmarked as speed-when-moving, not distance, because of the freezing-invisible tracking limitation. Latency is the least repeatable (R = 0.33); its composite α here (0.533) is lower than Section 11's 0.616 because Section 11 used listwise deletion of never-entered sessions while this analysis fills them as 20 min (n = 82), which adds a right-skewed extreme — latency reliability is sensitive to that choice.

### Two dissociable dimensions

Among-individual correlations (disattenuated / mixed-model cross-check) reveal two clusters:

- **Vertical-anxiety dimension:** pct_top ↔ latency strongly negatively correlated (among-individual r = −0.92 disattenuated, −0.70 mixed-model). More top time ↔ faster first entry. These are essentially one construct. (The disattenuated −0.92 over-corrects given latency's modest reliability; the mixed-model −0.70 is the more trustworthy estimate.)
- **Locomotor-activity dimension:** transitions ↔ velocity ↔ x_sd form a positive cluster (transitions×velocity r ≈ 0.71–0.75; transitions×x_sd ≈ 0.57–0.62; x_sd×velocity ≈ 0.17–0.24).
- **The two dimensions are largely independent:** pct_top correlates only 0.11–0.21 with the activity cluster. Vertical position and locomotor activity dissociate — directly matching blaser2012's finding that diving and other anxiety measures reflect separable mechanisms.

PCA on z-scored per-fish means confirms this: **two components have eigenvalues > 1 (Kaiser), explaining 70% of variance.**

| PC | Eigenvalue | % var | Interpretation |
|----|-----------|-------|----------------|
| PC1 | 2.17 | 43.5% | General activity: transitions (+.53), velocity (+.45), x_sd (+.41), pct_top (+.39), latency (−.44) |
| PC2 | 1.34 | 26.8% | Vertical anxiety vs. activity: pct_top (+.61), latency (−.54) against transitions (−.40), velocity (−.32), x_sd (−.28) |

PC1 is a general activity/exploration axis (everything loads positively, latency negatively). PC2 contrasts committed top-dwelling (high pct_top, low latency) against erratic shuttling/lateral activity.

### The strain difference lives on PC2 (and is large)

Welch t-tests of PC scores by strain:

- **PC1 (activity):** 5G M = +0.23, AB M = −0.19, t = 1.22, **p = 0.23, d = 0.28** — not significant.
- **PC2 (anxiety-vs-activity):** 5G M = −0.50, AB M = +0.41, t = −3.97, **p = 0.0002, d = −0.87** — large and highly significant.

**This is the most important result of the multivariate analysis.** The univariate pct_top strain comparison (Section 5) was non-significant with a small effect (d = 0.27, p = 0.17). Projecting onto the PC2 axis — which combines pct_top, latency, transitions, x_sd, and velocity in their natural correlational structure — yields a **large strain effect (d = 0.87)**. AB fish sit at the positive end (high pct_top, fast to top, few transitions, low lateral traversal = calm, committed top-dwellers); 5G fish sit at the negative end (lower pct_top, slower to top, more transitions, much wider lateral movement = active, erratic, wide-ranging). The strain difference is genuinely multivariate: it is diluted in any single measure and only emerges clearly in the composite. This is the quantitative synthesis of the qualitative picture built up across Sections 5, 10, 11, and 12.

### Implications

1. **pct_top is not the whole story.** It captures one of two dimensions. A complete behavioral description needs an activity axis (best indexed by transitions or velocity) alongside the vertical-anxiety axis (pct_top/latency).
2. **The strain difference is a profile, not a level.** 5G and AB differ most on the *combination* of being lower/slower in the water column while simultaneously more laterally active — not on mean top time alone.
3. **Next step (plan #7) is now well-posed.** Unsupervised clustering should be run on the two PC scores (or the trait set), with the expectation — per rajput2022 — that phenotypic clusters may cut across the strain boundary rather than reproduce it. *(Done — see Section 16.)*

---

---

## Section 16: Unsupervised Behavioral Phenotype Clustering

**Pipeline:** `src/10_clustering.py` → `output/clustering.txt`, `cluster_assignments.csv`, `figures/07a_clusters_pca.png`. Addresses analysis_plan.md item #7.

**Question.** Do natural behavioral clusters exist, and do they reproduce the 5G/AB boundary or cut across it (rajput2022; johnson2025)? Clustering used all five z-scored per-fish traits (pct_top, latency, transitions, x_sd, velocity) on the complete-case fish (n=82), with PC1/PC2 from Section 15 used only for visualization. k = 2, 3, 4 were tested with k-means and Ward hierarchical, scored by silhouette and bootstrap stability (mean ARI across 500 resamples).

### Cluster structure is weak — the data are closer to continuous than discrete

| k | KMeans silhouette | Ward silhouette | Stability ARI (mean ± sd) |
|---|------------------|-----------------|---------------------------|
| 2 | 0.261 | 0.305 | 0.530 ± 0.308 |
| 3 | 0.257 | 0.246 | 0.499 ± 0.230 |
| 4 | 0.239 | 0.238 | 0.481 ± 0.175 |

Silhouettes are modest at every k (≈0.24–0.31) and bootstrap stability is mediocre (best ARI 0.53 at k=2, with high variance). **There is no strongly separated cluster structure.** Behavioral variation in this dataset is better described as continuous (a spread along the two dimensions of Section 15) than as a set of discrete phenotypes. This is an honest negative on the "distinct anxiety types" framing — johnson2025's high/medium/low grouping was a post-hoc tertile split, not a demonstration of natural gaps, and the P1 data show no such gaps.

### To the extent any partition exists, it cuts across strain

| k | ARI(cluster, strain) | chi-square p | Notes |
|---|---------------------|--------------|-------|
| 2 | −0.002 | 0.71 | both clusters ≈50% AB |
| 3 | +0.042 | 0.004 | one AB-enriched cluster (89% AB, n=18) |
| 4 | +0.050 | 0.002 | one AB-enriched cluster (93% AB, n=15) |

The dominant k=2 split is **completely independent of strain** (ARI ≈ 0, p = 0.71). The two clusters are a *bold-active* group (n=56: high pct_top z=+0.43, low latency −0.48, elevated transitions/x_sd/velocity; PC1 = +0.75) and a *shy-inactive* group (n=26: pct_top −0.93, latency +1.04, low activity; PC1 = −1.61). This split runs along **PC1 (general activity/boldness)** — the axis that does *not* separate strains — which is exactly why the clusters do not track genotype. The figure makes this geometric: the cluster boundary is roughly orthogonal to the strain gradient (which lies along PC2).

At k = 3–4 a small AB-enriched subgroup does emerge (chi-square significant), consistent with the PC2 strain effect of Section 15, but it contributes little to overall structure (ARI still ≈ 0.05).

### Interpretation

The primary axis of individual variation that any clustering latches onto is the strain-independent boldness/activity continuum (PC1), directly supporting rajput2022's finding that exploratory phenotypes cross strain rather than reproduce it. The strain difference (PC2) is real and large (Section 15, d=0.87) but is a more subtle axis that does not generate separable clusters. Practically: 5G and AB are not two behavioral "types"; they are overlapping distributions whose *centroids* differ along PC2, embedded in a population whose largest behavioral variation (boldness/activity) is shared across both strains. Discrete-cluster or "anxiety subtype" interpretations of this dataset are not supported.
