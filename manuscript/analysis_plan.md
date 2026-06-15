# Analysis Plan — Future Directions

Analyses not yet completed, organized by priority. Literature sources are from `lit/` and `vault/`.

---

## 1. Velocity reliability

**Status:** ✅ DONE — `src/09_trait_structure.py`, analysis.md Section 15. Velocity single-session R = 0.514 (α = 0.828); x_sd R = 0.521 (α = 0.860) — both as repeatable as pct_top. Benchmarked vs thomson2020/baker2018; flagged the speed-when-moving caveat.

**What to do:** Run the same Cronbach's α and ICC battery on mean velocity that was run on pct_top (Sections 1, 3, 4 of analysis.md). Include inter-run correlations and strain comparison.

**Literature context:** thomson2020 found distance traveled R = 0.46 and stationary time R = 0.39 in AB zebrafish. baker2018 found swimming speed is the most repeatable locomotor measure across strains (R = 0.40–0.59). Direct benchmarks for comparison.

---

## 2. Repeatability conditional on temporal trend (Runs 2–6 α)

**Status:** ✅ DONE — `src/08_repeatability.py`, analysis.md Section 14. Runs 2–6 α = 0.876 (vs 0.855 all runs); single-measure R rises 0.50 → 0.59. Found the Run 1 effect is individual-level reordering (R_A ≈ R_C), not a mean trend.

**What to do:** Compute Cronbach's α, inter-run correlations, and ICC using only Runs 2–6, explicitly framing this as repeatability of trait anxiety in a familiar environment. Report alongside the all-runs estimate.

**Literature context:** biro2015 distinguishes R_A (agreement, includes mean-level change), R_C (consistency, removes mean trend), and R|time (conditional on time as covariate). The Run 1 → Run 2+ mean shift means the reported α = 0.855 conflates shared novelty-response with individual trait consistency. thomson2020 found week 1→2 repeatability (R = 0.25) much lower than subsequent weeks — a direct parallel that supports formally separating the acclimation session from the trait sessions.

---

## 3. Cross-trait correlation structure (behavioral syndrome analysis)

**Status:** ✅ DONE — `src/09_trait_structure.py`, analysis.md Section 15. Two dissociable dimensions (Kaiser PC1 activity 43.5%, PC2 vertical-anxiety 26.8%) via raw + repeatability-corrected (disattenuated and bivariate-mixed-model) correlations. Headline: the strain difference lives on PC2 with d=0.87 (vs d=0.27 for univariate pct_top).

**What to do:** Compute the cross-trait correlation matrix of per-fish means (pct_top_mean, latency_mean, transitions_mean, x_sd_mean). Test whether a single latent factor fits or whether dissociable dimensions emerge (e.g., a factor analysis or at minimum a table of pairwise r values).

**Literature context:** reale2007 frames animal personality as a multivariate structure; oswald2013 reports genetic and phenotypic correlations among bold-shy components in zebrafish. blaser2012 shows that diving response and scototaxis (dark preference) reflect dissociable mechanisms — supporting the expectation that multiple behavioral dimensions exist even within a single test. Knowing whether x_sd and transitions collapse onto a single "activity" factor vs. are independent of pct_top and latency would clarify the dimensionality of the strain difference.

---

## 4. Comparison to published repeatability benchmarks

**Status:** ✅ DONE — `src/08_repeatability.py`, analysis.md Section 14. Single-measure R ≈ 0.50 (all runs) computed with Spearman-Brown bridge to α; placed against bell2009 (0.37), thomson2020 (0.39), johnson2025 (0.29).

**What to do:** Add a paragraph (in the results or discussion) comparing the P1 α values to the animal personality repeatability literature. Key reference point: bell2009 meta-analysis of 750+ estimates found mean R ≈ 0.37 for lab studies, with publication bias pushing reported values upward. Note that α = 0.855 based on 6 sessions is not directly comparable to single-interval R without applying the Spearman-Brown correction (6-measure composite vs. single measure reliability).

**Literature context:** bell2009, nakagawa2010, baker2018, thomson2020. Placing the P1 findings in this context strengthens the argument that pct_top is an unusually reliable behavioral measure.

---

## 5. Erratic movement as an independent anxiety metric

**Status:** Not computed. Raw position data supports it.

**What to do:** Compute angular velocity or direction-change frequency from consecutive detected frames (angle between successive displacement vectors). Define an "erratic movement" bout (e.g., direction change > 90°) and compute rate per session. Run reliability and strain comparison.

**Literature context:** egan2009 defines erratic movement as sharp direction changes and rapid darting, shows it increases with alarm pheromone and decreases with fluoxetine — distinct from time-in-zone measures. A fish can be anxious at the top (erratic) or calm at the bottom (slow bottom-hugging), so this captures something pct_top misses. Relevant caution: the frame-differencing tracker only fires on movement, so direction-change counts will underestimate true erratic bouts during high-freeze sessions.

---

## 6. Batch as a source of behavioral variance

**Status:** ✅ DONE — `src/08_repeatability.py`, analysis.md Section 14. Nested variance decomposition: batch = 1.9%, among-individual = 49.8%, within-individual = 48.2%. Batch-adjusted R unchanged — batch variance is negligible.

**What to do:** Add batch as a random effect in the mixed models for pct_top reliability and strain comparisons. Report the proportion of variance attributable to batch vs. individual vs. residual. This is especially relevant given that the slot-adjustment analysis already found a within-batch spatial gradient.

**Literature context:** shishis2022 shows housing density and tank size meaningfully affect zebrafish behavior. The 13 batches ran across Jan–Mar 2014 under potentially varying conditions (room temperature, time of year, experimenter). Modeling batch variance is a standard control in repeated-measures animal behavior studies.

---

## 7. Unsupervised clustering of behavioral phenotypes

**Status:** Current ML work is supervised (strain as label). No unsupervised analysis has been run.

**Status:** ✅ DONE — `src/10_clustering.py`, analysis.md Section 16. Cluster structure is weak (silhouette ≈0.26, stability ARI ≈0.53); no discrete phenotypes. The dominant k=2 split (bold-active vs shy-inactive, along PC1) is strain-independent (ARI≈0, p=0.71) — clusters cut across strain, per rajput2022. A small AB-enriched subgroup appears at k=3–4 but contributes little. Discrete "anxiety subtype" interpretations are not supported.

**Literature context:** rajput2022 identified four behavioral clusters from 3D swim traces (bold, shy, wall-huggers, active explorers) that cut across strains. Given the large within-strain variance visible in the P1 data and the high LOOCV accuracy of the supervised classifier (~78%), it is plausible that natural phenotypic clusters exist that do not map cleanly onto genotype. johnson2025 took a similar approach, post-hoc grouping fish into high/medium/low anxiety categories that were stable across age.

---

## 8. Acclimation framing for Run 1 (literature parallel)

**Status:** ✅ DONE — `src/08_repeatability.py`, analysis.md Section 14. thomson2020 week-1→2 parallel made explicit; r(Run 1, trait composite) = 0.411 vs 0.584 within Runs 2–6 quantifies the acclimation pattern.

**What to do:** Add explicit discussion citing thomson2020 as a published parallel — they found week 1→2 repeatability was R = 0.25 and interpreted it as tank acclimation before stable trait behavior emerges. The P1 Run 1 pattern (r = 0.22–0.47 with later runs vs. r = 0.46–0.68 among Runs 2–6) is quantitatively similar. This reframes Run 1 not as a measurement anomaly but as a well-documented acclimation phenomenon, strengthening the argument for treating it separately.

---

## Dropped direction: sex differences

Sex of individual fish was not recorded in the P1 study and cannot be recovered. The literature (johnson2025, thomson2020, mustafa2019) finds substantial sex effects on both mean behavior and repeatability. This is a limitation of the dataset, not an addressable analysis. Note it as a limitation in the manuscript.
