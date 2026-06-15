# Supplemental / exploratory analyses

These scripts are **not part of the manuscript** (`manuscript/`, a reliability-only
paper). They are retained as a record of exploratory work. They run from the repo
root like the `src/` scripts (`python supplemental/<script>.py`) and read the same
`data/features/` inputs; their outputs go to `output/supplemental/`.

## Why these are demoted, not in the paper

The dataset contains two strains (5G inbred, AB wild type), but **every testing
batch was strain-pure**, so strain is completely confounded with batch (testing
date, room, handler, per-batch distortion calibration; verified in
`src/11_qc_checks.py` §0). Any "strain" result is therefore descriptive of the
batches, not attributable to genotype, and is not publishable as a strain effect.
All analyses below center on strain comparisons or on multivariate structure whose
chief use was the (confounded) strain story, so they are kept here rather than in
the manuscript.

## Contents

| Script | What it does | Status |
|--------|--------------|--------|
| `trait_structure.py` | Among-individual correlations + PCA dimensionality of the 5 measures; 5G/AB difference on the principal axes | Exploratory. Strain-confounded; the among-individual correlation is also only an approximation (the bivariate mixed model does not separate within- from among-individual covariance). |
| `clustering.py` | Unsupervised k-means/Ward clustering of fish on the 5 measures; cluster-vs-strain association | Exploratory. Weak cluster structure; strain-confounded. |
| `ml.py` | Supervised 5G-vs-AB classification (feature-importance framing) | Exploratory. The "strain" label is confounded with batch, so high accuracy may reflect batch/calibration signatures. |

## What the manuscript kept

The only reliability content extracted from this line of work is the **univariate
short-term reliability of the five measures** (manuscript Table 2), now produced by
`src/09_measure_reliability.py`, which also writes `data/features/trait_sessions.csv`
(the input these supplemental scripts read).
