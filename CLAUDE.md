# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R-based statistical analysis project for clinical research evaluating preoperative nodule marking in minimally invasive thoracic surgery (MITS). Contains propensity score and multivariable regression analyses comparing outcomes between marked and unmarked nodules.

## Repository Structure

```
nodule-marking/
├── R/                    # Analysis scripts
├── docs/                 # Paper drafts and tables
└── data/                 # Data files (not tracked in git)
```

## Running Scripts

The script expects the data file in `data/`. Run from the project root directory:

```bash
Rscript R/outcomes_table_analysis.R
```

Or in R/RStudio:
```r
source("R/outcomes_table_analysis.R")
```

## Required R Packages

- **Data/modeling**: tidyverse, broom, survey, Hmisc
- **Visualization**: ggplot2 (via tidyverse)

## Data

Data files are in `data/` (git-ignored). Scripts expect `data/mits_data_20251212_reformatted CLEAN.csv`.

Key variables:
- `marked` - treatment variable (preoperative nodule marking yes/no)
- `surgical_approach` - RATS or VATS
- Covariates: age, is_f, nod_size, dist_from_pleur_prox, dist_from_pleur_distal, fev_perc_pred, dlco_perc_pred, fvc_perc_pred, smoking_status, prior_thoracic_radiation, prior_thoracic_surgery, is_solid
- Outcomes: persistent_air_leak_binary, conversion_to_open_binary, length_of_stay

## Analysis Scripts

**R/sensitivity_supplementary_tables.R** - Supplementary Tables 2 (analyzed vs excluded) and 3 (sensitivity: doubly robust, overlap weights, PS without PFTs). Sources the main script, so it also regenerates Tables 3 & 4.

**R/outcomes_table_analysis.R** - Generates publication tables 3 & 4

Pipeline:
- Filters to RATS/VATS; derives margin binaries, `is_metastatic`, `is_benign`
- Propensity score model with 12 covariates (see Data section)
- Common support trimming (`in_common_support` flag)
- ATE weights and stabilized weights; stabilized weights used in all outcome models
- Covariate balance via SMD (unweighted vs IPTW-weighted)
- Per outcome, two adjustment methods: IPTW-weighted (primary) + multivariable regression (sensitivity)
- Descriptives from the full sample; adjusted effects from the trimmed sample
- Analyzes RATS, VATS, and pooled cohorts (pooled adjusts for surgical approach)

Table 3 outcomes (surgical): length of stay, persistent air leak, incision-to-close time, localization failure, conversion to open.
Table 4 outcomes (pathology): first resection diagnostic, margins clear, margin distance, resection mass, final margin status, specimen mass, metastatic, benign.

Outputs: `surgical_outcomes_table_wide.csv`/`.md`, `pathology_table_wide.csv`/`.md`, `pathology_table_results.csv`, `propensity_score_distribution.png`, `love_plot_smd.png`.

### Known caveats (do not "fix" silently)
- **Residual imbalance**: after IPTW, 5/12 covariates remain |SMD| ≥ 0.1 (max ≈ 0.15: FVC, prior radiation, female). IPTW resolves the largest baseline imbalances but worsens a few smaller ones; the "balanced count" is unchanged before/after by coincidence, not because weighting failed. Documented limitation, not a bug.
- **MV separation**: the multivariable sensitivity models emit many `glm.fit: fitted probabilities 0/1` warnings from separation on sparse binary outcomes (e.g. conversion, metastatic). Those MV odds ratios are unstable. IPTW is the primary estimate; MV is sensitivity only.
