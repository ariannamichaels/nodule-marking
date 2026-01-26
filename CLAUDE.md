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

Scripts expect data file in `data/`. Run from the project root directory:

```bash
Rscript R/IPTW_complete_analysis.R
Rscript R/multivariable_logistic_regression_RATS_VATS.R
```

Or in R/RStudio:
```r
source("R/IPTW_complete_analysis.R")
```

## Required R Packages

- **Data/modeling**: tidyverse, broom, survey, tableone, Hmisc
- **Diagnostics**: pROC, car, ResourceSelection, DescTools
- **Visualization**: forestplot

## Data

Data files are in `data/` (git-ignored). Scripts expect `data/mits_data_20251212_reformatted CLEAN.csv`.

Key variables:
- `marked` - treatment variable (preoperative nodule marking yes/no)
- `surgical_approach` - RATS or VATS
- Covariates: age, is_f, nod_size, dist_from_pleur_prox, dist_from_pleur_distal, fev_perc_pred, dlco_perc_pred, fvc_perc_pred, smoking_status, prior_thoracic_radiation, prior_thoracic_surgery, is_solid
- Outcomes: persistent_air_leak_binary, conversion_to_open_binary, length_of_stay

## Analysis Scripts

**R/IPTW_complete_analysis.R** - Inverse Probability of Treatment Weighting
- Propensity score model with 12 covariates
- ATE and stabilized weights calculation
- Covariate balance assessment via SMD
- Survey design object for weighted outcome analysis
- Outputs: iptw_dataset_with_weights.csv, iptw_smd_balance.csv, propensity_score_model_results.csv

**R/multivariable_logistic_regression_RATS_VATS.R** - Predictors of marking
- Separate RATS vs VATS analyses
- Log-transformed variables: nod_size_log, dist_deep_log, brock_log
- Diagnostics: VIF, Hosmer-Lemeshow, pseudo R-squared, ROC/AUC
- Outputs: *_marking_predictors_results.csv, marking_prediction_roc_curves.png, marking_predictors_forest_plot.png

**R/pathology_table_analysis.R** - Pathology outcomes (Table 4)
- Analyzes RATS, VATS, and pooled cohorts
- Binary outcomes: diagnostic first resection, margin clear, metastatic, benign
- Continuous outcomes: margin distance, resection mass, specimen mass
- Both IPTW-weighted (primary) and multivariable-adjusted (sensitivity) estimates
- Common support trimming for IPTW
- Outputs: pathology_table_results.csv
