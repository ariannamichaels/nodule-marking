# ==============================================================================
# SENSITIVITY ANALYSES AND SUPPLEMENTARY TABLES (REVIEWER ROUND 2)
# ==============================================================================
#
# PURPOSE:
#   Responds to reviewer comments on missing data, trimming, and residual
#   imbalance in the IPTW analysis.
#
# OUTPUTS:
#   - supp_table_2_included_excluded.csv / _generated.md
#       Baseline characteristics: analyzed vs excluded (missing covariates,
#       outside common support)
#   - supp_table_3_sensitivity.csv / _generated.md
#       All Table 3/4 outcomes under four specifications
#
# SPECIFICATIONS:
#   1. Primary: IPTW (stabilized weights), as in outcomes_table_analysis.R
#   2. Doubly robust: primary + covariates with residual |SMD| >= 0.1
#   3. Overlap weights (marked: 1 - PS; unmarked: PS) on the trimmed sample
#   4. PS model without PFT covariates (FEV1, DLCO, FVC), refit and retrimmed
#
# USAGE:
#   Rscript R/sensitivity_supplementary_tables.R
#   (Run from project root; sources outcomes_table_analysis.R, which also
#   regenerates the main tables)
#
# ==============================================================================

invisible(capture.output(source("R/outcomes_table_analysis.R")))

library(tidyverse)
library(survey)

# ==============================================================================
# 1. SUPPLEMENTARY TABLE 2: INCLUDED VS EXCLUDED
# ==============================================================================

df <- df %>%
  mutate(
    complete_cov = complete.cases(across(all_of(ps_covariates))),
    analysis_group = case_when(
      complete_cov & in_common_support ~ "Analyzed",
      complete_cov ~ "Outside common support",
      TRUE ~ "Missing covariates"
    ),
    is_rats = as.numeric(approach == "RATS"),
    marked_num = as.numeric(marked)
  )

supp2_vars <- c(
  marked_num = "Marked", is_rats = "RATS approach",
  age = "Age, years", is_f = "Female", nod_size = "Nodule size, mm",
  dist_from_pleur_prox = "Proximal pleural distance, mm",
  dist_from_pleur_distal = "Distal pleural distance, mm",
  is_solid = "Solid nodule", smoking_status = "Smoking history",
  prior_thoracic_radiation = "Prior thoracic radiation",
  prior_thoracic_surgery = "Prior thoracic surgery",
  fev_perc_pred = "FEV1, % predicted", dlco_perc_pred = "DLCO, % predicted",
  fvc_perc_pred = "FVC, % predicted",
  persistent_air_leak_binary = "Persistent air leak",
  conversion_to_open_binary = "Conversion to open",
  length_of_stay = "Length of stay, days"
)

# Three decimals below 0.1 so values near 0.05 (e.g. 0.045) are not rounded across it
fmt_p <- function(p) {
  if (p < 0.001) "<.001" else sub("^0", "", sprintf(if (p < 0.1) "%.3f" else "%.2f", p))
}

fmt_cell <- function(x, binary) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  if (binary) sprintf("%d/%d (%.0f%%)", sum(x), length(x), 100 * mean(x))
  else sprintf("%.1f (%.1f) [n=%d]", mean(x), sd(x), length(x))
}

groups_s2 <- c("Analyzed", "Missing covariates", "Outside common support")
supp2 <- map_dfr(names(supp2_vars), function(v) {
  x <- df[[v]]
  binary <- all(na.omit(x) %in% c(0, 1))
  inc <- df$analysis_group == "Analyzed"
  # Analyzed vs all excluded
  p <- if (binary) fisher.test(table(inc, x))$p.value else wilcox.test(x[inc], x[!inc])$p.value
  cells <- map_chr(groups_s2, ~ fmt_cell(x[df$analysis_group == .x], binary))
  tibble(Characteristic = supp2_vars[[v]],
         !!!set_names(as.list(cells), groups_s2),
         P = fmt_p(p))
})

n_s2 <- table(df$analysis_group)[groups_s2]
write.csv(supp2, "supp_table_2_included_excluded.csv", row.names = FALSE)

md_s2 <- c(
  "## Supplementary Table 2. Baseline Characteristics of Analyzed vs Excluded Cases",
  "",
  paste0("| **Characteristic** | **Analyzed (N=", n_s2[1], ")** | **Missing covariates (N=",
         n_s2[2], ")** | **Outside common support (N=", n_s2[3], ")** | ***P*** |"),
  "|---|---|---|---|---|",
  pmap_chr(supp2, function(...) paste0("| ", paste(c(...), collapse = " | "), " |")),
  "",
  paste0("**Note:** Continuous variables as mean (SD) [n non-missing]; binary as n/N (%). ",
         "P compares analyzed vs all excluded cases (Fisher's exact or Wilcoxon rank-sum).")
)
writeLines(md_s2, "supp_table_2_included_excluded_generated.md")

# ==============================================================================
# 2. ALTERNATIVE WEIGHTS
# ==============================================================================

# Overlap weights on the trimmed sample
df <- df %>%
  mutate(overlap_weight = case_when(
    !in_common_support ~ NA_real_,
    marked ~ 1 - propensity_score,
    TRUE ~ propensity_score
  ))

# PS model without PFT covariates (PFTs are the main source of missingness)
ps_cov_nopft <- setdiff(ps_covariates, c("fev_perc_pred", "dlco_perc_pred", "fvc_perc_pred"))
df_nopft <- df %>% filter(complete.cases(across(all_of(ps_cov_nopft))))
ps_model_nopft <- glm(reformulate(ps_cov_nopft, "marked"), data = df_nopft, family = binomial())
df_nopft$ps_nopft <- fitted(ps_model_nopft)
df_nopft <- df_nopft %>%
  filter(ps_nopft >= min(ps_nopft[marked]), ps_nopft <= max(ps_nopft[!marked])) %>%
  mutate(nopft_weight = ifelse(marked, mean(marked) / ps_nopft,
                               (1 - mean(marked)) / (1 - ps_nopft)))

cat("No-PFT PS model: N in common support =", nrow(df_nopft), "\n")

# Covariates with residual |SMD| >= 0.1 after primary IPTW (for doubly robust)
dr_covariates <- smd_results$Variable[abs(smd_results$SMD_Weighted) >= 0.1]
cat("Doubly robust adjustment covariates:", paste(dr_covariates, collapse = ", "), "\n")

# ==============================================================================
# 3. SUPPLEMENTARY TABLE 3: SENSITIVITY ANALYSES
# ==============================================================================

outcome_labels <- c(
  length_of_stay = "Length of stay, days",
  persistent_air_leak_binary = "Persistent air leak",
  incision_to_close_mins = "Incision to closing time, min",
  localization_failure_binary = "Localization failure",
  conversion_to_open_binary = "Conversion to open",
  diagnostic_first_resection_binary = "First resection diagnostic",
  margin_clear_first_resection_binary = "First resection margins clear",
  margin_dist_initial_mm_clean = "First resection margin distance, mm",
  resection_mass = "First resection mass, g",
  overall_margin_status_binary = "Final margin clear",
  specimen_mass = "Final mass tissue resected, g",
  is_metastatic = "Metastases resections",
  is_benign = "Benign lobectomies"
)

fit_weighted <- function(data, outcome, weight, pooled, adjust = NULL) {
  data <- data[!is.na(data[[outcome]]) & !is.na(data[[weight]]), ]
  if (length(unique(data[[outcome]])) < 2) return(NULL)
  binary <- all(data[[outcome]] %in% c(0, 1))
  rhs <- c("marked", if (pooled) "approach", adjust)
  design <- svydesign(ids = ~1, weights = reformulate(weight), data = data)
  model <- tryCatch(
    suppressWarnings(svyglm(reformulate(rhs, outcome), design = design,
                            family = if (binary) quasibinomial() else gaussian())),
    error = function(e) NULL)
  if (is.null(model)) return(NULL)
  est <- c(coef(model)["markedTRUE"], suppressWarnings(confint(model)["markedTRUE", ]))
  if (binary) est <- exp(est)
  list(est = est, p = summary(model)$coefficients["markedTRUE", 4], binary = binary)
}

fmt_effect <- function(r) {
  if (is.null(r)) return("n/a")
  p <- fmt_p(r$p)
  if (r$binary) sprintf("OR: %.2f (%.2f-%.2f); P=%s", r$est[1], r$est[2], r$est[3], p)
  else sprintf("β: %.1f (%.1f to %.1f); P=%s", r$est[1], r$est[2], r$est[3], p)
}

cohorts <- list(RATS = "RATS", VATS = "VATS", Pooled = c("RATS", "VATS"))
supp3 <- map_dfr(names(cohorts), function(cn) {
  pooled <- cn == "Pooled"
  d_main <- df %>% filter(approach %in% cohorts[[cn]], in_common_support)
  d_nopft <- df_nopft %>% filter(approach %in% cohorts[[cn]])
  map_dfr(names(outcome_labels), function(o) {
    tibble(
      Cohort = cn,
      Outcome = outcome_labels[[o]],
      `Primary IPTW` = fmt_effect(fit_weighted(d_main, o, "stabilized_weight", pooled)),
      `Doubly robust` = fmt_effect(fit_weighted(d_main, o, "stabilized_weight", pooled, dr_covariates)),
      `Overlap weights` = fmt_effect(fit_weighted(d_main, o, "overlap_weight", pooled)),
      `PS without PFTs` = fmt_effect(fit_weighted(d_nopft, o, "nopft_weight", pooled))
    )
  })
})

write.csv(supp3, "supp_table_3_sensitivity.csv", row.names = FALSE)

md_s3 <- c(
  "## Supplementary Table 3. Sensitivity Analyses of Adjusted Outcome Estimates",
  "",
  paste0("| **Cohort** | **Outcome** | **Primary IPTW (N=", sum(df$in_common_support, na.rm = TRUE),
         ")** | **Doubly robust** | **Overlap weights** | **PS without PFTs (N=", nrow(df_nopft), ")** |"),
  "|---|---|---|---|---|---|",
  pmap_chr(supp3, function(...) paste0("| ", paste(c(...), collapse = " | "), " |")),
  "",
  paste0("**Note:** OR = odds ratio (binary outcomes); β = regression coefficient (continuous outcomes); ",
         "marked vs unmarked. Doubly robust models add covariates with residual |SMD| ≥ 0.1 after IPTW (",
         paste(dr_covariates, collapse = ", "), "). Overlap weights use the common-support sample. ",
         "The PS model without PFTs omits FEV1, DLCO, and FVC % predicted and is refit and retrimmed. ",
         "Pooled models adjust for surgical approach. n/a: no events in one group.")
)
writeLines(md_s3, "supp_table_3_sensitivity_generated.md")

cat("Supplementary tables saved: supp_table_2_*, supp_table_3_*\n")
