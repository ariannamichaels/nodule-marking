# ==============================================================================
# SURGICAL AND PATHOLOGY OUTCOMES TABLE ANALYSIS
# ==============================================================================
#
# PURPOSE:
#   Generates publication-ready tables comparing outcomes between marked and
#   unmarked nodules in minimally invasive thoracic surgery (MITS).
#
# OUTPUTS:
#   - surgical_outcomes_table_wide.csv / .md  (Table 3: Surgical Outcomes)
#   - pathology_table_wide.csv / .md          (Table 4: Pathology Characteristics)
#   - propensity_score_distribution.png       (PS overlap diagnostic plot)
#   - pathology_table_results.csv             (detailed results for all analyses)
#
# METHODS:
#   1. Propensity score model fitted with 7 covariates (age, sex, nodule size,
#      pleural distance, consistency, smoking status, emphysema)
#   2. Common support trimming applied (excludes cases outside PS overlap region)
#   3. IPTW weights calculated for Average Treatment Effect (ATE) estimation
#   4. Descriptive statistics computed from FULL sample (all cases with outcome)
#   5. Adjusted effects computed from TRIMMED sample (common support only):
#      - Primary: IPTW-weighted regression
#      - Sensitivity: Multivariable-adjusted regression
#   6. Pooled analyses include surgical approach (RATS vs VATS) as covariate
#
# TABLE 3 OUTCOMES (Surgical):
#   - Length of stay, persistent air leak, incision-to-close time
#   - Localization failure, conversion to open
#
# TABLE 4 OUTCOMES (Pathology):
#   - First resection diagnostic, margins clear, margin distance, mass
#   - Final margin status, specimen mass
#   - Metastatic resections, benign lobectomies
#
# USAGE:
#   Rscript R/outcomes_table_analysis.R
#   (Run from project root; expects data in data/ directory)
#
# ==============================================================================

library(tidyverse)
library(broom)
library(survey)

cat(rep("=", 80), "\n", sep="")
cat("SURGICAL AND PATHOLOGY OUTCOMES TABLE ANALYSIS\n")
cat("Tables 3 & 4: Surgical and Pathology Outcomes by Marking Status\n")
cat(rep("=", 80), "\n\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

df <- read.csv("data/mits_data_20251212_reformatted CLEAN.csv", stringsAsFactors = FALSE)

# Standardize surgical approach and filter to RATS/VATS
df <- df %>%
  mutate(approach = toupper(trimws(surgical_approach))) %>%
  filter(approach %in% c("RATS", "VATS"))

cat("Dataset loaded:", nrow(df), "cases\n\n")

# ==============================================================================
# 2. DATA CLEANING - CREATE ANALYSIS VARIABLES
# ==============================================================================

cat(rep("=", 80), "\n", sep="")
cat("DATA CLEANING\n")
cat(rep("=", 80), "\n\n")

df <- df %>%
  mutate(
    # Ensure marked is logical
    marked = as.logical(marked),

    # Clean margin_clear_first_resection to binary
    margin_clear_first_resection_binary = case_when(
      tolower(margin_clear_first_resection) == "clear" ~ 1,
      tolower(margin_clear_first_resection) == "not clear" ~ 0,
      margin_clear_first_resection == "TRUE" ~ 1,
      margin_clear_first_resection == "FALSE" ~ 0,
      grepl("clear", tolower(margin_clear_first_resection)) ~ 1,
      TRUE ~ NA_real_
    ),

    # Clean overall_margin_status to binary
    overall_margin_status_binary = case_when(
      tolower(overall_margin_status) == "clear" ~ 1,
      tolower(overall_margin_status) == "not clear" ~ 0,
      overall_margin_status == "TRUE" ~ 1,
      overall_margin_status == "FALSE" ~ 0,
      grepl("clear", tolower(overall_margin_status)) ~ 1,
      TRUE ~ NA_real_
    ),

    # Clean margin_dist_initial_mm to numeric
    margin_dist_initial_mm_clean = suppressWarnings(as.numeric(margin_dist_initial_mm)),

    # Derive is_metastatic from final_pathology
    is_metastatic = as.numeric(grepl("metasta", tolower(final_pathology))),

    # Derive is_benign from final_pathology
    is_benign = as.numeric(grepl("benign|granuloma|hamartoma|inflammatory|fibrosis",
                                  tolower(final_pathology)) &
                            !grepl("adenocarcinoma|carcinoma|malignant|metasta",
                                   tolower(final_pathology)))
  )

cat("Variables created:\n")
cat("  - margin_clear_first_resection_binary\n")
cat("  - overall_margin_status_binary\n")
cat("  - margin_dist_initial_mm_clean\n")
cat("  - is_metastatic\n")
cat("  - is_benign\n\n")

# ==============================================================================
# 3. SAMPLE SIZE SUMMARY
# ==============================================================================

cat(rep("=", 80), "\n", sep="")
cat("SAMPLE SIZES\n")
cat(rep("=", 80), "\n\n")

sample_sizes <- df %>%
  group_by(approach, marked) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = c(approach, marked), values_from = n)

cat("By approach and marking status:\n")
print(sample_sizes)
cat("\n")

cat("RATS: Marked =", sum(df$approach == "RATS" & df$marked),
    ", Unmarked =", sum(df$approach == "RATS" & !df$marked), "\n")
cat("VATS: Marked =", sum(df$approach == "VATS" & df$marked),
    ", Unmarked =", sum(df$approach == "VATS" & !df$marked), "\n")
cat("Pooled: Marked =", sum(df$marked),
    ", Unmarked =", sum(!df$marked), "\n\n")

# ==============================================================================
# 4. PROPENSITY SCORE MODEL AND IPTW WEIGHTS
# ==============================================================================

cat(rep("=", 80), "\n", sep="")
cat("PROPENSITY SCORE MODEL\n")
cat(rep("=", 80), "\n\n")

# Define covariates (subset relevant to pathology outcomes)
ps_covariates <- c("age", "is_f", "nod_size", "dist_from_pleur_distal",
                   "is_solid", "smoking_status", "emphysema_on_ct")

# Create complete case dataset for PS model
df_ps <- df %>%
  select(marked, all_of(ps_covariates)) %>%
  na.omit()

cat("Complete cases for PS model:", nrow(df_ps), "/", nrow(df), "\n")

# Fit propensity score model
ps_formula <- as.formula(paste("marked ~", paste(ps_covariates, collapse = " + ")))
ps_model <- glm(ps_formula, data = df_ps, family = binomial(link = "logit"))

cat("PS model fitted successfully\n\n")

# Get propensity scores for full dataset
df$propensity_score <- predict(ps_model, newdata = df, type = "response")

# Common support trimming
ps_marked_min <- min(df$propensity_score[df$marked & !is.na(df$propensity_score)])
ps_unmarked_max <- max(df$propensity_score[!df$marked & !is.na(df$propensity_score)])

cat("Propensity Score Range:\n")
cat("  Marked: [", round(min(df$propensity_score[df$marked], na.rm=TRUE), 3),
    ", ", round(max(df$propensity_score[df$marked], na.rm=TRUE), 3), "]\n", sep="")
cat("  Unmarked: [", round(min(df$propensity_score[!df$marked], na.rm=TRUE), 3),
    ", ", round(max(df$propensity_score[!df$marked], na.rm=TRUE), 3), "]\n", sep="")
cat("  Common support: [", round(ps_marked_min, 3), ", ", round(ps_unmarked_max, 3), "]\n\n", sep="")

# Trim to common support
df <- df %>%
  mutate(
    in_common_support = propensity_score >= ps_marked_min &
                        propensity_score <= ps_unmarked_max &
                        !is.na(propensity_score)
  )

n_trimmed <- sum(!df$in_common_support, na.rm = TRUE)
cat("Cases trimmed for common support:", n_trimmed, "\n")
cat("Cases remaining:", sum(df$in_common_support, na.rm = TRUE), "\n\n")

# Calculate IPTW weights (ATE weights)
df <- df %>%
  mutate(
    iptw_weight = case_when(
      !in_common_support ~ NA_real_,
      marked ~ 1 / propensity_score,
      TRUE ~ 1 / (1 - propensity_score)
    )
  )

cat("IPTW Weight Summary (trimmed sample):\n")
cat("  Mean:", round(mean(df$iptw_weight, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(df$iptw_weight, na.rm = TRUE), 3), "\n")
cat("  Range: [", round(min(df$iptw_weight, na.rm = TRUE), 3), ", ",
    round(max(df$iptw_weight, na.rm = TRUE), 3), "]\n\n", sep="")

# Create propensity score histogram
cat("Creating propensity score distribution plot...\n")

png("propensity_score_distribution.png", width = 10, height = 6, units = "in", res = 300)

# Prepare data for plotting
ps_data <- df %>%
  filter(!is.na(propensity_score)) %>%
  mutate(Group = ifelse(marked, "Marked", "Unmarked"))

# Get sample sizes for legend labels
n_marked <- sum(ps_data$Group == "Marked")
n_unmarked <- sum(ps_data$Group == "Unmarked")

# Create legend labels with sample sizes
ps_data <- ps_data %>%
  mutate(Group = ifelse(marked,
                        paste0("Marked (n=", n_marked, ")"),
                        paste0("Unmarked (n=", n_unmarked, ")")))

# Create overlapping histograms
ggplot(ps_data, aes(x = propensity_score, fill = Group)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, alpha = 0.75, position = "identity", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = ps_marked_min, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = ps_unmarked_max, linetype = "dashed", color = "black", linewidth = 1) +
  annotate("text", x = ps_marked_min, y = Inf, label = paste0("Min marked\n(", round(ps_marked_min, 3), ")"),
           vjust = 1.5, hjust = -0.1, size = 3, color = "black") +
  annotate("text", x = ps_unmarked_max, y = Inf, label = paste0("Max unmarked\n(", round(ps_unmarked_max, 3), ")"),
           vjust = 1.5, hjust = 1.1, size = 3, color = "black") +
  scale_fill_manual(values = setNames(c("#E41A1C", "#377EB8"),
                                       c(paste0("Marked (n=", n_marked, ")"),
                                         paste0("Unmarked (n=", n_unmarked, ")")))) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Propensity Score Distribution by Treatment Group",
    x = "Propensity Score",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = c(0.9, 0.85),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = NA)
  )

dev.off()
cat("Saved: propensity_score_distribution.png\n\n")

# ==============================================================================
# 5. ANALYSIS FUNCTIONS
# ==============================================================================

# Function to analyze binary outcome
analyze_binary <- function(data, outcome_var, group_label, include_approach = FALSE) {

  # Full sample for descriptive stats (only filter on non-missing outcome)
  data_full <- data %>%
    filter(!is.na(.data[[outcome_var]]))

  # Trimmed sample for adjusted analyses (common support)
  data_trimmed <- data %>%
    filter(!is.na(.data[[outcome_var]]) & in_common_support)

  if (nrow(data_full) == 0) {
    return(list(
      group = group_label,
      outcome = outcome_var,
      n_marked = NA, pct_marked = NA,
      n_unmarked = NA, pct_unmarked = NA,
      p_unadjusted = NA,
      or_iptw = NA, ci_low_iptw = NA, ci_high_iptw = NA, p_iptw = NA,
      or_mv = NA, ci_low_mv = NA, ci_high_mv = NA, p_mv = NA
    ))
  }

  # Descriptive statistics from FULL sample (not trimmed)
  stats_marked <- data_full %>% filter(marked) %>%
    summarise(n = sum(.data[[outcome_var]], na.rm = TRUE),
              total = n(),
              pct = round(100 * n / total, 0))

  stats_unmarked <- data_full %>% filter(!marked) %>%
    summarise(n = sum(.data[[outcome_var]], na.rm = TRUE),
              total = n(),
              pct = round(100 * n / total, 0))

  # Unadjusted p-value from FULL sample (Fisher's exact test)
  tbl <- table(data_full$marked, data_full[[outcome_var]])
  p_unadj <- tryCatch({
    fisher.test(tbl)$p.value
  }, error = function(e) NA)

  # IPTW-weighted analysis from TRIMMED sample (include approach dummy for pooled)
  iptw_formula <- if (include_approach) {
    as.formula(paste(outcome_var, "~ marked + approach"))
  } else {
    as.formula(paste(outcome_var, "~ marked"))
  }

  iptw_result <- tryCatch({
    survey_design <- svydesign(ids = ~1, weights = ~iptw_weight, data = data_trimmed)
    model <- svyglm(iptw_formula, design = survey_design, family = quasibinomial())
    or <- exp(coef(model)["markedTRUE"])
    ci <- exp(confint(model)["markedTRUE", ])
    p <- summary(model)$coefficients["markedTRUE", "Pr(>|t|)"]
    list(or = or, ci_low = ci[1], ci_high = ci[2], p = p)
  }, error = function(e) list(or = NA, ci_low = NA, ci_high = NA, p = NA))

  # Multivariable adjusted analysis from TRIMMED sample (include approach dummy for pooled)
  mv_covars <- c("age", "is_f", "nod_size", "dist_from_pleur_distal", "is_solid", "smoking_status")
  if (include_approach) {
    mv_covars <- c(mv_covars, "approach")
  }
  mv_formula <- as.formula(paste(outcome_var, "~ marked +", paste(mv_covars, collapse = " + ")))

  mv_result <- tryCatch({
    model <- glm(mv_formula, data = data_trimmed, family = binomial())
    tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
    marked_row <- tidy_model %>% filter(term == "markedTRUE")
    list(or = marked_row$estimate, ci_low = marked_row$conf.low,
         ci_high = marked_row$conf.high, p = marked_row$p.value)
  }, error = function(e) list(or = NA, ci_low = NA, ci_high = NA, p = NA))

  return(list(
    group = group_label,
    outcome = outcome_var,
    n_marked = stats_marked$n, pct_marked = stats_marked$pct,
    total_marked = stats_marked$total,
    n_unmarked = stats_unmarked$n, pct_unmarked = stats_unmarked$pct,
    total_unmarked = stats_unmarked$total,
    p_unadjusted = p_unadj,
    or_iptw = iptw_result$or, ci_low_iptw = iptw_result$ci_low,
    ci_high_iptw = iptw_result$ci_high, p_iptw = iptw_result$p,
    or_mv = mv_result$or, ci_low_mv = mv_result$ci_low,
    ci_high_mv = mv_result$ci_high, p_mv = mv_result$p
  ))
}

# Function to analyze continuous outcome
analyze_continuous <- function(data, outcome_var, group_label, include_approach = FALSE) {

  # Full sample for descriptive stats (only filter on non-missing outcome)
  data_full <- data %>%
    filter(!is.na(.data[[outcome_var]]))

  # Trimmed sample for adjusted analyses (common support)
  data_trimmed <- data %>%
    filter(!is.na(.data[[outcome_var]]) & in_common_support)

  if (nrow(data_full) == 0) {
    return(list(
      group = group_label,
      outcome = outcome_var,
      mean_marked = NA, sd_marked = NA,
      mean_unmarked = NA, sd_unmarked = NA,
      p_unadjusted = NA,
      beta_iptw = NA, ci_low_iptw = NA, ci_high_iptw = NA, p_iptw = NA,
      beta_mv = NA, ci_low_mv = NA, ci_high_mv = NA, p_mv = NA
    ))
  }

  # Descriptive statistics from FULL sample (not trimmed)
  stats_marked <- data_full %>% filter(marked) %>%
    summarise(mean = mean(.data[[outcome_var]], na.rm = TRUE),
              sd = sd(.data[[outcome_var]], na.rm = TRUE),
              n = n())

  stats_unmarked <- data_full %>% filter(!marked) %>%
    summarise(mean = mean(.data[[outcome_var]], na.rm = TRUE),
              sd = sd(.data[[outcome_var]], na.rm = TRUE),
              n = n())

  # Unadjusted p-value from FULL sample (Wilcoxon test for robustness)
  p_unadj <- tryCatch({
    wilcox.test(data_full[[outcome_var]] ~ data_full$marked)$p.value
  }, error = function(e) NA)

  # IPTW-weighted analysis from TRIMMED sample (include approach dummy for pooled)
  iptw_formula <- if (include_approach) {
    as.formula(paste(outcome_var, "~ marked + approach"))
  } else {
    as.formula(paste(outcome_var, "~ marked"))
  }

  iptw_result <- tryCatch({
    survey_design <- svydesign(ids = ~1, weights = ~iptw_weight, data = data_trimmed)
    model <- svyglm(iptw_formula, design = survey_design)
    beta <- coef(model)["markedTRUE"]
    ci <- confint(model)["markedTRUE", ]
    p <- summary(model)$coefficients["markedTRUE", "Pr(>|t|)"]
    list(beta = beta, ci_low = ci[1], ci_high = ci[2], p = p)
  }, error = function(e) list(beta = NA, ci_low = NA, ci_high = NA, p = NA))

  # Multivariable adjusted analysis from TRIMMED sample (include approach dummy for pooled)
  mv_covars <- c("age", "is_f", "nod_size", "dist_from_pleur_distal", "is_solid", "smoking_status")
  if (include_approach) {
    mv_covars <- c(mv_covars, "approach")
  }
  mv_formula <- as.formula(paste(outcome_var, "~ marked +", paste(mv_covars, collapse = " + ")))

  mv_result <- tryCatch({
    model <- lm(mv_formula, data = data_trimmed)
    tidy_model <- tidy(model, conf.int = TRUE)
    marked_row <- tidy_model %>% filter(term == "markedTRUE")
    list(beta = marked_row$estimate, ci_low = marked_row$conf.low,
         ci_high = marked_row$conf.high, p = marked_row$p.value)
  }, error = function(e) list(beta = NA, ci_low = NA, ci_high = NA, p = NA))

  return(list(
    group = group_label,
    outcome = outcome_var,
    mean_marked = stats_marked$mean, sd_marked = stats_marked$sd,
    n_marked = stats_marked$n,
    mean_unmarked = stats_unmarked$mean, sd_unmarked = stats_unmarked$sd,
    n_unmarked = stats_unmarked$n,
    p_unadjusted = p_unadj,
    beta_iptw = iptw_result$beta, ci_low_iptw = iptw_result$ci_low,
    ci_high_iptw = iptw_result$ci_high, p_iptw = iptw_result$p,
    beta_mv = mv_result$beta, ci_low_mv = mv_result$ci_low,
    ci_high_mv = mv_result$ci_high, p_mv = mv_result$p
  ))
}

# ==============================================================================
# 6. HELPER FUNCTIONS FOR FORMATTING
# ==============================================================================

format_p_table <- function(p) {
  if (is.na(p)) return("n/a")
  if (p < 0.001) return("<.001")
  if (p > 0.99) return(">.9")
  paste0(".", substr(sprintf("%.2f", p), 3, 4))
}

format_or_table <- function(or, ci_low, ci_high) {
  if (is.na(or)) return("OR: [not estimable]")
  paste0("OR: ", round(or, 2), " (", round(ci_low, 2), "-", round(ci_high, 2), ")")
}

format_beta_table <- function(beta, ci_low, ci_high) {
  if (is.na(beta)) return("--")
  paste0("β: ", round(beta, 1), " (", round(ci_low, 1), " to ", round(ci_high, 1), ")")
}

format_binary_stat <- function(n, pct) {
  if (is.na(n)) return("--")
  paste0(n, " (", pct, ")")
}

format_continuous_stat <- function(mean_val, sd_val) {
  if (is.na(mean_val)) return("--")
  paste0(round(mean_val, 1), " (", round(sd_val, 1), ")")
}

# ==============================================================================
# 7. RUN TABLE 3 ANALYSES (SURGICAL OUTCOMES)
# ==============================================================================

cat(rep("=", 80), "\n", sep="")
cat("TABLE 3: SURGICAL CHARACTERISTICS & OUTCOMES ANALYSIS\n")
cat(rep("=", 80), "\n\n")

# Define Table 3 outcomes
table3_binary_outcomes <- c(
  "persistent_air_leak_binary",
  "localization_failure_binary",
  "conversion_to_open_binary"
)

table3_continuous_outcomes <- c(
  "length_of_stay",
  "incision_to_close_mins"
)

# Store Table 3 results
table3_results <- list()

# Define analysis groups (used for both Table 3 and Table 4)
groups <- list(
  list(filter_expr = quote(approach == "RATS"), label = "RATS", include_approach = FALSE),
  list(filter_expr = quote(approach == "VATS"), label = "VATS", include_approach = FALSE),
  list(filter_expr = quote(TRUE), label = "Pooled", include_approach = TRUE)
)

# Analyze Table 3 outcomes by group
for (grp in groups) {
  cat("\n", rep("-", 60), "\n", sep="")
  cat("Table 3 - ", grp$label, " ANALYSIS\n", sep="")
  if (grp$include_approach) {
    cat("(adjusted for surgical approach: RATS vs VATS)\n")
  }
  cat(rep("-", 60), "\n\n")

  data_subset <- df %>% filter(eval(grp$filter_expr))

  # Binary outcomes
  for (outcome in table3_binary_outcomes) {
    result <- analyze_binary(data_subset, outcome, grp$label, include_approach = grp$include_approach)
    table3_results[[paste(grp$label, outcome, sep = "_")]] <- result

    cat(outcome, ":\n")
    cat("  Marked: ", result$n_marked, "/", result$total_marked,
        " (", result$pct_marked, "%)\n", sep="")
    cat("  Unmarked: ", result$n_unmarked, "/", result$total_unmarked,
        " (", result$pct_unmarked, "%)\n", sep="")
    cat("  Unadjusted p:", format.pval(result$p_unadjusted, digits = 2), "\n")
    cat("  IPTW OR:", round(result$or_iptw, 2),
        " (", round(result$ci_low_iptw, 2), "-", round(result$ci_high_iptw, 2), ")",
        " p=", format.pval(result$p_iptw, digits = 2), "\n\n", sep="")
  }

  # Continuous outcomes
  for (outcome in table3_continuous_outcomes) {
    result <- analyze_continuous(data_subset, outcome, grp$label, include_approach = grp$include_approach)
    table3_results[[paste(grp$label, outcome, sep = "_")]] <- result

    cat(outcome, ":\n")
    cat("  Marked: ", round(result$mean_marked, 1), " (", round(result$sd_marked, 1),
        ") n=", result$n_marked, "\n", sep="")
    cat("  Unmarked: ", round(result$mean_unmarked, 1), " (", round(result$sd_unmarked, 1),
        ") n=", result$n_unmarked, "\n", sep="")
    cat("  Unadjusted p:", format.pval(result$p_unadjusted, digits = 2), "\n")
    cat("  IPTW beta:", round(result$beta_iptw, 2),
        " (", round(result$ci_low_iptw, 2), " to ", round(result$ci_high_iptw, 2), ")",
        " p=", format.pval(result$p_iptw, digits = 2), "\n\n", sep="")
  }
}

# ==============================================================================
# 8. EXPORT TABLE 3 (SURGICAL OUTCOMES)
# ==============================================================================

# Build Table 3 wide format
cat("\nCreating Table 3 wide-format table...\n")

table3_outcome_order <- c(
  "length_of_stay",
  "persistent_air_leak_binary",
  "incision_to_close_mins",
  "localization_failure_binary",
  "conversion_to_open_binary"
)

table3_outcome_labels <- c(
  "length_of_stay" = "Mean length of stay, days (SD)",
  "persistent_air_leak_binary" = "Persistent air leak (%)",
  "incision_to_close_mins" = "Incision to closing time, min (SD)",
  "localization_failure_binary" = "Localization failure (%)",
  "conversion_to_open_binary" = "Conversion to open (%)"
)

table3_wide <- data.frame(Characteristic = character(), stringsAsFactors = FALSE)

for (outcome in table3_outcome_order) {
  is_binary <- outcome %in% table3_binary_outcomes

  row_data <- list(Characteristic = table3_outcome_labels[outcome])

  for (grp in c("RATS", "VATS", "Pooled")) {
    result <- table3_results[[paste(grp, outcome, sep = "_")]]

    if (is_binary) {
      row_data[[paste0(grp, "_Marked")]] <- format_binary_stat(result$n_marked, result$pct_marked)
      row_data[[paste0(grp, "_Unmarked")]] <- format_binary_stat(result$n_unmarked, result$pct_unmarked)
      row_data[[paste0(grp, "_P")]] <- format_p_table(result$p_unadjusted)
      row_data[[paste0(grp, "_Adjusted_Effect")]] <- format_or_table(result$or_iptw, result$ci_low_iptw, result$ci_high_iptw)
      row_data[[paste0(grp, "_Adj_P")]] <- format_p_table(result$p_iptw)
    } else {
      row_data[[paste0(grp, "_Marked")]] <- format_continuous_stat(result$mean_marked, result$sd_marked)
      row_data[[paste0(grp, "_Unmarked")]] <- format_continuous_stat(result$mean_unmarked, result$sd_unmarked)
      row_data[[paste0(grp, "_P")]] <- format_p_table(result$p_unadjusted)
      row_data[[paste0(grp, "_Adjusted_Effect")]] <- format_beta_table(result$beta_iptw, result$ci_low_iptw, result$ci_high_iptw)
      row_data[[paste0(grp, "_Adj_P")]] <- format_p_table(result$p_iptw)
    }
  }

  table3_wide <- bind_rows(table3_wide, as.data.frame(row_data, stringsAsFactors = FALSE))
}

# Rename columns
colnames(table3_wide) <- c(
  "Characteristic",
  "RATS Marked (N=38)", "RATS Unmarked (N=45)", "RATS P", "RATS Adjusted Effect (95% CI)", "RATS Adj. P",
  "VATS Marked (N=33)", "VATS Unmarked (N=75)", "VATS P", "VATS Adjusted Effect (95% CI)", "VATS Adj. P",
  "Pooled Marked (N=71)", "Pooled Unmarked (N=120)", "Pooled P", "Pooled Adjusted Effect (95% CI)", "Pooled Adj. P"
)

write.csv(table3_wide, "surgical_outcomes_table_wide.csv", row.names = FALSE)
cat("Table 3 CSV saved to: surgical_outcomes_table_wide.csv\n")

# Generate Table 3 markdown
md3_header <- paste0(
  "| **Characteristic** | **RATS Marked (N=38)** | **RATS Unmarked (N=45)** | **RATS *P*** | **RATS Adjusted Effect (95% CI)** | **RATS Adj. *P*** | ",
  "**VATS Marked (N=33)** | **VATS Unmarked (N=75)** | **VATS *P*** | **VATS Adjusted Effect (95% CI)** | **VATS Adj. *P*** | ",
  "**Pooled Marked (N=71)** | **Pooled Unmarked (N=120)** | **Pooled *P*** | **Pooled Adjusted Effect (95% CI)** | **Pooled Adj. *P*** |"
)

md3_separator <- "|-------------------|----------------------|------------------------|------------|------------------------------|-----------------|----------------------|------------------------|------------|------------------------------|-----------------|------------------------|---------------------------|--------------|-------------------------------------|-------------------|"

md3_lines <- c("## Table 3. Surgical Characteristics & Outcomes", "", md3_header, md3_separator)

for (i in 1:nrow(table3_wide)) {
  row <- table3_wide[i, ]
  md_row <- paste0("| ", paste(row, collapse = " | "), " |")
  md3_lines <- c(md3_lines, md_row)
}

md3_lines <- c(md3_lines, "",
               "**Note:** Adjusted effects are from IPTW-weighted regression models. OR = Odds Ratio for binary outcomes; β = regression coefficient for continuous outcomes. Pooled analysis adjusted for surgical approach (RATS vs VATS). Common support trimming applied.")

writeLines(md3_lines, "surgical_outcomes_table_wide_generated.md")
cat("Table 3 markdown saved to: surgical_outcomes_table_wide_generated.md\n\n")

# ==============================================================================
# 9. RUN TABLE 4 ANALYSES (PATHOLOGY OUTCOMES)
# ==============================================================================

cat(rep("=", 80), "\n", sep="")
cat("TABLE 4: PATHOLOGY OUTCOMES ANALYSIS\n")
cat(rep("=", 80), "\n\n")

# Define pathology outcomes
binary_outcomes <- c(
  "diagnostic_first_resection_binary",
  "margin_clear_first_resection_binary",
  "overall_margin_status_binary",
  "is_metastatic",
  "is_benign"
)

continuous_outcomes <- c(
  "margin_dist_initial_mm_clean",
  "mass_initial_wedges_g",
  "specimen_mass"
)

# Store pathology results
all_results <- list()

# Analyze pathology outcomes by group
for (grp in groups) {
  cat("\n", rep("-", 80), "\n", sep="")
  cat(grp$label, " ANALYSIS\n", sep="")
  if (grp$include_approach) {
    cat("(adjusted for surgical approach: RATS vs VATS)\n")
  }
  cat(rep("-", 80), "\n\n")

  data_subset <- df %>% filter(eval(grp$filter_expr))

  cat("Sample: Marked =", sum(data_subset$marked & data_subset$in_common_support, na.rm=TRUE),
      ", Unmarked =", sum(!data_subset$marked & data_subset$in_common_support, na.rm=TRUE), "\n\n")

  # Binary outcomes
  for (outcome in binary_outcomes) {
    result <- analyze_binary(data_subset, outcome, grp$label, include_approach = grp$include_approach)
    all_results[[paste(grp$label, outcome, sep = "_")]] <- result

    cat(outcome, ":\n")
    cat("  Marked: ", result$n_marked, "/", result$total_marked,
        " (", result$pct_marked, "%)\n", sep="")
    cat("  Unmarked: ", result$n_unmarked, "/", result$total_unmarked,
        " (", result$pct_unmarked, "%)\n", sep="")
    cat("  Unadjusted p:", format.pval(result$p_unadjusted, digits = 2), "\n")
    cat("  IPTW OR:", round(result$or_iptw, 2),
        " (", round(result$ci_low_iptw, 2), "-", round(result$ci_high_iptw, 2), ")",
        " p=", format.pval(result$p_iptw, digits = 2), "\n", sep="")
    cat("  MV OR:", round(result$or_mv, 2),
        " (", round(result$ci_low_mv, 2), "-", round(result$ci_high_mv, 2), ")",
        " p=", format.pval(result$p_mv, digits = 2), "\n\n", sep="")
  }

  # Continuous outcomes
  for (outcome in continuous_outcomes) {
    result <- analyze_continuous(data_subset, outcome, grp$label, include_approach = grp$include_approach)
    all_results[[paste(grp$label, outcome, sep = "_")]] <- result

    cat(outcome, ":\n")
    cat("  Marked: ", round(result$mean_marked, 1), " (", round(result$sd_marked, 1),
        ") n=", result$n_marked, "\n", sep="")
    cat("  Unmarked: ", round(result$mean_unmarked, 1), " (", round(result$sd_unmarked, 1),
        ") n=", result$n_unmarked, "\n", sep="")
    cat("  Unadjusted p:", format.pval(result$p_unadjusted, digits = 2), "\n")
    cat("  IPTW beta:", round(result$beta_iptw, 2),
        " (", round(result$ci_low_iptw, 2), " to ", round(result$ci_high_iptw, 2), ")",
        " p=", format.pval(result$p_iptw, digits = 2), "\n", sep="")
    cat("  MV beta:", round(result$beta_mv, 2),
        " (", round(result$ci_low_mv, 2), " to ", round(result$ci_high_mv, 2), ")",
        " p=", format.pval(result$p_mv, digits = 2), "\n\n", sep="")
  }
}

# ==============================================================================
# 10. EXPORT TABLE 4 (PATHOLOGY)
# ==============================================================================

cat(rep("=", 80), "\n", sep="")
cat("EXPORTING TABLE 4: PATHOLOGY CHARACTERISTICS\n")
cat(rep("=", 80), "\n\n")

# Convert results list to dataframe
results_df <- bind_rows(lapply(all_results, as.data.frame))

write.csv(results_df, "pathology_table_results.csv", row.names = FALSE)
cat("Results saved to: pathology_table_results.csv\n\n")

cat("Creating Table 4 wide-format table...\n")

# Define outcome order and labels matching the markdown table
outcome_order <- c(
  "diagnostic_first_resection_binary",
  "margin_clear_first_resection_binary",
  "margin_dist_initial_mm_clean",
  "mass_initial_wedges_g",
  "overall_margin_status_binary",
  "specimen_mass",
  "is_metastatic",
  "is_benign"
)

outcome_labels_table <- c(

"diagnostic_first_resection_binary" = "First resection diagnostic (%)",
  "margin_clear_first_resection_binary" = "First resection margins clear (%)",
  "margin_dist_initial_mm_clean" = "Mean first resection margin distance, mm (SD)",
  "mass_initial_wedges_g" = "Mean first resection mass, g (SD)",
  "overall_margin_status_binary" = "Final margin clear (%)",
  "specimen_mass" = "Final mass tissue resected, g (SD)",
  "is_metastatic" = "Metastases resections (%)",
  "is_benign" = "Benign lobectomies (%)"
)

# Build wide table
wide_table <- data.frame(Characteristic = character(), stringsAsFactors = FALSE)

for (outcome in outcome_order) {
  is_binary <- outcome %in% binary_outcomes

  row_data <- list(Characteristic = outcome_labels_table[outcome])

  for (grp in c("RATS", "VATS", "Pooled")) {
    result <- all_results[[paste(grp, outcome, sep = "_")]]

    if (is_binary) {
      row_data[[paste0(grp, "_Marked")]] <- format_binary_stat(result$n_marked, result$pct_marked)
      row_data[[paste0(grp, "_Unmarked")]] <- format_binary_stat(result$n_unmarked, result$pct_unmarked)
      row_data[[paste0(grp, "_P")]] <- format_p_table(result$p_unadjusted)
      row_data[[paste0(grp, "_Adjusted_Effect")]] <- format_or_table(result$or_iptw, result$ci_low_iptw, result$ci_high_iptw)
      row_data[[paste0(grp, "_Adj_P")]] <- format_p_table(result$p_iptw)
    } else {
      row_data[[paste0(grp, "_Marked")]] <- format_continuous_stat(result$mean_marked, result$sd_marked)
      row_data[[paste0(grp, "_Unmarked")]] <- format_continuous_stat(result$mean_unmarked, result$sd_unmarked)
      row_data[[paste0(grp, "_P")]] <- format_p_table(result$p_unadjusted)
      row_data[[paste0(grp, "_Adjusted_Effect")]] <- format_beta_table(result$beta_iptw, result$ci_low_iptw, result$ci_high_iptw)
      row_data[[paste0(grp, "_Adj_P")]] <- format_p_table(result$p_iptw)
    }
  }

  wide_table <- bind_rows(wide_table, as.data.frame(row_data, stringsAsFactors = FALSE))
}

# Rename columns to match markdown table headers
colnames(wide_table) <- c(
  "Characteristic",
  "RATS Marked (N=38)", "RATS Unmarked (N=45)", "RATS P", "RATS Adjusted Effect (95% CI)", "RATS Adj. P",
  "VATS Marked (N=33)", "VATS Unmarked (N=75)", "VATS P", "VATS Adjusted Effect (95% CI)", "VATS Adj. P",
  "Pooled Marked (N=71)", "Pooled Unmarked (N=120)", "Pooled P", "Pooled Adjusted Effect (95% CI)", "Pooled Adj. P"
)

write.csv(wide_table, "pathology_table_wide.csv", row.names = FALSE)
cat("Wide-format table saved to: pathology_table_wide.csv\n\n")

# Also output as markdown
cat("Generating markdown table...\n")

md_header <- paste0(
  "| **Characteristic** | **RATS Marked (N=38)** | **RATS Unmarked (N=45)** | **RATS *P*** | **RATS Adjusted Effect (95% CI)** | **RATS Adj. *P*** | ",
  "**VATS Marked (N=33)** | **VATS Unmarked (N=75)** | **VATS *P*** | **VATS Adjusted Effect (95% CI)** | **VATS Adj. *P*** | ",
  "**Pooled Marked (N=71)** | **Pooled Unmarked (N=120)** | **Pooled *P*** | **Pooled Adjusted Effect (95% CI)** | **Pooled Adj. *P*** |"
)

md_separator <- "|-------------------|----------------------|------------------------|------------|------------------------------|-----------------|----------------------|------------------------|------------|------------------------------|-----------------|------------------------|---------------------------|--------------|-------------------------------------|-------------------|"

md_lines <- c("## Table 4. Pathology Characteristics", "", md_header, md_separator)

for (i in 1:nrow(wide_table)) {
  row <- wide_table[i, ]
  md_row <- paste0("| ", paste(row, collapse = " | "), " |")
  md_lines <- c(md_lines, md_row)
}

md_lines <- c(md_lines, "",
              "**Note:** Adjusted effects are from IPTW-weighted regression models. OR = Odds Ratio for binary outcomes; β = regression coefficient for continuous outcomes. Common support trimming applied.")

writeLines(md_lines, "pathology_table_wide_generated.md")
cat("Markdown table saved to: pathology_table_wide_generated.md\n\n")

# ==============================================================================
# 11. METHODOLOGY NOTES
# ==============================================================================

cat("\n\n", rep("=", 80), "\n", sep="")
cat("METHODOLOGY NOTES\n")
cat(rep("=", 80), "\n\n")

cat("IPTW vs Multivariable Regression:\n\n")

cat("IPTW (Inverse Probability of Treatment Weighting):\n")
cat("  - Creates pseudo-population where treatment is independent of covariates\n")
cat("  - Separates design (PS model) from outcome analysis\n")
cat("  - Preferred when: good PS discrimination, adequate overlap\n")
cat("  - Used as PRIMARY adjustment method in this analysis\n\n")

cat("Multivariable Regression:\n")
cat("  - Directly adjusts for covariates in outcome model\n")
cat("  - More efficient with small samples\n")
cat("  - Provided as SENSITIVITY analysis\n\n")

cat("When IPTW and MV estimates diverge substantially, this may indicate:\n")
cat("  - Model misspecification\n")
cat("  - Effect modification by covariates\n")
cat("  - Limited overlap in PS distributions\n\n")

cat("Covariates adjusted for:\n")
cat("  age, sex, nodule size, pleural distance, nodule consistency, smoking status\n\n")

cat(rep("=", 80), "\n", sep="")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n")
