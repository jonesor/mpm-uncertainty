# 99: create publication-ready tables for manuscript and appendices.

# libraries ----
source("code/setup.R")
setup_packages(c("dplyr", "readr", "stringr", "tibble", "officer", "flextable"))


# helpers ----
table_dir <- "docs/tables"
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

write_table_docx <- function(df, title, path) {
  ft <- flextable::flextable(df)
  ft <- flextable::autofit(ft)
  ft <- flextable::align(ft, align = "left", part = "all")

  if ("Parameter" %in% names(df)) {
    vals <- df$Parameter

    idx <- which(vals == "Pop.~growth~rate~(italic(lambda))")
    if (length(idx) > 0) {
      ft <- flextable::compose(
        ft,
        i = idx,
        j = "Parameter",
        value = flextable::as_paragraph(
          "Pop. growth rate (", "λ", ")"
        )
      )
    }

    idx <- which(vals == "Damping~ratio~(italic(rho))")
    if (length(idx) > 0) {
      ft <- flextable::compose(
        ft,
        i = idx,
        j = "Parameter",
        value = flextable::as_paragraph(
          "Damping ratio (", "ρ", ")"
        )
      )
    }

    idx <- which(vals == "Repro.~value~(Veg.)")
    if (length(idx) > 0) {
      ft <- flextable::compose(
        ft,
        i = idx,
        j = "Parameter",
        value = flextable::as_paragraph("Reproductive value (vegetative stage)")
      )
    }

    idx <- which(vals == "Generation~time~(italic(T))")
    if (length(idx) > 0) {
      ft <- flextable::compose(
        ft,
        i = idx,
        j = "Parameter",
        value = flextable::as_paragraph(
          "Generation time (", flextable::as_i("T"), ")"
        )
      )
    }

    idx <- which(vals == "Life~expectancy~(italic(l[0]))")
    if (length(idx) > 0) {
      ft <- flextable::compose(
        ft,
        i = idx,
        j = "Parameter",
        value = flextable::as_paragraph(
          "Life expectancy (", flextable::as_i("l"), flextable::as_sub("0"), ")"
        )
      )
    }
  }

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "Normal")
  doc <- flextable::body_add_flextable(doc, value = ft)
  print(doc, target = path)
}

clean_par_label <- function(x) {
  x %>%
    stringr::str_replace_all("~", " ") %>%
    stringr::str_replace_all("italic\\(([^)]+)\\)", "\\1") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()
}

fmt_mean_sd <- function(x, s, digits = 2) {
  paste0(formatC(x, digits = digits, format = "f"), " ± ",
         formatC(s, digits = digits, format = "f"))
}

fmt_n_prop <- function(n, p, digits = 1) {
  paste0(n, " (", formatC(100 * p, digits = digits, format = "f"), "%)")
}


# Table 1: variance components ----
var_path <- "data/derived/analysis_cache/case1_variance_ratios.csv"
if (file.exists(var_path)) {
  table1 <- readr::read_csv(var_path, show_col_types = FALSE) %>%
    transmute(
      Parameter = parameter,
      `Sampling:other variance ratio` = signif(ratio, 3),
      `95% interval` = sprintf("[%.2f, %.2f]", low95, upp95)
    )

  write_table_docx(
    table1,
    "Table 1. Sampling-to-other variance ratios for parameters derived from MPMs, with 95% intervals from the variance-component model.",
    file.path(table_dir, "Table1_variance_components.docx")
  )
}


# Table S1-2: single MPM derived parameters ----
s1_2_path <- "data/derived/analysis_cache/fig1_derived_param_summary.csv"
if (file.exists(s1_2_path)) {
  table_s1_2 <- readr::read_csv(s1_2_path, show_col_types = FALSE) %>%
    mutate(par = clean_par_label(par)) %>%
    transmute(
      Parameter = par,
      Median = med,
      `2.5%` = low,
      `97.5%` = upp,
      `Point estimate` = value
    ) %>%
    mutate(across(where(is.numeric), ~ signif(.x, 3)))

  write_table_docx(
    table_s1_2,
    "Table S1-2. Point estimates and draws from the sampling distribution for derived parameters from a single MPM.",
    file.path(table_dir, "Table_S1-2_single_mpm_derived.docx")
  )
}


# Table S1-3 fallback ----
s1_3_path <- "data/derived/analysis_cache/fig1_mean_mpm_param_summary.csv"
if (file.exists(s1_3_path)) {
  table_s1_3 <- readr::read_csv(s1_3_path, show_col_types = FALSE) %>%
    mutate(par = clean_par_label(par)) %>%
    transmute(
      Parameter = par,
      Median = med,
      `2.5%` = low,
      `97.5%` = upp,
      `Point estimate` = value
    ) %>%
    mutate(across(where(is.numeric), ~ signif(.x, 3)))

  write_table_docx(
    table_s1_3,
    "Table S1-3. Point estimates and draws from the sampling distribution for derived parameters from a mean MPM.",
    file.path(table_dir, "Table_S1-3_mean_mpm_derived.docx")
  )
} else {
  fallback_note <- tibble::tibble(
    Note = "Source data unavailable: data/derived/analysis_cache/fig1_mean_mpm_param_summary.csv"
  )
  write_table_docx(
    fallback_note,
    "Table S1-3. Point estimates and draws from the sampling distribution for derived parameters from a mean MPM.",
    file.path(table_dir, "Table_S1-3_mean_mpm_derived.docx")
  )
}


# Table S1-7: included versus excluded screened study comparison ----
study_sel_num_path <- "data/derived/analysis_cache/study_selection_comparison_summary_numeric.csv"
study_sel_life_path <- "data/derived/analysis_cache/study_selection_comparison_summary_life_form.csv"
if (file.exists(study_sel_num_path) && file.exists(study_sel_life_path)) {
  sel_num <- readr::read_csv(study_sel_num_path, show_col_types = FALSE)
  sel_life <- readr::read_csv(study_sel_life_path, show_col_types = FALSE)

  included_num <- sel_num %>% filter(selection_group == "Included")
  excluded_num <- sel_num %>% filter(selection_group == "Excluded")

  life_levels <- c(
    "Herbaceous perennial",
    "Succulent",
    "Shrub",
    "Tree",
    "Palm"
  )

  life_tbl <- tibble::tibble(dominant_life_form = life_levels) %>%
    left_join(
      sel_life %>%
        filter(selection_group == "Included") %>%
        select(dominant_life_form, n_studies, prop_studies) %>%
        rename(
          included_n = n_studies,
          included_prop = prop_studies
        ),
      by = "dominant_life_form"
    ) %>%
    left_join(
      sel_life %>%
        filter(selection_group == "Excluded") %>%
        select(dominant_life_form, n_studies, prop_studies) %>%
        rename(
          excluded_n = n_studies,
          excluded_prop = prop_studies
        ),
      by = "dominant_life_form"
    ) %>%
    mutate(
      included_n = coalesce(included_n, 0),
      included_prop = coalesce(included_prop, 0),
      excluded_n = coalesce(excluded_n, 0),
      excluded_prop = coalesce(excluded_prop, 0)
    )

  table_s1_7 <- tibble::tibble(
    Characteristic = c(
      "Studies, n",
      "Mean matrix dimension",
      "Mean stage-specific survival",
      "Proportion of survival estimates equal to 0",
      "Proportion of survival estimates equal to 1",
      "Proportion of survival estimates at the boundary (0 or 1)",
      paste("Dominant life form:", life_tbl$dominant_life_form)
    ),
    Included = c(
      as.character(included_num$n_studies),
      fmt_mean_sd(
        included_num$mean_matrix_dimension_mean,
        included_num$mean_matrix_dimension_sd
      ),
      fmt_mean_sd(
        included_num$mean_stage_survival_mean,
        included_num$mean_stage_survival_sd
      ),
      fmt_mean_sd(
        included_num$prop_boundary_zero_mean,
        included_num$prop_boundary_zero_sd
      ),
      fmt_mean_sd(
        included_num$prop_boundary_one_mean,
        included_num$prop_boundary_one_sd
      ),
      fmt_mean_sd(
        included_num$prop_boundary_any_mean,
        included_num$prop_boundary_any_sd
      ),
      mapply(
        fmt_n_prop,
        life_tbl$included_n,
        life_tbl$included_prop
      )
    ),
    Excluded = c(
      as.character(excluded_num$n_studies),
      fmt_mean_sd(
        excluded_num$mean_matrix_dimension_mean,
        excluded_num$mean_matrix_dimension_sd
      ),
      fmt_mean_sd(
        excluded_num$mean_stage_survival_mean,
        excluded_num$mean_stage_survival_sd
      ),
      fmt_mean_sd(
        excluded_num$prop_boundary_zero_mean,
        excluded_num$prop_boundary_zero_sd
      ),
      fmt_mean_sd(
        excluded_num$prop_boundary_one_mean,
        excluded_num$prop_boundary_one_sd
      ),
      fmt_mean_sd(
        excluded_num$prop_boundary_any_mean,
        excluded_num$prop_boundary_any_sd
      ),
      mapply(
        fmt_n_prop,
        life_tbl$excluded_n,
        life_tbl$excluded_prop
      )
    )
  )

  write_table_docx(
    table_s1_7,
    paste(
      "Table S1-7. Comparison of screened target studies that were included",
      "versus excluded under the current reproducible COMPADRE-based screen.",
      "Numeric entries are study-level means ± SD; life-form entries are",
      "numbers of studies with percentages in parentheses."
    ),
    file.path(table_dir, "Table_S1-7_included_vs_excluded_screened_studies.docx")
  )
}


# Table S12: boundary estimates and small sample sizes ----
boundary_group_path <- "data/derived/analysis_cache/boundary_smallN_group_summary.csv"
if (file.exists(boundary_group_path)) {
  boundary_groups <- readr::read_csv(boundary_group_path, show_col_types = FALSE) %>%
    mutate(
      group_type = factor(group_type, levels = c("Life form", "Matrix dimension")),
      group_value = case_when(
        group_type == "Matrix dimension" & group_value == "3-5" ~ "3-5 stages",
        group_type == "Matrix dimension" & group_value == "6+" ~ "6+ stages",
        TRUE ~ group_value
      )
    ) %>%
    arrange(group_type, group_value)

  table_s12 <- boundary_groups %>%
    transmute(
      Group = group_type,
      Category = group_value,
      `Population entries, n` = n_populations,
      `Any boundary survival estimate (%)` = round(pct_any_boundary, 1),
      `Any stage with N < 20 (%)` = round(pct_any_small_n20, 1),
      `Any stage with N < 50 (%)` = round(pct_any_small_n50, 1),
      `Mean minimum stage N` = round(mean_min_stage_n, 1)
    )

  write_table_docx(
    table_s12,
    paste(
      "Table S12. Distribution of boundary survival estimates and small",
      "stage-specific sample sizes across analysis 1 population-level MPM",
      "entries, grouped by the COMPADRE life-form categories available for",
      "these studies and by broad matrix-dimension class. Percentages give",
      "the share of population entries with at least one boundary survival",
      "estimate or at least one stage with sample size below the stated",
      "threshold."
    ),
    file.path(table_dir, "Table_S12_boundary_smallN_patterns.docx")
  )
}


# Table S13: displacement models ----
disp_model_path <- "data/derived/analysis_cache/case1_displacement_boundary_models.csv"
if (file.exists(disp_model_path)) {
  table_s13 <- readr::read_csv(disp_model_path, show_col_types = FALSE) %>%
    mutate(
      Response = dplyr::case_when(
        response == "abs_disp_S" ~ "Shape displacement |point - median|",
        response == "abs_disp_log_L" ~ "Life expectancy displacement |log10(point) - log10(median)|",
        TRUE ~ response
      ),
      Predictor = dplyr::case_when(
        predictor == "n_boundary_stages" ~ "Number of boundary survival stages",
        predictor == "prop_individuals_boundary" ~ "Proportion of individuals in boundary-survival stages",
        TRUE ~ predictor
      ),
      Slope = sprintf("%.3f", slope),
      `95% CI` = sprintf("[%.3f, %.3f]", conf_low, conf_upp),
      `P value` = format.pval(p_value, digits = 2, eps = 1e-4),
      `R²` = sprintf("%.3f", r_squared)
    ) %>%
    transmute(
      Response,
      Predictor,
      `Population entries, n` = n,
      Slope,
      `95% CI`,
      `P value`,
      `R²`
    )

  write_table_docx(
    table_s13,
    paste(
      "Table S13. Linear models relating displacement magnitude in analysis 1",
      "to boundary-survival predictors. Displacement is defined as the absolute",
      "difference between the point estimate and the median of the sampling",
      "distribution; life expectancy displacement is modelled on the log10 scale.",
      "Positive slopes indicate greater point-estimate displacement in matrices",
      "with more boundary survival structure."
    ),
    file.path(table_dir, "Table_S13_displacement_boundary_models.docx")
  )
}


# Table S14: quasi-stable threshold sensitivity ----
qsd_sens_path <- "data/derived/analysis_cache/case1_qsd_threshold_sensitivity_summary.csv"
if (file.exists(qsd_sens_path)) {
  table_s14 <- readr::read_csv(qsd_sens_path, show_col_types = FALSE) %>%
    transmute(
      Comparison = comparison,
      `Population entries, n` = n,
      `Correlation of point estimates` = sprintf("%.3f", cor_shape_pt),
      `Median |ΔS| (point estimate)` = sprintf("%.3f", median_abs_diff_pt),
      `95th percentile |ΔS| (point estimate)` = sprintf("%.3f", p95_abs_diff_pt),
      `Max |ΔS| (point estimate)` = sprintf("%.3f", max_abs_diff_pt),
      `Median |ΔS| (sampling median)` = sprintf("%.3f", median_abs_diff_med),
      `Sign changes in S` = n_sign_change_pt,
      `Median Δq` = median_delta_q,
      `Range of Δq` = range_delta_q
    )

  write_table_docx(
    table_s14,
    paste(
      "Table S14. Sensitivity of analysis 1 shape estimates to the",
      "quasi-stable-distribution truncation threshold used when calculating",
      "the shape metric. Comparisons are against the 1% threshold used in",
      "the main analysis. Small differences and near-unit correlations show",
      "that the substantive shape results are stable to plausible threshold",
      "choices."
    ),
    file.path(table_dir, "Table_S14_qsd_threshold_sensitivity.docx")
  )
}
