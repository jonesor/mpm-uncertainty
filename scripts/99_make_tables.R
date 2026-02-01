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


# Table 1: variance components ----
var_path <- "data/derived/analysis_cache/case1_variance_ratios.csv"
if (file.exists(var_path)) {
  table1 <- readr::read_csv(var_path, show_col_types = FALSE) %>%
    transmute(
      Parameter = parameter,
      `Variance ratio (sampling / point)` = signif(ratio, 3)
    )

  write_table_docx(
    table1,
    "Table 1. Variance components in parameters derived from MPMs.",
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


# Table S1-3 placeholder ----
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
  placeholder <- tibble::tibble(
    Note = "Source data not yet generated: data/derived/analysis_cache/fig1_mean_mpm_param_summary.csv"
  )
  write_table_docx(
    placeholder,
    "Table S1-3. Point estimates and draws from the sampling distribution for derived parameters from a mean MPM.",
    file.path(table_dir, "Table_S1-3_mean_mpm_derived.docx")
  )
}
