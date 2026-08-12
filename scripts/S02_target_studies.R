# S02: identify target studies and compare included versus excluded study sets.


# libraries ----
source("code/setup.R")
setup_packages(c("tidyverse", "Rcompadre"))
source("code/functions.R")


# helpers ----
norm_key <- function(x) {
  x %>%
    coalesce("") %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]", "")
}

first_author <- function(x) {
  str_split(coalesce(x, ""), ";|,| and ", simplify = TRUE)[, 1] %>%
    str_trim()
}


# load compadre data ----
compadre <- load_compadre(corrected = TRUE)


# subset COMPADRE to screened target studies ----
target_compadre <- compadre %>%
  mutate(
    YearPublication = suppressWarnings(as.integer(YearPublication)),
    MatrixDimension = suppressWarnings(as.integer(MatrixDimension))
  ) %>%
  filter(!is.na(YearPublication), YearPublication >= 2010) %>%
  filter(
    MatrixSplit == "Divided",
    MatrixFec == "Yes",
    is.na(MatrixTreatment) | MatrixTreatment == "Unmanipulated",
    MatrixDimension > 2,
    MatrixCaptivity == "W",
    ProjectionInterval == "1",
    OrganismType %in% c(
      "Herbaceous perennial",
      "Succulent",
      "Shrub",
      "Tree",
      "Palm"
    )
  )

studies_check <- target_compadre %>%
  group_by(Authors, Journal, YearPublication, DOI_ISBN) %>%
  summarize(
    AdditionalSource = collapse_fn(AdditionalSource),
    n_spp = length(unique(SpeciesAccepted)),
    n_pop = length(unique(paste(SpeciesAccepted, MatrixPopulation))),
    n_mat = n(),
    .groups = "drop"
  ) %>%
  arrange(YearPublication, Authors)

studies_sources <- target_compadre %>%
  filter(!is.na(AdditionalSource), AdditionalSource != "") %>%
  group_by(Authors, YearPublication, Journal, DOI_ISBN) %>%
  summarize(
    SpeciesAuthor = paste(sort(unique(SpeciesAuthor)), collapse = "; "),
    AdditionalSource = collapse_fn(AdditionalSource),
    .groups = "drop"
  ) %>%
  mutate(
    DOI_ISBN = case_when(
      Authors == "Villellas; Morris; Garcia" & is.na(DOI_ISBN) ~ "10.1111/j.1600-0587.2012.07425.x",
      Authors == "Navarro; Galeano; Bernal" & is.na(DOI_ISBN) ~ "10.1177/19400829110040010",
      TRUE ~ DOI_ISBN
    )
  ) %>%
  arrange(YearPublication, Authors)


# compare included and excluded screened studies ----
included_studies <- read_csv(
  "data/derived/studies/_data_sources.csv",
  col_names = c(
    "Authors",
    "YearPublication",
    "Journal",
    "DOI_ISBN",
    "SpeciesAccepted"
  ),
  trim_ws = TRUE,
  show_col_types = FALSE
) %>%
  mutate(
    YearPublication = suppressWarnings(as.integer(YearPublication)),
    doi_norm = norm_key(DOI_ISBN),
    authors_norm = norm_key(Authors),
    journal_norm = norm_key(Journal),
    first_author_norm = norm_key(first_author(Authors))
  )

included_studies_2010plus <- included_studies %>%
  filter(!is.na(YearPublication), YearPublication >= 2010)

target_studies_current <- studies_check %>%
  mutate(
    doi_norm = norm_key(DOI_ISBN),
    authors_norm = norm_key(Authors),
    journal_norm = norm_key(Journal),
    first_author_norm = norm_key(first_author(Authors)),
    match_doi = doi_norm != "" & doi_norm %in% included_studies_2010plus$doi_norm,
    match_author_journal_year = paste(
      authors_norm, journal_norm, YearPublication, sep = "|"
    ) %in% paste(
      included_studies_2010plus$authors_norm,
      included_studies_2010plus$journal_norm,
      included_studies_2010plus$YearPublication,
      sep = "|"
    ),
    match_first_author_year = paste(
      first_author_norm, YearPublication, sep = "|"
    ) %in% paste(
      included_studies_2010plus$first_author_norm,
      included_studies_2010plus$YearPublication,
      sep = "|"
    ),
    included_current_screen = match_doi |
      match_author_journal_year |
      match_first_author_year
  )

study_level_stats <- target_compadre %>%
  as_tibble() %>%
  transmute(
    Authors,
    Journal,
    YearPublication,
    DOI_ISBN,
    OrganismType,
    MatrixDimension,
    surv = map(mat, ~ colSums(matU(.x)))
  ) %>%
  mutate(
    surv_mean = map_dbl(surv, ~ mean(.x, na.rm = TRUE)),
    surv_boundary_zero = map_dbl(surv, ~ mean(.x == 0, na.rm = TRUE)),
    surv_boundary_one = map_dbl(surv, ~ mean(.x == 1, na.rm = TRUE)),
    surv_boundary_any = map_dbl(surv, ~ mean(.x %in% c(0, 1), na.rm = TRUE))
  ) %>%
  select(-surv) %>%
  group_by(Authors, Journal, YearPublication, DOI_ISBN) %>%
  summarize(
    dominant_life_form = names(sort(table(OrganismType), decreasing = TRUE))[1],
    mean_matrix_dimension = mean(MatrixDimension, na.rm = TRUE),
    mean_stage_survival = mean(surv_mean, na.rm = TRUE),
    prop_boundary_zero = mean(surv_boundary_zero, na.rm = TRUE),
    prop_boundary_one = mean(surv_boundary_one, na.rm = TRUE),
    prop_boundary_any = mean(surv_boundary_any, na.rm = TRUE),
    .groups = "drop"
  )

study_selection_comparison <- target_studies_current %>%
  left_join(
    study_level_stats,
    by = c("Authors", "Journal", "YearPublication", "DOI_ISBN")
  ) %>%
  mutate(selection_group = if_else(included_current_screen, "Included", "Excluded"))

selection_summary_numeric <- study_selection_comparison %>%
  group_by(selection_group) %>%
  summarize(
    n_studies = n(),
    mean_matrix_dimension_mean = mean(mean_matrix_dimension, na.rm = TRUE),
    mean_matrix_dimension_sd = sd(mean_matrix_dimension, na.rm = TRUE),
    mean_stage_survival_mean = mean(mean_stage_survival, na.rm = TRUE),
    mean_stage_survival_sd = sd(mean_stage_survival, na.rm = TRUE),
    prop_boundary_zero_mean = mean(prop_boundary_zero, na.rm = TRUE),
    prop_boundary_zero_sd = sd(prop_boundary_zero, na.rm = TRUE),
    prop_boundary_one_mean = mean(prop_boundary_one, na.rm = TRUE),
    prop_boundary_one_sd = sd(prop_boundary_one, na.rm = TRUE),
    prop_boundary_any_mean = mean(prop_boundary_any, na.rm = TRUE),
    prop_boundary_any_sd = sd(prop_boundary_any, na.rm = TRUE),
    .groups = "drop"
  )

selection_summary_life_form <- study_selection_comparison %>%
  count(selection_group, dominant_life_form, name = "n_studies") %>%
  group_by(selection_group) %>%
  mutate(prop_studies = n_studies / sum(n_studies)) %>%
  ungroup()

included_not_in_current_screen <- included_studies_2010plus %>%
  mutate(
    matched_current_screen = doi_norm %in% target_studies_current$doi_norm |
      paste(authors_norm, journal_norm, YearPublication, sep = "|") %in%
      paste(
        target_studies_current$authors_norm,
        target_studies_current$journal_norm,
        target_studies_current$YearPublication,
        sep = "|"
      ) |
      paste(first_author_norm, YearPublication, sep = "|") %in%
      paste(
        target_studies_current$first_author_norm,
        target_studies_current$YearPublication,
        sep = "|"
      )
  ) %>%
  filter(!matched_current_screen) %>%
  select(Authors, YearPublication, Journal, DOI_ISBN, SpeciesAccepted)


# compare final 31-study dataset to the full COMPADRE plant database ----
included_exact_keys <- included_studies %>%
  mutate(
    exact_author_journal_year = paste(
      authors_norm, journal_norm, YearPublication, sep = "|"
    )
  ) %>%
  select(doi_norm, exact_author_journal_year) %>%
  distinct()

included_exact_compadre <- compadre %>%
  mutate(
    YearPublication = suppressWarnings(as.integer(YearPublication)),
    doi_norm = norm_key(DOI_ISBN),
    authors_norm = norm_key(Authors),
    journal_norm = norm_key(Journal),
    exact_author_journal_year = paste(
      authors_norm, journal_norm, YearPublication, sep = "|"
    )
  ) %>%
  filter(
    (doi_norm != "" & doi_norm %in% included_exact_keys$doi_norm) |
      exact_author_journal_year %in% included_exact_keys$exact_author_journal_year
  )

summarise_life_form_distribution <- function(df, unit = c("study", "species")) {
  unit <- match.arg(unit)
  out <- df

  if (unit == "study") {
    out <- out %>%
      group_by(Authors, Journal, YearPublication, DOI_ISBN) %>%
      summarize(
        OrganismType = names(sort(table(OrganismType), decreasing = TRUE))[1],
        .groups = "drop"
      )
  } else {
    out <- out %>%
      group_by(SpeciesAuthor) %>%
      summarize(
        OrganismType = names(sort(table(OrganismType), decreasing = TRUE))[1],
        .groups = "drop"
      )
  }

  out %>%
    count(OrganismType, name = "n") %>%
    mutate(prop = n / sum(n)) %>%
    arrange(desc(n), OrganismType)
}

life_form_comp_study <- full_join(
  summarise_life_form_distribution(included_exact_compadre, "study") %>%
    rename(
      included_n = n,
      included_prop = prop
    ),
  summarise_life_form_distribution(compadre %>% as_tibble(), "study") %>%
    rename(
      full_n = n,
      full_prop = prop
    ),
  by = "OrganismType"
) %>%
  mutate(
    included_n = coalesce(included_n, 0L),
    included_prop = coalesce(included_prop, 0),
    full_n = coalesce(full_n, 0L),
    full_prop = coalesce(full_prop, 0),
    prop_diff = included_prop - full_prop
  ) %>%
  arrange(desc(included_prop))

life_form_comp_species <- full_join(
  summarise_life_form_distribution(included_exact_compadre, "species") %>%
    rename(
      included_n = n,
      included_prop = prop
    ),
  summarise_life_form_distribution(compadre %>% as_tibble(), "species") %>%
    rename(
      full_n = n,
      full_prop = prop
    ),
  by = "OrganismType"
) %>%
  mutate(
    included_n = coalesce(included_n, 0L),
    included_prop = coalesce(included_prop, 0),
    full_n = coalesce(full_n, 0L),
    full_prop = coalesce(full_prop, 0),
    prop_diff = included_prop - full_prop
  ) %>%
  arrange(desc(included_prop))


# write to file ----
if (!dir.exists("data/derived/studies")) dir.create("data/derived/studies", recursive = TRUE)
if (!dir.exists("data/derived/analysis_cache")) dir.create("data/derived/analysis_cache", recursive = TRUE)
if (!dir.exists("docs/tables")) dir.create("docs/tables", recursive = TRUE)

write_csv(studies_check, "data/derived/studies/target_studies.csv")
write.csv(studies_check, "studies_check.csv", row.names = FALSE)
write_csv(studies_sources, "docs/tables/studies_additional_sources.csv")

write_csv(
  study_selection_comparison,
  "data/derived/analysis_cache/study_selection_comparison.csv"
)
write_csv(
  selection_summary_numeric,
  "data/derived/analysis_cache/study_selection_comparison_summary_numeric.csv"
)
write_csv(
  selection_summary_life_form,
  "data/derived/analysis_cache/study_selection_comparison_summary_life_form.csv"
)
write_csv(
  included_not_in_current_screen,
  "data/derived/analysis_cache/study_selection_included_not_in_current_screen.csv"
)
write_csv(
  life_form_comp_study,
  "data/derived/analysis_cache/final_dataset_vs_compadre_life_form_study.csv"
)
write_csv(
  life_form_comp_species,
  "data/derived/analysis_cache/final_dataset_vs_compadre_life_form_species.csv"
)
