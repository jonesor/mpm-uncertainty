# Check concordance between project target studies and latest public COMPADRE release ----

source(here::here("code", "setup.R"))

setup_packages(c("Rcompadre", "dplyr", "readr", "stringr", "tibble"))

# Normalize free-text keys so metadata differences (case/punctuation) do not block matching.
norm_key <- function(x) {
  x |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]", "")
}

first_author <- function(x) {
  stringr::str_split(x, ";|,| and ", simplify = TRUE)[, 1] |>
    stringr::str_trim()
}

compadre_latest <- Rcompadre::cdb_fetch("compadre") |>
  as_tibble() |>
  dplyr::mutate(
    doi_norm = norm_key(dplyr::coalesce(DOI_ISBN, "")),
    authors_norm = norm_key(dplyr::coalesce(Authors, "")),
    first_author_norm = norm_key(first_author(dplyr::coalesce(Authors, ""))),
    journal_norm = norm_key(dplyr::coalesce(Journal, "")),
    year_int = suppressWarnings(as.integer(YearPublication))
  )

target <- readr::read_csv(
  here::here("data", "derived", "studies", "target_studies.csv"),
  show_col_types = FALSE
) |>
  dplyr::mutate(
    doi_norm = norm_key(dplyr::coalesce(DOI_ISBN, "")),
    authors_norm = norm_key(dplyr::coalesce(Authors, "")),
    first_author_norm = norm_key(first_author(dplyr::coalesce(Authors, ""))),
    journal_norm = norm_key(dplyr::coalesce(Journal, "")),
    year_int = suppressWarnings(as.integer(YearPublication))
  )

key_doi <- compadre_latest |>
  dplyr::filter(doi_norm != "") |>
  dplyr::distinct(doi_norm) |>
  dplyr::pull(doi_norm)

key_auth_jour_year <- compadre_latest |>
  dplyr::distinct(authors_norm, journal_norm, year_int) |>
  dplyr::mutate(key = paste(authors_norm, journal_norm, year_int, sep = "|")) |>
  dplyr::pull(key)

key_firstauth_year <- compadre_latest |>
  dplyr::distinct(first_author_norm, year_int) |>
  dplyr::mutate(key = paste(first_author_norm, year_int, sep = "|")) |>
  dplyr::pull(key)

concordance <- target |>
  dplyr::mutate(
    match_doi = doi_norm != "" & doi_norm %in% key_doi,
    match_author_journal_year = paste(authors_norm, journal_norm, year_int, sep = "|") %in% key_auth_jour_year,
    match_first_author_year = paste(first_author_norm, year_int, sep = "|") %in% key_firstauth_year,
    matched_latest = match_doi | match_author_journal_year | match_first_author_year
  )

# Manual linkage for known study-level records where additional source material
# is attached to a DOI-linked study entry in the analysis dataset.
manual_support_link_doi <- c(
  "101890es13000941"
)

concordance <- concordance |>
  dplyr::mutate(
    match_manual_support_link = doi_norm %in% manual_support_link_doi &
      !is.na(AdditionalSource) &
      AdditionalSource != "",
    matched_concordance = matched_latest | match_manual_support_link
  ) |>
  dplyr::select(
    Authors,
    Journal,
    YearPublication,
    DOI_ISBN,
    AdditionalSource,
    match_doi,
    match_author_journal_year,
    match_first_author_year,
    matched_latest,
    match_manual_support_link,
    matched_concordance
  )

readr::write_csv(
  concordance,
  here::here("data", "derived", "studies", "compadre_latest_concordance.csv")
)

summary_all <- concordance |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_matched = sum(matched_concordance, na.rm = TRUE),
    n_unmatched = n_total - n_matched
  )

summary_with_sources <- concordance |>
  dplyr::filter(!is.na(AdditionalSource), AdditionalSource != "") |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_matched = sum(matched_concordance, na.rm = TRUE),
    n_unmatched = n_total - n_matched
  )

message("Wrote: data/derived/studies/compadre_latest_concordance.csv")
message(
  "Target studies matched under concordance check: ",
  summary_all$n_matched,
  "/",
  summary_all$n_total,
  " (unmatched ",
  summary_all$n_unmatched,
  ")"
)
message(
  "Target studies with AdditionalSource matched under concordance check: ",
  summary_with_sources$n_matched,
  "/",
  summary_with_sources$n_total,
  " (unmatched ",
  summary_with_sources$n_unmatched,
  ")"
)
