# S02: subset COMPADRE to target studies and save derived study list.

# libraries ----
source("code/setup.R")
setup_packages(c("tidyverse", "Rcompadre"))
source("code/functions.R")


# load compadre data ----
compadre <- load_compadre(corrected = TRUE)


# subset COMPADRE to studies of interest ----
studies_check <- compadre %>%
  filter(YearPublication >= 2010) %>%
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
  ) %>%
  group_by(Authors, Journal, YearPublication, DOI_ISBN) %>%
  summarize(
    AdditionalSource = collapse_fn(AdditionalSource),
    n_spp = length(unique(SpeciesAccepted)),
    n_pop = length(unique(paste(SpeciesAccepted, MatrixPopulation))),
    n_mat = n(),
    .groups = "drop"
  ) %>%
  arrange(YearPublication, Authors)

studies_sources <- compadre %>%
  filter(YearPublication >= 2010) %>%
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
    ),
    !is.na(AdditionalSource),
    AdditionalSource != ""
  ) %>%
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


# write to file ----
if (!dir.exists("data/derived/studies")) dir.create("data/derived/studies", recursive = TRUE)
write_csv(studies_check, "data/derived/studies/target_studies.csv")
write.csv(studies_check, "studies_check.csv", row.names = FALSE)

if (!dir.exists("docs/tables")) dir.create("docs/tables", recursive = TRUE)
write_csv(studies_sources, "docs/tables/studies_additional_sources.csv")
