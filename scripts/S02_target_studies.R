# S02: subset COMPADRE to target studies and save derived study list.

# libraries ----
source("code/setup.R")
setup_packages(c("tidyverse", "Rcompadre"))
source("code/functions.R")


# load compadre data ----
compadre <- cdb_fetch("data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData")


# subset COMPADRE to studies of interest ----
studies_check <- compadre %>%
  filter(YearPublication >= 2010) %>% 
  filter(MatrixSplit == "Divided",
         MatrixFec == "Yes",
         is.na(MatrixTreatment) | MatrixTreatment == "Unmanipulated",
         MatrixDimension > 2,
         MatrixCaptivity == "W",
         ProjectionInterval == "1",
         OrganismType %in% c("Herbaceous perennial",
                             "Succulent",
                             "Shrub",
                             "Tree",
                             "Palm")) %>% 
  group_by(Authors, Journal, YearPublication, DOI_ISBN) %>% 
  summarize(AdditionalSource = collapse_fn(AdditionalSource),
            n_spp = length(unique(SpeciesAccepted)),
            n_pop = length(unique(paste(SpeciesAccepted, MatrixPopulation))),
            n_mat = n(),
            .groups = "drop") %>% 
  arrange(YearPublication, Authors)


# write to file ----
if (!dir.exists("data/derived/studies")) dir.create("data/derived/studies", recursive = TRUE)
write_csv(studies_check, "data/derived/studies/target_studies.csv")
write.csv(studies_check, "studies_check.csv", row.names = FALSE)

