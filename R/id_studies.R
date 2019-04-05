
### libraries
library(tidyverse)
library(Rcompadre)


### load compadre data
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X.RData")


## fix Lemke et al. treatment (should be Unmanipulated)
lemke_fix <- which(compadre$SpeciesAuthor == "Trollius_europaeus")
compadre$MatrixTreatment[lemke_fix] <- "Unmanipulated"


### function to collapse AdditionalSource column
collapse_fn <- function(x) {
  ifelse(all(is.na(x)),
         NA_character_,
         paste(unique(x[!is.na(x)]), collapse = "; "))
}


## subset COMPADRE to studies of interest
studies_check <- compadre %>%
  filter(YearPublication >= 2010) %>% 
  filter(MatrixSplit == "Divided",
         MatrixFec == "Yes",
         is.na(MatrixTreatment) | MatrixTreatment == "Unmanipulated",
         MatrixCaptivity == "W",
         AnnualPeriodicity == "1",
         OrganismType %in% c("Herbaceous perennial",
                             "Annual",
                             "Succulent",
                             "Shrub",
                             "Tree",
                             "Palm")) %>% 
  as_tibble() %>% 
  group_by(Authors, Journal, YearPublication, DOI.ISBN) %>% 
  summarize(AdditionalSource = collapse_fn(AdditionalSource),
            n_spp = length(unique(SpeciesAccepted)),
            n_pop = length(unique(paste(SpeciesAccepted, MatrixPopulation))),
            n_mat = n()) %>% 
  ungroup() %>% 
  arrange(YearPublication, Authors)


### write to file
write.csv("studies_check.csv", row.names = FALSE)



# ### median number of years of temporal replication
# comp_rep <- compadre %>% 
#   filter(MatrixTreatment == "Unmanipulated",
#          MatrixCaptivity == "W",
#          MatrixComposite == "Individual",
#          AnnualPeriodicity == "1") %>% 
#   as_tibble() %>% 
#   group_by(SpeciesAuthor, MatrixPopulation) %>% 
#   summarize(n = length(unique(paste(MatrixStartYear, MatrixEndYear)))) %>% 
#   ungroup() %>% 
#   filter(n > 1)
# 
# table(comp_rep$n)
# median(comp_rep$n)

