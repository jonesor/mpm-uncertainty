
### libraries
library(tidyverse)
library(Rcompadre)
library(popdemo)


### load COMPADRE
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


## possible columns to collapse on
col_spp <- c("id_stage", "SpeciesAuthor", "AnnualPeriodicity")
col_pop <- c(col_spp, "MatrixPopulation", "MatrixTreatment")


### subset
compadre_sub <- compadre %>% 
  cdb_flag() %>% 
  filter(check_NA_U == FALSE, check_zero_U == FALSE, check_NA_A == FALSE) %>% 
  filter(MatrixComposite != "Seasonal") %>% 
  mutate(id_stage = cdb_id_stages(.)) %>% 
  mutate(id1 = cdb_id(., col_spp)) %>% 
  mutate(id2 = cdb_id(., col_pop))


# full db
comp_tot <- compadre_sub %>% 
  cdb_unnest() %>% 
  mutate(ergodic = map_lgl(matA, isErgodic)) %>% 
  mutate(surv = map(matU, colSums)) %>% 
  mutate(ns1 = map_int(surv, ~ length(.x[.x >= 1])))

# collpase by population/treatment
comp_pop <- compadre_sub %>% 
  cdb_collapse("id2") %>% 
  cdb_unnest() %>% 
  mutate(ergodic = map_lgl(matA, isErgodic)) %>% 
  mutate(surv = map(matU, colSums)) %>% 
  mutate(ns1 = map_int(surv, ~ length(.x[.x >= 1])))

# collapse by species
comp_spp <- compadre_sub %>% 
  cdb_collapse("id1") %>% 
  cdb_unnest() %>% 
  mutate(ergodic = map_lgl(matA, isErgodic)) %>% 
  mutate(surv = map(matU, colSums)) %>% 
  mutate(ns1 = map_int(surv, ~ length(.x[.x >= 1])))

# number with x stage-specific survival >= 1
table(comp_tot$ns1)
table(comp_pop$ns1)
table(comp_spp$ns1)

# proportion MPMs with 1+ stage-specific survival >= 1
sum(table(comp_tot$ns1)[-1]) / sum(table(comp_tot$ns1))
sum(table(comp_pop$ns1)[-1]) / sum(table(comp_pop$ns1))
sum(table(comp_spp$ns1)[-1]) / sum(table(comp_spp$ns1))

# proportion MPMs with 2+ stage-specific survival >= 1
sum(table(comp_tot$ns1)[-(1:2)]) / sum(table(comp_tot$ns1))
sum(table(comp_pop$ns1)[-(1:2)]) / sum(table(comp_pop$ns1))
sum(table(comp_spp$ns1)[-(1:2)]) / sum(table(comp_spp$ns1))


