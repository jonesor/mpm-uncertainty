
### libraries
library(tidyverse)
library(Rcompadre)
library(Rage)
library(popbio)
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




## number of reproductive stages
out <- comp_spp %>% 
  mutate(nrep = map_int(matF, ~ length(which(colSums(.x) > 0)))) %>% 
  filter(nrep > 0)

sum(out$nrep <= 3) / nrow(out)




out_shape <- comp_spp %>% 
  filter(check_NA_F == FALSE) %>% 
  filter(SurvivalIssue < 1.01) %>% 
  mutate(has_active = mpm_has_active(.)) %>% 
  filter(has_active == TRUE) %>% 
  mutate(rep_stages = map(matF, ~ colSums(.x) > 0)) %>% 
  mutate(start = mpm_first_active(.)) %>% 
  mutate(perennial = map2_lgl(matU, start, ~ mpm_to_lx(.x, .y, N = 3)[4] > 0)) %>% 
  mutate(any_rep = map_lgl(matF, ~ any(.x > 0))) %>% 
  filter(perennial == TRUE, any_rep == TRUE) %>% 
  mutate(rep_prop1 = pmap(list(matU, start, rep_stages), repro_prop_start)) %>% 
  mutate(rep_na = map_lgl(rep_prop1, ~ any(is.nan(.x)))) %>% 
  filter(rep_na == FALSE) %>% 
  mutate(lx = map2(matU, rep_prop1, lx_from_mature)) %>%
  mutate(lx_n = map_int(lx, length)) %>% 
  mutate(q = map2_int(matU, rep_prop1, qsd, nmax = 1e5)) %>%
  filter(!is.na(q)) %>% 
  mutate(lxs = map2(lx, q, lx_submax)) %>% 
  mutate(lxs_min = map_dbl(lxs, min)) %>% 
  mutate(l0_pt = map_dbl(lx, sum)) %>% 
  mutate(l0_pt_int = as.integer(round(l0_pt, 0))) %>% 
  mutate(lxs_n = map_int(lxs, length)) %>% 
  mutate(nrep = map_int(matF, ~ length(which(colSums(.x) > 0))))# %>%
  #filter(nrep > 1)

table(out_shape$q)[1:20]

sum(out_shape$q < out_shape$l0_pt) / nrow(out_shape)



sum(out_shape$lxs_min > 0.1) / nrow(out_shape)



%>%
  filter(lxs_n >= 3) %>% 
  mutate(shape_pt = map2_dbl(lxs, q+1, shape_surv2)) %>% 
  mutate(shape_l0_pt = map2_dbl(lx, l0_pt_int, ~ 1 + log(.x[.y]))) %>% 
  as_tibble() %>% 
  # filter(lx_n > 2) %>%
  filter(!is.na(shape_pt)) %>%
  # mutate(shape_pt = shape_l0_pt) %>%
  mutate(id_shape = fct_reorder(fct_drop(id), shape_pt)) %>% 
  mutate(id_l0 = fct_reorder(fct_drop(id), l0_pt))




