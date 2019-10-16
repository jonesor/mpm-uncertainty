
### libraries
library(tidyverse)
library(popbio)
library(popdemo)
library(Rcompadre)
library(Rage)
source("R/functions.R")


### set options for rstan library
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


### load compadre data
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


### load study-specific sampling distribution files
sd_files <- paste0("analysis/", list.files("analysis"))
sd_files <- sd_files[grep("/sd_", sd_files)]


### bind sampling distributions into single tibble
mpm_draws <- cdb_bind_rows(lapply(sd_files, rdata_load)) %>% 
  mutate(id = as.factor(1:n())) %>% 
  cdb_unnest() %>% 
  mutate(any_repro = map_lgl(matF, ~ any(.x > 0))) %>% 
  filter(any_repro == TRUE) %>% # make sure some repro
  mutate(matU = map(matU, scale_U)) %>% 
  mutate(matA = pmap(list(matU, matF, matC), ~ ..1 + ..2 + ..3)) %>% 
  mutate(start = map_int(mat, Rcompadre::mpm_first_active)) %>% 
  mutate(rep_stages = map(matF, ~ colSums(.x) > 0)) %>% 
  mutate(exclude_stages = map(MatrixClassOrganized,
                              ~ ifelse(.x == "active", FALSE, TRUE)))

# summary counts
length(unique(mpm_draws$Authors))
length(unique(mpm_draws$MatrixPopulation))
length(unique(mpm_draws$SpeciesAccepted))

# n_pops by species
mpm_draws %>% 
  as_tibble() %>% 
  count(SpeciesAuthor) %>% 
  arrange(desc(n))


### point estimates for parameters of interest
pt_shape <- mpm_draws %>% 
  mutate(rep_prop1 = pmap(list(matU, start, rep_stages), Rage::mature_distrib)) %>% 
  mutate(lx4 = map2_dbl(matU, rep_prop1,
                        ~ Rage::mpm_to_lx(.x, .y, lx_crit = -1, xmax = 3)[4])) %>% 
  filter(lx4 > 0) %>% # check perennial
  mutate(q = map2_int(matU, rep_prop1, Rage::qsd_converge, conv = 0.01, N = 1e5)) %>% 
  filter(q >= 3) %>%   # make sure at least 3 time steps
  mutate(lx = pmap(list(matU, rep_prop1, q),
                   ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3), lx_crit = -1)) %>% 
  mutate(L_pt = map2_dbl(matU, rep_prop1, Rage::life_expect)) %>% 
  mutate(lx_min = map_dbl(lx, min)) %>% 
  mutate(S_pt = map_dbl(lx, Rage::shape_surv)) %>% 
  as_tibble() %>% 
  mutate(id_L = fct_reorder(fct_drop(id), L_pt)) %>%
  mutate(id_S = fct_reorder(fct_drop(id), S_pt))
  

pt_other <- mpm_draws %>% 
  mutate(loglam_pt = map_dbl(matA, ~ log(popbio::lambda(.x)))) %>% 
  mutate(damp_pt = map_dbl(matA, popbio::damping.ratio)) %>% 
  mutate(gen_pt = map2_dbl(matU, matF, Rage::gen_time)) %>%
  mutate(pmature_pt = pmap_dbl(list(matU, matF, start), Rage::mature_prob)) %>% 
  mutate(growth_pt = pmap_dbl(list(matU, exclude_stages),
                              ~ Rage::vr_growth(..1, exclude = ..2))) %>% 
  mutate(elast_pt = pmap_dbl(list(matU, matF, exclude_stages),
                             ~ perturb_cust(..1, ..2, exclude = ..3,
                                            type = "elasticity")$progr)) %>% 
  as_tibble() %>% 
  mutate(id_loglam = fct_reorder(fct_drop(id), loglam_pt)) %>%
  mutate(id_damp = fct_reorder(fct_drop(id), damp_pt)) %>%
  mutate(id_gen = fct_reorder(fct_drop(id), gen_pt)) %>%
  mutate(id_pmature = fct_reorder(fct_drop(id), pmature_pt)) %>%
  mutate(id_growth = fct_reorder(fct_drop(id), growth_pt)) %>%
  mutate(id_elast = fct_reorder(fct_drop(id), elast_pt))



### sampling distributions for derived parameters
sd_shape <- pt_shape %>%
  select(id, SpeciesAuthor, MatrixPopulation, simU, simF, q) %>%
  unnest() %>%
  left_join(select(pt_shape, id, start, rep_stages)) %>%
  mutate(rep_prop1 = pmap(list(simU, start, rep_stages), Rage::mature_distrib)) %>%
  mutate(lx = pmap(list(simU, rep_prop1, q),
                   ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3), lx_crit = -1)) %>%
  mutate(L = map2_dbl(simU, rep_prop1, Rage::life_expect)) %>%
  mutate(S = map_dbl(lx, Rage::shape_surv))

sd_other <- pt_other %>%
  select(id, SpeciesAuthor, MatrixPopulation, simU, simF) %>%
  unnest() %>%
  left_join(select(pt_other, id, exclude_stages), by = "id") %>%
  mutate(simA = pmap(list(simU, simF), ~ ..1 + ..2)) %>%
  left_join(select(as_tibble(pt_other), id, start, rep_stages), by = "id") %>%
  mutate(loglam = map_dbl(simA, ~ log(popbio::lambda(.x)))) %>%
  mutate(damp = map_dbl(simA, popbio::damping.ratio)) %>%
  mutate(gen = map2_dbl(simU, simF, Rage::gen_time)) %>%
  mutate(pmature = pmap_dbl(list(simU, simF, start), Rage::mature_prob)) %>%
  mutate(growth = pmap_dbl(list(simU, exclude_stages),
                           ~ Rage::vr_growth(..1, exclude = ..2))) %>%
  mutate(elast = pmap_dbl(list(simU, simF, exclude_stages),
                          ~ perturb_cust(..1, ..2, exclude = ..3,
                                         type = "elasticity")$progr))






### write to file
pt_shape_out <- pt_shape %>% 
  dplyr::select(id, SpeciesAuthor, MatrixPopulation, q,
                ends_with("_pt"), starts_with("id_"))

pt_other_out <- pt_other %>% 
  dplyr::select(id, SpeciesAuthor, MatrixPopulation,
                ends_with("_pt"), starts_with("id_"))

sd_shape_out <- sd_shape %>% 
  left_join(dplyr::select(pt_shape_out, id, starts_with("id_"))) %>% 
  dplyr::select(id, SpeciesAuthor, MatrixPopulation, start, L, S, starts_with("id_"))

sd_other_out <- sd_other %>% 
  left_join(dplyr::select(pt_other_out, id, starts_with("id_"))) %>% 
  dplyr::select(id, SpeciesAuthor, MatrixPopulation, start,
                loglam, damp, gen, pmature, starts_with("id_"))

levels()


# sd_shape <- sd_shape %>%
#   select(which(sapply(sd_shape, class) != "list"))
# 
# sd_other <- sd_other %>%
#   select(which(sapply(sd_other, class) != "list"))
# 


save(pt_shape, file = "analysis/full_pt_shape.RData")
save(pt_other, file = "analysis/full_pt_other.RData")

save(sd_shape, file = "analysis/full_sd_shape.RData")
save(sd_other, file = "analysis/full_sd_other.RData")



### load sampling distributions
load(file = "analysis/full_sd_shape.RData")
load(file = "analysis/full_sd_other.RData")

