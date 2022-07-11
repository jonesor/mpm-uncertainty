
### libraries
library(tidyverse)
library(popbio)
library(popdemo)
library(Rcompadre)
library(Rage)
library(rstan)
library(loo)
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


### point estimates for parameters of interest
pt_shape <- mpm_draws %>% 
  as_tibble() %>% 
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
  unnest(cols = c(simU, simF)) %>%
  left_join(select(pt_shape, id, start, rep_stages)) %>%
  mutate(rep_prop1 = pmap(list(simU, start, rep_stages), Rage::mature_distrib)) %>%
  mutate(lx = pmap(list(simU, rep_prop1, q),
                   ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3), lx_crit = -1)) %>%
  mutate(L = map2_dbl(simU, rep_prop1, Rage::life_expect)) %>%
  mutate(S = map_dbl(lx, Rage::shape_surv)) %>%
  left_join(select(pt_shape, id, id_L, id_S, ends_with("pt")), by = "id")

# sd_other <- pt_other %>%
#   select(id, SpeciesAuthor, MatrixPopulation, simU, simF) %>%
#   unnest() %>%
#   left_join(select(pt_other, id, exclude_stages), by = "id") %>%
#   mutate(simA = pmap(list(simU, simF), ~ ..1 + ..2)) %>%
#   left_join(select(as_tibble(pt_other), id, start, rep_stages), by = "id") %>%
#   mutate(loglam = map_dbl(simA, ~ log(popbio::lambda(.x)))) %>%
#   mutate(damp = map_dbl(simA, popbio::damping.ratio)) %>%
#   mutate(gen = map2_dbl(simU, simF, Rage::gen_time)) %>%
#   mutate(pmature = pmap_dbl(list(simU, simF, start), Rage::mature_prob)) %>%
#   mutate(growth = pmap_dbl(list(simU, exclude_stages),
#                            ~ Rage::vr_growth(..1, exclude = ..2))) %>%
#   mutate(elast = pmap_dbl(list(simU, simF, exclude_stages),
#                           ~ perturb_cust(..1, ..2, exclude = ..3,
#                                          type = "elasticity")$progr)) %>%
#   left_join(select(pt_other, id, starts_with("id_"), ends_with("pt")), by = "id")


# ### write to file
# sd_shape <- sd_shape %>%
#   select(which(sapply(sd_shape, class) != "list"))
# 
# sd_other <- sd_other %>%
#   select(which(sapply(sd_other, class) != "list"))
# 
# save(sd_shape, file = "analysis/full_sd_shape.RData")
# save(sd_other, file = "analysis/full_sd_other.RData")


sd_shape
sd_full <- full_join()


### load sampling distributions
load(file = "analysis/full_sd_shape.RData")
load(file = "analysis/full_sd_other.RData")




### prep df for variance component analysis
df_shape <- sd_shape %>% 
  mutate(log_L = log10(L)) %>% 
  group_by(SpeciesAuthor, MatrixPopulation) %>% 
  summarize(S_med = quantile(S, 0.500),
            S_mean = mean(S),
            S_se = sd(S),
            S_low = quantile(S, 0.025),
            S_upp = quantile(S, 0.975),
            L_med = quantile(L, 0.500),
            L_mean = mean(L),
            L_se = sd(L),
            log_L_med = quantile(log_L, 0.500),
            log_L_mean = mean(log_L),
            log_L_se = sd(log_L),
            L_low = quantile(L, 0.025),
            L_upp = quantile(L, 0.975),
            log_L_low = quantile(log_L, 0.025),
            log_L_upp = quantile(log_L, 0.975)) %>% 
  ungroup() %>% 
  left_join(pt_shape) %>% 
  mutate(spp_int = as.integer(as.factor(SpeciesAuthor)))

df_other <- sd_other %>% 
  group_by(SpeciesAuthor, MatrixPopulation) %>% 
  summarize(loglam_mean = mean(loglam),
            loglam_se = sd(loglam),
            damp_mean = mean(log10(damp)),
            damp_se = sd(log10(damp)),
            gen_mean = mean(log10(gen)),
            gen_se = sd(log10(gen)),
            pmature_mean = mean(logit(pmature)),
            pmature_se = sd(logit(pmature)),
            growth_mean = mean(logit(growth)),
            growth_se = sd(logit(growth)),
            elast_mean = mean(elast),
            elast_se = sd(elast)) %>% 
  ungroup() %>% 
  left_join(pt_other) %>% 
  mutate(damp_pt = log10(damp_pt)) %>% 
  mutate(gen_pt = log10(gen_pt)) %>% 
  mutate(pmature_pt = logit(pmature_pt)) %>% 
  mutate(growth_pt = logit(growth_pt))



### Variance components analysis
stan_varcomp <- stan_model("stan/varcomp.stan")



dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$log_l0_mean,
                 y_se = df_shape$log_l0_se,
                 y_pt = log10(df_shape$l0_pt))

dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$shape_mean,
                 y_se = df_shape$shape_se,
                 y_pt = df_shape$shape_pt)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$loglam_mean,
                 y_se = df_other$loglam_se,
                 y_pt = df_other$loglam_pt)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$damp_mean,
                 y_se = df_other$damp_se,
                 y_pt = df_other$damp_pt)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$gen_mean,
                 y_se = df_other$gen_se,
                 y_pt = df_other$gen_pt)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$growth_mean,
                 y_se = df_other$growth_se,
                 y_pt = df_other$growth_pt)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$pmature_mean,
                 y_se = df_other$pmature_se,
                 y_pt = df_other$pmature_pt)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$elast_mean,
                 y_se = df_other$elast_se,
                 y_pt = df_other$elast_pt)

# fit stan model
stan_fit_varcomp <- sampling(
  stan_varcomp,
  data = dat_stan,
  warmup = 3000,
  iter = 4000,
  thin = 2,
  chains = 2,
  control = list(adapt_delta = 0.95, stepsize  = 0.05, max_treedepth = 12)
)

pvar_w <- rstan_extract(stan_fit_varcomp, "pvar_w")
quantile(pvar_w, c(0.025, 0.500, 0.975))




df_theta <- posterior_vec(stan_fit_varcomp, x = df_shape$id_shape, "theta") %>% 
  mutate(x = fct_reorder(x, med))

ggplot(df_theta, aes(x = x)) +
  # geom_hline(yintercept = 0, alpha = 0.5) +
  geom_point(aes(y = med)) +
  geom_errorbar(aes(ymin = low95, ymax = upp95)) +
  coord_flip()


var_a_pt <- var(dat_stan$y_pt)
var_a <- rstan_extract(stan_fit_varcomp, "var_a")
quantile(var_a_pt / var_a, c(0.025, 0.500, 0.975))


var(df_shape$shape_pt) / var(df_shape$shape_mean)
var(log10(df_shape$l0_pt)) / var(df_shape$log_l0_mean)

var(df_other$loglam_pt) / var(df_other$loglam_mean)
var(df_other$damp_pt) / var(df_other$damp_mean)
var(df_other$gen_pt) / var(df_other$gen_mean)
var(df_growth$growth_pt) / var(df_growth$growth_mean)



