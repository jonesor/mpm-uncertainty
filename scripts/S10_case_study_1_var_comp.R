
### libraries
source("code/setup.R")
setup_packages(c("tidyverse", "popbio", "popdemo", "Rcompadre", "Rage",
                 "rstan", "loo"))
setup_rstan()
source("code/functions.R")
seed <- 12345
set.seed(seed)


### set options for rstan library
# handled in setup_rstan()


### load compadre data
compadre <- cdb_fetch("data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData")


### load study-specific sampling distribution files
sd_files <- paste0("data/derived/analysis_cache/", list.files("data/derived/analysis_cache"))
sd_files <- sd_files[grep("/sd_", sd_files)]


### bind sampling distributions into single tibble
mpm_draws <- cdb_bind_rows(map(sd_files, rdata_load)) %>% 
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
  mutate(L_pt = map2_dbl(matU, rep_prop1, life_expect)) %>% 
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
  mutate(L = map2_dbl(simU, rep_prop1, life_expect)) %>%
  mutate(S = map_dbl(lx, Rage::shape_surv)) %>%
  left_join(select(pt_shape, id, id_L, id_S, ends_with("pt")), by = "id")



### load sampling distributions
load(file = "data/derived/analysis_cache/full_sd_shape.RData")
load(file = "data/derived/analysis_cache/full_sd_other.RData")
if (!exists("sd_other") && exists("sd_other_out")) {
  sd_other <- sd_other_out
}




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
            pmature_se = sd(logit(pmature))) %>% 
  ungroup() %>% 
  left_join(pt_other) %>% 
  mutate(damp_pt = log10(damp_pt)) %>% 
  mutate(gen_pt = log10(gen_pt)) %>% 
  mutate(pmature_pt = logit(pmature_pt))



### Variance components analysis
stan_varcomp <- stan_model("models/varcomp.stan")



dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$log_L_mean,
                 y_se = df_shape$log_L_se,
                 y_pt = log10(df_shape$L_pt))

dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$S_mean,
                 y_se = df_shape$S_se,
                 y_pt = df_shape$S_pt)

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
                 y_mean = df_other$pmature_mean,
                 y_se = df_other$pmature_se,
                 y_pt = df_other$pmature_pt)

dat_stan$y_se <- pmax(dat_stan$y_se, 1e-6)

theta_x <- if (dat_stan$N == nrow(df_other)) {
  df_other$SpeciesAuthor
} else {
  df_shape$SpeciesAuthor
}

# fit stan model
stan_fit_varcomp <- sampling(
  stan_varcomp,
  data = dat_stan,
  warmup = 3000,
  iter = 4000,
  thin = 2,
  chains = 2,
  control = list(adapt_delta = 0.95, stepsize  = 0.05, max_treedepth = 12),
  seed = seed
)

pvar_w <- rstan_extract(stan_fit_varcomp, "pvar_w")
quantile(pvar_w, c(0.025, 0.500, 0.975))




df_theta <- posterior_vec(stan_fit_varcomp, x = theta_x, "theta") %>% 
  mutate(x = fct_reorder(x, med))

ggplot(df_theta, aes(x = x)) +
  geom_point(aes(y = med)) +
  geom_errorbar(aes(ymin = low95, ymax = upp95)) +
  coord_flip()


var_a_pt <- var(dat_stan$y_pt)
var_a <- rstan_extract(stan_fit_varcomp, "var_a")
quantile(var_a_pt / var_a, c(0.025, 0.500, 0.975))


var(df_shape$S_pt) / var(df_shape$S_mean)
var(log10(df_shape$L_pt)) / var(df_shape$log_L_mean)

var(df_other$loglam_pt) / var(df_other$loglam_mean)
var(df_other$damp_pt) / var(df_other$damp_mean)
var(df_other$gen_pt) / var(df_other$gen_mean)
var(df_other$pmature_pt) / var(df_other$pmature_mean)
