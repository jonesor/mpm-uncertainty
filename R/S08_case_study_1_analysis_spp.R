

### libraries
library(tidyverse)
library(ggridges)
library(cowplot)
library(Rcompadre)
library(Rage)
library(popbio)
library(popdemo)
library(gridExtra)
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
sd_files <- sd_files[grep("analysis/sds_", sd_files)]



### bind sampling distributions into single tibble
mpm_draws <- cdb_bind_rows(lapply(sd_files, rdata_load2)) %>% 
  mutate(id = as.factor(1:n())) %>% 
  cdb_unnest() %>% 
  mutate(matU = map(matU, scale_U)) %>% 
  mutate(matA = pmap(list(matU, matF, matC), ~ ..1 + ..2 + ..3)) %>% 
  mutate(start = map_int(mat, mpm_first_active)) %>% 
  mutate(rep_stages = map(matF, ~ colSums(.x) > 0))

### point estimates for parameters of interest
pt_shape <- mpm_draws %>% 
  mutate(perennial = map2_lgl(matU, start, ~ mpm_to_lx(.x, .y, N = 3)[4] > 0)) %>% 
  mutate(any_rep = map_lgl(matF, ~ any(.x > 0))) %>% 
  filter(perennial == TRUE, any_rep == TRUE) %>% 
  mutate(first_rep = map_int(rep_stages, ~ min(which(.x)))) %>% 
  mutate(rep_prop1 = pmap(list(matU, start, rep_stages), repro_prop_start)) %>% 
  mutate(lx = map2(matU, rep_prop1, lx_from_mature)) %>% 
  mutate(lx_min = map_dbl(lx, min)) %>% 
  mutate(lx_n = map_int(lx, length)) %>% 
  mutate(q = map2_int(matU, rep_prop1, qsd)) %>% 
  mutate(lxs = map2(lx, q, lx_submax)) %>% 
  mutate(lxs_min = map_dbl(lxs, min)) %>% 
  mutate(l0_pt = map_dbl(lx, sum)) %>% 
  mutate(l0_pt_int = as.integer(round(l0_pt, 0))) %>% 
  mutate(shape_pt = map2_dbl(lxs, q+1, shape_surv2)) %>% 
  mutate(shape_l0_pt = map2_dbl(lx, l0_pt_int, ~ 1 + log(.x[.y]))) %>% 
  as_tibble() %>% 
  filter(lx_n > 2) %>%
  # filter(!is.na(shape_pt)) %>%
  mutate(shape_pt = shape_l0_pt) %>% 
  mutate(id_shape = fct_reorder(fct_drop(id), shape_pt)) %>% 
  mutate(id_l0 = fct_reorder(fct_drop(id), l0_pt))

pt_other <- mpm_draws %>% 
  mutate(any_rep = map_lgl(matF, ~ any(.x > 0))) %>% 
  filter(any_rep == TRUE) %>% 
  mutate(loglam_pt = map_dbl(matA, ~ log(lambda(.x)))) %>% 
  mutate(damp_pt = map_dbl(matA, damping.ratio)) %>% 
  mutate(R0_pt = map2_dbl(matU, matF, R0)) %>%
  mutate(gen_pt = log(R0_pt) / loglam_pt) %>% 
  mutate(w = map(matA, stable.stage)) %>%
  mutate(growth_pt = map2_dbl(matU, w, ~ vr_growth(.x, weights_col = .y))) %>% 
  as_tibble() %>% 
  mutate(id_loglam = fct_reorder(fct_drop(id), loglam_pt)) %>% 
  mutate(id_damp = fct_reorder(fct_drop(id), damp_pt)) %>% 
  mutate(id_gen = fct_reorder(fct_drop(id), gen_pt)) %>% 
  mutate(id_growth = fct_reorder(fct_drop(id), growth_pt))


### point estimate of shape vs l0
ggplot(pt_shape) +
  geom_point(aes(l0_pt, shape_pt, color = lxs_min), size = 3) +
  scale_x_log10() +
  scale_color_gradient(low = "navyblue", high = "orange")


# ### plot individual hazard trajectories
# out <- df %>% 
#   as_tibble() %>% 
#   mutate(id = 1:n()) %>% 
#   filter(id == sample(id, 1)) %>% 
#   mutate(x = map(lxs, ~ seq_along(.x) - 1)) %>% 
#   mutate(hx = map(lxs, lx_to_hx)) %>% 
#   select(SpeciesAuthor, MatrixPopulation, id, q, lxs_min, x, lxs, shape_pt, hx) %>% 
#   unnest()
#   
# ggplot(out, aes(x, hx)) +
#   geom_line() +
#   ggtitle(paste(out$SpeciesAuthor[1], round(out$shape_pt, 3), sep = "; "))


## sampling distributions for derived parameters
sd_shape <- pt_shape %>%
  select(id, SpeciesAuthor, MatrixPopulation, simU, simF, q) %>%
  unnest() %>%
  left_join(select(as_tibble(pt_shape), id, start, rep_stages)) %>%
  mutate(rep_prop1 = pmap(list(simU, start, rep_stages), repro_prop_start)) %>%
  mutate(lx = map2(simU, rep_prop1, lx_from_mature)) %>%
  mutate(lx_min = map_dbl(lx, min)) %>%
  mutate(lx_n = map_int(lx, length)) %>%
  mutate(lxs = map2(lx, q, lx_submax)) %>%
  mutate(lxs_min = map_dbl(lxs, min)) %>%
  mutate(l0 = map_dbl(lx, sum)) %>%
  mutate(l0_int = as.integer(round(l0, 0))) %>%
  mutate(shape = map2_dbl(lxs, q+1, shape_surv2)) %>%
  mutate(shape_l0 = map2_dbl(lx, l0_int, ~ 1 + log(.x[.y]))) %>%
  mutate(shape = shape_l0) %>% 
  left_join(select(pt_shape, id, id_l0, id_shape, ends_with("pt")), by = "id")

# sd_other <- pt_other %>%
#   select(id, SpeciesAuthor, MatrixPopulation, simU, simF) %>%
#   unnest() %>%
#   mutate(simA = pmap(list(simU, simF), ~ ..1 + ..2)) %>% 
#   left_join(select(as_tibble(pt_other), id, start, rep_stages), by = "id") %>% 
#   mutate(loglam = map_dbl(simA, ~ log(lambda(.x)))) %>% 
#   mutate(damp = map_dbl(simA, damping.ratio)) %>% 
#   mutate(R0 = map2_dbl(simU, simF, R0)) %>%
#   mutate(gen = log(R0) / loglam) %>% 
#   mutate(w = map(simA, stable.stage)) %>%
#   mutate(growth = map2_dbl(simU, w, ~ vr_growth(.x, weights_col = .y))) %>% 
#   left_join(select(pt_other, id, starts_with("id_"), ends_with("pt")), by = "id")
# 
# ## write to file
# sd_shape <- sd_shape %>%
#   select(which(sapply(sd_shape, class) != "list"))
# 
# sd_other <- sd_other %>%
#   select(which(sapply(sd_other, class) != "list"))
# 
# save(sd_shape, file = "analysis/full_sd_spp_shape.RData")
# save(sd_other, file = "analysis/full_sd_spp_other.RData")

load(file = "analysis/full_sd_spp_shape.RData")
load(file = "analysis/full_sd_spp_other.RData")




### plot sampling distributions vs. point estimate for shape and l0
tt <- theme(panel.grid = element_blank(),
            axis.title = element_text(size = 12.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA))

p1 <- ggplot(sd_shape, aes(y = id_shape)) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  geom_density_ridges(aes(x = shape), rel_min_height = 0.01,
                      scale = 2.5, fill = "#9ebcda", size = 0.4) +
  geom_point(data = pt_shape, aes(x = shape_pt), size = 0.9) +
  # coord_flip(xlim = c(-0.3, 0.2)) +
  coord_flip(xlim = c(-0.5, 0.5)) +
  labs(y = expression(paste("Population (ranked by ", italic(S), ")")),
       x = expression(paste("Mortality trajectory shape (", italic(S), ")"))) +
  tt

p2 <- ggplot(sd_shape, aes(y = id_l0)) +
  geom_density_ridges(aes(x = l0), rel_min_height = 0.01,
                      scale = 2.5, fill = "#9ebcda", size = 0.4) +
  geom_point(data = pt_shape, aes(x = l0_pt), size = 0.9) +
  scale_x_log10() +
  coord_flip() +
  labs(y = expression(paste("Population (ranked by ", italic(L[alpha]), ")")),
       x = expression(paste("Mature life expectancy (", italic(L[alpha]), ")"))) +
  tt

# arrange plot panels
g <- rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")

# print to screen
dev.off()
quartz(height = 5.5, width = 5.5, dpi = 120)
grid.arrange(g)

# save png
# ggsave("img/sds_shape_spp.png", g, height = 5.5, width = 5.5, units = "in", dpi = 300)




### plot sampling distributions vs. point estimate for other parameters
tt <- theme(panel.grid = element_blank(),
            axis.title = element_text(size = 13),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA))

p1 <- ggplot(sd_other, aes(y = id_loglam)) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  geom_density_ridges(aes(x = loglam), rel_min_height = 0.01,
                      scale = 2.5, fill = "#9ebcda", size = 0.4) +
  geom_point(data = pt_other, aes(x = loglam_pt), size = 0.9) +
  scale_x_continuous(breaks = seq(-0.4, 0.6, 0.2)) +
  coord_flip(xlim = c(-0.4, 0.7)) +
  labs(y = expression(paste("Population (ranked by ", log~italic(lambda), ")")),
       x = expression(paste("Population growth (", log~italic(lambda), ")"))) +
  tt

# p2 <- ggplot(sd_other, aes(y = id_damp)) +
#   geom_density_ridges(aes(x = damp), rel_min_height = 0.01,
#                       scale = 3, fill = "#9ebcda", size = 0.4) +
#   geom_point(data = pt_other, aes(x = damp_pt), size = 0.9) +
#   scale_x_log10() +
#   coord_flip(xlim = c(1, 12)) +
#   labs(y = expression(paste("Population (ranked by ", italic(rho), ")")),
#        x = expression(paste("Damping ratio (", italic(rho), ")"))) +
#   tt

p3 <- ggplot(sd_other, aes(y = id_growth)) +
  geom_density_ridges(aes(x = growth), rel_min_height = 0.01,
                      scale = 2.5, fill = "#9ebcda", size = 0.4) +
  geom_point(data = pt_other, aes(x = growth_pt), size = 0.9) +
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(y = expression(paste("Population (ranked by ", italic(gamma), ")")),
       x = expression(paste("Individual growth (", italic(gamma), ")"))) +
  tt

p4 <- ggplot(sd_other, aes(y = id_gen)) +
  geom_density_ridges(aes(x = gen), rel_min_height = 0.01,
                      scale = 2.5, fill = "#9ebcda", size = 0.4) +
  geom_point(data = pt_other, aes(x = gen_pt), size = 0.9) +
  scale_x_log10() +
  coord_flip(xlim = c(1, 350)) +
  labs(y = expression(paste("Population (ranked by ", italic(T), ")")),
       x = expression(paste("Generation time (", italic(T), ")"))) +
  tt

# arrange plot panels
g <- rbind(ggplotGrob(p1), ggplotGrob(p3), ggplotGrob(p4), size = "last")

# print to screen
dev.off()
quartz(height = 7, width = 5, dpi = 140)
grid.arrange(g)

# save png
# ggsave("img/sd_other_spp.png", g, height = 7, width = 5, units = "in", dpi = 300)



### prep df for shape vs. pace analysis
df_shape <- sd_shape %>% 
  mutate(log_l0 = log10(l0)) %>% 
  group_by(SpeciesAuthor, MatrixPopulation) %>% 
  summarize(shape_med = quantile(shape, 0.500),
            shape_mean = mean(shape),
            shape_se = sd(shape),
            shape_low = quantile(shape, 0.025),
            shape_upp = quantile(shape, 0.975),
            l0_med = quantile(l0, 0.500),
            l0_mean = mean(l0),
            l0_se = sd(l0),
            log_l0_med = quantile(log_l0, 0.500),
            log_l0_mean = mean(log_l0),
            log_l0_se = sd(log_l0),
            l0_low = quantile(l0, 0.025),
            l0_upp = quantile(l0, 0.975),
            log_l0_low = quantile(log_l0, 0.025),
            log_l0_upp = quantile(log_l0, 0.975)) %>% 
  ungroup() %>% 
  left_join(pt_shape) %>% 
  mutate(spp_int = as.integer(as.factor(SpeciesAuthor)))




### variance components
vw <- mean(df_shape$l0_se^2)
va <- var(df_shape$l0_mean)
vw / (va + vw)

vw <- mean(df_shape$shape_se^2)
va <- var(df_shape$shape_mean)
vw / (va + vw)

# ggplot(df_shape, aes(x = l0_pt, y = shape_pt)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   scale_x_log10()
# 
# ggplot(df_shape, aes(x = l0_mean, y = shape_mean, col = lxs_min)) +
#   geom_point(size = 2.5) +
#   geom_smooth(method = "lm") +
#   scale_x_log10() +
#   scale_color_gradient(low = "navyblue", high = "orange")
# 
# ggplot(df_shape) +
#   geom_segment(aes(x = l0_pt, y = shape_pt, xend = l0_mean, yend = shape_mean),
#                size = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
#   geom_point(aes(l0_pt, shape_pt)) +
#   scale_x_log10()



### Model relationship between l0 and shape, assuming no sampling uncertainty

# compile stan models
stan_regress_hier <- stan_model("stan/regress2.stan")
stan_regress_hier_error <- stan_model("stan/regress_error2.stan")

# arrange data for stan
x_cent <- mean(log10(df_shape$l0_pt))
dat_stan <- list(N = nrow(df_shape),
                 N_spp = length(unique(df_shape$SpeciesAuthor)),
                 spp = df_shape$spp_int,
                 x = log10(df_shape$l0_pt) - x_cent,
                 y = df_shape$shape_pt)

# fit stan model
stan_fit <- sampling(
  stan_regress_hier,
  data = dat_stan,
  warmup = 2000,
  iter = 4000,
  thin = 2,
  chains = 2
)

# model diagnostics
# library(shinystan)
# launch_shinystan(stan_fit)

# extract posterior samples for intercept and slope
mu_alpha <- rstan_extract(stan_fit, "mu_alpha")
mu_beta <- rstan_extract(stan_fit, "mu_beta")

# extract posterior samples for best fit line and 95% credible interval
pred_x <- seq(min(dat_stan$x), max(dat_stan$x), length.out = 50)

pred_reg <- tibble(mu_alpha, mu_beta, pred_x = list(pred_x)) %>% 
  mutate(pred = pmap(list(mu_alpha, mu_beta, pred_x), ~ ..1 + ..2 * ..3)) %>% 
  dplyr::select(pred_x, pred) %>% 
  unnest() %>% 
  mutate(pred_x = 10^(pred_x + x_cent)) %>% 
  group_by(pred_x) %>% 
  summarize(pred_med = quantile(pred, 0.500),
            pred_low = quantile(pred, 0.025),
            pred_upp = quantile(pred, 0.975))



### Model relationship between l0 and shape, with sampling uncertainty
# arrange data for stan
x_cent_error <- mean(df_shape$log_l0_mean)
dat_stan <- list(N = nrow(df_shape),
                 N_spp = length(unique(df_shape$spp_int)),
                 spp = df_shape$spp_int,
                 x_mean = df_shape$log_l0_mean - x_cent_error,
                 x_se = df_shape$log_l0_se,
                 y_mean = df_shape$shape_mean,
                 y_se = df_shape$shape_se)

# fit stan model
stan_fit_error <- sampling(
  stan_regress_hier_error,
  data = dat_stan,
  warmup = 3000,
  iter = 4000,
  thin = 2,
  chains = 2,
  control = list(adapt_delta = 0.95, stepsize  = 0.05, max_treedepth = 12)
)

# model diagnostics
# shinystan::launch_shinystan(stan_fit_error)

# posterior samples for intercept and slope
mu_alpha_error <- rstan_extract(stan_fit_error, "mu_alpha")
mu_beta_error <- rstan_extract(stan_fit_error, "mu_beta")

quantile(mu_beta, c(0.025, 0.500, 0.975))
quantile(mu_beta_error, c(0.025, 0.500, 0.975))

df_beta_error <- data.frame(mu_beta_error)
pred_x_error <- seq(min(df_shape$log_l0_low - x_cent_error),
                    max(df_shape$log_l0_upp - x_cent_error),
                    length.out = 50)

# posterior samples for best fit line and 95% credible interval
pred_error <- tibble(mu_alpha_error, mu_beta_error, pred_x = list(pred_x_error)) %>% 
  mutate(pred = pmap(list(mu_alpha_error, mu_beta_error, pred_x), ~ ..1 + ..2 * ..3)) %>% 
  dplyr::select(pred_x, pred) %>% 
  unnest() %>% 
  mutate(pred_x = 10^(pred_x + x_cent_error)) %>% 
  group_by(pred_x) %>% 
  summarize(pred_med = quantile(pred, 0.500),
            pred_low = quantile(pred, 0.025),
            pred_upp = quantile(pred, 0.975))



### prepare plot data
# left panel
lev <- c("Model of point estimates", "Model with sampling uncertainty")

# fit line
pred_full <- rbind(
  mutate(pred_reg, model = lev[1]),
  mutate(pred_error, model = lev[2])
) %>% mutate(model = factor(model, levels = lev))

# points and error bars
bars_reg <- df_shape %>% 
  select(SpeciesAuthor, MatrixPopulation,
         l0_pt, shape_pt,
         l0_low, l0_upp,
         shape_low, shape_upp) %>% 
  mutate(model = lev[1],
         l0_med = l0_pt, shape_med = shape_pt,
         l0_low = NA_real_, l0_upp = NA_real_,
         shape_low = NA_real_, shape_upp = NA_real_)

bars_err <- df_shape %>% 
  select(SpeciesAuthor, MatrixPopulation,
         l0_med, l0_low, l0_upp,
         shape_med, shape_low, shape_upp) %>% 
  mutate(model = lev[2], l0_pt = NA_real_, shape_pt = NA_real_)

bars_full <- bind_rows(bars_reg, bars_err) %>% 
  mutate(model = factor(model, levels = lev))

df_alpha <- bind_rows(
  tibble(alpha = mu_alpha, model = lev[1]),
  tibble(alpha = mu_alpha_error, model = lev[2])
) %>% mutate(model = factor(model, levels = lev))

df_beta <- bind_rows(
  tibble(beta = mu_beta, model = lev[1]),
  tibble(beta = mu_beta_error, model = lev[2])
) %>% mutate(model = factor(model, levels = lev))


### plot
tt <- theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 11.5),
        axis.ticks = element_line(size = 0.4))

p1 <- ggplot(pred_full) +
  geom_point(data = bars_full, aes(x = l0_pt, y = shape_pt), size = 1.3) +
  geom_linerange(data = bars_full, aes(x = l0_med, ymin = shape_low, ymax = shape_upp), size = 0.3, alpha = 0.6) +
  geom_errorbarh(data = bars_full, aes(y = shape_med, xmin = l0_low, xmax = l0_upp), size = 0.3, alpha = 0.6) +
  geom_line(aes(x = pred_x, y = pred_med), col = "darkblue") +
  geom_ribbon(aes(x = pred_x, ymin = pred_low, ymax = pred_upp), fill = "darkblue", alpha = 0.2) +
  scale_x_log10() +
  coord_cartesian(ylim = c(-0.3, 0.2)) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = expression(paste("Mature life expectancy (", italic(L[alpha]), ")")),
       y = expression(paste("Shape of mortality trajectory (", italic(S), ")"))) +
  tt

p2 <- ggplot(df_beta, aes(x = beta)) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  geom_density(fill = "darkred", alpha = 0.4, size = 0) +
  coord_cartesian(xlim = c(-0.12, 0.12)) +
  facet_wrap(~ model, ncol = 1) +
  labs(x = expression(paste("Slope coefficient (", italic(beta), ")")), y = "Posterior density") +
  tt

# combine both plots
p <- plot_grid(p1, p2, labels = c("A", "B"), rel_widths = c(1.08, 1), nrow = 1)

# print to screen
dev.off()
quartz(height = 4.5, width = 6.25, dpi = 160)
print(p)

# save to png
# ggsave2("img/shape_spp.png", p, height = 4.5, width = 6.25)



### posterior summary
# intercept
median(mu_alpha)
median(mu_alpha_error)

var(mu_alpha)
var(mu_alpha_error)

# slope
median(mu_beta)
median(mu_beta_error)

length(mu_beta[mu_beta > 0]) / length(mu_beta)
length(mu_beta_error[mu_beta_error > 0]) / length(mu_beta_error)




### variance components model
stan_varcomp <- stan_model("stan/varcomp.stan")

df_other <- sd_other %>% 
  group_by(SpeciesAuthor, MatrixPopulation) %>% 
  summarize(loglam_mean = mean(loglam),
            loglam_se = sd(loglam),
            damp_mean = mean(damp),
            damp_se = sd(damp),
            gen_mean = mean(log(gen)),
            gen_se = sd(log(gen)),
            growth_mean = mean(logit(growth)),
            growth_se = sd(logit(growth))) %>% 
  ungroup() %>% 
  left_join(pt_other)

dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$log_l0_mean,
                 y_se = df_shape$log_l0_se)

dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$shape_mean,
                 y_se = df_shape$shape_se)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$loglam_mean,
                 y_se = df_other$loglam_se)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$damp_mean,
                 y_se = df_other$damp_se)

dat_stan <- list(N = nrow(df_other),
                 y_mean = df_other$gen_mean,
                 y_se = df_other$gen_se)

df_growth <- filter(df_other, !is.na(growth_se))

dat_stan <- list(N = nrow(df_growth),
                 y_mean = df_growth$growth_mean,
                 y_se = df_growth$growth_se)


dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$shape_mean,
                 y_se = df_shape$shape_se)

dat_stan <- list(N = nrow(df_shape),
                 y_mean = df_shape$log_l0_mean,
                 y_se = df_shape$log_l0_se)


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

var_a_pt <- var(dat_stan$y_mean)
var_a <- rstan_extract(stan_fit_varcomp, "var_a")
quantile(var_a_pt / var_a, c(0.025, 0.500, 0.975))






### examine hazard trajectories for select populations

pt <- pt_shape %>% 
  filter(shape_pt > 0.1) %>% 
  mutate(lx = map(lx, ~ .x[1:11])) %>%
  mutate(hx = map(lx, lx_to_px)) %>%
  mutate(x = map(lx, ~ seq_along(.x) - 1)) %>% 
  select(SpeciesAuthor, MatrixPopulation, x, hx) %>% 
  unnest() #%>% 
# filter(!is.na(hxs))

sdist <- sd_shape %>% 
  filter(shape_pt > 0.1) %>% 
  mutate(lx = map(lx, ~ .x[1:11])) %>%
  mutate(hx = map(lx, lx_to_px)) %>%
  mutate(x = map(lx, ~ seq_along(.x) - 1)) %>% 
  select(SpeciesAuthor, MatrixPopulation, x, hx) %>% 
  group_by(SpeciesAuthor, MatrixPopulation) %>% 
  mutate(rep = 1:n()) %>% 
  slice(sample(rep, 50)) %>% 
  ungroup() %>% 
  unnest()# %>% 
# filter(!is.na(hxs))

pt$matF[[3]]
pt$matU[[3]] %>% colSums()
pt$Authors[2]

lapply(pt$matU, colSums)


ggplot(sdist, aes(x, hx)) +
  geom_line(aes(group = rep), alpha = 0.4, size = 0.3) +
  geom_line(data = pt, col = "darkred", size = 1.2) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  # scale_y_log10() +
  # coord_cartesian(ylim = c(0, 10)) +
  facet_wrap(~ SpeciesAuthor, ncol = 1)



