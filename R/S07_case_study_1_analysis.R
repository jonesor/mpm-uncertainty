
### libraries
library(tidyverse)
library(popbio)
library(popdemo)
library(Rcompadre)
library(Rage)
library(ggridges)
library(cowplot)
library(gridExtra)
library(rstan)
library(loo)
source("R/functions.R")


### set options for rstan library
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


### load compadre data
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X_Corrected.RData")


### load sampling distributions
load(file = "analysis/full_pt_shape.RData")
load(file = "analysis/full_pt_other.RData")

load(file = "analysis/full_sd_shape.RData")
load(file = "analysis/full_sd_other.RData")



### plot sampling distributions vs. point estimate for shape and l0
tt <- theme(panel.grid = element_blank(),
            axis.title = element_text(size = 12.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA))

p1 <- ggplot(sd_shape_out, aes(y = id_S)) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  geom_density_ridges(aes(x = S), rel_min_height = 1e-2,
                      scale = 3, fill = "#9ebcda", size = 0.4) +
  geom_point(data = pt_shape_out, aes(x = S_pt), size = 0.9) +
  annotate("text", x = Inf, y = 2.6, label = "A", vjust = 1.5, size = 5, fontface = "bold") +
  coord_flip(xlim = c(-0.3, 0.2)) +
  labs(y = expression(paste("Population (ranked by ", italic(S), ")")),
       x = expression(paste("Mortality trajectory shape (", italic(S), ")"))) +
  tt

p2 <- ggplot(sd_shape_out, aes(y = id_L)) +
  geom_density_ridges(aes(x = L), rel_min_height = 1e-2,
                      scale = 3, fill = "#9ebcda", size = 0.4) +
  geom_point(data = pt_shape_out, aes(x = L_pt), size = 0.9) +
  annotate("text", x = Inf, y = 2.6, label = "B", vjust = 1.5, size = 5, fontface = "bold") +
  scale_x_log10(limits = c(1.2, 1500)) +
  coord_flip() +
  labs(y = expression(paste("Population (ranked by ", italic(L), ")")),
       x = expression(paste("Mature life expectancy (", italic(L), ")"))) +
  tt

# arrange plot panels
g <- rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")

# print to screen
dev.off()
quartz(height = 5.5, width = 5.5, dpi = 160)
grid.arrange(g)

# save png
# ggsave("img/raw/Fig_2.png", g, height = 5.5, width = 5.5, units = "in", dpi = 300)



### plot sampling distributions vs. point estimate for other parameters
tt <- theme(panel.grid = element_blank(),
            axis.title = element_text(size = 10),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 8.5, angle = 90, hjust = 0.5),
            axis.ticks.x = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(size = 0.5, fill = NA))

pt_size <- 0.7
pt_shp <- 19
rdg_scale <- 7
rdg_size <- 0.2

p1 <- ggplot(sd_other_out, aes(y = id_loglam)) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  geom_density_ridges(aes(x = loglam), rel_min_height = 0.01,
                      scale = rdg_scale, fill = "#9ebcda", size = rdg_size) +
  geom_point(data = pt_other_out, aes(x = loglam_pt), size = pt_size, shape = pt_shp) +
  annotate("text", x = Inf, y = 4, label = "A", vjust = 1.5, size = 4, fontface = "bold") +
  scale_x_continuous(breaks = seq(-0.4, 0.6, 0.2)) +
  coord_flip(xlim = c(-0.4, 0.7)) +
  labs(y = expression(paste("Population (ranked by ", log~italic(lambda), ")")),
       x = expression(paste("Population growth (", log~italic(lambda), ")"))) +
  tt

p2 <- ggplot(sd_other_out, aes(y = id_damp)) +
  geom_density_ridges(aes(x = damp), rel_min_height = 0.01,
                      scale = rdg_scale, fill = "#9ebcda", size = 0.3) +
  geom_point(data = pt_other_out, aes(x = damp_pt), size = pt_size, shape = pt_shp) +
  annotate("text", x = Inf, y = 4, label = "B", vjust = 1.5, size = 4, fontface = "bold") +
  scale_x_log10() +
  coord_flip() +
  labs(y = expression(paste("Population (ranked by ", italic(rho), ")")),
       x = expression(paste("Damping ratio (", italic(rho), ")"))) +
  tt

p3 <- ggplot(sd_other_out, aes(y = id_pmature)) +
  geom_density_ridges(aes(x = pmature), rel_min_height = 0.01,
                      scale = rdg_scale, fill = "#9ebcda", size = rdg_size) +
  geom_point(data = pt_other_out, aes(x = pmature_pt), size = pt_size, shape = pt_shp) +
  annotate("text", x = Inf, y = 4, label = "C", vjust = 1.5, size = 4, fontface = "bold") +
  coord_flip() +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(y = expression(paste("Population (ranked by Pr[Maturity])")),
       x = expression(paste("Pr[Maturity]"))) +
  tt

p4 <- ggplot(sd_other_out, aes(y = id_gen)) +
  geom_density_ridges(aes(x = gen), rel_min_height = 0.01,
                      scale = rdg_scale, fill = "#9ebcda", size = rdg_size) +
  geom_point(data = pt_other_out, aes(x = gen_pt), size = pt_size, shape = pt_shp) +
  annotate("text", x = Inf, y = 4, label = "D", vjust = 1.5, size = 4, fontface = "bold") +
  scale_x_log10() +
  coord_flip() +
  labs(y = expression(paste("Population (ranked by ", italic(T), ")")),
       x = expression(paste("Generation time (", italic(T), ")"))) +
  tt

p5 <- ggplot(sd_other_out, aes(y = id_growth)) +
  geom_density_ridges(aes(x = growth), rel_min_height = 0.01,
                      scale = rdg_scale, fill = "#9ebcda", size = rdg_size) +
  geom_point(data = pt_other_out, aes(x = growth_pt), size = pt_size, shape = pt_shp) +
  annotate("text", x = Inf, y = 4, label = "E", vjust = 1.5, size = 4, fontface = "bold") +
  coord_flip() +
  labs(y = expression(paste("Population (ranked by ", italic(gamma), ")")),
       x = expression(paste("Pr[Growth|Survival] (", italic(gamma), ")"))) +
  tt

p6 <- ggplot(sd_other_out, aes(y = id_elast)) +
  geom_density_ridges(aes(x = elast), rel_min_height = 0.01,
                      scale = rdg_scale, fill = "#9ebcda", size = rdg_size) +
  geom_point(data = pt_other_out, aes(x = elast_pt), size = pt_size, shape = pt_shp) +
  annotate("text", x = Inf, y = 4, label = "F", vjust = 1.5, size = 4, fontface = "bold") +
  coord_flip() +
  labs(y = expression(paste("Population (ranked by ", italic(E[pi]), ")")),
       x = expression(paste("Elasticity of progression (", italic(E[pi]), ")"))) +
  tt

# arrange plot panels
g1 <- rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last")
g2 <- rbind(ggplotGrob(p4), ggplotGrob(p5), ggplotGrob(p6), size = "last")
g <- cbind(g1, g2, size = "last")

# print to screen
dev.off()
quartz(height = 6, width = 6.5, dpi = 150)
grid.arrange(g)

# save png
# ggsave("img/sd_other_out.png", g, height = 6, width = 6.5, units = "in", dpi = 300)






### prep df for shape vs. pace analysis
df_shape <- sd_shape_out %>% 
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
  left_join(pt_shape_out) %>% 
  mutate(spp_int = as.integer(as.factor(SpeciesAuthor)))



### plot pace vs. shape, point estimates vs posterior means
ggplot(df_shape) +
  geom_segment(aes(x = L_pt, y = S_pt, xend = L_mean, yend = S_mean),
               size = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_point(aes(L_pt, S_pt)) +
  scale_x_log10()




### model relationship between l0 and shape, assuming no sampling uncertainty

# compile stan models
stan_regress_hier <- stan_model("stan/regress2.stan")
stan_regress_hier_error <- stan_model("stan/regress_hier_error.stan")

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
  control = list(adapt_delta = 0.9, stepsize  = 0.1, max_treedepth = 12)
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
  labs(x = expression(paste("Mature life expectancy (", italic(L), ")")),
       y = expression(paste("Shape of mortality trajectory (", italic(S), ")"))) +
  tt

p2 <- ggplot(df_beta, aes(x = beta)) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  geom_density(fill = "darkred", alpha = 0.4, size = 0) +
  coord_cartesian(xlim = c(-0.07, 0.07)) +
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
# ggsave2("img/shape.png", p, height = 4.5, width = 6.25)





### posterior summary
median(mu_beta)
median(mu_beta_error)

var(mu_beta_error) / var(mu_beta)

length(mu_beta[mu_beta > 0]) / length(mu_beta)
length(mu_beta_error[mu_beta_error > 0]) / length(mu_beta_error)








### Variance components analysis
stan_varcomp <- stan_model("stan/varcomp.stan")



df_other <- sd_other_out %>% 
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
  left_join(pt_other_out) %>% 
  mutate(damp_pt = log10(damp_pt)) %>% 
  mutate(gen_pt = log10(gen_pt)) %>% 
  mutate(pmature_pt = logit(pmature_pt)) %>% 
  mutate(growth_pt = logit(growth_pt))





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



