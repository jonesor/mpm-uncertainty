# S12: case study 2 climate analysis and figures (uses PRISM inputs from S11).

# libraries ----
source("code/setup.R")
setup_packages(c(
  "tidyverse", "cowplot", "Rcompadre", "Rage", "popbio",
  "gridExtra", "rstan", "loo", "patchwork"
))
setup_rstan()
source("code/functions.R")
# Plot style helpers (theme_mpm(), mpm_colors()) from code/functions.R
set_mpm_plot_defaults()
seed <- 5654
set.seed(seed)
cols <- mpm_colors()


# set options for rstan library ----
# handled in setup_rstan()


# load compadre data ----
compadre <- cdb_fetch("data/raw/compadre/COMPADRE_v.X.X.X_Corrected.RData")


# find long time-series for climate analysis ----
# first subset to wild, unmanipulated populations with 1-yr periodicity
comp_sub <- compadre %>%
  filter(
    MatrixComposite == "Individual",
    MatrixTreatment == "Unmanipulated",
    ProjectionInterval == "1",
    MatrixCaptivity == "W"
  )

# find populations with time-series >= 5 years
comp_time_series <- comp_sub %>%
  as_tibble() %>%
  filter(!is.na(Lon) & !is.na(Lat)) %>%
  group_by(SpeciesAuthor) %>%
  mutate(n_year = length(unique(MatrixStartYear))) %>%
  ungroup() %>%
  filter(n_year >= 5) %>%
  group_by(SpeciesAuthor, MatrixPopulation) %>%
  # note Eryngium_alpinum MatrixPopulation "PRC" has two coords
  summarize(
    Lon = unique(Lon)[1],
    Lat = unique(Lat)[1],
    n_year = unique(n_year)
  ) %>%
  ungroup()

# write lat/lon to file for climate analysis
# write.csv(comp_time_series, "species_coords.csv", row.names = FALSE)


# Model climate for Silene_spaldingii ----
# load Ellis et al (2012) data
ellis_data <- read.table("data/raw/ellis_2012/Transition_Matrices.txt",
  sep = "\t",
  header = TRUE, stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(matA = map(Mx, string_to_mat)) %>%
  mutate(matU = map(Tmx, string_to_mat)) %>%
  mutate(matF = map2(matA, matU, ~ .x - .y)) %>%
  mutate(N = map(Nx, nx_to_vec))

# stage-specific sample sizes for Silene
silene_n <- ellis_data %>%
  filter(SPP == "SISP") %>%
  mutate(
    MatrixStartYear = YR,
    SpeciesAuthor = "Silene_spaldingii"
  ) %>%
  dplyr::select(MatrixStartYear, SpeciesAuthor, N)

# load weather data
wx <- read_csv("data/derived/climate/species_clim_prism.csv") %>%
  filter(SpeciesAuthor == "Silene_spaldingii")

wx_spring <- wx %>%
  filter(Month %in% 2:4) %>%
  group_by(Year) %>%
  summarize(tmp = mean(tmp), ppt = sum(ppt)) %>%
  mutate(MatrixStartYear = Year)

# derive fecundity transition, and underlying counts
silene <- comp_sub %>%
  filter(SpeciesAuthor == "Silene_spaldingii") %>%
  left_join(silene_n, by = c("MatrixStartYear", "SpeciesAuthor")) %>%
  cdb_unnest() %>%
  as_tibble() %>%
  mutate(fecund = map_dbl(matF, sum)) %>%
  mutate(n_repro = map_int(N, ~ .x[3])) %>%
  mutate(n_offsp = round(n_repro * fecund, 0)) %>%
  mutate(fec_low = qgamma(0.025, 1 + n_offsp, n_repro)) %>%
  mutate(fec_upp = qgamma(0.975, 1 + n_offsp, n_repro)) %>%
  mutate(fec_sim = map2(n_offsp, n_repro, ~ rgamma(2000, 1 + .x, .y))) %>%
  left_join(wx_spring, by = "MatrixStartYear") %>%
  mutate(
    tmp = as.numeric(scale(tmp)),
    ppt = as.numeric(scale(ppt))
  )


# plot Spring temp vs. fecundity ----
p_scatter <- ggplot(silene, aes(tmp, fecund)) +
  geom_point(color = cols$dark) +
  geom_linerange(aes(ymin = fec_low, ymax = fec_upp), color = cols$dark) +
  geom_smooth(method = "lm", color = cols$accent, fill = cols$accent, alpha = 0.25) +
  scale_y_log10()

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
ggsave("figures/Analysis2_spring_temperature_recruitment_scatter.png", p_scatter, height = 4, width = 4.5, units = "in", dpi = 300)


# Simple regression of fecundity vs. spring temperature ----

# compile stan models
stan_regress <- stan_model("models/regress.stan")
stan_regress_err <- stan_model("models/regress_log_err.stan")

# sequence of x values for prediction line
xpred <- seq(min(silene$tmp), max(silene$tmp), length.out = 100)

# arrange data for stan
dat_reg <- dat_raw <- list(
  N = nrow(silene),
  P = length(xpred),
  x = silene$tmp,
  xpred = xpred,
  y = log(silene$fecund)
)

dat_reg$y <- log(silene$fecund)
dat_raw$y <- silene$n_offsp
dat_raw$offset <- silene$n_repro


# fit stan models
fit_reg <- stanfn(stan_regress, data = dat_reg, iter = 5000, seed = seed)
fit_err <- stanfn(stan_regress_err, data = dat_raw, iter = 5000, seed = seed)


# extract posterior samples for intercept and slope
beta_reg <- rstan::extract(fit_reg, "beta")$beta
beta_err <- rstan::extract(fit_err, "beta")$beta

# compare posterior median beta
median(beta_reg)
median(beta_err)
(median(beta_err) - median(beta_reg)) / median(beta_reg)

# compare uncertainty in beta
var_beta_reg <- var(beta_reg)
var_beta_err <- var(beta_err)
var_beta_err / var_beta_reg # var(beta_err) is ~12% higher

lev <- c("Point estimates", "Sampling uncertainty")
model_cols <- c(
  "Point estimates" = cols$light,
  "Sampling uncertainty" = cols$dark
)
model_alpha <- c(
  "Point estimates" = 1.0,
  "Sampling uncertainty" = 0.72
)

# plot beta by model type
df_beta <- tibble(reg = beta_reg, err = beta_err) %>%
  gather(model, val) %>%
  mutate(model = recode(model,
    reg = "Point estimates",
    err = "Sampling uncertainty"
  )) %>%
  mutate(model = factor(model, levels = lev)) %>%
  group_by(model) %>%
  summarize(
    med = quantile(val, 0.500),
    low80 = quantile(val, 0.10),
    upp80 = quantile(val, 0.90),
    low95 = quantile(val, 0.025),
    upp95 = quantile(val, 0.975)
  )

p_beta <- ggplot(df_beta, aes(x = model)) +
  geom_point(aes(y = med, color = model), size = 2.5) +
  geom_linerange(aes(ymin = low80, ymax = upp80, color = model), linewidth = 1.5) +
  geom_linerange(aes(ymin = low95, ymax = upp95, color = model)) +
  geom_hline(yintercept = 0, alpha = 0.5, linetype = 2) +
  coord_flip() +
  scale_color_manual(values = model_cols) +
  guides(color = "none")

ggsave("figures/Analysis2_spring_temperature_beta_summary.png", p_beta, height = 3.5, width = 4, units = "in", dpi = 300)

if (!dir.exists("data/derived/analysis_cache")) {
  dir.create("data/derived/analysis_cache",
    recursive = TRUE
  )
}
beta_summary <- df_beta %>%
  select(model, med, low95, upp95)
write_csv(
  beta_summary %>% mutate(across(where(is.numeric), ~ signif(.x, 3))),
  "data/derived/analysis_cache/case2_spring_beta_summary.csv"
)

# plot fit lines ----
# uses lev/model_cols defined above

# fit line
pred_full <- rbind(
  mutate(posterior_vec(fit_reg, xpred, "pred"), model = lev[1]),
  mutate(posterior_vec(fit_err, xpred, "pred"), model = lev[2])
) %>% mutate(model = factor(model, levels = lev))

# points and error bars
year_err <- silene %>%
  dplyr::select(fecund, fec_low, fec_upp) %>%
  mutate(x = silene$tmp) %>%
  mutate(model = lev[2])

year_full <- year_err %>%
  mutate(model = lev[1], fec_low = NA_real_, fec_upp = NA_real_) %>%
  rbind(year_err) %>%
  mutate(model = factor(model, levels = lev))

# plot
tt <- theme_mpm() +
  theme(
    text = element_text(size = 11.5),
    axis.ticks = element_line(linewidth = 0.4),
    plot.margin = margin(2, 2, 2, 2)
  )

p1 <- ggplot(pred_full, aes(x = x)) +
  geom_ribbon(aes(ymin = low95, ymax = upp95, fill = model), alpha = 0.20, color = NA) +
  geom_line(aes(y = med, color = model, alpha = model), linewidth = 0.8) +
  geom_point(
    data = year_full,
    aes(y = fecund, color = model, group = model),
    size = 1,
    alpha = 0.85,
    position = position_dodge(width = 0.07)
  ) +
  geom_linerange(
    data = year_full,
    aes(ymin = fec_low, ymax = fec_upp, color = model, group = model),
    alpha = 0.85,
    position = position_dodge(width = 0.07)
  ) +
  scale_y_log10(breaks = 10^(-2:0), labels = c("0.01", "0.1", "1")) +
  scale_color_manual(values = model_cols) +
  scale_fill_manual(values = model_cols) +
  scale_alpha_manual(values = model_alpha) +
  labs(x = "Spring temperature (Feb-Apr)", y = "Recruitment") +
  tt +
  guides(color = "none", fill = "none", alpha = "none")

ggsave("figures/Analysis2_recruitment_model_fits.png", p1, height = 4.5, width = 3.5, units = "in", dpi = 300)


# Stage-specific survival vs spring temperature ----
stage_surv <- silene %>%
  transmute(MatrixStartYear, tmp, N, matU) %>%
  mutate(
    stage_tbl = map2(N, matU, ~ {
      surv_prob <- colSums(.y)
      ssd <- tryCatch(
        {
          w <- popbio::stable.stage(.y)
          as.numeric(w / sum(w))
        },
        error = function(e) rep(NA_real_, length(.x))
      )
      tibble(
        stage = seq_along(.x),
        n_stage = as.integer(round(.x)),
        n_surv = as.integer(round(.x * surv_prob)),
        ssd_w = pmax(ssd, 1e-6)
      )
    })
  ) %>%
  select(-N, -matU) %>%
  unnest(stage_tbl) %>%
  mutate(
    n_fail = pmax(n_stage - n_surv, 0L),
    surv = n_surv / pmax(n_stage, 1L)
  ) %>%
  filter(n_stage > 0)

stage_beta <- stage_surv %>%
  group_by(stage) %>%
  group_modify(~ {
    if (nrow(.x) < 5 || n_distinct(.x$tmp) < 2 || sum(.x$n_surv) == 0 || sum(.x$n_fail) == 0) {
      return(tibble(
        model = c("Point estimates", "Sampling uncertainty"),
        beta = NA_real_,
        low95 = NA_real_,
        upp95 = NA_real_,
        n_year = nrow(.x)
      ))
    }

    surv_clamped <- pmin(pmax(.x$surv, 1e-6), 1 - 1e-6)
    fit_pt <- glm(surv_clamped ~ tmp, weights = ssd_w, family = quasibinomial(), data = .x)
    beta_pt <- unname(stats::coef(fit_pt)["tmp"])
    se_pt <- sqrt(stats::vcov(fit_pt)["tmp", "tmp"])

    fit_bin <- glm(cbind(n_surv, n_fail) ~ tmp, family = binomial(), data = .x)
    beta_bin <- unname(stats::coef(fit_bin)["tmp"])
    se_bin <- sqrt(stats::vcov(fit_bin)["tmp", "tmp"])

    tibble(
      model = c("Point estimates", "Sampling uncertainty"),
      beta = c(beta_pt, beta_bin),
      low95 = c(beta_pt - 1.96 * se_pt, beta_bin - 1.96 * se_bin),
      upp95 = c(beta_pt + 1.96 * se_pt, beta_bin + 1.96 * se_bin),
      n_year = nrow(.x)
    )
  }) %>%
  ungroup() %>%
  mutate(
    model = factor(model, levels = lev),
    stage_label = paste0("Stage ", stage)
  )

stage_weight_sensitivity <- stage_surv %>%
  group_by(stage) %>%
  group_modify(~ {
    if (nrow(.x) < 5 || n_distinct(.x$tmp) < 2 || sum(.x$n_surv) == 0 || sum(.x$n_fail) == 0) {
      return(tibble(
        weighting = c("No weights", "SSD weights"),
        beta = NA_real_,
        low95 = NA_real_,
        upp95 = NA_real_,
        n_year = nrow(.x)
      ))
    }
    surv_clamped <- pmin(pmax(.x$surv, 1e-6), 1 - 1e-6)
    fit_unw <- glm(surv_clamped ~ tmp, family = quasibinomial(), data = .x)
    fit_ssd <- glm(surv_clamped ~ tmp, weights = ssd_w, family = quasibinomial(), data = .x)
    beta_unw <- unname(stats::coef(fit_unw)["tmp"])
    se_unw <- sqrt(stats::vcov(fit_unw)["tmp", "tmp"])
    beta_ssd <- unname(stats::coef(fit_ssd)["tmp"])
    se_ssd <- sqrt(stats::vcov(fit_ssd)["tmp", "tmp"])
    tibble(
      weighting = c("No weights", "SSD weights"),
      beta = c(beta_unw, beta_ssd),
      low95 = c(beta_unw - 1.96 * se_unw, beta_ssd - 1.96 * se_ssd),
      upp95 = c(beta_unw + 1.96 * se_unw, beta_ssd + 1.96 * se_ssd),
      n_year = nrow(.x)
    )
  }) %>%
  ungroup()
write_csv(
  stage_weight_sensitivity %>% mutate(across(where(is.numeric), ~ signif(.x, 3))),
  "data/derived/analysis_cache/case2_stage_survival_weight_sensitivity.csv"
)

stage_levels <- stage_surv %>%
  distinct(stage) %>%
  arrange(stage) %>%
  mutate(stage_label = paste0("Stage ", stage)) %>%
  pull(stage_label)

# For non-estimable stages, show descriptive markers only (no uncertainty bars).
stage_placeholders <- stage_surv %>%
  group_by(stage) %>%
  summarize(
    mean_surv = mean(surv, na.rm = TRUE),
    all_survive = all(n_fail == 0),
    .groups = "drop"
  ) %>%
  mutate(
    stage_label = paste0("Stage ", stage),
    x = case_when(
      stage == 1 ~ mean_surv,
      stage == 4 & all_survive ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(x))

p3 <- ggplot(stage_beta %>% filter(!is.na(beta)),
  aes(x = beta, y = stage_label, color = model)
) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5, color = cols$accent) +
  geom_linerange(
    aes(xmin = low95, xmax = upp95),
    position = position_dodge(width = 0.6),
    linewidth = 0.7
  ) +
  geom_point(position = position_dodge(width = 0.6), size = 1.6) +
  geom_point(
    data = stage_placeholders,
    aes(x = x, y = stage_label),
    inherit.aes = FALSE,
    shape = 1,
    size = 2.2,
    stroke = 0.8,
    color = cols$light
  ) +
  scale_color_manual(values = model_cols, drop = FALSE) +
  scale_y_discrete(limits = rev(stage_levels), drop = FALSE) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(
    x = expression(paste("Temperature coefficient (logit scale, ", italic(beta), ")")),
    y = "Stage-specific survival"
  ) +
  tt +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.key.width = unit(10, "pt")
  )

ggsave("figures/Analysis2_stage_specific_survival_beta_summary.png", p3, height = 4.5, width = 3.5, units = "in", dpi = 300)


# Moving beta model ----

# compile stan models
mod_null_reg <- stan_model("models/null.stan")
mod_null_err <- stan_model("models/null_err.stan")
mod_gprc_reg <- stan_model("models/movbeta_gprc.stan")
mod_gprc_err <- stan_model("models/movbeta_gprc_err.stan")

# focal years for Silene series
focal_yrs <- seq(
  min(silene_n$MatrixStartYear) - 1,
  max(silene_n$MatrixStartYear) + 1
)

# arrange weather data for moving-beta model
wx_mb <- wx %>%
  filter(Year %in% focal_yrs) %>%
  group_by(Month) %>%
  mutate(tmp = as.numeric(scale(tmp)), ppt = as.numeric(scale(ppt))) %>%
  ungroup() %>%
  mutate(date = as.Date(paste(Year, Month, "01", sep = "-")))


# prepare data for stan ----
year <- silene$MatrixEndYear
y <- log(silene$fecund)
N <- length(y)
K <- 24
month_start <- "07"


# assemble matrix of climate data ----
X <- map(seq_along(y), ~ {
  yr_focal <- year[.x]
  date_origin <- as.Date(paste(yr_focal, month_start, "01", sep = "-"))
  dates_focal <- sort(seq(date_origin, by = "-1 month", length.out = K))
  filter(wx_mb, date %in% dates_focal)$tmp
}) %>%
  map(~ tibble_row(!!!setNames(as.list(.x), paste0("m", seq_len(K))))) %>%
  list_rbind() %>%
  as.matrix()


# arrange data for stan ----
dat_reg <- list(N = N, K = K, X = X, y = y)
dat_err <- list(N = N, K = K, X = X, y = silene$n_offsp, offset = silene$n_repro)


# fit stan models ----
fit_null_reg <- stanfn(mod_null_reg, data = dat_reg, seed = seed)
fit_null_err <- stanfn(mod_null_err, data = dat_err, seed = seed)
fit_gprc_reg <- stanfn(mod_gprc_reg, data = dat_reg, control = ctrl2, seed = seed)
fit_gprc_err <- stanfn(mod_gprc_err, data = dat_err, control = ctrl2, seed = seed)


# diagnostics ----


# measures of fit ----
rbind(
  summarize_fit(fit_null_reg, "null"),
  summarize_fit(fit_gprc_reg, "moving-beta"),
  summarize_fit(fit_null_err, "null (err)"),
  summarize_fit(fit_gprc_err, "moving-beta (err)")
)


# plot moving-beta coefficients ----
lev <- c("Point estimates", "Sampling uncertainty")

gprc_betas <- rbind(
  summarize_beta(fit_gprc_reg, "Point estimates"),
  summarize_beta(fit_gprc_err, "Sampling uncertainty")
) %>%
  mutate(model = factor(model, levels = lev)) %>%
  mutate(
    lag_err = lag + if_else(model == "Point estimates", -0.18, 0.18)
  )


p2 <- ggplot(gprc_betas, aes(x = lag)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_line(aes(y = beta_med, color = model, alpha = model), linewidth = 0.8) +
  geom_linerange(
    aes(x = lag_err, ymin = beta_low95, ymax = beta_upp95, color = model),
    alpha = 0.85
  ) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  scale_y_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  scale_color_manual(values = model_cols) +
  scale_alpha_manual(values = model_alpha) +
  labs(
    x = "Months before survey",
    y = expression(paste("Temperature coefficient (", italic(b), ")"))
  ) +
  tt +
  theme(plot.margin = margin(2, 2, 2, 2)) +
  guides(color = "none", alpha = "none")

ggsave("figures/Analysis2_monthly_lag_coefficients.png", p2, height = 4.5, width = 3.5, units = "in", dpi = 300)

write_csv(
  gprc_betas %>% mutate(across(where(is.numeric), ~ signif(.x, 3))),
  "data/derived/analysis_cache/case2_gprc_beta_summary.csv"
)


# combine both climate plots
p <- p2 / p1 / p3 +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.margin = margin(2, 2, 2, 2))
  )
ggsave("figures/Figure_5_analysis2_climate_effects_recruitment.png", p, height = 195, width = 90, units = "mm", dpi = 300)
