# S16: simulation analysis of climate effects on fecundity under finite sampling.


# libraries and setup ----
source("code/setup.R")
setup_packages(c(
  "tidyverse", "here", "mpmsim", "popbio", "rstan", "patchwork"
))
setup_rstan()
source("code/functions.R")
# Plot style helpers (theme_mpm(), mpm_colors()) from code/functions.R
set_mpm_plot_defaults()

seed <- 5654
set.seed(seed)
cols <- mpm_colors()


# simulation settings ----
run_mode <- "production"

valid_modes <- c("quick", "intermediate", "production")
if (!run_mode %in% valid_modes) {
  stop("run_mode must be one of: ", paste(valid_modes, collapse = ", "), call. = FALSE)
}

if (run_mode == "quick") {
  beta_vals <- c(0.3)
  n_vals <- c(20)
  n_reps <- 5L
  stan_iter <- 1200L
  stan_warmup <- 600L
  stan_chains <- 2L
} else if (run_mode == "intermediate") {
  beta_vals <- c(0.1, 0.3, 0.5)
  n_vals <- c(10, 20, 50, 100)
  n_reps <- 20L
  stan_iter <- 1600L
  stan_warmup <- 800L
  stan_chains <- 2L
} else {
  beta_vals <- c(0.1, 0.3, 0.5)
  n_vals <- c(10, 20, 50, 100)
  n_reps <- 200L
  stan_iter <- 3000L
  stan_warmup <- 1500L
  stan_chains <- 4L
}

t_years <- 15L
sigma_process_true <- 0.2
target_lambda <- 1.0
lambda_tol <- 0.02
fecundity_base <- 1
stan_control <- if (run_mode == "production") ctrl3 else ctrl2

if (!dir.exists("data/derived/analysis_cache")) {
  dir.create("data/derived/analysis_cache", recursive = TRUE)
}

if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}


# baseline MPM ----
generate_baseline_lefko <- function(seed,
                                    target_lambda = 1,
                                    lambda_tol = 0.02,
                                    n_stages = 4L,
                                    archetype = 2L,
                                    fecundity = 1) {
  set.seed(seed)

  repeat {
    mats <- mpmsim::rand_lefko_mpm(
      n_stages = n_stages,
      fecundity = fecundity,
      archetype = archetype,
      split = TRUE
    )

    lam <- popbio::lambda(mats$mat_A)

    if (is.finite(lam) &&
      abs(lam - target_lambda) <= lambda_tol &&
      mats$mat_F[1, n_stages] > 0) {
      out_tbl <- bind_rows(
        as_tibble(as.data.frame(as.table(mats$mat_U))) %>%
          transmute(matrix = "U", to_stage = as.integer(Var1), from_stage = as.integer(Var2), value = Freq),
        as_tibble(as.data.frame(as.table(mats$mat_F))) %>%
          transmute(matrix = "F", to_stage = as.integer(Var1), from_stage = as.integer(Var2), value = Freq),
        as_tibble(as.data.frame(as.table(mats$mat_A))) %>%
          transmute(matrix = "A", to_stage = as.integer(Var1), from_stage = as.integer(Var2), value = Freq)
      )

      return(list(
        mat_u = mats$mat_U,
        mat_f = mats$mat_F,
        mat_a = mats$mat_A,
        lambda = lam,
        matrix_long = out_tbl
      ))
    }
  }
}

baseline <- generate_baseline_lefko(
  seed = seed,
  target_lambda = target_lambda,
  lambda_tol = lambda_tol,
  fecundity = fecundity_base
)

u_base <- baseline$mat_u
f_base <- baseline$mat_f
write_csv(
  baseline$matrix_long,
  "data/derived/analysis_cache/analysis3_sim_baseline_mpm.csv"
)


# design grid ----
design <- tidyr::expand_grid(
  beta_true = beta_vals,
  n_stage = n_vals,
  replicate = seq_len(n_reps)
) %>%
  mutate(
    sim_id = row_number(),
    seed = seed + sim_id
  )

write_csv(
  design,
  "data/derived/analysis_cache/analysis3_sim_design_grid.csv"
)


# helper functions ----
simulate_one_dataset <- function(beta_true, n_stage, seed, t_years,
                                 f_base, sigma_process_true) {
  set.seed(seed)

  x_t <- rnorm(t_years, mean = 0, sd = 1)
  z_t <- rnorm(t_years, mean = 0, sd = sigma_process_true)
  log_f_true <- log(f_base[1, 4]) + beta_true * x_t + z_t
  f_true_14 <- exp(log_f_true)
  recruits <- rpois(t_years, lambda = n_stage * f_true_14)
  fhat_14 <- recruits / n_stage

  tibble(
    year = seq_len(t_years),
    x_t = x_t,
    z_t = z_t,
    f_true_14 = f_true_14,
    recruits = recruits,
    n_repro = n_stage,
    fhat_14 = fhat_14
  )
}

fit_gaussian_point <- function(dat) {
  # Match the generating model's log-scale effect without logging observed zeros.
  fit <- glm(
    fhat_14 ~ x_t,
    family = gaussian(link = "log"),
    start = c(log(mean(dat$fhat_14)), 0),
    data = dat
  )
  if (!isTRUE(fit$converged)) stop("Gaussian log-link model did not converge.")
  est <- unname(coef(fit)[["x_t"]])
  se <- coef(summary(fit))["x_t", "Std. Error"]
  int <- est + qt(c(0.025, 0.975), df = df.residual(fit)) * se

  tibble(
    model = "Gaussian point estimate",
    beta_est = est,
    beta_lo = int[1],
    beta_hi = int[2],
    fit_ok = TRUE
  )
}

fit_poisson_point <- function(dat) {
  fit <- glm(
    recruits ~ x_t + offset(log(n_repro)),
    family = poisson(),
    data = dat
  )
  est <- unname(coef(fit)[["x_t"]])
  se <- coef(summary(fit))["x_t", "Std. Error"]

  tibble(
    model = "Poisson point estimate",
    beta_est = est,
    beta_lo = est + qnorm(0.025) * se,
    beta_hi = est + qnorm(0.975) * se,
    fit_ok = TRUE
  )
}

safe_fit <- function(expr, model_label) {
  warns <- character()
  result <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(result, "error")) {
    return(tibble(
      model = model_label,
      beta_est = NA_real_,
      beta_lo = NA_real_,
      beta_hi = NA_real_,
      fit_ok = FALSE,
      warning_text = paste(c(conditionMessage(result), warns), collapse = " | "),
      rhat_high = NA_real_,
      n_eff_low = NA_real_,
      mcse_high = NA_real_,
      n_diverg = NA_real_
    ))
  }

  if (!"warning_text" %in% names(result)) {
    result$warning_text <- if (length(warns) == 0) NA_character_ else paste(unique(warns), collapse = " | ")
  }

  return(result)
}

fit_poisson_uncertainty <- function(dat, stan_mod, seed) {
  stan_dat <- list(
    T = nrow(dat),
    recruits = dat$recruits,
    x = dat$x_t,
    n_repro = as.vector(dat$n_repro)
  )

  fit <- rstan::sampling(
    object = stan_mod,
    data = stan_dat,
    chains = stan_chains,
    iter = stan_iter,
    warmup = stan_warmup,
    seed = seed,
    refresh = 0,
    control = stan_control
  )

  beta_draws <- rstan::extract(fit, pars = "beta")$beta
  diag_tbl <- stan_diagnostics(fit)

  tibble(
    model = "Poisson uncertainty-aware",
    beta_est = median(beta_draws),
    beta_lo = quantile(beta_draws, 0.025),
    beta_hi = quantile(beta_draws, 0.975),
    fit_ok = TRUE,
    warning_text = NA_character_,
    rhat_high = diag_tbl$rhat_high,
    n_eff_low = diag_tbl$n_eff_low,
    mcse_high = diag_tbl$mcse_high,
    n_diverg = diag_tbl$n_diverg
  )
}

summarise_fit <- function(fit_tbl, beta_true, n_stage, replicate, sim_id) {
  fit_tbl %>%
    mutate(
      beta_true = beta_true,
      n_stage = n_stage,
      replicate = replicate,
      sim_id = sim_id,
      bias = beta_est - beta_true,
      sq_error = (beta_est - beta_true)^2,
      covered = beta_lo <= beta_true & beta_hi >= beta_true,
      ci_width = beta_hi - beta_lo
    )
}


# compile Stan model ----
stan_mod <- rstan::stan_model("models/poisson_climate_latent.stan")


# example dataset ----
example_dat <- simulate_one_dataset(
  beta_true = design$beta_true[[1]],
  n_stage = design$n_stage[[1]],
  seed = design$seed[[1]],
  t_years = t_years,
  f_base = f_base,
  sigma_process_true = sigma_process_true
)

write_csv(
  example_dat,
  "data/derived/analysis_cache/analysis3_sim_example_timeseries.csv"
)


# run simulation ----
results <- purrr::pmap_dfr(
  design,
  function(beta_true, n_stage, replicate, sim_id, seed) {
    dat <- simulate_one_dataset(
      beta_true = beta_true,
      n_stage = n_stage,
      seed = seed,
      t_years = t_years,
      f_base = f_base,
      sigma_process_true = sigma_process_true
    )

    fits <- bind_rows(
      safe_fit(fit_gaussian_point(dat), "Gaussian point estimate"),
      safe_fit(fit_poisson_point(dat), "Poisson point estimate"),
      safe_fit(fit_poisson_uncertainty(dat, stan_mod = stan_mod, seed = seed), "Poisson uncertainty-aware")
    )

    summarise_fit(
      fit_tbl = fits,
      beta_true = beta_true,
      n_stage = n_stage,
      replicate = replicate,
      sim_id = sim_id
    )
  }
)

write_csv(
  results,
  "data/derived/analysis_cache/analysis3_sim_replicate_summary.csv"
)

diagnostics <- results %>%
  group_by(model, beta_true, n_stage) %>%
  summarise(
    n_total = dplyr::n(),
    n_fit_ok = sum(fit_ok, na.rm = TRUE),
    n_warn = sum(!is.na(warning_text) & warning_text != "", na.rm = TRUE),
    n_rhat = sum(replace_na(rhat_high, 0) > 0),
    n_neff = sum(replace_na(n_eff_low, 0) > 0),
    n_mcse = sum(replace_na(mcse_high, 0) > 0),
    n_divergent = sum(replace_na(n_diverg, 0) > 0),
    .groups = "drop"
  )

write_csv(
  diagnostics,
  "data/derived/analysis_cache/analysis3_sim_diagnostics.csv"
)


# summary outputs ----
cell_summary <- results %>%
  group_by(model, beta_true, n_stage) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    median_bias = median(bias, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    coverage = mean(covered, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    n_ok = sum(fit_ok, na.rm = TRUE),
    n_warn = sum(!is.na(warning_text) & warning_text != "", na.rm = TRUE),
    n_divergent = sum(replace_na(n_diverg, 0) > 0),
    n_total = dplyr::n(),
    .groups = "drop"
  )

write_csv(
  cell_summary,
  "data/derived/analysis_cache/analysis3_sim_cell_summary.csv"
)

errorbar_summary <- bind_rows(
  results %>%
    group_by(model, beta_true, n_stage) %>%
    summarise(
      metric = "Bias",
      value = mean(bias, na.rm = TRUE),
      se = stats::sd(bias, na.rm = TRUE) / sqrt(sum(!is.na(bias))),
      ymin = value - stats::qt(0.975, df = sum(!is.na(bias)) - 1) * se,
      ymax = value + stats::qt(0.975, df = sum(!is.na(bias)) - 1) * se,
      .groups = "drop"
    ),
  results %>%
    group_by(model, beta_true, n_stage) %>%
    summarise(
      metric = "Coverage",
      value = mean(covered, na.rm = TRUE),
      n_cover = sum(covered, na.rm = TRUE),
      n_total = sum(!is.na(covered)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      ci = list(stats::binom.test(n_cover, n_total)$conf.int),
      ymin = ci[[1]],
      ymax = ci[[2]]
    ) %>%
    ungroup() %>%
    select(-ci, -n_cover, -n_total)
)

write_csv(
  errorbar_summary,
  "data/derived/analysis_cache/analysis3_sim_errorbar_summary.csv"
)


# Figure 6 simulation performance summary ----
model_cols <- c(
  "Gaussian point estimate" = cols$point,
  "Poisson point estimate" = "#2C7FB8",
  "Poisson uncertainty-aware" = cols$sampling
)

plot_dat <- bind_rows(
  cell_summary %>%
    transmute(model, beta_true, n_stage, metric = "Bias", value = mean_bias, ref = 0),
  cell_summary %>%
    transmute(model, beta_true, n_stage, metric = "Coverage", value = coverage, ref = 0.95)
) %>%
  mutate(
    metric = factor(metric, levels = c("Bias", "Coverage")),
    beta_lab = factor(
      paste0("beta == ", beta_true),
      levels = paste0("beta == ", sort(unique(beta_true)))
    ),
    n_stage_f = factor(n_stage, levels = sort(unique(n_stage)))
  )

ref_dat <- plot_dat %>%
  distinct(metric, beta_lab, ref) %>%
  filter(!is.na(ref))

p_fig6 <- ggplot(plot_dat, aes(x = n_stage_f, y = value, color = model, group = model)) +
  geom_hline(
    data = ref_dat,
    aes(yintercept = ref),
    inherit.aes = FALSE,
    linetype = 2,
    color = "grey70",
    linewidth = 0.35
  ) +
  geom_point() +
  geom_line(linewidth = 0.5) +
  facet_grid(
    metric ~ beta_lab,
    scales = "free_y",
    switch = "y",
    labeller = labeller(beta_lab = label_parsed, metric = label_value)
  ) +
  scale_color_manual(values = model_cols) +
  labs(
    x = "Sample size per stage per year",
    y = NULL,
    color = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(size = 9.5),
    strip.placement = "outside",
    axis.text.x = element_text(size = 8.5),
    panel.spacing = unit(4, "mm"),
    plot.margin = margin(6, 8, 4, 6)
  )

ggsave(
  "figures/Figure_6_analysis3_simulation_performance.png",
  p_fig6,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 300
)

plot_dat_err <- errorbar_summary %>%
  mutate(
    metric = factor(metric, levels = c("Bias", "Coverage")),
    beta_lab = factor(
      paste0("beta == ", beta_true),
      levels = paste0("beta == ", sort(unique(beta_true)))
    ),
    n_stage_f = factor(n_stage, levels = sort(unique(n_stage)))
  )

ref_dat_err <- plot_dat_err %>%
  mutate(ref = ifelse(metric == "Bias", 0, 0.95)) %>%
  distinct(metric, beta_lab, ref)

p_fig6_err <- ggplot(plot_dat_err, aes(x = n_stage_f, y = value, color = model, group = model)) +
  geom_hline(
    data = ref_dat_err,
    aes(yintercept = ref),
    inherit.aes = FALSE,
    linetype = 2,
    color = "grey70",
    linewidth = 0.35
  ) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.08, linewidth = 0.35) +
  geom_point() +
  geom_line(linewidth = 0.5) +
  facet_grid(
    metric ~ beta_lab,
    scales = "free_y",
    switch = "y",
    labeller = labeller(beta_lab = label_parsed, metric = label_value)
  ) +
  scale_color_manual(values = model_cols) +
  labs(
    x = "Sample size per stage per year",
    y = NULL,
    color = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(size = 9.5),
    strip.placement = "outside",
    axis.text.x = element_text(size = 8.5),
    panel.spacing = unit(4, "mm"),
    plot.margin = margin(6, 8, 4, 6)
  )

ggsave(
  "figures/Figure_6_analysis3_simulation_performance_errorbars.png",
  p_fig6_err,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 300
)
