# S14: analysis 3 climate analysis for Astragalus scaphoides across sites.


# libraries ----
source("code/setup.R")
setup_packages(c(
  "tidyverse", "Rcompadre", "rstan", "patchwork", "popdemo"
))
setup_rstan()
source("code/functions.R")
# Plot style helpers (theme_mpm(), mpm_colors()) from code/functions.R
set_mpm_plot_defaults()
seed <- 5654
set.seed(seed)
cols <- mpm_colors()


# analysis settings ----
site_key <- tibble(
  SpeciesAuthor = "Astragalus_scaphoides_2",
  MatrixPopulation = c("Haynes Creek", "Sheep Corral Gulch", "McDevitt Creek"),
  ellis_spp = "ASSC",
  ellis_pop = c("ASSC_haynes", "ASSC_sheep", "ASSC_mcdevi")
)

site_component_panels <- list()


# load data ----
compadre <- load_compadre(corrected = TRUE) %>%
  mutate(
    MatrixStartYear = suppressWarnings(as.integer(MatrixStartYear)),
    MatrixEndYear = suppressWarnings(as.integer(MatrixEndYear))
  )
ellis_data <- read.table(
  "data/raw/ellis_2012/Transition_Matrices.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(matA = map(Mx, string_to_mat)) %>%
  mutate(matU = map(Tmx, string_to_mat)) %>%
  mutate(matF = map2(matA, matU, ~ .x - .y)) %>%
  mutate(N = map(Nx, nx_to_vec))

wx <- read_csv("data/derived/climate/species_clim_prism.csv")


# check that PRISM extraction exists for all target sites ----
clim_sites <- wx %>%
  distinct(SpeciesAuthor, MatrixPopulation)
missing_sites <- anti_join(site_key, clim_sites, by = c("SpeciesAuthor", "MatrixPopulation"))
if (nrow(missing_sites) > 0) {
  stop(
    paste0(
      "Missing PRISM climate rows for: ",
      paste(paste(missing_sites$SpeciesAuthor, missing_sites$MatrixPopulation, sep = " / "), collapse = "; "),
      ". Run scripts/S11_case_study_2_prism.R first."
    )
  )
}

if (!dir.exists("data/derived/analysis_cache")) {
  dir.create("data/derived/analysis_cache", recursive = TRUE)
}
if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

# compile stan models once ----
stan_regress <- stan_model("models/regress.stan")
stan_regress_err <- stan_model("models/regress_log_err.stan")
mod_null_reg <- stan_model("models/null.stan")
mod_null_err <- stan_model("models/null_err.stan")
mod_gprc_reg <- stan_model("models/movbeta_gprc.stan")
mod_gprc_err <- stan_model("models/movbeta_gprc_err.stan")


# helpers ----
slugify_site <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

analyze_site <- function(SpeciesAuthor, MatrixPopulation, ellis_spp, ellis_pop) {
  spp_target <- SpeciesAuthor
  pop_target <- MatrixPopulation
  site_slug <- slugify_site(pop_target)

  comp_sub <- compadre %>%
    filter(
      MatrixComposite == "Individual",
      MatrixTreatment == "Unmanipulated",
      ProjectionInterval == "1",
      MatrixCaptivity == "W",
      SpeciesAuthor == spp_target,
      MatrixPopulation == pop_target
    )

  astr_n <- ellis_data %>%
    filter(SPP == ellis_spp, POP == ellis_pop) %>%
    transmute(MatrixStartYear = YR, SpeciesAuthor = spp_target, N)

  wx_spring <- wx %>%
    filter(SpeciesAuthor == spp_target, MatrixPopulation == pop_target) %>%
    filter(Month %in% 2:4) %>%
    group_by(Year) %>%
    summarize(tmp = mean(tmp), ppt = sum(ppt), .groups = "drop") %>%
    mutate(MatrixStartYear = Year)

  astr <- comp_sub %>%
    left_join(astr_n, by = c("MatrixStartYear", "SpeciesAuthor")) %>%
    cdb_unnest() %>%
    as_tibble() %>%
    mutate(rep_idx = map(matF, ~ which(colSums(.x) > 0))) %>%
    mutate(n_repro = map2_int(N, rep_idx, ~ sum(.x[.y]))) %>%
    mutate(fecund = map_dbl(matF, sum)) %>%
    mutate(n_offsp = round(n_repro * fecund, 0)) %>%
    mutate(fec_low = qgamma(0.025, 1 + n_offsp, pmax(n_repro, 1))) %>%
    mutate(fec_upp = qgamma(0.975, 1 + n_offsp, pmax(n_repro, 1))) %>%
    mutate(repro_surv = map2_dbl(matU, rep_idx, ~ sum(colSums(.x)[.y]))) %>%
    mutate(n_surv = round(n_repro * repro_surv, 0)) %>%
    mutate(surv = n_surv / pmax(n_repro, 1)) %>%
    mutate(surv_low = qbeta(0.025, 1 + n_surv, 1 + pmax(n_repro - n_surv, 0))) %>%
    mutate(surv_upp = qbeta(0.975, 1 + n_surv, 1 + pmax(n_repro - n_surv, 0))) %>%
    mutate(ergodic = map_lgl(matA, ~ tryCatch(popdemo::isErgodic(.x), error = function(e) FALSE))) %>%
    left_join(wx_spring, by = "MatrixStartYear") %>%
    mutate(
      tmp = as.numeric(scale(tmp)),
      ppt = as.numeric(scale(ppt))
    ) %>%
    arrange(MatrixStartYear)

  astr_fit <- astr %>%
    filter(n_repro > 0) %>%
    filter(!is.na(tmp), !is.na(n_offsp))

  coverage_tbl <- astr %>%
    transmute(
      SpeciesAuthor,
      MatrixPopulation,
      MatrixStartYear,
      ergodic,
      n_repro,
      n_offsp,
      n_surv,
      include_in_fit = n_repro > 0 & !is.na(tmp) & !is.na(n_offsp)
    )

  write_csv(
    coverage_tbl,
    paste0("data/derived/analysis_cache/case3_astragalus_coverage_", site_slug, ".csv")
  )

  if (nrow(astr_fit) < 5) {
    return(tibble(
      SpeciesAuthor = spp_target,
      MatrixPopulation = pop_target,
      site_slug = site_slug,
      years_total = nrow(coverage_tbl),
      years_fit = nrow(astr_fit),
      years_excluded = nrow(coverage_tbl) - nrow(astr_fit),
      beta_med_point = NA_real_,
      beta_med_sampling = NA_real_,
      status = "skipped_too_few_years"
    ))
  }

  p_scatter <- ggplot(astr_fit, aes(tmp, fecund)) +
    geom_point(color = cols$dark) +
    geom_linerange(aes(ymin = fec_low, ymax = fec_upp), color = cols$dark) +
    geom_smooth(method = "lm", color = cols$accent, fill = cols$accent, alpha = 0.25) +
    scale_y_log10() +
    labs(
      x = "Spring temperature (Feb-Apr; z-scored)",
      y = "Recruitment",
      title = pop_target
    ) +
    theme_mpm()


  xpred <- seq(min(astr_fit$tmp), max(astr_fit$tmp), length.out = 100)

  dat_reg <- list(
    N = nrow(astr_fit),
    P = length(xpred),
    x = astr_fit$tmp,
    xpred = xpred,
    y = log(astr_fit$fecund)
  )

  dat_raw <- list(
    N = nrow(astr_fit),
    P = length(xpred),
    x = astr_fit$tmp,
    xpred = xpred,
    y = astr_fit$n_offsp,
    offset = astr_fit$n_repro
  )

  fit_reg <- stanfn(stan_regress, data = dat_reg, control = ctrl3, iter = 5000, seed = seed)
  fit_err <- stanfn(stan_regress_err, data = dat_raw, control = ctrl3, iter = 5000, seed = seed)

  dat_surv_reg <- list(
    N = nrow(astr_fit),
    P = length(xpred),
    x = astr_fit$tmp,
    xpred = xpred,
    y = log(pmax(astr_fit$surv, 1e-6))
  )

  dat_surv_raw <- list(
    N = nrow(astr_fit),
    P = length(xpred),
    x = astr_fit$tmp,
    xpred = xpred,
    y = astr_fit$n_surv,
    offset = astr_fit$n_repro
  )

  fit_surv_reg <- stanfn(stan_regress, data = dat_surv_reg, control = ctrl3, iter = 5000, seed = seed)
  fit_surv_err <- stanfn(stan_regress_err, data = dat_surv_raw, control = ctrl3, iter = 5000, seed = seed)

  diag_tbl <- bind_rows(
    stan_diagnostics(fit_reg) %>% mutate(model = "point_estimate"),
    stan_diagnostics(fit_err) %>% mutate(model = "sampling_uncertainty"),
    stan_diagnostics(fit_surv_reg) %>% mutate(model = "point_estimate_survival"),
    stan_diagnostics(fit_surv_err) %>% mutate(model = "sampling_uncertainty_survival")
  ) %>%
    mutate(SpeciesAuthor = spp_target, MatrixPopulation = pop_target) %>%
    relocate(SpeciesAuthor, MatrixPopulation, model)

  write_csv(
    diag_tbl,
    paste0("data/derived/analysis_cache/case3_astragalus_stan_diagnostics_", site_slug, ".csv")
  )

  beta_reg <- rstan::extract(fit_reg, "beta")$beta
  beta_err <- rstan::extract(fit_err, "beta")$beta

  lev <- c("Model of point estimates", "Model with sampling uncertainty")
  model_cols <- c(
    "Model of point estimates" = cols$point,
    "Model with sampling uncertainty" = cols$sampling
  )
  model_alpha <- c(
    "Model of point estimates" = 1.0,
    "Model with sampling uncertainty" = 0.72
  )

  df_beta <- tibble(reg = beta_reg, err = beta_err) %>%
    pivot_longer(cols = everything(), names_to = "model", values_to = "val") %>%
    mutate(model = recode(model,
      reg = "Model of point estimates",
      err = "Model with sampling uncertainty"
    )) %>%
    mutate(model = factor(model, levels = lev)) %>%
    group_by(model) %>%
    summarize(
      med = quantile(val, 0.500),
      low80 = quantile(val, 0.10),
      upp80 = quantile(val, 0.90),
      low95 = quantile(val, 0.025),
      upp95 = quantile(val, 0.975),
      .groups = "drop"
    )

  p_beta <- ggplot(df_beta, aes(x = model)) +
    geom_point(aes(y = med, color = model), size = 2.5) +
    geom_linerange(aes(ymin = low80, ymax = upp80, color = model), linewidth = 1.5) +
    geom_linerange(aes(ymin = low95, ymax = upp95, color = model)) +
    geom_hline(yintercept = 0, alpha = 0.5, linetype = 2) +
    coord_flip() +
    scale_color_manual(values = model_cols) +
    labs(
      x = NULL,
      y = expression(paste("Spring-temperature coefficient (", italic(beta), ")")),
      title = pop_target
    ) +
    theme_mpm() +
    guides(color = "none")


  write_csv(
    df_beta %>%
      mutate(
        SpeciesAuthor = spp_target,
        MatrixPopulation = pop_target,
        site_slug = site_slug,
        across(where(is.numeric), ~ signif(.x, 3))
      ) %>%
      relocate(SpeciesAuthor, MatrixPopulation, site_slug, model),
    paste0("data/derived/analysis_cache/case3_astragalus_spring_beta_summary_", site_slug, ".csv")
  )

  tt <- theme_mpm() +
    theme(
      text = element_text(size = 11.5),
      axis.ticks = element_line(linewidth = 0.4),
      plot.margin = margin(2, 2, 2, 2)
    )

  pred_full <- bind_rows(
    mutate(posterior_vec(fit_reg, xpred, "pred"), model = lev[1]),
    mutate(posterior_vec(fit_err, xpred, "pred"), model = lev[2])
  ) %>% mutate(model = factor(model, levels = lev))

  year_err <- astr_fit %>%
    select(fecund, fec_low, fec_upp, tmp) %>%
    rename(x = tmp) %>%
    mutate(model = lev[2])

  year_full <- bind_rows(
    year_err %>% mutate(model = lev[1], fec_low = NA_real_, fec_upp = NA_real_),
    year_err
  ) %>% mutate(model = factor(model, levels = lev))

  p_fit <- ggplot(pred_full, aes(x = x)) +
    geom_ribbon(aes(ymin = low95, ymax = upp95, fill = model), alpha = 0.18, color = NA) +
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
    scale_y_log10() +
    scale_color_manual(values = model_cols) +
    scale_fill_manual(values = model_cols) +
    scale_alpha_manual(values = model_alpha) +
    labs(
      x = "Spring temperature (Feb-Apr; z-scored)",
      y = "Recruitment",
      title = pop_target
    ) +
    tt +
    guides(color = "none", fill = "none", alpha = "none")


  pred_surv <- bind_rows(
    mutate(posterior_vec(fit_surv_reg, xpred, "pred", exp = TRUE), model = lev[1]),
    mutate(posterior_vec(fit_surv_err, xpred, "pred", exp = TRUE), model = lev[2])
  ) %>%
    mutate(model = factor(model, levels = lev))

  year_surv_err <- astr_fit %>%
    select(surv, surv_low, surv_upp, tmp) %>%
    rename(x = tmp) %>%
    mutate(model = lev[2])

  year_surv_pts <- bind_rows(
    year_surv_err %>% mutate(model = lev[1]),
    year_surv_err
  ) %>%
    mutate(model = factor(model, levels = lev))

  p_surv <- ggplot(pred_surv, aes(x = x)) +
    geom_line(aes(y = med, color = model), linewidth = 0.7) +
    geom_ribbon(aes(ymin = low95, ymax = upp95, fill = model), alpha = 0.20, color = NA) +
    geom_point(data = year_surv_pts, aes(y = surv, color = model), size = 1) +
    geom_linerange(data = year_surv_err, aes(ymin = surv_low, ymax = surv_upp, color = model)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(values = model_cols) +
    scale_fill_manual(values = model_cols) +
    facet_wrap(~model, ncol = 1) +
    labs(
      x = "Spring temperature (Feb-Apr; z-scored)",
      y = "Survival",
      title = pop_target
    ) +
    tt +
    guides(color = "none", fill = "none")


  # stage-specific survival diagnostics (GLM, logit link) ----
  stage_surv <- astr_fit %>%
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
      surv = n_surv / pmax(n_stage, 1L),
      surv_low = qbeta(0.025, 1 + n_surv, 1 + n_fail),
      surv_upp = qbeta(0.975, 1 + n_surv, 1 + n_fail),
      stage_label = paste0("Stage ", stage)
    ) %>%
    filter(n_stage > 0)

  n_stages <- astr_fit %>%
    mutate(n_stage_class = map_int(N, length)) %>%
    pull(n_stage_class) %>%
    max(na.rm = TRUE)

  stage_beta <- stage_surv %>%
    group_by(stage, stage_label) %>%
    group_modify(~ {
      if (nrow(.x) < 5 || n_distinct(.x$tmp) < 2 || sum(.x$n_surv) == 0 || sum(.x$n_fail) == 0) {
        return(tibble(
          model = c("Model of point estimates", "Model with sampling uncertainty"),
          beta = NA_real_, low95 = NA_real_, upp95 = NA_real_, n_year = nrow(.x)
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
        model = c("Model of point estimates", "Model with sampling uncertainty"),
        beta = c(beta_pt, beta_bin),
        low95 = c(beta_pt - 1.96 * se_pt, beta_bin - 1.96 * se_bin),
        upp95 = c(beta_pt + 1.96 * se_pt, beta_bin + 1.96 * se_bin),
        n_year = nrow(.x)
      )
    }) %>%
    ungroup() %>%
    mutate(model = factor(model, levels = c("Model of point estimates", "Model with sampling uncertainty"))) %>%
    arrange(stage, model)

  stage_levels <- paste0("Stage ", seq_len(n_stages))

  stage_placeholders <- stage_surv %>%
    group_by(stage, stage_label) %>%
    summarize(
      mean_surv = mean(surv, na.rm = TRUE),
      all_survive = all(n_fail == 0),
      .groups = "drop"
    ) %>%
    right_join(
      tibble(stage = seq_len(n_stages), stage_label = paste0("Stage ", seq_len(n_stages))),
      by = c("stage", "stage_label")
    ) %>%
    filter(!(stage %in% unique(stage_beta$stage[!is.na(stage_beta$beta)]))) %>%
    mutate(
      x = case_when(
        isTRUE(all_survive) ~ 1,
        TRUE ~ mean_surv
      )
    )

  stage_pred <- stage_surv %>%
    group_by(stage, stage_label) %>%
    group_modify(~ {
      if (nrow(.x) < 5 || n_distinct(.x$tmp) < 2 || sum(.x$n_surv) == 0 || sum(.x$n_fail) == 0) {
        return(tibble())
      }
      surv_clamped <- pmin(pmax(.x$surv, 1e-6), 1 - 1e-6)
      fit_pt <- glm(surv_clamped ~ tmp, weights = ssd_w, family = quasibinomial(), data = .x)
      xseq <- seq(min(.x$tmp), max(.x$tmp), length.out = 100)
      pred_pt <- predict(fit_pt, newdata = tibble(tmp = xseq), type = "link", se.fit = TRUE)

      fit_bin <- glm(cbind(n_surv, n_fail) ~ tmp, family = binomial(), data = .x)
      pred <- predict(fit_bin, newdata = tibble(tmp = xseq), type = "link", se.fit = TRUE)
      bind_rows(
        tibble(
          model = "Model of point estimates",
          tmp = xseq,
          med = plogis(pred_pt$fit),
          low95 = plogis(pred_pt$fit - 1.96 * pred_pt$se.fit),
          upp95 = plogis(pred_pt$fit + 1.96 * pred_pt$se.fit)
        ),
        tibble(
          model = "Model with sampling uncertainty",
          tmp = xseq,
          med = plogis(pred$fit),
          low95 = plogis(pred$fit - 1.96 * pred$se.fit),
          upp95 = plogis(pred$fit + 1.96 * pred$se.fit)
        )
      )
    }) %>%
    ungroup() %>%
    mutate(model = factor(model, levels = c("Model of point estimates", "Model with sampling uncertainty")))

  stage_surv_plot <- stage_surv %>%
    crossing(model = factor(c("Model of point estimates", "Model with sampling uncertainty"),
      levels = c("Model of point estimates", "Model with sampling uncertainty")
    ))

  stage_weight_sensitivity <- stage_surv %>%
    group_by(stage, stage_label) %>%
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
    ungroup() %>%
    mutate(SpeciesAuthor = spp_target, MatrixPopulation = pop_target, site_slug = site_slug)
  write_csv(
    stage_weight_sensitivity %>% mutate(across(where(is.numeric), ~ signif(.x, 3))),
    paste0("data/derived/analysis_cache/case3_astragalus_stage_weight_sensitivity_", site_slug, ".csv")
  )

  p_surv_stage_beta <- ggplot(
    stage_beta %>% filter(!is.na(beta)),
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
    scale_color_manual(
      values = model_cols,
      drop = FALSE,
      labels = c("Point estimates", "Sampling uncertainty")
    ) +
    scale_y_discrete(limits = rev(stage_levels), drop = FALSE) +
    labs(
      x = expression(paste("Temperature coefficient (logit scale, ", italic(beta), ")")),
      y = NULL,
      title = pop_target
    ) +
    tt +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(2, 16, 2, 2)
    )


  p_surv_stage_curves <- ggplot(stage_pred, aes(x = tmp, y = med, color = model, fill = model)) +
    geom_ribbon(aes(ymin = low95, ymax = upp95), alpha = 0.20, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(data = stage_surv_plot, aes(y = surv, color = model), size = 0.8) +
    geom_linerange(
      data = stage_surv_plot,
      aes(x = tmp, ymin = surv_low, ymax = surv_upp, color = model),
      inherit.aes = FALSE
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(values = model_cols) +
    scale_fill_manual(values = model_cols) +
    facet_grid(model ~ stage_label) +
    labs(
      x = "Spring temperature (Feb-Apr; z-scored)",
      y = "Stage-specific survival",
      title = pop_target
    ) +
    tt +
    guides(color = "none", fill = "none")


  # moving-beta model (Figure 5 style) ----
  focal_yrs <- seq(min(astr$MatrixStartYear) - 1, max(astr$MatrixStartYear) + 1)
  wx_mb <- wx %>%
    filter(SpeciesAuthor == spp_target, MatrixPopulation == pop_target, Year %in% focal_yrs) %>%
    group_by(Month) %>%
    mutate(tmp = as.numeric(scale(tmp)), ppt = as.numeric(scale(ppt))) %>%
    ungroup() %>%
    mutate(date = as.Date(paste(Year, Month, "01", sep = "-")))

  year <- astr_fit$MatrixEndYear
  y <- log(astr_fit$fecund)
  N <- length(y)
  K <- 24
  month_start <- "07"

  X <- map(seq_along(y), ~ {
    yr_focal <- year[.x]
    date_origin <- as.Date(paste(yr_focal, month_start, "01", sep = "-"))
    dates_focal <- sort(seq(date_origin, by = "-1 month", length.out = K))
    filter(wx_mb, date %in% dates_focal)$tmp
  }) %>%
    map(~ tibble_row(!!!setNames(as.list(.x), paste0("m", seq_len(K))))) %>%
    list_rbind() %>%
    as.matrix()

  dat_mb_reg <- list(N = N, K = K, X = X, y = y)
  dat_mb_err <- list(N = N, K = K, X = X, y = astr_fit$n_offsp, offset = astr_fit$n_repro)

  fit_null_reg <- stanfn(mod_null_reg, data = dat_mb_reg, seed = seed)
  fit_null_err <- stanfn(mod_null_err, data = dat_mb_err, seed = seed)
  fit_gprc_reg <- stanfn(mod_gprc_reg, data = dat_mb_reg, control = ctrl2, seed = seed)
  fit_gprc_err <- stanfn(mod_gprc_err, data = dat_mb_err, control = ctrl2, seed = seed)

  gprc_betas <- bind_rows(
    summarize_beta(fit_gprc_reg, "Model of point estimates"),
    summarize_beta(fit_gprc_err, "Model with sampling uncertainty")
  ) %>%
    mutate(
      model = factor(model, levels = c("Model of point estimates", "Model with sampling uncertainty")),
      lag_err = lag + if_else(model == "Model of point estimates", -0.18, 0.18),
      SpeciesAuthor = spp_target,
      MatrixPopulation = pop_target,
      site_slug = site_slug
    ) %>%
    relocate(SpeciesAuthor, MatrixPopulation, site_slug, model)

  write_csv(
    gprc_betas %>% mutate(across(where(is.numeric), ~ signif(.x, 3))),
    paste0("data/derived/analysis_cache/case3_astragalus_gprc_beta_summary_", site_slug, ".csv")
  )

  p2 <- ggplot(gprc_betas, aes(x = lag)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5, color = cols$accent) +
    geom_linerange(
      aes(x = lag_err, ymin = beta_low95, ymax = beta_upp95, color = model),
      alpha = 0.85
    ) +
    geom_line(aes(y = beta_med, color = model, alpha = model), linewidth = 0.8) +
    scale_x_continuous(breaks = seq(0, 24, 6)) +
    scale_color_manual(values = model_cols) +
    scale_alpha_manual(values = model_alpha) +
    labs(
      x = "Months before survey",
      y = expression(paste("Temperature coefficient (", italic(b), ")"))
    ) +
    tt +
    theme(plot.margin = margin(2, 2, 2, 2)) +
    guides(color = "none", alpha = "none")


  site_component_panels[[site_slug]] <<- list(
    p2 = p2 + labs(title = pop_target) + theme(plot.title = element_text(hjust = 0.5)),
    p_fit = p_fit + labs(title = NULL),
    p3 = p_surv_stage_beta + labs(title = NULL)
  )

  tibble(
    SpeciesAuthor = spp_target,
    MatrixPopulation = pop_target,
    site_slug = site_slug,
    years_total = nrow(coverage_tbl),
    years_fit = nrow(astr_fit),
    years_excluded = nrow(coverage_tbl) - nrow(astr_fit),
    beta_med_point = df_beta %>% filter(model == "Model of point estimates") %>% pull(med),
    beta_med_sampling = df_beta %>% filter(model == "Model with sampling uncertainty") %>% pull(med),
    status = "ok"
  )
}


site_results <- purrr::pmap_dfr(site_key, analyze_site)
write_csv(site_results, "data/derived/analysis_cache/case3_astragalus_site_summary.csv")

all_beta <- list.files(
  "data/derived/analysis_cache",
  pattern = "^case3_astragalus_spring_beta_summary_.*\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(read_csv, show_col_types = FALSE)
write_csv(all_beta, "data/derived/analysis_cache/case3_astragalus_spring_beta_summary.csv")

all_diag <- list.files(
  "data/derived/analysis_cache",
  pattern = "^case3_astragalus_stan_diagnostics_.*\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(read_csv, show_col_types = FALSE)
write_csv(all_diag, "data/derived/analysis_cache/case3_astragalus_stan_diagnostics.csv")

all_gprc <- list.files(
  "data/derived/analysis_cache",
  pattern = "^case3_astragalus_gprc_beta_summary_.*\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(read_csv, show_col_types = FALSE)
write_csv(all_gprc, "data/derived/analysis_cache/case3_astragalus_gprc_beta_summary.csv")

all_stage_w <- list.files(
  "data/derived/analysis_cache",
  pattern = "^case3_astragalus_stage_weight_sensitivity_.*\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(read_csv, show_col_types = FALSE)
write_csv(all_stage_w, "data/derived/analysis_cache/case3_astragalus_stage_weight_sensitivity.csv")

# Combine site-level analysis 3 figure panels into one multisite figure ----
panel_order <- c("haynes_creek", "mcdevitt_creek", "sheep_corral_gulch")
valid_sites <- panel_order[panel_order %in% names(site_component_panels)]

if (length(valid_sites) > 0) {
  p_multi <- wrap_plots(
    c(
      map(valid_sites, ~ site_component_panels[[.x]]$p2),
      map(valid_sites, ~ site_component_panels[[.x]]$p_fit),
      map(valid_sites, ~ site_component_panels[[.x]]$p3)
    ),
    ncol = 3
  ) +
    plot_annotation(tag_levels = "A")

  ggsave(
    "figures/Figure_6_analysis3_climate_effects_multisite.png",
    p_multi,
    height = 225,
    width = 270,
    units = "mm",
    dpi = 300
  )
}
