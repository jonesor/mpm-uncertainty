# S09: diagnose stage-specific survival issues in COMPADRE matrices.

# libraries ----
source("code/setup.R")
setup_packages(c("tidyverse", "Rcompadre", "Rage", "popbio", "popdemo"))
source("code/functions.R")
set_mpm_plot_defaults()
cols <- mpm_colors()


# load COMPADRE ----
compadre <- load_compadre(corrected = TRUE)


# possible columns to collapse on ----
col_spp <- c("id_stage", "SpeciesAuthor", "ProjectionInterval")
col_pop <- c(col_spp, "MatrixPopulation", "MatrixTreatment")


# subset ----
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


# number of reproductive stages ----
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
  mutate(perennial = map2_lgl(matU, start, ~ mpm_to_lx(.x, .y, xmax = 3)[4] > 0)) %>%
  mutate(any_rep = map_lgl(matF, ~ any(.x > 0))) %>%
  filter(perennial == TRUE, any_rep == TRUE) %>%
  mutate(rep_prop1 = pmap(list(matU, start, rep_stages), repro_prop_start)) %>%
  mutate(rep_na = map_lgl(rep_prop1, ~ any(is.nan(.x)))) %>%
  filter(rep_na == FALSE) %>%
  mutate(lx = map2(matU, rep_prop1, lx_from_mature)) %>%
  mutate(lx_n = map_int(lx, length)) %>%
  mutate(q = map2_int(matU, rep_prop1, ~ qsd_safe(.x, .y))) %>%
  filter(!is.na(q)) %>%
  mutate(lxs = map2(lx, q, lx_submax)) %>%
  mutate(lxs_min = map_dbl(lxs, min)) %>%
  mutate(l0_pt = map_dbl(lx, sum)) %>%
  mutate(l0_pt_int = as.integer(round(l0_pt, 0))) %>%
  mutate(lxs_n = map_int(lxs, length)) %>%
  mutate(nrep = map_int(matF, ~ length(which(colSums(.x) > 0))))

table(out_shape$q)[1:20]

sum(out_shape$q < out_shape$l0_pt) / nrow(out_shape)


sum(out_shape$lxs_min > 0.1) / nrow(out_shape)

out_shape <- out_shape %>%
  filter(lxs_n >= 3) %>%
  mutate(shape_pt = map2_dbl(lxs, q + 1, shape_surv2)) %>%
  mutate(shape_l0_pt = map2_dbl(lx, l0_pt_int, ~ 1 + log(.x[.y]))) %>%
  as_tibble() %>%
  filter(!is.na(shape_pt)) %>%
  mutate(id = row_number()) %>%
  mutate(id_shape = fct_reorder(fct_drop(as.factor(id)), shape_pt)) %>%
  mutate(id_l0 = fct_reorder(fct_drop(as.factor(id)), l0_pt))


# Analysis 1 displacement diagnostics ----
sd_files <- paste0("data/derived/analysis_cache/", list.files("data/derived/analysis_cache"))
sd_files <- sd_files[grep("/sd_", sd_files)]

mpm_draws <- cdb_bind_rows(map(sd_files, rdata_load)) %>%
  mutate(id = as.factor(dplyr::row_number())) %>%
  cdb_unnest() %>%
  mutate(any_repro = map_lgl(matF, ~ any(.x > 0))) %>%
  filter(any_repro == TRUE) %>%
  mutate(matU = map(matU, scale_U)) %>%
  mutate(matA = pmap(list(matU, matF, matC), ~ ..1 + ..2 + ..3)) %>%
  mutate(start = map_int(mat, Rcompadre::mpm_first_active)) %>%
  mutate(rep_stages = map(matF, ~ colSums(.x) > 0))

pt_shape <- mpm_draws %>%
  mutate(rep_prop1 = pmap(list(matU, start, rep_stages), Rage::mature_distrib)) %>%
  mutate(lx4 = map2_dbl(
    matU, rep_prop1,
    ~ Rage::mpm_to_lx(.x, .y, lx_crit = -1, xmax = 3)[4]
  )) %>%
  filter(lx4 > 0) %>%
  mutate(q = map2_int(matU, rep_prop1, ~ Rage::qsd_converge(.x, .y, conv = 0.01, N = 1e5))) %>%
  filter(q >= 3) %>%
  mutate(lx = pmap(list(matU, rep_prop1, q),
    ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3),
    lx_crit = -1
  )) %>%
  mutate(
    L_pt = map2_dbl(matU, rep_prop1, life_expect),
    S_pt = map_dbl(lx, Rage::shape_surv),
    stage_survival = map(matU, colSums),
    stage_N = map(N, as.numeric)
  ) %>%
  mutate(boundary_tbl = map2(stage_survival, stage_N, function(surv, n_stage) {
    if (length(n_stage) == 0 || length(n_stage) != length(surv)) {
      n_stage <- rep(NA_real_, length(surv))
    }

    tibble(
      stage_index = seq_along(surv),
      survival = as.numeric(surv),
      n_stage = as.numeric(n_stage),
      boundary_any = survival %in% c(0, 1)
    )
  })) %>%
  mutate(
    n_boundary_stages = map_int(boundary_tbl, ~ sum(.x$boundary_any, na.rm = TRUE)),
    prop_boundary_stages = map_dbl(boundary_tbl, ~ mean(.x$boundary_any, na.rm = TRUE)),
    prop_individuals_boundary = map_dbl(boundary_tbl, function(x) {
      if (all(is.na(x$n_stage)) || sum(x$n_stage, na.rm = TRUE) == 0) {
        return(NA_real_)
      }
      sum(x$n_stage[x$boundary_any], na.rm = TRUE) / sum(x$n_stage, na.rm = TRUE)
    })
  ) %>%
  as_tibble()

sd_shape <- pt_shape %>%
  select(id, SpeciesAuthor, MatrixPopulation, simU, simF, q) %>%
  unnest(cols = c("simU", "simF")) %>%
  left_join(select(pt_shape, id, start, rep_stages), by = "id") %>%
  mutate(rep_prop1 = pmap(list(simU, start, rep_stages), Rage::mature_distrib)) %>%
  mutate(lx = pmap(list(simU, rep_prop1, q),
    ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3),
    lx_crit = -1
  )) %>%
  mutate(
    L = map2_dbl(simU, rep_prop1, life_expect),
    S = map_dbl(lx, Rage::shape_surv)
  )

disp_summary <- sd_shape %>%
  group_by(id) %>%
  summarize(
    L_med = median(L),
    L_low = quantile(L, 0.025),
    L_upp = quantile(L, 0.975),
    S_med = median(S),
    S_low = quantile(S, 0.025),
    S_upp = quantile(S, 0.975),
    .groups = "drop"
  ) %>%
  left_join(
    pt_shape %>%
      select(
        id, SpeciesAuthor, SpeciesAccepted, MatrixPopulation, OrganismType,
        MatrixDimension, L_pt, S_pt, n_boundary_stages, prop_boundary_stages,
        prop_individuals_boundary
      ),
    by = "id"
  ) %>%
  mutate(
    disp_L = L_pt - L_med,
    disp_log_L = log10(L_pt) - log10(L_med),
    disp_S = S_pt - S_med,
    abs_disp_L = abs(disp_L),
    abs_disp_log_L = abs(disp_log_L),
    abs_disp_S = abs(disp_S)
  )

model_tbl <- tribble(
  ~response, ~predictor,
  "abs_disp_S", "n_boundary_stages",
  "abs_disp_S", "prop_individuals_boundary",
  "abs_disp_log_L", "n_boundary_stages",
  "abs_disp_log_L", "prop_individuals_boundary"
) %>%
  mutate(
    fit = map2(response, predictor, ~ lm(stats::as.formula(paste(.x, "~", .y)), data = disp_summary)),
    fit_sum = map(fit, summary),
    conf = map(fit, confint)
  ) %>%
  mutate(
    slope = map2_dbl(fit, predictor, ~ coef(.x)[[.y]]),
    conf_low = map2_dbl(conf, predictor, ~ .x[.y, 1]),
    conf_upp = map2_dbl(conf, predictor, ~ .x[.y, 2]),
    p_value = map2_dbl(fit_sum, predictor, ~ coef(.x)[.y, "Pr(>|t|)"]),
    r_squared = map_dbl(fit_sum, "r.squared"),
    n = map_int(fit, stats::nobs)
  ) %>%
  select(response, predictor, n, slope, conf_low, conf_upp, p_value, r_squared)

write_csv(
  disp_summary %>% mutate(across(where(is.numeric), ~ signif(.x, 4))),
  "data/derived/analysis_cache/case1_displacement_boundary_summary.csv"
)
write_csv(
  model_tbl %>% mutate(across(where(is.numeric), ~ signif(.x, 4))),
  "data/derived/analysis_cache/case1_displacement_boundary_models.csv"
)

plot_df <- disp_summary %>%
  select(
    SpeciesAccepted, MatrixPopulation, OrganismType,
    n_boundary_stages, prop_individuals_boundary,
    abs_disp_S, abs_disp_log_L
  ) %>%
  pivot_longer(
    cols = c(abs_disp_S, abs_disp_log_L),
    names_to = "response",
    values_to = "displacement"
  ) %>%
  pivot_longer(
    cols = c(n_boundary_stages, prop_individuals_boundary),
    names_to = "predictor",
    values_to = "predictor_value"
  ) %>%
  mutate(
    response = recode(
      response,
      abs_disp_S = "Shape displacement |point - median|",
      abs_disp_log_L = "Life expectancy displacement |log10(point) - log10(median)|"
    ),
    predictor = recode(
      predictor,
      n_boundary_stages = "Number of boundary survival stages",
      prop_individuals_boundary = "Proportion of individuals in boundary-survival stages"
    )
  )

p_disp <- ggplot(plot_df, aes(predictor_value, displacement)) +
  geom_point(color = cols$accent, alpha = 0.8, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, color = cols$dark, fill = cols$fill, linewidth = 0.7) +
  facet_grid(response ~ predictor, scales = "free_x") +
  labs(
    x = NULL,
    y = "Displacement magnitude",
    title = NULL
  ) +
  theme_mpm() +
  theme(
    strip.text = element_text(size = 9.5),
    panel.spacing = unit(6, "pt"),
    plot.margin = margin(4, 6, 4, 4)
  )

ggsave(
  "figures/Figure_S3_displacement_boundary_diagnostics.png",
  p_disp,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 300
)
