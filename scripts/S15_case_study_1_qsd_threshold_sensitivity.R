# S15: sensitivity of analysis 1 shape estimates to the quasi-stable threshold.

# libraries ----
source("code/setup.R")
setup_packages(c(
  "tidyverse", "Rcompadre", "Rage", "patchwork", "viridisLite"
))
source("code/functions.R")
# Plot style helpers (theme_mpm(), mpm_colors()) from code/functions.R
set_mpm_plot_defaults()
cols <- mpm_colors()


# load study-specific sampling distribution files ----
sd_files <- file.path(
  "data/derived/analysis_cache",
  list.files("data/derived/analysis_cache", pattern = "^sd_.*\\.RData$")
)


# bind sampling distributions into a single tibble ----
mpm_draws <- cdb_bind_rows(map(sd_files, rdata_load)) %>%
  mutate(id = as.factor(dplyr::row_number())) %>%
  cdb_unnest() %>%
  mutate(any_repro = map_lgl(matF, ~ any(.x > 0))) %>%
  filter(any_repro) %>%
  mutate(matU = map(matU, scale_U)) %>%
  mutate(start = map_int(mat, Rcompadre::mpm_first_active)) %>%
  mutate(rep_stages = map(matF, ~ colSums(.x) > 0)) %>%
  as_tibble()


# calculate point-estimate and bootstrap summaries across thresholds ----
threshold_tbl <- tibble(
  threshold = c(0.005, 0.01, 0.05),
  threshold_label = c("0.5%", "1%", "5%")
)

base_pt <- mpm_draws %>%
  mutate(rep_prop1 = pmap(list(matU, start, rep_stages), Rage::mature_distrib)) %>%
  mutate(lx4 = map2_dbl(
    matU, rep_prop1,
    ~ Rage::mpm_to_lx(.x, .y, lx_crit = -1, xmax = 3)[4]
  )) %>%
  filter(lx4 > 0)

threshold_results <- vector("list", length = nrow(threshold_tbl))

for (i in seq_len(nrow(threshold_tbl))) {
  threshold_value <- threshold_tbl$threshold[[i]]
  threshold_label <- threshold_tbl$threshold_label[[i]]

  pt_i <- base_pt %>%
    mutate(
      q = map2_int(
        matU, rep_prop1,
        ~ Rage::qsd_converge(.x, .y, conv = threshold_value, N = 1e5)
      )
    ) %>%
    filter(q >= 3) %>%
    mutate(
      lx = pmap(
        list(matU, rep_prop1, q),
        ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3),
        lx_crit = -1
      ),
      shape_pt = map_dbl(lx, Rage::shape_surv)
    ) %>%
    transmute(
      id, SpeciesAuthor, MatrixPopulation,
      threshold = threshold_value,
      threshold_label,
      q,
      shape_pt
    )

  sd_i <- pt_i %>%
    select(id, threshold, threshold_label, q) %>%
    left_join(
      base_pt %>%
        select(id, SpeciesAuthor, MatrixPopulation, simU, start, rep_stages),
      by = "id"
    ) %>%
    unnest(cols = "simU") %>%
    mutate(rep_prop1 = pmap(list(simU, start, rep_stages), Rage::mature_distrib)) %>%
    mutate(
      lx = pmap(
        list(simU, rep_prop1, q),
        ~ Rage::mpm_to_lx(..1, ..2, xmax = ..3),
        lx_crit = -1
      ),
      shape = map_dbl(lx, Rage::shape_surv)
    ) %>%
    group_by(id, threshold, threshold_label) %>%
    summarize(
      shape_med = quantile(shape, 0.5),
      shape_low = quantile(shape, 0.025),
      shape_upp = quantile(shape, 0.975),
      shape_sd = sd(shape),
      .groups = "drop"
    )

  threshold_results[[i]] <- pt_i %>%
    left_join(sd_i, by = c("id", "threshold", "threshold_label"))
}

shape_threshold_long <- bind_rows(threshold_results)


# keep only entries estimable under all three thresholds ----
common_ids <- shape_threshold_long %>%
  count(id, name = "n_thresholds") %>%
  filter(n_thresholds == nrow(threshold_tbl)) %>%
  pull(id)

shape_threshold_long <- shape_threshold_long %>%
  filter(id %in% common_ids)

shape_threshold_wide <- shape_threshold_long %>%
  select(
    id, SpeciesAuthor, MatrixPopulation, threshold_label,
    q, shape_pt, shape_med, shape_low, shape_upp, shape_sd
  ) %>%
  pivot_wider(
    names_from = threshold_label,
    values_from = c(q, shape_pt, shape_med, shape_low, shape_upp, shape_sd),
    names_glue = "{.value}_{threshold_label}"
  ) %>%
  mutate(
    abs_diff_pt_0p5_vs_1 = abs(`shape_pt_0.5%` - `shape_pt_1%`),
    abs_diff_pt_5_vs_1 = abs(`shape_pt_5%` - `shape_pt_1%`),
    abs_diff_med_0p5_vs_1 = abs(`shape_med_0.5%` - `shape_med_1%`),
    abs_diff_med_5_vs_1 = abs(`shape_med_5%` - `shape_med_1%`),
    sign_change_pt_0p5_vs_1 = sign(`shape_pt_0.5%`) != sign(`shape_pt_1%`),
    sign_change_pt_5_vs_1 = sign(`shape_pt_5%`) != sign(`shape_pt_1%`),
    sign_change_med_0p5_vs_1 = sign(`shape_med_0.5%`) != sign(`shape_med_1%`),
    sign_change_med_5_vs_1 = sign(`shape_med_5%`) != sign(`shape_med_1%`)
  )

summary_threshold_cmp <- tibble(
  comparison = c("0.5% vs 1%", "5% vs 1%"),
  n = nrow(shape_threshold_wide),
  cor_shape_pt = c(
    cor(shape_threshold_wide$`shape_pt_0.5%`, shape_threshold_wide$`shape_pt_1%`),
    cor(shape_threshold_wide$`shape_pt_5%`, shape_threshold_wide$`shape_pt_1%`)
  ),
  cor_shape_med = c(
    cor(shape_threshold_wide$`shape_med_0.5%`, shape_threshold_wide$`shape_med_1%`),
    cor(shape_threshold_wide$`shape_med_5%`, shape_threshold_wide$`shape_med_1%`)
  ),
  median_abs_diff_pt = c(
    median(shape_threshold_wide$abs_diff_pt_0p5_vs_1),
    median(shape_threshold_wide$abs_diff_pt_5_vs_1)
  ),
  p95_abs_diff_pt = c(
    quantile(shape_threshold_wide$abs_diff_pt_0p5_vs_1, 0.95),
    quantile(shape_threshold_wide$abs_diff_pt_5_vs_1, 0.95)
  ),
  max_abs_diff_pt = c(
    max(shape_threshold_wide$abs_diff_pt_0p5_vs_1),
    max(shape_threshold_wide$abs_diff_pt_5_vs_1)
  ),
  median_abs_diff_med = c(
    median(shape_threshold_wide$abs_diff_med_0p5_vs_1),
    median(shape_threshold_wide$abs_diff_med_5_vs_1)
  ),
  p95_abs_diff_med = c(
    quantile(shape_threshold_wide$abs_diff_med_0p5_vs_1, 0.95),
    quantile(shape_threshold_wide$abs_diff_med_5_vs_1, 0.95)
  ),
  max_abs_diff_med = c(
    max(shape_threshold_wide$abs_diff_med_0p5_vs_1),
    max(shape_threshold_wide$abs_diff_med_5_vs_1)
  ),
  n_sign_change_pt = c(
    sum(shape_threshold_wide$sign_change_pt_0p5_vs_1),
    sum(shape_threshold_wide$sign_change_pt_5_vs_1)
  ),
  n_sign_change_med = c(
    sum(shape_threshold_wide$sign_change_med_0p5_vs_1),
    sum(shape_threshold_wide$sign_change_med_5_vs_1)
  ),
  median_delta_q = c(
    median(shape_threshold_wide$`q_0.5%` - shape_threshold_wide$`q_1%`),
    median(shape_threshold_wide$`q_5%` - shape_threshold_wide$`q_1%`)
  ),
  range_delta_q = c(
    paste0(
      min(shape_threshold_wide$`q_0.5%` - shape_threshold_wide$`q_1%`), " to ",
      max(shape_threshold_wide$`q_0.5%` - shape_threshold_wide$`q_1%`)
    ),
    paste0(
      min(shape_threshold_wide$`q_5%` - shape_threshold_wide$`q_1%`), " to ",
      max(shape_threshold_wide$`q_5%` - shape_threshold_wide$`q_1%`)
    )
  )
)


# write outputs ----
if (!dir.exists("data/derived/analysis_cache")) {
  dir.create("data/derived/analysis_cache", recursive = TRUE)
}

write_csv(
  shape_threshold_long,
  "data/derived/analysis_cache/case1_qsd_threshold_sensitivity_long.csv"
)
write_csv(
  shape_threshold_wide,
  "data/derived/analysis_cache/case1_qsd_threshold_sensitivity_by_population.csv"
)
write_csv(
  summary_threshold_cmp,
  "data/derived/analysis_cache/case1_qsd_threshold_sensitivity_summary.csv"
)


# figure ----
plot_dat <- bind_rows(
  shape_threshold_wide %>%
    transmute(
      comparison = "0.5% threshold versus 1% threshold",
      baseline = `shape_pt_1%`,
      alternative = `shape_pt_0.5%`
    ),
  shape_threshold_wide %>%
    transmute(
      comparison = "5% threshold versus 1% threshold",
      baseline = `shape_pt_1%`,
      alternative = `shape_pt_5%`
    )
)

p <- ggplot(plot_dat, aes(x = baseline, y = alternative)) +
  geom_abline(
    intercept = 0, slope = 1,
    color = cols$mid, linetype = 2, linewidth = 0.5
  ) +
  geom_point(color = cols$accent, alpha = 0.7, size = 1.6) +
  facet_wrap(~ comparison, ncol = 2) +
  labs(
    x = "Shape at 1% threshold",
    y = "Shape at alternative threshold"
  ) +
  theme_mpm() +
  theme(
    strip.text = element_text(size = 10.5),
    plot.margin = margin(3, 3, 3, 3)
  )

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
ggsave(
  "figures/Figure_S4_shape_threshold_sensitivity.png",
  p,
  height = 90,
  width = 180,
  units = "mm",
  dpi = 300
)
