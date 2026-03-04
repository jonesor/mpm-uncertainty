# Render main manuscript and supplement ----
#
# Usage (from project root):
#   source('docs/manuscript/render_manuscripts.R')
# This script also runs:
#   source("scripts/prepare_zenodo_bundle.R")
#   prepare_zenodo_bundle(skip_render = TRUE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render manuscripts.")
}

main_rmd <- "docs/manuscript/manuscript_main.Rmd"
supp_rmd <- "docs/manuscript/manuscript_supplement.Rmd"

if (!file.exists(main_rmd) || !file.exists(supp_rmd)) {
  stop("Expected manuscript files were not found in docs/manuscript/.")
}

main_out <- rmarkdown::render(
  main_rmd,
  output_format = "bookdown::word_document2",
  output_file = "sampling_uncertainty_mpm_main.docx",
  quiet = TRUE
)
supp_out <- rmarkdown::render(
  supp_rmd,
  output_format = "bookdown::word_document2",
  output_file = "sampling_uncertainty_mpm_supplement.docx",
  quiet = TRUE
)

cat("Rendered manuscript files:\n")
cat("- ", main_out, "\n", sep = "")
cat("- ", supp_out, "\n", sep = "")

# Copy manuscript figure PNGs to docs/manuscript with manuscript-numbered names.
figure_map <- c(
  "figures/Figure_1_transition_rates_single_mpm.png" = "docs/manuscript/Figure_1_transition_rates_single_mpm.png",
  "figures/Figure_2_derived_parameters_single_mpm.png" = "docs/manuscript/Figure_2_derived_parameters_single_mpm.png",
  "figures/Figure_3_analysis1_point_vs_sampling_distributions.png" = "docs/manuscript/Figure_3_analysis1_point_vs_sampling_distributions.png",
  "figures/Figure_4_analysis1_life_expectancy_shape_relationship.png" = "docs/manuscript/Figure_4_analysis1_life_expectancy_shape_relationship.png",
  "figures/Figure_5_analysis2_climate_effects_recruitment.png" = "docs/manuscript/Figure_5_analysis2_climate_effects_recruitment.png",
  "figures/Figure_S1_boundary_estimate_diagnostic.png" = "docs/manuscript/Figure_S1_boundary_estimate_diagnostic.png",
  "figures/Figure_S2_boundary_survivorship_illustration.png" = "docs/manuscript/Figure_S2_boundary_survivorship_illustration.png"
)

for (src in names(figure_map)) {
  dst <- figure_map[[src]]
  if (!file.exists(src)) {
    warning("Missing figure source file: ", src)
    next
  }
  file.copy(src, dst, overwrite = TRUE)
}

cat("\nCopied manuscript figure PNGs:\n")
for (dst in unname(figure_map)) {
  if (file.exists(dst)) cat("- ", dst, "\n", sep = "")
}

# Build Zenodo bundle after rendering (skip re-render to avoid duplication).
source("scripts/prepare_zenodo_bundle.R")
bundle_dir <- prepare_zenodo_bundle(skip_render = TRUE)
cat("\nZenodo bundle prepared at:\n")
cat("- ", bundle_dir, "\n", sep = "")
