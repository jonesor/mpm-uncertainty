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

# Build Zenodo bundle after rendering (skip re-render to avoid duplication).
source("scripts/prepare_zenodo_bundle.R")
bundle_dir <- prepare_zenodo_bundle(skip_render = TRUE)
cat("\nZenodo bundle prepared at:\n")
cat("- ", bundle_dir, "\n", sep = "")
