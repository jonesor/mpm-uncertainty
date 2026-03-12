# Validate reproducibility outputs from a Zenodo bundle ----
#
# Usage from bundle/project root:
#   source("scripts/validate_reproduction.R")
#   validate_reproduction(run_pipeline = FALSE)
#
# CLI:
#   Rscript --vanilla scripts/validate_reproduction.R
#   Rscript --vanilla scripts/validate_reproduction.R --run

validate_reproduction <- function(run_pipeline = FALSE, check_manifest = TRUE) {
  if (!file.exists("RUN_REPRODUCTION.R")) {
    stop("RUN_REPRODUCTION.R not found. Run this script from the bundle root.")
  }

  if (run_pipeline) {
    message("Running full reproduction pipeline via RUN_REPRODUCTION.R ...")
    source("RUN_REPRODUCTION.R")
  }

  source("code/setup.R")
  active_raw <- get_compadre_path(corrected = FALSE)
  active_corrected <- get_compadre_path(corrected = TRUE)

  expected_bundle_files <- c(
    "docs/manuscript/sampling_uncertainty_mpm_main.docx",
    "docs/manuscript/sampling_uncertainty_mpm_supplement.docx",
    "figures/Figure_1_transition_rates_single_mpm.png",
    "figures/Figure_2_derived_parameters_single_mpm.png",
    "figures/Figure_3_analysis1_point_vs_sampling_distributions.png",
    "figures/Figure_4_analysis1_life_expectancy_shape_relationship.png",
    "figures/Figure_5_analysis2_climate_effects_recruitment.png",
    "figures/Figure_6_analysis3_climate_effects_multisite.png",
    active_raw,
    active_corrected,
    "data/raw/ellis_2012/Transition_Matrices.txt",
    "data/derived/climate/species_coords.csv",
    "data/derived/climate/species_clim_prism.csv",
    "data/derived/studies/_data_sources.csv"
  )

  expected_generated_files <- c(
    "docs/manuscript/sampling_uncertainty_mpm_main.docx",
    "docs/manuscript/sampling_uncertainty_mpm_supplement.docx",
    "figures/Figure_1_transition_rates_single_mpm.png",
    "figures/Figure_2_derived_parameters_single_mpm.png",
    "figures/Figure_3_analysis1_point_vs_sampling_distributions.png",
    "figures/Figure_4_analysis1_life_expectancy_shape_relationship.png",
    "figures/Figure_5_analysis2_climate_effects_recruitment.png",
    "figures/Figure_6_analysis3_climate_effects_multisite.png",
    "data/derived/analysis_cache/case1_variance_ratios.csv",
    "data/derived/analysis_cache/case2_spring_beta_summary.csv",
    "data/derived/analysis_cache/case3_astragalus_site_summary.csv"
  )

  expected_files <- if (isTRUE(run_pipeline)) {
    expected_generated_files
  } else {
    expected_bundle_files
  }

  missing <- expected_files[!file.exists(expected_files)]
  if (length(missing) > 0) {
    message("Missing expected outputs:")
    for (f in missing) message("- ", f)
  } else {
    message("All expected output files are present.")
  }

  manifest_ok <- NA
  if (check_manifest && file.exists("MANIFEST.csv")) {
    if (!requireNamespace("readr", quietly = TRUE)) {
      warning("Package 'readr' is not installed. Skipping MANIFEST.csv check.")
    } else {
      manifest <- readr::read_csv("MANIFEST.csv", show_col_types = FALSE)
      has_cols <- all(c("path", "md5") %in% names(manifest))
      if (!has_cols) {
        warning("MANIFEST.csv is missing required columns ('path', 'md5').")
      } else {
        present <- manifest$path[file.exists(manifest$path)]
        calc <- unname(tools::md5sum(present))
        ref <- manifest$md5[match(present, manifest$path)]
        manifest_ok <- all(calc == ref)
        if (isTRUE(manifest_ok)) {
          message("MANIFEST.csv checksum check passed for present files.")
        } else {
          warning("MANIFEST.csv checksum check failed for one or more files.")
        }
      }
    }
  } else if (check_manifest) {
    warning("MANIFEST.csv not found. Skipping checksum validation.")
  }

  status <- list(
    outputs_missing_n = length(missing),
    outputs_missing = missing,
    manifest_ok = manifest_ok
  )

  if (length(missing) == 0) {
    message("Validation PASS.")
  } else {
    message("Validation FAIL.")
  }

  invisible(status)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  run_pipeline <- "--run" %in% args
  validate_reproduction(run_pipeline = run_pipeline, check_manifest = TRUE)
}
