# Prepare a reproducible Zenodo deposition bundle ----
#
# This script creates a timestamped snapshot in
# docs/manuscript/zenodo/staging/ that includes code, data, models,
# generated outputs, and a manifest. It can also create an anonymised
# review bundle that excludes manuscript-facing files.
#
# RStudio usage:
#   source('scripts/prepare_zenodo_bundle.R')
#   prepare_zenodo_bundle(skip_render = FALSE)
#   prepare_zenodo_bundle(skip_render = FALSE, bundle_type = "anonymized_review")
#
# CLI usage:
#   Rscript --vanilla scripts/prepare_zenodo_bundle.R
#   Rscript --vanilla scripts/prepare_zenodo_bundle.R --skip-render

prepare_zenodo_bundle <- function(
  skip_render = FALSE,
  clean_old_staging = TRUE,
  bundle_type = c("standard", "anonymized_review")
) {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required.")
  }

  bundle_type <- match.arg(bundle_type)

  # Make behaviour predictable in RStudio even if the working directory changed.
  setwd(here::here())

  source("code/setup.R")
  setup_packages(c("fs", "here", "readr", "tibble", "dplyr", "purrr"))

  staging_root <- here::here("docs", "manuscript", "zenodo", "staging")
  archive_root <- here::here("docs", "manuscript", "zenodo", "archive")

  fs::dir_create(staging_root, recurse = TRUE)
  fs::dir_create(archive_root, recurse = TRUE)

  if (clean_old_staging) {
    old_bundles <- fs::dir_ls(
      staging_root,
      regexp = "mpm-uncertainty_(zenodo|review)_",
      type = "directory"
    )
    if (length(old_bundles) > 0) {
      purrr::walk(old_bundles, ~unlink(.x, recursive = TRUE, force = TRUE))
    }
  }

  bundle_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  bundle_prefix <- if (bundle_type == "standard") {
    "mpm-uncertainty_zenodo_"
  } else {
    "mpm-uncertainty_anonymised_review_"
  }
  bundle_name <- paste0(bundle_prefix, bundle_timestamp)
  bundle_dir <- fs::path(staging_root, bundle_name)

  fs::dir_create(bundle_dir, recurse = TRUE)

  if (!skip_render && bundle_type == "standard") {
    source("docs/manuscript/render_manuscripts.R")
  }

  copy_paths <- c(
    "README.md",
    "code",
    "models/README.md",
    "models/categorical_logit_hier_nc.stan",
    "models/categorical_logit_hier_nc_pop.stan",
    "models/categorical_logit_hier_nc_spp.stan",
    "models/movbeta_gprc.stan",
    "models/movbeta_gprc_err.stan",
    "models/null.stan",
    "models/null_err.stan",
    "models/poisson_hier.stan",
    "models/regress.stan",
    "models/regress2.stan",
    "models/regress_error.stan",
    "models/regress_error2.stan",
    "models/regress_hier.stan",
    "models/regress_hier_error.stan",
    "models/regress_log_err.stan",
    "models/varcomp.stan",
    "data/README.md",
    "data/metadata",
    "data/raw/README.md",
    "data/raw/compadre/README.md",
    "data/raw/compadre/COMPADRE_v6.26.3.0.RData",
    "data/raw/compadre/COMPADRE_v6.26.3.0_Corrected.RData",
    "data/raw/ellis_2012",
    "data/derived/README.md",
    "data/derived/studies/README.md",
    "data/derived/studies/_data_sources.csv",
    "data/derived/studies/aschero_U.RData",
    "data/derived/studies/andrieu_n.csv",
    "data/derived/studies/andrello_n.csv",
    "data/derived/studies/arroyo_n.csv",
    "data/derived/studies/auestad_n.csv",
    "data/derived/studies/crone_n.csv",
    "data/derived/studies/csergo_n.csv",
    "data/derived/studies/dias_n.csv",
    "data/derived/studies/dostalek_n.csv",
    "data/derived/studies/ehrlen_n.csv",
    "data/derived/studies/eriksson_n.csv",
    "data/derived/studies/ferrer_n.csv",
    "data/derived/studies/flores_n.csv",
    "data/derived/studies/jacquemyns_n.csv",
    "data/derived/studies/jolls_n.csv",
    "data/derived/studies/keller_n.csv",
    "data/derived/studies/kiviniemi_n.csv",
    "data/derived/studies/law_n.csv",
    "data/derived/studies/lazaro_n.csv",
    "data/derived/studies/lemke_n.csv",
    "data/derived/studies/lopez_n.csv",
    "data/derived/studies/martin_n.csv",
    "data/derived/studies/martinez_n.csv",
    "data/derived/studies/noel_n.csv",
    "data/derived/studies/plank_n.csv",
    "data/derived/studies/portela_n.csv",
    "data/derived/studies/raghu_n.csv",
    "data/derived/studies/satterthwaite_n.csv",
    "data/derived/studies/scanga_n.csv",
    "data/derived/studies/shryock_n.csv",
    "data/derived/studies/target_studies.csv",
    "data/derived/studies/toledo_n.csv",
    "data/derived/studies/torres_n.csv",
    "data/derived/climate/README.md",
    "data/derived/climate/species_coords.csv",
    "data/derived/climate/species_clim_prism.csv",
    "figures/README.md",
    "figures/Figure_1_transition_rates_single_mpm.png",
    "figures/Figure_2_derived_parameters_single_mpm.png",
    "figures/Figure_3_analysis1_point_vs_sampling_distributions.png",
    "figures/Figure_4_analysis1_life_expectancy_shape_relationship.png",
    "figures/Figure_5_analysis2_climate_effects_recruitment.png",
    "figures/Figure_6_analysis3_climate_effects_multisite.png",
    "figures/Figure_S1_boundary_estimate_diagnostic.png",
    "figures/Figure_S2_boundary_survivorship_illustration.png",
    "data/metadata/sources.md",
    ".lintr"
  )

  if (bundle_type == "standard") {
    copy_paths <- c(
      copy_paths,
      "docs/manuscript/manuscript_main.Rmd",
      "docs/manuscript/manuscript_supplement.Rmd",
      "docs/manuscript/render_manuscripts.R",
      "docs/manuscript/styles",
      "docs/tables"
    )
  }

  script_paths <- c(
    "scripts/00_check_setup.R",
    "scripts/99_make_tables.R",
    "scripts/README.md",
    "scripts/S01_compadre_correct.R",
    "scripts/S02_target_studies.R",
    "scripts/S03_sampling_distrib_single_mpm.R",
    "scripts/S04_case_study_1_studies.R",
    "scripts/S05_case_study_1_studies_spp.R",
    "scripts/S06_case_study_1_derived.R",
    "scripts/S07_case_study_1_analysis.R",
    "scripts/S08_case_study_1_analysis_spp.R",
    "scripts/S09_case_study_1_surv_issue.R",
    "scripts/S10_case_study_1_var_comp.R",
    "scripts/S11_case_study_2_prism.R",
    "scripts/S12_case_study_2.R",
    "scripts/S13_supplement.R",
    "scripts/S14_case_study_3_astragalus.R",
    "scripts/prepare_zenodo_bundle.R",
    "scripts/validate_reproduction.R"
  )

  copy_paths <- c(copy_paths, script_paths)

  copy_item <- function(path) {
    if (!file.exists(path)) {
      warning("Skipping missing path: ", path)
      return(invisible(NULL))
    }

    target <- fs::path(bundle_dir, path)
    fs::dir_create(fs::path_dir(target), recurse = TRUE)

    if (fs::is_dir(path)) {
      fs::dir_copy(path, target, overwrite = TRUE)
    } else {
      fs::file_copy(path, target, overwrite = TRUE)
    }
  }

  purrr::walk(copy_paths, copy_item)

  run_script <- c(
    if (bundle_type == "standard") {
      "# Reproduce analyses and manuscripts from this Zenodo bundle ----"
    } else {
      "# Reproduce analyses from this anonymised review bundle ----"
    },
    "source('scripts/00_check_setup.R')",
    "source('scripts/S01_compadre_correct.R')",
    "source('scripts/S02_target_studies.R')",
    "source('scripts/S03_sampling_distrib_single_mpm.R')",
    "source('scripts/S04_case_study_1_studies.R')",
    "source('scripts/S05_case_study_1_studies_spp.R')",
    "source('scripts/S06_case_study_1_derived.R')",
    "source('scripts/S07_case_study_1_analysis.R')",
    "source('scripts/S08_case_study_1_analysis_spp.R')",
    "source('scripts/S09_case_study_1_surv_issue.R')",
    "source('scripts/S10_case_study_1_var_comp.R')",
    "# Optional PRISM download/re-extraction step:",
    "# source('scripts/S11_case_study_2_prism.R')",
    "source('scripts/S12_case_study_2.R')",
    "source('scripts/S13_supplement.R')",
    "source('scripts/99_make_tables.R')"
  )
  if (bundle_type == "standard") {
    run_script <- c(
      run_script,
      "rmarkdown::render(",
      "  'docs/manuscript/manuscript_main.Rmd',",
      "  output_format = 'bookdown::word_document2',",
      "  output_file = 'sampling_uncertainty_mpm_main.docx',",
      "  quiet = TRUE",
      ")",
      "rmarkdown::render(",
      "  'docs/manuscript/manuscript_supplement.Rmd',",
      "  output_format = 'bookdown::word_document2',",
      "  output_file = 'sampling_uncertainty_mpm_supplement.docx',",
      "  quiet = TRUE",
      ")"
    )
  }
  readr::write_lines(run_script, fs::path(bundle_dir, "RUN_REPRODUCTION.R"))

  all_files <- fs::dir_ls(bundle_dir, recurse = TRUE, type = "file")
  relative_paths <- fs::path_rel(all_files, start = bundle_dir)
  file_info <- fs::file_info(all_files)
  checksums <- unname(tools::md5sum(all_files))

  manifest <- tibble::tibble(
    path = as.character(relative_paths),
    size_bytes = as.numeric(file_info$size),
    modified_utc = format(
      as.POSIXct(file_info$modification_time, tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    ),
    md5 = checksums
  ) |>
    dplyr::arrange(path)

  readr::write_csv(manifest, fs::path(bundle_dir, "MANIFEST.csv"))

  manifest_md <- c(
    if (bundle_type == "standard") {
      "# Zenodo bundle manifest"
    } else {
      "# Anonymised review bundle manifest"
    },
    "",
    paste0("- Bundle: `", bundle_name, "`"),
    paste0("- Created (UTC): ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("- File count: ", nrow(manifest)),
    "",
    "## Reproduction",
    "",
    "From the bundle root, run:",
    "",
    "```r",
    "source('RUN_REPRODUCTION.R')",
    "```",
    "",
    "## Integrity",
    "",
    "- `MANIFEST.csv` records path, file size, modification timestamp, and MD5 checksum for each file.",
    "- Recompute checksums and compare to `MANIFEST.csv` before archiving/deposition."
  )
  readr::write_lines(manifest_md, fs::path(bundle_dir, "MANIFEST.md"))

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(staging_root)
  zip_file <- fs::path(archive_root, paste0(bundle_name, ".zip"))
  if (file.exists(zip_file)) {
    fs::file_delete(zip_file)
  }
  zip_status <- system2(
    command = "zip",
    args = c("-r", "-q", "-9", "-X", zip_file, bundle_name)
  )
  if (!identical(zip_status, 0L)) {
    stop("Failed to create Zenodo archive zip.")
  }

  cat("Zenodo bundle prepared at:\n", bundle_dir, "\n", sep = "")
  cat("Zenodo archive zip written to:\n", zip_file, "\n", sep = "")
  invisible(bundle_dir)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  skip_render <- "--skip-render" %in% args
  bundle_type <- if ("--anonymized-review" %in% args) {
    "anonymized_review"
  } else {
    "standard"
  }
  prepare_zenodo_bundle(skip_render = skip_render, bundle_type = bundle_type)
}
