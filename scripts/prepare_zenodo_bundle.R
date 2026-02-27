# Prepare a reproducible Zenodo deposition bundle ----
#
# This script creates a timestamped snapshot in
# docs/manuscript/zenodo/staging/ that includes code, data, models,
# generated outputs, manuscript files, and a manifest.
#
# RStudio usage:
#   source('scripts/prepare_zenodo_bundle.R')
#   prepare_zenodo_bundle(skip_render = FALSE)
#
# CLI usage:
#   Rscript --vanilla scripts/prepare_zenodo_bundle.R
#   Rscript --vanilla scripts/prepare_zenodo_bundle.R --skip-render

prepare_zenodo_bundle <- function(skip_render = FALSE, clean_old_staging = TRUE) {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required.")
  }

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
      regexp = "mpm-uncertainty_zenodo_",
      type = "directory"
    )
    if (length(old_bundles) > 0) {
      purrr::walk(old_bundles, ~unlink(.x, recursive = TRUE, force = TRUE))
    }
  }

  bundle_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  bundle_name <- paste0("mpm-uncertainty_zenodo_", bundle_timestamp)
  bundle_dir <- fs::path(staging_root, bundle_name)

  fs::dir_create(bundle_dir, recurse = TRUE)

  if (!skip_render) {
    source("docs/manuscript/render_manuscripts.R")
  }

  copy_paths <- c(
    "README.md",
    "code",
    "scripts",
    "models",
    "data",
    "figures",
    "docs/manuscript/manuscript_main.Rmd",
    "docs/manuscript/sampling_uncertainty_mpm_main.docx",
    "docs/manuscript/manuscript_supplement.Rmd",
    "docs/manuscript/sampling_uncertainty_mpm_supplement.docx",
    "docs/manuscript/render_manuscripts.R",
    "docs/manuscript/styles",
    "docs/tables",
    "docs/analysis_summary.md",
    "data/metadata/sources.md",
    ".lintr"
  )

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
    "# Reproduce analyses and manuscripts from this Zenodo bundle ----",
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
    "source('scripts/S11_case_study_2_prism.R')",
    "source('scripts/S12_case_study_2.R')",
    "source('scripts/S13_supplement.R')",
    "source('scripts/99_make_tables.R')",
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
    "# Zenodo bundle manifest",
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
  prepare_zenodo_bundle(skip_render = skip_render)
}
