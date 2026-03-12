# Check case study 2 outputs using the latest local COMPADRE snapshot ----
#
# Usage:
# 1) Put latest COMPADRE file in data/raw/compadre/, e.g.
#    COMPADRE_6.25.8.0_Aug_14_2025.RData
# 2) source("scripts/check_case2_with_latest_compadre.R")
#
# This script:
# - backs up current active COMPADRE files
# - activates the latest local COMPADRE snapshot
# - reruns S01 + case-study-2 scripts (S02, S12)
# - compares key case-study-2 outputs to baseline
# - restores original active COMPADRE files

source("code/setup.R")
setup_packages(c("readr", "dplyr", "stringr", "purrr", "fs", "tibble", "tools"))

raw_dir <- "data/raw/compadre"
active_raw <- file.path(raw_dir, "COMPADRE_v.X.X.X.RData")
active_corrected <- file.path(raw_dir, "COMPADRE_v.X.X.X_Corrected.RData")

# Optional manual override. Set this to a file path if you want to force
# which snapshot is used for the comparison.
latest_file <- "data/raw/compadre/COMPADRE_v.X.X.X_backup_20260310_125942.RData"

if (is.null(latest_file)) {
  latest_candidates <- fs::dir_ls(raw_dir, regexp = "COMPADRE_[^v].*\\.RData$", type = "file")
  if (length(latest_candidates) == 0) {
    stop(
      "No latest COMPADRE snapshot found in data/raw/compadre (expected COMPADRE_*.RData). ",
      "Set `latest_file` manually in scripts/check_case2_with_latest_compadre.R."
    )
  }
  latest_file <- latest_candidates[length(latest_candidates)]
} else {
  if (!file.exists(latest_file)) {
    stop("Configured latest_file does not exist: ", latest_file)
  }
}
message("Using latest snapshot: ", latest_file)

if (!file.exists(active_raw) || !file.exists(active_corrected)) {
  stop("Active COMPADRE files not found (expected COMPADRE_v.X.X.X*.RData).")
}

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
backup_raw <- file.path(raw_dir, paste0("COMPADRE_v.X.X.X_pre_case2_check_", ts, ".RData"))
backup_corrected <- file.path(raw_dir, paste0("COMPADRE_v.X.X.X_Corrected_pre_case2_check_", ts, ".RData"))
file.copy(active_raw, backup_raw, overwrite = TRUE)
file.copy(active_corrected, backup_corrected, overwrite = TRUE)

on.exit({
  file.copy(backup_raw, active_raw, overwrite = TRUE)
  file.copy(backup_corrected, active_corrected, overwrite = TRUE)
  message("Restored active COMPADRE files from backups.")
}, add = TRUE)

compare_dir <- file.path("/tmp", paste0("case2_compare_", ts))
fs::dir_create(compare_dir, recurse = TRUE)
baseline_dir <- file.path(compare_dir, "baseline")
latest_dir <- file.path(compare_dir, "latest")
fs::dir_create(baseline_dir, recurse = TRUE)
fs::dir_create(latest_dir, recurse = TRUE)

key_outputs <- c(
  "data/derived/analysis_cache/case2_gprc_beta_summary.csv",
  "data/derived/analysis_cache/case2_spring_beta_summary.csv",
  "data/derived/analysis_cache/case2_stage_survival_weight_sensitivity.csv",
  "figures/Figure_5_analysis2_climate_effects_recruitment.png",
  "data/derived/studies/target_studies.csv"
)

for (f in key_outputs) {
  if (file.exists(f)) {
    fs::dir_create(file.path(baseline_dir, dirname(f)), recurse = TRUE)
    file.copy(f, file.path(baseline_dir, f), overwrite = TRUE)
  }
}

file.copy(latest_file, active_raw, overwrite = TRUE)

scripts_to_run <- c(
  "scripts/S01_compadre_correct.R",
  "scripts/S02_target_studies.R",
  "scripts/S12_case_study_2.R"
)

log_dir <- file.path(compare_dir, "logs")
fs::dir_create(log_dir, recurse = TRUE)

for (s in scripts_to_run) {
  log_file <- file.path(log_dir, paste0(basename(tools::file_path_sans_ext(s)), ".log"))
  message("Running ", s)
  out <- system2("Rscript", c("--vanilla", s), stdout = TRUE, stderr = TRUE)
  readr::write_lines(out, log_file)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("Script failed: ", s, " (see ", log_file, ")")
  }
}

for (f in key_outputs) {
  if (file.exists(f)) {
    fs::dir_create(file.path(latest_dir, dirname(f)), recurse = TRUE)
    file.copy(f, file.path(latest_dir, f), overwrite = TRUE)
  }
}

cmp_rows <- purrr::map_dfr(key_outputs, function(f) {
  old_path <- file.path(baseline_dir, f)
  new_path <- file.path(latest_dir, f)
  if (!file.exists(old_path) || !file.exists(new_path)) {
    return(tibble::tibble(file = f, status = "missing"))
  }

  if (grepl("\\.png$", f)) {
    return(tibble::tibble(
      file = f,
      md5_old = unname(tools::md5sum(old_path)),
      md5_new = unname(tools::md5sum(new_path)),
      identical = identical(unname(tools::md5sum(old_path)), unname(tools::md5sum(new_path))),
      status = "ok"
    ))
  }

  old <- readr::read_csv(old_path, show_col_types = FALSE)
  new <- readr::read_csv(new_path, show_col_types = FALSE)
  tibble::tibble(
    file = f,
    rows_old = nrow(old),
    rows_new = nrow(new),
    cols_old = ncol(old),
    cols_new = ncol(new),
    names_identical = identical(names(old), names(new)),
    md5_old = unname(tools::md5sum(old_path)),
    md5_new = unname(tools::md5sum(new_path)),
    identical = identical(unname(tools::md5sum(old_path)), unname(tools::md5sum(new_path))),
    status = "ok"
  )
})

report_file <- file.path(compare_dir, "case2_latest_compadre_check.csv")
readr::write_csv(cmp_rows, report_file)

message("Case study 2 COMPADRE check complete.")
message("Compare dir: ", compare_dir)
message("Summary report: ", report_file)
