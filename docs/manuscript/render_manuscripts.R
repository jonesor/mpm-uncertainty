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
  "figures/Figure_6_analysis3_climate_effects_multisite.png" = "docs/manuscript/Figure_6_analysis3_climate_effects_multisite.png",
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

# Refresh submission bundle artifacts (cover letters are preserved as-is).
submission_dir <- "docs/manuscript/submission"
if (!dir.exists(submission_dir)) dir.create(submission_dir, recursive = TRUE, showWarnings = FALSE)

main_doc <- "docs/manuscript/sampling_uncertainty_mpm_main.docx"
supp_doc <- "docs/manuscript/sampling_uncertainty_mpm_supplement.docx"
combined_doc <- file.path(submission_dir, "sampling_uncertainty_mpm_combined.docx")

file.copy(main_doc, file.path(submission_dir, basename(main_doc)), overwrite = TRUE)
file.copy(supp_doc, file.path(submission_dir, basename(supp_doc)), overwrite = TRUE)

# Build combined manuscript+supplement DOCX by concatenating Rmd sources and rendering.
strip_yaml <- function(lines) {
  if (length(lines) < 3 || trimws(lines[1]) != "---") {
    return(lines)
  }
  end_idx <- which(trimws(lines[-1]) == "---")[1] + 1
  if (is.na(end_idx)) {
    return(lines)
  }
  if (end_idx >= length(lines)) {
    return(character(0))
  }
  lines[(end_idx + 1):length(lines)]
}

prefix_r_chunk_labels <- function(lines, prefix = "supp-") {
  starts <- grepl("^```\\{r[^}]*\\}$", lines)
  if (!any(starts)) return(lines)
  lines[starts] <- vapply(
    lines[starts],
    FUN.VALUE = character(1),
    function(x) {
      inside <- sub("^```\\{r\\s*", "", x)
      inside <- sub("\\}$", "", inside)
      comma_pos <- regexpr(",", inside, fixed = TRUE)

      if (comma_pos[1] > 0) {
        label <- trimws(substr(inside, 1, comma_pos[1] - 1))
        rest <- substr(inside, comma_pos[1], nchar(inside))
      } else {
        label <- trimws(inside)
        rest <- ""
      }

      if (identical(label, "")) return(x)
      paste0("```{r ", prefix, label, rest, "}")
    }
  )
  lines
}

main_lines <- readLines(main_rmd, warn = FALSE, encoding = "UTF-8")
supp_lines <- readLines(supp_rmd, warn = FALSE, encoding = "UTF-8")
supp_body <- strip_yaml(supp_lines)
supp_body <- prefix_r_chunk_labels(supp_body, prefix = "supp-")

combined_rmd <- file.path("docs", "manuscript", "manuscript_combined.Rmd")
combined_lines <- c(
  main_lines,
  "",
  "\\newpage",
  "",
  supp_body
)
writeLines(combined_lines, combined_rmd, useBytes = TRUE)

combined_out <- tryCatch(
  rmarkdown::render(
    input = combined_rmd,
    output_format = "bookdown::word_document2",
    output_file = basename(combined_doc),
    output_dir = submission_dir,
    quiet = TRUE
  ),
  error = function(e) {
    warning("Could not build combined DOCX from combined Rmd: ", conditionMessage(e))
    NA_character_
  }
)

if (is.na(combined_out) || !file.exists(combined_doc)) {
  warning("Combined DOCX was not produced at: ", combined_doc)
}

zip_path <- file.path(
  "docs", "manuscript", "zenodo", "archive",
  paste0(basename(bundle_dir), ".zip")
)
if (file.exists(zip_path)) {
  file.copy(zip_path, file.path(submission_dir, basename(zip_path)), overwrite = TRUE)
}

# Keep only the latest Zenodo zip in submission to avoid stale large artifacts.
submission_zips <- list.files(
  submission_dir,
  pattern = "^mpm-uncertainty_zenodo_.*\\.zip$",
  full.names = TRUE
)
latest_zip <- file.path(submission_dir, basename(zip_path))
stale_zips <- setdiff(submission_zips, latest_zip)
if (length(stale_zips) > 0) {
  invisible(file.remove(stale_zips))
}

# Remove legacy combined-Rmd copy in submission if present.
legacy_combined_rmd <- file.path(submission_dir, "sampling_uncertainty_mpm_combined.Rmd")
if (file.exists(legacy_combined_rmd)) {
  file.remove(legacy_combined_rmd)
}

main_figs <- c(
  "docs/manuscript/Figure_1_transition_rates_single_mpm.png",
  "docs/manuscript/Figure_2_derived_parameters_single_mpm.png",
  "docs/manuscript/Figure_3_analysis1_point_vs_sampling_distributions.png",
  "docs/manuscript/Figure_4_analysis1_life_expectancy_shape_relationship.png",
  "docs/manuscript/Figure_5_analysis2_climate_effects_recruitment.png",
  "docs/manuscript/Figure_6_analysis3_climate_effects_multisite.png"
)
for (fig in main_figs) {
  if (file.exists(fig)) file.copy(fig, file.path(submission_dir, basename(fig)), overwrite = TRUE)
}

desc_path <- file.path(submission_dir, "zenodo_zip_description.txt")
desc_lines <- c(
  "Zenodo archive description",
  "==========================",
  "",
  "Archive file",
  "------------",
  paste0("- Filename: ", basename(zip_path)),
  "- Purpose: Reproducible snapshot of data, code, models, and manuscript sources/outputs used in the paper.",
  "",
  "What the zip contains",
  "---------------------",
  "Top-level project folders included in the archive:",
  "- apps/                  Shiny app source code.",
  "- code/                  Reusable R helper functions.",
  "- data/                  Raw, derived, and metadata files used by the analyses.",
  "- docs/                  Manuscript source files and analysis documentation.",
  "- figures/               Figure outputs generated by scripts.",
  "- models/                Stan model files.",
  "- scripts/               End-to-end analysis scripts (S01-S14).",
  "",
  "Top-level files included in the archive:",
  "- README.md              Repository overview and run order.",
  "- .gitignore             Ignore rules used in the project.",
  "- .lintr                 Lint configuration.",
  "- mpm-uncertainty.Rproj  RStudio project file.",
  "",
  "Reproducibility notes",
  "---------------------",
  "1) Open the project in R/RStudio from the extracted archive root.",
  "2) Run scripts in order (S01 to S14) for full regeneration of analysis outputs.",
  "3) Render manuscript files via:",
  "   source(\"docs/manuscript/render_manuscripts.R\")",
  "4) This command renders the main manuscript and supplement and rebuilds a Zenodo staging bundle.",
  "",
  "Data provenance",
  "---------------",
  "- Matrix data from COMPADRE and associated study-level sources.",
  "- Climate covariates from PRISM (downloaded by project scripts where required).",
  "- Provenance and source notes are documented in data/metadata/."
)
writeLines(desc_lines, desc_path, useBytes = TRUE)

cat("\nSubmission bundle refreshed at:\n")
cat("- ", submission_dir, "\n", sep = "")
