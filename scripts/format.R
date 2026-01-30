# Format R code using styler

args <- commandArgs(trailingOnly = TRUE)
paths <- if (length(args) > 0) args else c("scripts", "code")

if (!requireNamespace("styler", quietly = TRUE)) {
  stop(
    "Missing package: styler. Install with install.packages(\"styler\").",
    call. = FALSE
  )
}

styler::style_dir(path = paths, recursive = TRUE)
